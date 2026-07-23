#!/usr/bin/env python3
"""Reject changes outside the documentation-only readability contract.

Run this script from either repository checkout (or pass ``--repo``). It uses
only the Python standard library except for PyYAML when YAML files are present.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote


FORBIDDEN_SUFFIXES = {
    ".h5ad", ".h5", ".hdf5", ".loom", ".mtx", ".rds", ".rdata",
    ".ckpt", ".pt", ".pth", ".pkl", ".pickle", ".joblib", ".npy",
    ".npz", ".pyc",
}
NOTEBOOK_AREAS = (
    "notebooks/benchmarks/",
    "notebooks/public_validation/",
    "notebooks/revision_finalization/",
)
PRESENTATION_PYTHON_FILES = {
    "scripts/assemble_gse211713_manuscript_outputs.py",
    "scripts/assemble_gse249479_manuscript_outputs.py",
    "scripts/assemble_reviewer2_repairs.py",
    "scripts/execute_gse249479_validation.py",
    "scripts/execute_pancreas_validation.py",
    "scripts/finalize_revision_evidence.py",
    "scripts/gse249479_memory_safe.py",
    "scripts/gse249479_qc_common.py",
    "scripts/gse249479_scgeo_common.py",
    "scripts/pancreas_validation_common.py",
}


def git(repo: Path, *args: str, check: bool = True) -> str:
    result = subprocess.run(
        ["git", *args], cwd=repo, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if check and result.returncode:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    return result.stdout


def changed_paths(repo: Path, base: str) -> dict[str, str]:
    changes: dict[str, str] = {}
    for line in git(repo, "diff", "--name-status", "--find-renames", base).splitlines():
        fields = line.split("\t")
        status = fields[0]
        if status.startswith("R"):
            changes[fields[1]] = "D"
            changes[fields[2]] = "A"
        else:
            changes[fields[-1]] = status[0]
    for path in git(repo, "ls-files", "--others", "--exclude-standard").splitlines():
        changes[path] = "A"
    return changes


def show_json(repo: Path, revision: str, path: str) -> dict:
    raw = git(repo, "show", f"{revision}:{path}")
    return json.loads(raw)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def notebook_signature(nb: dict) -> dict:
    cells = nb.get("cells", [])
    code = [cell for cell in cells if cell.get("cell_type") == "code"]
    return {
        "cell_count": len(cells),
        "cell_types": [cell.get("cell_type") for cell in cells],
        "code_sources": [cell.get("source", []) for cell in code],
        "code_outputs": [cell.get("outputs", []) for cell in code],
        "execution_counts": [cell.get("execution_count") for cell in code],
    }


def validate_notebook_change(repo: Path, base: str, rel: str, errors: list[str]) -> None:
    before = notebook_signature(show_json(repo, base, rel))
    after = notebook_signature(read_json(repo / rel))
    for key in ("cell_count", "cell_types", "code_sources", "code_outputs", "execution_counts"):
        if before[key] != after[key]:
            errors.append(f"{rel}: notebook {key} changed")


class _StringNeutralizer(ast.NodeTransformer):
    def visit_Constant(self, node: ast.Constant) -> ast.AST:  # noqa: N802
        if isinstance(node.value, str):
            return ast.copy_location(ast.Constant(value="<display-text>"), node)
        return node


def validate_python_display_change(repo: Path, base: str, rel: str, errors: list[str]) -> None:
    """Allow explanatory string/comment edits while rejecting structural Python changes."""
    try:
        before = ast.parse(git(repo, "show", f"{base}:{rel}"), filename=rel)
        after = ast.parse((repo / rel).read_text(encoding="utf-8"), filename=rel)
    except Exception as exc:  # noqa: BLE001
        errors.append(f"{rel}: Python parse failed ({exc})")
        return
    before = _StringNeutralizer().visit(before)
    after = _StringNeutralizer().visit(after)
    if ast.dump(before, include_attributes=False) != ast.dump(after, include_attributes=False):
        errors.append(f"{rel}: change is not limited to comments or string display text")


def validate_all_notebooks(repo: Path, errors: list[str]) -> int:
    count = 0
    for rel in git(repo, "ls-files", "*.ipynb").splitlines():
        count += 1
        try:
            nb = read_json(repo / rel)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{rel}: invalid notebook JSON ({exc})")
            continue
        for index, cell in enumerate(nb.get("cells", [])):
            if cell.get("cell_type") != "code":
                continue
            if cell.get("outputs", []):
                errors.append(f"{rel}: code cell {index} contains outputs")
            if cell.get("execution_count") is not None:
                errors.append(f"{rel}: code cell {index} has an execution count")
    return count


def validate_structured_files(repo: Path, errors: list[str]) -> tuple[int, int]:
    json_count = 0
    yaml_count = 0
    tracked = git(repo, "ls-files").splitlines()
    for rel in tracked:
        path = repo / rel
        if rel.endswith((".json", ".ipynb")):
            json_count += 1
            try:
                json.loads(path.read_text(encoding="utf-8"))
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{rel}: JSON parse failed ({exc})")
        elif rel.endswith((".yaml", ".yml")):
            yaml_count += 1
            try:
                import yaml  # type: ignore
                yaml.safe_load(path.read_text(encoding="utf-8"))
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{rel}: YAML parse failed ({exc})")
    return json_count, yaml_count


LINK_RE = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)|!\[[^\]]*\]\(([^)]+)\)")


def validate_markdown_links(repo: Path, errors: list[str]) -> int:
    checked = 0
    for rel in git(repo, "ls-files", "*.md").splitlines():
        path = repo / rel
        text = path.read_text(encoding="utf-8")
        for match in LINK_RE.finditer(text):
            raw = (match.group(1) or match.group(2)).strip()
            if raw.startswith("<") and raw.endswith(">"):
                raw = raw[1:-1]
            target = raw.split()[0].split("#", 1)[0]
            if not target or re.match(r"^[a-z][a-z0-9+.-]*:", target, re.I):
                continue
            checked += 1
            resolved = (path.parent / unquote(target)).resolve()
            try:
                resolved.relative_to(repo.resolve())
            except ValueError:
                errors.append(f"{rel}: local link leaves repository: {raw}")
                continue
            if not resolved.exists():
                errors.append(f"{rel}: broken local link: {raw}")
    return checked


def is_package_repo(repo: Path) -> bool:
    return (repo / "scgeo").is_dir() and (repo / "tests").is_dir()


def allowed_package_path(rel: str) -> bool:
    return rel == "README.md" or rel.startswith("docs/") or rel in {
        "CONTRIBUTING.md", "CITATION.cff"
    }


def allowed_notebook_path(rel: str, status: str) -> bool:
    if rel == "README.md" or rel.startswith("docs/"):
        return True
    if rel == "scripts/validate_readability_only.py":
        return True
    if rel in PRESENTATION_PYTHON_FILES:
        return True
    if rel.endswith("/README.md"):
        return True
    if rel.endswith(".ipynb") and rel.startswith(NOTEBOOK_AREAS):
        return True
    if status == "D" and ("/__pycache__/" in rel or rel.endswith(".pyc")):
        return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--base", default="main")
    args = parser.parse_args()
    repo = args.repo.resolve()
    base = git(repo, "rev-parse", "--verify", args.base).strip()
    mode = "package" if is_package_repo(repo) else "notebooks"
    changes = changed_paths(repo, args.base)
    errors: list[str] = []

    for rel, status in sorted(changes.items()):
        allowed = allowed_package_path(rel) if mode == "package" else allowed_notebook_path(rel, status)
        if not allowed:
            errors.append(f"{rel}: {status} is outside the readability-only allowlist")
        if rel.startswith(("configs/", "tests/")):
            errors.append(f"{rel}: config or test changed")
        if mode == "package" and rel.startswith("scgeo/"):
            errors.append(f"{rel}: package implementation changed")
        if rel.startswith("results/") or Path(rel).suffix.lower() in FORBIDDEN_SUFFIXES:
            if not (status == "D" and rel.endswith(".pyc")):
                errors.append(f"{rel}: numerical/model/data/cache artifact changed or added")
        if mode == "notebooks" and rel.endswith(".ipynb") and status != "D":
            if status == "A":
                errors.append(f"{rel}: adding notebooks is not allowed in this pass")
            else:
                validate_notebook_change(repo, args.base, rel, errors)
        if mode == "notebooks" and rel in PRESENTATION_PYTHON_FILES and status != "D":
            validate_python_display_change(repo, args.base, rel, errors)

    tracked_cache = [
        rel for rel in git(repo, "ls-files").splitlines()
        if "/__pycache__/" in rel or rel.endswith(".pyc")
    ]
    for rel in tracked_cache:
        if changes.get(rel) != "D":
            errors.append(f"{rel}: tracked Python cache remains")

    notebook_count = validate_all_notebooks(repo, errors) if mode == "notebooks" else 0
    json_count, yaml_count = validate_structured_files(repo, errors) if mode == "notebooks" else (0, 0)
    link_count = validate_markdown_links(repo, errors)

    report = {
        "repository": str(repo),
        "mode": mode,
        "base": base,
        "changed_paths": len(changes),
        "tracked_notebooks_checked": notebook_count,
        "json_files_checked": json_count,
        "yaml_files_checked": yaml_count,
        "local_markdown_links_checked": link_count,
        "errors": errors,
        "status": "passed" if not errors else "failed",
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
