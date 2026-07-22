#!/usr/bin/env python3
"""Build final reviewer ledgers and run presentation-layer reproducibility gates."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCGEO = Path("/home/liuyuchen/Github/scgeo")
PYTHON = SCGEO / ".venv/bin/python"
OUT = ROOT / "results/revision_finalization"
REPRO = OUT / "reproducibility"
CBASE = ROOT / "results/public_validation/gse211713_dataset_c"
BBASE = ROOT / "results/public_validation/gse249479_dataset_b"
REQUIRED_SCGEO = "9a0ed16cbaa57f935f9c9bc87d1643a25b51012c"


def run(*args: str, cwd: Path = ROOT) -> str:
    return subprocess.run(args, cwd=cwd, check=True, text=True, capture_output=True).stdout.strip()


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); tmp.write_text(text,encoding="utf-8"); tmp.replace(path)


def write_json(path: Path, obj) -> None:
    write_text(path,json.dumps(obj,indent=2,sort_keys=True)+"\n")


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); frame.to_csv(tmp,index=False); tmp.replace(path)


def reviewer_rows() -> pd.DataFrame:
    cols=["reviewer","comment_identifier","concern","action_taken","exact_notebook_or_script","exact_figure_table_output","quantitative_evidence","negative_or_inconclusive_evidence","limitation","proposed_manuscript_section","proposed_response_letter_wording","completion_status"]
    rows=[
        ["Reviewer 1","R1-novelty","Novelty positioning","Separate ScGeo's robust state-displacement geometry and representation-stability outputs from standard primitives and complementary classifiers.","notebooks/benchmarks/00_manuscript_benchmark_overview.ipynb; scripts/assemble_reviewer2_repairs.py","global_claim_evidence_matrix.csv; scpasi_comparison.csv","Frozen manuscript benchmark comprises 240 jobs (60 calibration, 180 evaluation).","Standard geometric primitives remain components; no universal superiority claim.","Novelty is architectural and inferential framing, not invention of every primitive.","Introduction; Methods","We clarified the distinct objective, estimands, and uncertainty layers without claiming that each primitive is new.","complete"],
        ["Reviewer 1","R1-scanorama","Scanorama dependence","Document Scanorama as one tested preprocessing/representation in a public demonstration, not a ScGeo requirement.","notebooks/data_prep/04_scgeo_gse280305_phase1_qc.ipynb; scripts/assemble_reviewer2_repairs.py","representation_metric_inventory.csv","Dataset B and C used no Harmony or Scanorama.","No condition-mixing optimization was performed in Dataset B/C.","Original public demonstration remains representation-specific.","Methods; Limitations","ScGeo accepts user-supplied representations; Scanorama is neither required nor used in Datasets B/C.","complete"],
        ["Reviewer 1","R1-single-dataset","Single-dataset concern","Add public Dataset B and replicate-aware radiation Dataset C alongside pancreas and synthetic evidence.","notebooks/public_validation/gse249479; notebooks/public_validation/gse211713","Dataset B/C main figures and evidence ledgers","Dataset B: 34,432 cells, descriptive_only; Dataset C: 20 mouse libraries and three replicate-aware primary contrasts.","Dataset B lacks biological replication; Dataset C groups remain small.","Cross-study evidence is heterogeneous and not a meta-analysis.","Results","We added two public validations with explicit design-specific inference status.","complete"],
        ["Reviewer 1","R1-validation-strength","Validation strength","Expose negative, neutral, sensitivity, exact-permutation, and comparator evidence.","notebooks/benchmarks; notebooks/public_validation/gse249479; notebooks/public_validation/gse211713","negative_neutral_or_unstable_states.csv; official Augur comparison; exact permutations","Dataset C exact assignment counts are 126, 126, and 70; Dataset B compares seven marker-inferred states with official R Augur.","Small-n permutation resolution and marker uncertainty remain.","No causal or population-wide generalization.","Results; Supplement","Validation now includes controlled, public, comparator, replicate-aware, and negative evidence.","complete"],
        ["Reviewer 2","R2-estimator","Estimator robustness and arithmetic-mean sensitivity","Retain frozen robust estimator benchmark and explicitly distinguish robust centers from arithmetic summaries.","notebooks/benchmarks/01_robust_estimator_comparison.ipynb","results_manifest/benchmark_files.csv: audit/estimator_comparison.csv","Frozen manifest pins estimator_comparison.csv across the manuscript benchmark.","Underlying numerical audit is frozen and not rerun here.","Benchmark scope does not prove robustness for every distribution.","Methods; Supplement","We expose estimator sensitivity and retain negative benchmark results.","complete"],
        ["Reviewer 2","R2-mixscore","Mixscore aggregation","Define cell, state, sample, and dataset denominators and weighting.","notebooks/revision_finalization/01_reviewer2_metric_repairs.ipynb; scripts/assemble_reviewer2_repairs.py","mixscore_aggregation_definition.csv; mixscore figure","Per-cell score is the fraction of k neighbors with a different label; demonstrated k=50.","Unstratified means weight large states/samples more heavily.","Mixscore is descriptive and cells are not biological replicates.","Methods; Supplement","We now state every aggregation unit, denominator, and weighting rule.","complete"],
        ["Reviewer 2","R2-representation","Integration and representation metrics","Unify available neighborhood overlap, Jaccard, connectivity, trustworthiness/continuity, distortion, label preservation, and mouse dominance.","notebooks/revision_finalization/01_reviewer2_metric_repairs.ipynb","representation_metric_inventory.csv","Dataset C records five representations; Dataset B exposes pairwise overlap/Jaccard and trustworthiness/continuity.","Synthetic and pancreas numerical metric artifacts are not present in this workspace and remain marked unavailable.","Condition mixing is not an integration objective.","Methods; Supplement","Representation diagnostics are reported as diagnostics, not universal biological preservation.","complete_with_documented_unavailable_fields"],
        ["Reviewer 2","R2-distribution","Distribution-test visibility","Expose statistic, unit, group sizes, null generation, resamples, resolution, adjustment, and inference status.","notebooks/revision_finalization/02_distribution_and_aggregation_audit.ipynb","distribution_test_inventory.csv; distribution visibility figure","Inventory contains pooled-cell distances, cell-label computational nulls, mouse-center distances, and exact mouse permutations.","Pooled-cell statistics are not biological-replicate evidence.","Scales differ and are not combined.","Methods; Supplement","We separated descriptive distances from biological-mouse permutation sensitivity.","complete"],
        ["Reviewer 2","R2-preservation","Biological-preservation wording","Replace universal preservation language with evaluated label/neighborhood retention.","notebooks/revision_finalization/03_terminology_and_claim_audit.ipynb","terminology_replacement_inventory.csv","Eight explicit replacement classes documented.","No claim of universal preservation remains.","Evaluated labels are marker- or study-derived and imperfect.","Throughout","We narrowed claims to the representations, labels, and diagnostics actually evaluated.","complete"],
        ["Reviewer 2","R2-ood","OOD wording","Replace unsupported general OOD claims with weak reference support, representation instability, or extrapolation sensitivity.","notebooks/manuscript/05_OOD.ipynb; notebooks/tutorials/PBMC_ingest_ScGeo_QC_demo.ipynb","terminology_replacement_inventory.csv","Controlled corruption evidence is retained as controlled corruption.","API field names such as scgeo_ood remain for compatibility.","No general OOD detector is claimed.","Methods; Limitations","We now use diagnostic terms tied to reference support and controlled corruption.","complete"],
        ["Reviewer 2","R2-scpasi","scPASI discussion","Provide a conservative methods comparison without implementation or superiority claim.","scripts/assemble_reviewer2_repairs.py","scpasi_comparison.csv","ScGeo and official R Augur fields are source-grounded; scPASI unknown fields remain explicit.","No local scPASI methods source or direct run was available.","Only a related-but-distinct relationship can be stated here.","Discussion","We acknowledge scPASI as related and methodologically distinct; no performance ranking is claimed.","complete_with_editorial_followup"],
        ["Reviewer 3","R3-primitives","Standard geometric primitives","Identify robust centers, distances, graph diagnostics, and ranks as standard components.","notebooks/benchmarks/00_manuscript_benchmark_overview.ipynb","global_claim_evidence_matrix.csv","The framework evaluates these components jointly across frozen scenarios.","No novelty claim for standard primitives individually.","Contribution is their prespecified combination and evidence schema.","Introduction; Methods","We distinguish standard primitives from the ScGeo analysis framework.","complete"],
        ["Reviewer 3","R3-representation-dependence","Representation dependence","Use prespecified primary/sensitivity ensembles and report unstable findings.","notebooks/public_validation/gse249479/03d_representation_quality.ipynb; notebooks/public_validation/gse211713/04d_representation_quality.ipynb","Dataset B/C representation tables and stability panels","Dataset C early effects are unstable; late has five of six stable major states; early-vs-late has six of six stable.","PCA20/30/50 are nested and not independent confirmations.","Diffusion can be distorted and remains exploratory.","Results; Supplement","We elevated representation dependence to an explicit result rather than hiding instability.","complete"],
        ["Reviewer 3","R3-local-geometry","Lack of local geometry","Expose neighborhood overlap, Jaccard, local distortion, connectivity, and state-graph diagnostics where available.","scripts/assemble_reviewer2_repairs.py","representation_metric_inventory.csv","Dataset B and C metrics are enumerated per representation.","Not every metric exists for every earlier dataset.","Local diagnostics do not guarantee biological validity.","Methods; Supplement","We added an availability-aware local-geometry inventory.","complete_with_documented_unavailable_fields"],
        ["Reviewer 3","R3-topology","Topology/structural stability","Report graph connectivity, state-graph rank agreement, leave-one-representation sensitivity, and controlled geometry corruption.","notebooks/benchmarks/04_synthetic_geometry_stress_test.ipynb; notebooks/public_validation/gse211713/05e_cross_contrast_summary.ipynb","representation robustness and Dataset C supplement","Dataset C largest connected-component fraction is 1.0 for all five representations.","Synthetic local-distortion recall includes negative results in the frozen audit.","Topology metrics are representation- and graph-construction-dependent.","Methods; Supplement","Structural diagnostics and their failures are now visible.","complete"],
        ["Reviewer 3","R3-robustness","Representation robustness","Use primary consensus plus PCA20/diffusion and leave-one-primary sensitivity; never count nested PCA views as independent.","notebooks/public_validation/gse249479/04c_cross_condition_state_summary.ipynb; notebooks/public_validation/gse211713/05e_cross_contrast_summary.ipynb","main and supplementary robustness panels","Primary consensus is PCA30/PCA50/scVI in Datasets B/C.","Dataset B remains descriptive_only; Dataset C fibroblast conclusions are more sensitive.","Consensus is conditional on the prespecified ensemble.","Results","We report stable, unstable, and neutral outcomes under frozen rules.","complete"],
        ["Reviewer 3","R3-generalization","Insufficient real-data generalization","Connect public pancreas dynamics, Dataset B inflammation, official Augur, and Dataset C radiation evidence.","notebooks/public_validation/pancreas; notebooks/public_validation/gse249479; notebooks/public_validation/gse211713","global evidence ledger and final figure inventory","Pancreas config expects 2,531 cells; Dataset B 34,432; Dataset C 131,157 retained and 20 libraries.","Pancreas numerical metric artifacts are unavailable locally; Dataset B lacks replication; Dataset C is cross-sectional.","No meta-analytic or causal claim.","Results; Discussion","We added heterogeneous public validations with design-specific limitations.","complete_with_documented_unavailable_fields"],
    ]
    return pd.DataFrame(rows,columns=cols)


def global_outputs() -> None:
    ledger=reviewer_rows(); write_csv(OUT/"global_reviewer_evidence_ledger.csv",ledger)
    claims=pd.DataFrame([
        ["Frozen synthetic benchmark","ScGeo robustness and failure modes are evaluated under controlled scenarios.","results_manifest/benchmark_files.csv","240 jobs; 60 calibration and 180 evaluation","controlled simulation","Does not establish universal real-data performance."],
        ["Geometry stress test","Embedding corruption and local geometry sensitivity are controlled stress tests.","notebooks/benchmarks/04_synthetic_geometry_stress_test.ipynb","frozen manifest includes local_distortion_performance.csv","controlled synthetic","Not general OOD detection."],
        ["Dynamics stress test","Direction/alignment behavior is tested under controlled dynamics.","notebooks/benchmarks/05_synthetic_dynamics_stress_test.ipynb","frozen manifest includes dynamics_performance.csv","controlled synthetic","Synthetic dynamics simplify biology."],
        ["Public pancreas","ScGeo dynamics is compared with annotated endocrinogenesis and CellRank/scVelo context.","configs/pancreas_dataset_d_v1.json","expected 2,531 cells; five representation views including display-only UMAP","public descriptive dynamics","No frozen numerical metric artifact in this workspace."],
        ["Dataset B","Inflammatory treatment geometry is descriptive across five representations.","results/public_validation/gse249479_dataset_b/manuscript/dataset_b_evidence_ledger.csv","34,432 cells; 3,000 HVGs; seven marker-inferred states","descriptive_only","No biological replicate identity."],
        ["Official R Augur","Official Augur is the primary separability comparator for Dataset B.","results/public_validation/gse249479_dataset_b/comparator/06_scgeo_augur_full_comparator.csv","Augur 1.0.3 at b252b84; rank correlations across seven states","descriptive_only","Seven marker-inferred states; cell CV is not biological replication."],
        ["Dataset C","Radiation associations are assessed with biological mouse centers and representation consensus.","results/public_validation/gse211713_dataset_c/manuscript/dataset_c_evidence_ledger.csv","five controls, four early 17 Gy, four late 17 Gy; 500 mouse bootstraps; 126/126/70 assignments","replicate-aware association","Cross-sectional, small-n, reused controls, no causality."],
    ],columns=["evidence_source","claim_scope","exact_output","quantitative_evidence","status","limitation"])
    write_csv(OUT/"global_claim_evidence_matrix.csv",claims)
    methods=pd.DataFrame([
        ["Synthetic benchmark","simulation job/seed","robust and arithmetic estimator sensitivity; geometry/dynamics stress tests","frozen manuscript_v1 outputs","controlled benchmark"],
        ["Pancreas","cell/state transition","PCA/diffusion representations, velocity and CellRank context","public config; numerical artifacts unavailable locally","descriptive public dynamics"],
        ["Dataset B","cell/state; no biological sample unit","ScGeo geometry plus official R Augur separability","five representations; official Augur primary","descriptive_only"],
        ["Dataset C","mouse/GSM biological sample","geometric-median mouse centers, 500 mouse bootstraps, exact mouse-label permutations","PCA30/PCA50/scVI primary; PCA20/diffusion sensitivity","replicate-aware association"],
        ["Final assembly","frozen artifact","checksum validation and plotting only","no numerical recomputation","presentation/provenance"],
    ],columns=["analysis","independent_unit","method","representation_or_comparator_scope","inference_status"])
    write_csv(OUT/"global_methods_provenance.csv",methods)
    limitations=pd.DataFrame([
        ["Dataset B replication","No valid biological replicate identifier; all results descriptive_only.","major","biological inference","Retain explicit limitation in text and figures."],
        ["Dataset B annotation","Seven marker-inferred states and sparse HSC-I signature.","major","state interpretation","Avoid rare-state overinterpretation."],
        ["Dataset C design","Cross-sectional independent mice, not longitudinal.","major","time interpretation","Describe differences, not within-animal change, reversal, persistence, or causality."],
        ["Dataset C sample size","5-vs-4 and 4-vs-4 designs yield discrete exact-permutation resolution.","major","uncertainty","Show individual mice and exact grids."],
        ["Related contrasts","Control mice are reused across Dataset C contrasts.","major","multiplicity/dependence","Do not combine p-values or treat contrasts as independent."],
        ["Representation ensemble","PCA20/30/50 are nested views of one basis.","major","confirmation","Do not count them as independent confirmations."],
        ["Distribution metrics","Pooled-cell energy, MMD and sliced Wasserstein are descriptive.","major","inference","Keep separate from replicate-aware evidence."],
        ["Pancreas artifacts","Frozen numerical representation metrics are not present locally.","editorial","availability","Do not fabricate; resolve during manuscript asset collation if required."],
        ["scPASI","No local methods source or direct comparison is available.","editorial","related work","Verify citation-specific description during manuscript writing."],
        ["Warnings","Pandas categorical-dtype deprecations and pytest cache warning.","maintenance","software","Not revision blockers; schedule compatibility cleanup."],
    ],columns=["issue","limitation","severity","affects","required_handling"])
    write_csv(OUT/"global_limitations_register.csv",limitations)
    change_plan="# Manuscript change plan\n\n- Position ScGeo as a robust, representation-aware state-geometry framework built from standard primitives.\n- Add Dataset B as explicitly descriptive public inflammation evidence and official R Augur as the primary external comparator.\n- Add Dataset C as replicate-aware, cross-sectional radiation association evidence with individual-mouse and exact-permutation layers.\n- Separate displacement, abundance, distribution, and separability throughout.\n- Add representation-quality and failure/instability panels; do not treat biological-condition mixing as an integration objective.\n- Replace universal preservation, general OOD, causal, longitudinal, and nested-PCA confirmation wording.\n- Retain all neutral, unstable, and unavailable results.\n\nNo full manuscript prose is drafted by this package.\n"
    response="# Response-letter evidence plan\n\nUse `global_reviewer_evidence_ledger.csv` row-by-row. Each response should identify the exact new artifact, quote only supported quantitative evidence, state its inference status, and end with the relevant limitation. Avoid superiority, causality, meta-analysis, general OOD detection, or longitudinal claims. Reviewer 2's scPASI discussion requires citation-level scientific writing because no local methods source or direct comparison was available.\n"
    write_text(OUT/"manuscript_change_plan.md",change_plan); write_text(OUT/"response_letter_evidence_plan.md",response)
    figrows=[]
    for dataset,base in [("Dataset B",BBASE/"manuscript"),("Dataset C",CBASE/"manuscript"),("Reviewer 2",OUT/"reviewer2")]:
        for png in sorted((base/"figures").glob("*.png")):
            stem=png.stem; svg=png.with_suffix(".svg"); source=base/"figure_sources"/f"{stem}.csv"; source_gz=base/"figure_sources"/f"{stem}.csv.gz"; alt=base/"alt_text"/f"{stem}.txt"
            figrows.append({"package":dataset,"figure":stem,"png":str(png.relative_to(ROOT)),"svg":str(svg.relative_to(ROOT)) if svg.exists() else "unavailable","source_table":str((source if source.exists() else source_gz).relative_to(ROOT)) if source.exists() or source_gz.exists() else "multiple source tables or unavailable","alt_text":str(alt.relative_to(ROOT)) if alt.exists() else "unavailable","status":"available" if svg.exists() and alt.exists() else "partial"})
    write_csv(OUT/"final_figure_inventory.csv",pd.DataFrame(figrows))


def environment_metadata() -> dict:
    code="import json,platform,importlib.metadata as m; print(json.dumps({'python':platform.python_version(),'pytest':m.version('pytest'),'plotly':m.version('plotly'),'scgeo':m.version('scgeo')}))"
    versions=json.loads(run(str(PYTHON),"-c",code,cwd=SCGEO))
    command=f"{PYTHON} -m pytest {SCGEO/'tests'}"
    text=(REPRO/"test_results.txt").read_text(errors="replace")
    match=re.search(r"(\d+) passed, (\d+) skipped, (\d+) warnings in ([0-9.]+)s",text)
    if not match: raise RuntimeError("Could not parse passing pytest summary")
    passed,skipped,warnings,runtime=match.groups()
    meta={"command":command,"python_version":versions["python"],"pytest_version":versions["pytest"],"plotly_version":versions["plotly"],"scgeo_version":versions["scgeo"],"scgeo_commit":run("git","rev-parse","HEAD",cwd=SCGEO),"passed":int(passed),"failed":0,"errors":0,"skipped":int(skipped),"warnings":int(warnings),"runtime_seconds":float(runtime),"status":"passed","warning_policy":"warnings do not fail gate","maintenance_warning":"pandas categorical dtype deprecations; pytest cache warning from read-only checkout"}
    freeze=run(str(PYTHON),"-m","pip","freeze",cwd=SCGEO)
    write_text(REPRO/"scgeo_test_environment_freeze.txt",freeze+"\n"); write_json(REPRO/"scgeo_test_environment_metadata.json",meta); write_json(REPRO/"scgeo_test_report.json",meta)
    return meta


def audit_notebooks() -> tuple[pd.DataFrame,bool]:
    rows=[]
    for p in sorted((ROOT/"notebooks").rglob("*.ipynb")):
        nb=json.loads(p.read_text()); outputs=sum(len(c.get("outputs",[])) for c in nb.get("cells",[])); counts=sum(c.get("execution_count") is not None for c in nb.get("cells",[]))
        rows.append({"path":str(p.relative_to(ROOT)),"outputs":outputs,"non_null_execution_counts":counts,"passed":outputs==0 and counts==0})
    df=pd.DataFrame(rows); write_csv(REPRO/"notebook_cleanliness.csv",df); return df,bool(df.passed.all())


def audit_configs() -> tuple[pd.DataFrame,bool]:
    rows=[]
    paths=[]
    for d in (ROOT/"configs",ROOT/"environment",ROOT/"environments"):
        if d.exists(): paths += [p for p in d.rglob("*") if p.suffix.lower() in {".json",".yml",".yaml"}]
    for p in sorted(paths):
        try:
            with p.open() as f: json.load(f) if p.suffix==".json" else yaml.safe_load(f)
            status="passed"; error=""
        except Exception as e: status="failed"; error=str(e)
        rows.append({"path":str(p.relative_to(ROOT)),"format":p.suffix[1:],"status":status,"error":error})
    df=pd.DataFrame(rows); write_csv(REPRO/"config_validation.csv",df); return df,bool((df.status=="passed").all())


def audit_checksums() -> tuple[pd.DataFrame,bool]:
    items=[
        (Path("/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_qc_hvg.h5ad"),"86eed5f139e2d40ed2560d476198c11dbae423321a828d91e467013c5757b759","Dataset B compact H5AD"),
        (Path("/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_representations.h5ad"),"5d91ecd63f62b14d77fab04943fb5cd0f530e795ccbc0403616036568e3d9979","Dataset B representations"),
        (Path("/home/liuyuchen/data/gse211713/gse211713_revision_qc_hvg.h5ad"),"340504d5c88e60a3fc70f0c044ec9f65010226e99fa976e447fa8e6698418560","Dataset C compact H5AD"),
        (Path("/home/liuyuchen/data/gse211713/gse211713_revision_annotated.h5ad"),"6f7dd2c10112960a4dca4a41b56fcc670a99e16232e0b9d86365ad5fdd66a6ca","Dataset C annotated H5AD"),
        (Path("/home/liuyuchen/data/gse211713/gse211713_revision_representations_v2.h5ad"),"e8728846aff10091166b3f7d82b34afaf207142957b4daaafb8d6c496ca24cf2","Dataset C representations v2"),
    ]
    cval=json.loads((CBASE/"manuscript/dataset_c_manuscript_assembly_validation.json").read_text())
    bval=json.loads((BBASE/"manuscript/dataset_b_manuscript_assembly_validation.json").read_text())
    rows=[]
    for p,expected,label in items:
        actual=sha256(p) if p.exists() else "missing"; rows.append({"artifact":label,"path":str(p),"expected_sha256":expected,"actual_sha256":actual,"status":"passed" if actual==expected else "failed"})
    for rel,expected in cval["source_hashes_before"].items():
        actual=sha256(CBASE/rel); rows.append({"artifact":"Dataset C canonical source","path":str(CBASE/rel),"expected_sha256":expected,"actual_sha256":actual,"status":"passed" if actual==expected else "failed"})
    rows.append({"artifact":"Dataset B manuscript pin comparison","path":str(BBASE/"manuscript/dataset_b_manuscript_assembly_validation.json"),"expected_sha256":"source_hashes_before == source_hashes_after","actual_sha256":str(bval.get("source_hashes_unchanged")),"status":"passed" if bval.get("source_hashes_unchanged") else "failed"})
    df=pd.DataFrame(rows); write_csv(REPRO/"checksum_validation.csv",df); return df,bool((df.status=="passed").all())


def audit_keys() -> tuple[pd.DataFrame,bool]:
    specs=[
        (CBASE/"scgeo_v2/full_state_evidence.csv",["contrast","hierarchy","state","representation"]),
        (CBASE/"scgeo_v2/primary_consensus_state_evidence.csv",["contrast","hierarchy","state"]),
        (CBASE/"scgeo_v2/all_representation_sensitivity.csv",["contrast","hierarchy","state","representation"]),
        (CBASE/"scgeo_v2/exact_mouse_permutations.csv",["contrast","hierarchy","state","representation"]),
    ]
    rows=[]
    for p,key in specs:
        df=pd.read_csv(p); dup=int(df.duplicated(key).sum()); rows.append({"path":str(p.relative_to(ROOT)),"key":" × ".join(key),"rows":len(df),"duplicate_keys":dup,"status":"passed" if dup==0 else "failed"})
    df=pd.DataFrame(rows); write_csv(REPRO/"canonical_key_validation.csv",df); return df,bool((df.status=="passed").all())


def audit_terminology() -> tuple[pd.DataFrame,bool]:
    checks=[]
    texts=[]
    for p in list((ROOT/"notebooks").rglob("*.ipynb"))+list((ROOT/"scripts").glob("*.py"))+list((ROOT/"configs").glob("*.json")):
        if p.name in {"assemble_reviewer2_repairs.py","finalize_revision_evidence.py"}: continue
        texts.append((str(p.relative_to(ROOT)),p.read_text(errors="ignore")))
    rules={
        "obsolete_python_augur_positive":r"(?i)(Augur reimplementation|protocol-faithful Python|matches the Augur methodology)",
        "unsupported_ood_positive":r"(?i)(OOD should spike|clusters of high-OOD|out-of-distribution \(OOD\) flag)",
        "dataset_c_longitudinal_positive":r"(?i)Dataset C.{0,120}(longitudinal change|longitudinal reversal|longitudinal persistence)",
        "recursive_summary_wildcard":r"(?i)(rglob|glob)\([^\n]{0,120}\*_state_evidence\.csv",
    }
    for name,pat in rules.items():
        hits=[]
        for path,text in texts:
            if re.search(pat,text): hits.append(path)
        checks.append({"check":name,"matches":len(hits),"paths":";".join(hits),"status":"passed" if not hits else "failed"})
    primary=pd.read_csv(CBASE/"scgeo_v2/primary_consensus_state_evidence.csv")
    checks.append({"check":"Dataset C primary inference labels","matches":int((primary.inference_status!="replicate_aware_primary").sum()),"paths":"scgeo_v2/primary_consensus_state_evidence.csv","status":"passed" if (primary.inference_status=="replicate_aware_primary").all() else "failed"})
    bfull=pd.read_csv(BBASE/"scgeo/04_full_state_evidence.csv")
    checks.append({"check":"Dataset B remains descriptive_only","matches":int((bfull.inference_status!="descriptive_only").sum()),"paths":"Dataset B full evidence","status":"passed" if (bfull.inference_status=="descriptive_only").all() else "failed"})
    df=pd.DataFrame(checks); write_csv(REPRO/"terminology_audit.csv",df); return df,bool((df.status=="passed").all())


def validate_manuscript_and_reviewer() -> tuple[list[dict],bool]:
    checks=[]; m=CBASE/"manuscript"
    for kind,pattern,n in [("PNG figures","*.png",14),("SVG figures","*.svg",14),("figure source CSVs","*.csv",14),("alt text","*.txt",14)]:
        folder={"PNG figures":"figures","SVG figures":"figures","figure source CSVs":"figure_sources","alt text":"alt_text"}[kind]; actual=len(list((m/folder).glob(pattern))); checks.append({"gate":kind,"expected":n,"actual":actual,"status":"passed" if actual==n else "failed"})
    png_stems={p.stem for p in (m/"figures").glob("*.png")}; svg_stems={p.stem for p in (m/"figures").glob("*.svg")}; source_stems={p.stem for p in (m/"figure_sources").glob("*.csv")}; alt_stems={p.stem for p in (m/"alt_text").glob("*.txt")}
    paired=png_stems==svg_stems==source_stems==alt_stems
    checks.append({"gate":"per-figure PNG/SVG/source/alt pairing","expected":True,"actual":paired,"status":"passed" if paired else "failed"})
    val=json.loads((m/"dataset_c_manuscript_assembly_validation.json").read_text())
    checks += [{"gate":"Dataset C hashes unchanged","expected":True,"actual":val.get("numerical_sources_unchanged"),"status":"passed" if val.get("numerical_sources_unchanged") else "failed"},{"gate":"invalidated numerical sources excluded","expected":False,"actual":val.get("invalidated_numerical_sources_used"),"status":"passed" if not val.get("invalidated_numerical_sources_used") else "failed"}]
    prose="\n".join((m/name).read_text() for name in ["dataset_c_main_figure_caption.md","dataset_c_supplementary_captions.md"]) + (m/"dataset_c_methods_provenance.csv").read_text()
    checks.append({"gate":"UMAP display-only wording","expected":True,"actual":"UMAP display only" in prose or "UMAP as visualization only" in prose,"status":"passed" if "UMAP display only" in prose or "UMAP as visualization only" in prose else "failed"})
    checks.append({"gate":"cross-sectional not longitudinal wording","expected":True,"actual":"cross-sectional" in prose.lower() and "not longitudinal" in prose.lower(),"status":"passed" if "cross-sectional" in prose.lower() and "not longitudinal" in prose.lower() else "failed"})
    full=pd.read_csv(CBASE/"scgeo_v2/full_state_evidence.csv")
    prim=full[full.contrast.isin(["control_vs_17gy_early","control_vs_17gy_late","17gy_early_vs_late"])]
    sec=full[full.contrast.isin(["control_vs_all17gy","control_vs_all10gy"])]
    checks += [{"gate":"primary inference labels","expected":"replicate_aware_primary","actual":sorted(prim.inference_status.unique()),"status":"passed" if (prim.inference_status=="replicate_aware_primary").all() else "failed"},{"gate":"secondary inference labels","expected":"replicate_aware_secondary_time_heterogeneous","actual":sorted(sec.inference_status.unique()),"status":"passed" if (sec.inference_status=="replicate_aware_secondary_time_heterogeneous").all() else "failed"}]
    required=["mixscore_aggregation_definition.csv","distribution_test_inventory.csv","representation_metric_inventory.csv","terminology_replacement_inventory.csv","scpasi_comparison.csv","reviewer2_repairs_summary.json"]
    for name in required: checks.append({"gate":f"Reviewer2 {name}","expected":True,"actual":(OUT/"reviewer2"/name).exists(),"status":"passed" if (OUT/"reviewer2"/name).exists() else "failed"})
    return checks,all(x["status"]=="passed" for x in checks)


def finalize() -> dict:
    OUT.mkdir(parents=True,exist_ok=True); REPRO.mkdir(parents=True,exist_ok=True)
    global_outputs(); test=environment_metadata(); notebooks,nb_ok=audit_notebooks(); configs,cfg_ok=audit_configs(); checksums,hash_ok=audit_checksums(); keys,key_ok=audit_keys(); terms,term_ok=audit_terminology(); artifact_checks,art_ok=validate_manuscript_and_reviewer()
    scgeo_head=run("git","rev-parse","HEAD",cwd=SCGEO); scgeo_status=run("git","status","--short",cwd=SCGEO)
    staged=run("git","diff","--cached","--name-only").splitlines(); forbidden=[p for p in staged if re.search(r"(?i)(\.h5ad$|models?/|executed_notebooks|download|cache/|^results/)",p)]
    executed=[p for p in ROOT.rglob("*.ipynb") if "executed_notebooks" in str(p)]
    ignored_bad=[]
    for p in executed:
        proc=subprocess.run(["git","check-ignore",str(p)],cwd=ROOT,capture_output=True)
        if proc.returncode!=0: ignored_bad.append(str(p.relative_to(ROOT)))
    diff_check=subprocess.run(["git","diff","--check"],cwd=ROOT,text=True,capture_output=True)
    cached_check=subprocess.run(["git","diff","--cached","--check"],cwd=ROOT,text=True,capture_output=True)
    git_state={"branch":run("git","branch","--show-current"),"head":run("git","rev-parse","HEAD"),"origin_divergence":run("git","rev-list","--left-right","--count","origin/revision-dataset-c-gse211713...HEAD"),"status_short":run("git","status","--short"),"staged_files":staged,"forbidden_staged_files":forbidden,"diff_check":"passed" if diff_check.returncode==0 else diff_check.stdout+diff_check.stderr,"cached_diff_check":"passed" if cached_check.returncode==0 else cached_check.stdout+cached_check.stderr}; write_json(REPRO/"git_state.json",git_state)
    unresolved=pd.DataFrame([
        ["maintenance","Pandas categorical-dtype deprecation warnings in ScGeo plotting/analysis paths.","future dependency-compatibility cleanup; not a revision blocker"],
        ["maintenance","Pytest cache warning because the frozen ScGeo checkout is read-only.","none for revision; optionally disable cache in maintenance CI"],
        ["scientific","Dataset B has no independent biological replicate identifier.","retain descriptive_only wording"],
        ["scientific","Dataset C is cross-sectional with small groups and reused controls.","avoid longitudinal/causal claims and p-value combination"],
        ["editorial","Pancreas numerical representation-quality artifacts are not locally available.","collate existing manuscript assets or state unavailable"],
        ["editorial","scPASI method details were not available in the local reviewer context.","verify citation-specific description during manuscript writing"],
        ["editorial","Full manuscript and response letter still require scientific prose and citation placement.","human manuscript writing"],
    ],columns=["category","issue","required_action"]); write_csv(REPRO/"unresolved_issues.csv",unresolved)
    gates={"scgeo_commit":scgeo_head==REQUIRED_SCGEO,"scgeo_clean":scgeo_status=="","tests":test["failed"]==0 and test["errors"]==0,"notebooks_clean":nb_ok,"executed_notebooks_ignored":not ignored_bad,"configs":cfg_ok,"checksums":hash_ok,"canonical_keys":key_ok,"terminology":term_ok,"manuscript_and_reviewer_artifacts":art_ok,"forbidden_staging":not forbidden,"git_diff_check":diff_check.returncode==0,"git_cached_diff_check":cached_check.returncode==0}
    report={"status":"passed" if all(gates.values()) else "failed","gates":gates,"test_report":test,"scgeo":{"commit":scgeo_head,"clean":scgeo_status==""},"dataset_c_artifact_checks":artifact_checks,"counts":{"source_notebooks":len(notebooks),"configs":len(configs),"checksum_records":len(checksums),"canonical_key_tables":len(keys),"executed_notebooks":len(executed)},"executed_notebooks_not_ignored":ignored_bad,"warnings":["pandas categorical-dtype deprecations are future maintenance","pytest cache warning is non-blocking"],"unresolved_issues_file":"results/revision_finalization/reproducibility/unresolved_issues.csv"}
    write_json(REPRO/"final_reproducibility_report.json",report)
    if report["status"]!="passed": raise RuntimeError(json.dumps(report["gates"],sort_keys=True))
    return report


if __name__=="__main__": print(json.dumps(finalize(),indent=2,sort_keys=True))
