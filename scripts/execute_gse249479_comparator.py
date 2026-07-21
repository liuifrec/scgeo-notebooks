#!/usr/bin/env python3
"""Execute GSE249479 Phase 3C comparator notebooks in fresh processes."""

from __future__ import annotations

import argparse, hashlib, json, os, subprocess, sys, time, traceback
from pathlib import Path
import nbformat
from nbclient import NotebookClient

ROOT=Path(__file__).resolve().parents[1]; OUTPUT=ROOT/"results/public_validation/gse249479_dataset_b"
NOTEBOOKS=["notebooks/public_validation/gse249479/05_augur_comparator.ipynb","notebooks/public_validation/gse249479/06_scgeo_augur_comparison.ipynb"]

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def output_bytes(nb):
    total=0
    for c in nb.cells:
        for o in c.get("outputs",[]):
            if "text" in o: total+=len("".join(o["text"]) if isinstance(o["text"],list) else str(o["text"]))
            for v in o.get("data",{}).values(): total+=len("".join(v) if isinstance(v,list) else str(v))
    return total

def worker(rel,timeout):
    source=ROOT/rel; executed=OUTPUT/"executed_notebooks"/rel; report=OUTPUT/"execution"/f"{source.stem}_execution_report.json"; before=sha256(source); started=time.perf_counter()
    try:
        nb=nbformat.read(source,as_version=4); embedded=output_bytes(nb)
        if embedded or any(c.get("execution_count") is not None for c in nb.cells): raise RuntimeError("Source notebook is not clean-kernel/output-free")
        NotebookClient(nb,timeout=timeout,kernel_name="scgeo_pre",resources={"metadata":{"path":str(ROOT)}}).execute()
        after=sha256(source)
        if after!=before: raise RuntimeError("Source notebook changed")
        executed.parent.mkdir(parents=True,exist_ok=True); tmp=executed.with_name(f".{executed.name}.tmp.{os.getpid()}"); nbformat.write(nb,tmp); os.replace(tmp,executed)
        row={"notebook":rel,"kernel":"scgeo_pre","status":"passed","inference_status":"descriptive_only","runtime_seconds":time.perf_counter()-started,"source_sha256_before":before,"source_sha256_after":after,"source_output_bytes":embedded,"executed_notebook":str(executed.relative_to(ROOT)),"executed_notebook_sha256":sha256(executed),"executed_output_bytes":output_bytes(nb),"worker_pid":os.getpid()}; code=0
    except Exception as exc:
        row={"notebook":rel,"kernel":"scgeo_pre","status":"failed","inference_status":"descriptive_only","runtime_seconds":time.perf_counter()-started,"source_sha256_before":before,"source_sha256_after":sha256(source),"error":repr(exc),"traceback":traceback.format_exc(),"worker_pid":os.getpid()}; code=1
    report.parent.mkdir(parents=True,exist_ok=True); tmp=report.with_name(f".{report.name}.tmp.{os.getpid()}"); tmp.write_text(json.dumps(row,indent=2,sort_keys=True)+"\n"); os.replace(tmp,report); print(json.dumps(row,indent=2),flush=True); return code

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--worker-notebook"); ap.add_argument("--timeout",type=int,default=14400); a=ap.parse_args()
    if a.worker_notebook: return worker(a.worker_notebook,a.timeout)
    if subprocess.run(["git","-C",str(ROOT),"check-ignore","-q",str((OUTPUT/"executed_notebooks").relative_to(ROOT))]).returncode: raise RuntimeError("Executed notebooks must be ignored")
    env=os.environ.copy(); env.update({"SCGEO_GSE249479_COMPACT_H5AD":"/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_qc_hvg.h5ad","SCGEO_GSE249479_REPRESENTATION_H5AD":"/home/liuyuchen/hsc_memory_nature_2026/results/gse249479_revision_representations.h5ad","SCGEO_GSE249479_OUTPUT_DIR":str(OUTPUT),"SCGEO_SOURCE_REPO":"/home/liuyuchen/Github/scgeo","PYTHONPATH":"/home/liuyuchen/Github/scgeo:"+str(ROOT/"scripts"),"NUMBA_CACHE_DIR":str(OUTPUT/"_numba_cache_augur"),"MPLCONFIGDIR":str(OUTPUT/"_matplotlib_cache_augur"),"PYTHONHASHSEED":"0","OMP_NUM_THREADS":"1","MKL_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1"})
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True,exist_ok=True); Path(env["MPLCONFIGDIR"]).mkdir(parents=True,exist_ok=True)
    rows=[]; overall=0
    for rel in NOTEBOOKS:
        print(f"[gse249479-phase3c] fresh worker: {rel}",flush=True)
        r=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--worker-notebook",rel,"--timeout",str(a.timeout)],cwd=ROOT,env=env,check=False)
        rows.append(json.loads((OUTPUT/"execution"/f"{Path(rel).stem}_execution_report.json").read_text()))
        if r.returncode: overall=r.returncode; break
    final={"workflow":"gse249479_phase3c_augur_comparator","status":"passed" if overall==0 else "failed","inference_status":"descriptive_only","fresh_worker_per_notebook":True,"notebooks":rows}
    (OUTPUT/"execution"/"gse249479_phase3c_execution_report.json").write_text(json.dumps(final,indent=2,sort_keys=True)+"\n"); return overall

if __name__=="__main__": raise SystemExit(main())
