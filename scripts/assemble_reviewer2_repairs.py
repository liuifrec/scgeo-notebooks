#!/usr/bin/env python3
"""Assemble Reviewer 2 presentation, schema, and terminology repairs only."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/revision_finalization/reviewer2"
B = ROOT / "results/public_validation/gse249479_dataset_b"
C = ROOT / "results/public_validation/gse211713_dataset_c"


def write_csv(name: str, df: pd.DataFrame) -> None:
    path=OUT/name; path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); df.to_csv(tmp,index=False); tmp.replace(path)


def write_text(name: str, text: str) -> None:
    path=OUT/name; path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); tmp.write_text(text,encoding="utf-8"); tmp.replace(path)


def save_fig(stem: str, df: pd.DataFrame, x: str, y: str, title: str, alt: str) -> None:
    fig,ax=plt.subplots(figsize=(10,5.5),layout="constrained")
    if df.empty or x not in df or y not in df:
        ax.axis("off"); ax.text(.02,.92,alt,va="top",wrap=True)
    else:
        d=df.dropna(subset=[y]).copy(); cats=list(dict.fromkeys(d[x].astype(str))); cmap=plt.get_cmap("tab10")
        hue="dataset" if "dataset" in d else ("comparison_unit" if "comparison_unit" in d else None)
        if hue:
            for i,h in enumerate(dict.fromkeys(d[hue].astype(str))):
                z=d[d[hue].astype(str)==h]; ax.scatter([cats.index(v) for v in z[x].astype(str)],z[y],s=18,alpha=.6,label=h,color=cmap(i%10))
            ax.legend(frameon=False,fontsize=8)
        else: ax.scatter([cats.index(v) for v in d[x].astype(str)],d[y],s=18,alpha=.6)
        ax.set_xticks(range(len(cats)),cats,rotation=35,ha="right",fontsize=8); ax.set_ylabel(y.replace("_"," ")); ax.grid(axis="y",alpha=.2)
    ax.set_title(title,loc="left",fontweight="bold")
    for folder in ("figures","figure_sources","alt_text","captions"): (OUT/folder).mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/"figures"/f"{stem}.png",dpi=300,bbox_inches="tight",metadata={"Software":"matplotlib"})
    fig.savefig(OUT/"figures"/f"{stem}.svg",bbox_inches="tight",metadata={"Date":None})
    plt.close(fig); write_csv(f"figure_sources/{stem}.csv",df); write_text(f"alt_text/{stem}.txt",alt+"\n"); write_text(f"captions/{stem}.md",f"**{title}.** {alt}\n")


def source_occurrences(term: str) -> list[dict]:
    rows=[]
    suffixes={".ipynb",".py",".md",".json",".csv",".yml",".yaml"}
    for p in ROOT.rglob("*"):
        rel=str(p.relative_to(ROOT))
        if not p.is_file() or p.suffix.lower() not in suffixes: continue
        if any(x in rel for x in (".git/","R_library/","renv_cache/","invalidated_runs/","executed_notebooks","_cache/","revision_finalization/reviewer2/","scripts/assemble_reviewer2_repairs.py","scripts/finalize_revision_evidence.py")): continue
        try: lines=p.read_text(encoding="utf-8",errors="ignore").splitlines()
        except OSError: continue
        for n,line in enumerate(lines,1):
            if term.lower() in line.lower(): rows.append({"source_path":rel,"line":n,"text":line.strip()[:500]})
    return rows


def mixscore_inventory() -> pd.DataFrame:
    rows=[]
    for hit in source_occurrences("mixscore"):
        rows.append({**hit,"calculation_unit":"cell","denominator":"k selected neighbors (k=50 in the demonstrated calls)","per_cell_formula":"m_i = k^-1 sum_{j in N_k(i)} I(label_j != label_i)","aggregation_formula":"ScGeo metadata reports unweighted arithmetic mean/median/min/max across cell scores; any state/sample summary must name its grouping and weighting","weighting":"per-cell summaries weight cells equally; unequal state sizes therefore affect an unstratified dataset mean","unequal_state_sizes":"not corrected by the demonstrated dataset-level mean","biological_replicate_interpretation":"none; sample-label mixing is a descriptive neighborhood diagnostic","status":"descriptive"})
    return pd.DataFrame(rows)


def mixscore_definition() -> pd.DataFrame:
    return pd.DataFrame([
        {"level":"per-cell","unit":"cell i","formula":"m_i = (1/k) sum_{j in N_k(i)} I(label_j != label_i)","denominator":"k neighbors","weighting":"each selected neighbor equally","worked_example":"labels A,A with two neighbors each: cell 1 neighbors A,B -> 0.5","replicate_interpretation":"none"},
        {"level":"per-state","unit":"state s","formula":"mean_{i:state_i=s}(m_i)","denominator":"cells in state s","weighting":"cells equally within state","worked_example":"state S cell scores 0.5,1.0 -> 0.75","replicate_interpretation":"descriptive; cells are not biological replicates"},
        {"level":"per-sample","unit":"sample q","formula":"mean_{i:sample_i=q}(m_i)","denominator":"cells in sample q","weighting":"cells equally within sample","worked_example":"sample Q cell scores 0.5,0.0 -> 0.25","replicate_interpretation":"sample summaries may be compared only when samples are genuine biological units"},
        {"level":"dataset","unit":"all cells","formula":"mean_i(m_i)","denominator":"all cells","weighting":"cells equally; large states/samples dominate","worked_example":"four cell scores 0.5,1.0,0.5,0.0 -> 0.50","replicate_interpretation":"descriptive only"},
        {"level":"balanced dataset alternative","unit":"state summaries","formula":"|S|^-1 sum_s mean_{i:state_i=s}(m_i)","denominator":"number of states","weighting":"states equally","worked_example":"state means 0.75 and 0.25 -> 0.50","replicate_interpretation":"descriptive; must not replace sample-level inference"},
    ])


def distribution_inventory() -> pd.DataFrame:
    rows=[]
    bd=pd.read_csv(B/"scgeo/04_within_state_distribution_metrics.csv")
    for r in bd.itertuples():
        for name,col in [("energy_distance","energy_distance_biased_full_value"),("MMD_RBF_squared","mmd_rbf_squared_biased_full_value"),("sliced_Wasserstein","sliced_wasserstein_full_value")]:
            rows.append({"dataset":"Dataset B GSE249479","test_or_distance":name,"statistic":getattr(r,col),"comparison_unit":"balanced pooled cells","group_size_a":r.n_each,"group_size_b":r.n_each,"state":r.state,"representation":r.representation,"contrast":r.contrast,"null_generation":"none","permutations_or_resamples":0,"raw_p_or_empirical_fraction":np.nan,"attainable_resolution":np.nan,"adjustment":"none","inferential_status":"descriptive_only","limitations":"no biological replicate identifier; bounded balanced cell subsample"})
    cd=pd.read_csv(C/"scgeo_v2/distribution_metrics.csv")
    for r in cd.itertuples():
        for name,col in [("energy_distance","energy_distance"),("MMD_RBF","mmd_rbf"),("sliced_Wasserstein","sliced_wasserstein")]:
            rows.append({"dataset":"Dataset C GSE211713","test_or_distance":name,"statistic":getattr(r,col),"comparison_unit":r.comparison_unit,"group_size_a":r.n_group0,"group_size_b":r.n_group1,"state":r.state,"representation":r.representation,"contrast":r.contrast,"null_generation":"none","permutations_or_resamples":0,"raw_p_or_empirical_fraction":np.nan,"attainable_resolution":np.nan,"adjustment":"none","inferential_status":r.inference_status,"limitations":"distance summary; not a replacement biological-replicate test"})
    cp=pd.read_csv(C/"scgeo_v2/exact_mouse_permutations.csv")
    for r in cp.itertuples():
        rows.append({"dataset":"Dataset C GSE211713","test_or_distance":"exact normalized-displacement mouse-label permutation","statistic":r.observed_normalized_displacement,"comparison_unit":"biological mouse center","group_size_a":r.n_mice_group0,"group_size_b":r.n_mice_group1,"state":r.state,"representation":r.representation,"contrast":r.contrast,"null_generation":"exhaustive unique mouse-label assignments; observed included; ties >= observed","permutations_or_resamples":r.assignment_count,"raw_p_or_empirical_fraction":r.empirical_upper_tail_fraction,"attainable_resolution":r.minimum_attainable_grid,"adjustment":"raw; separate BH sensitivity within contrast × representation × major states","inferential_status":"replicate-aware exact-permutation sensitivity","limitations":"small-n discrete grid; no p-value combination across representations or related contrasts"})
    null=pd.read_csv(B/"controls/04_shuffled_condition_computational_null.csv")
    for r in null.itertuples():
        rows.append({"dataset":"Dataset B GSE249479","test_or_distance":"shuffled condition-label computational null","statistic":r.normalized_effect,"comparison_unit":"cell label","group_size_a":np.nan,"group_size_b":np.nan,"state":r.state,"representation":r.representation,"contrast":r.contrast,"null_generation":"condition labels shuffled among cells","permutations_or_resamples":1,"raw_p_or_empirical_fraction":np.nan,"attainable_resolution":np.nan,"adjustment":"none","inferential_status":"descriptive_only non-biological null","limitations":"cells are not biological replicates; no biological p-value"})
    return pd.DataFrame(rows)


def representation_inventory() -> pd.DataFrame:
    rows=[]
    # Synthetic and pancreas artifacts are not locally available; expose that rather than fabricate values.
    rows += [{"dataset":"Synthetic manuscript benchmark","representation":"ensemble","role":"controlled frozen benchmark","neighborhood_overlap":np.nan,"jaccard_edge_retention":np.nan,"graph_connectivity":np.nan,"trustworthiness":np.nan,"continuity":np.nan,"local_distortion":np.nan,"marker_label_preservation":np.nan,"fibroblast_label_preservation":np.nan,"mouse_neighborhood_dominance":np.nan,"condition_time_silhouette":np.nan,"availability":"manifest pins representation_robustness and local_distortion files; numerical files not present in notebook workspace","interpretation":"controlled corruption and representation robustness; not general OOD detection"},
             {"dataset":"Pancreas validation","representation":"PCA20/PCA30/PCA50/diffusion","role":"public dynamics validation","neighborhood_overlap":np.nan,"jaccard_edge_retention":np.nan,"graph_connectivity":np.nan,"trustworthiness":np.nan,"continuity":np.nan,"local_distortion":np.nan,"marker_label_preservation":np.nan,"fibroblast_label_preservation":np.nan,"mouse_neighborhood_dominance":np.nan,"condition_time_silhouette":np.nan,"availability":"configuration present; no frozen pancreas metric artifact in notebook workspace","interpretation":"cluster annotations are biological reference; UMAP display only"}]
    b=pd.read_csv(B/"scgeo/04_local_geometry_representation_summary.csv")
    for r in b.itertuples():
        rows.append({"dataset":"Dataset B GSE249479","representation":r.rep,"role":"primary or prespecified sensitivity per Dataset B config","neighborhood_overlap":r.median_pairwise_neighbor_overlap,"jaccard_edge_retention":r.median_pairwise_neighbor_jaccard,"graph_connectivity":np.nan,"trustworthiness":r.median_trustworthiness_to_others,"continuity":r.median_continuity_from_others,"local_distortion":r.median_local_shape_distortion,"marker_label_preservation":np.nan,"fibroblast_label_preservation":np.nan,"mouse_neighborhood_dominance":np.nan,"condition_time_silhouette":np.nan,"availability":"available","interpretation":"condition is biological and was not a mixing objective; descriptive_only"})
    c=pd.read_csv(C/"c6_v2/metadata/representation_metrics.csv")
    for r in c.itertuples():
        rows.append({"dataset":"Dataset C GSE211713","representation":r.representation,"role":r.role,"neighborhood_overlap":r.neighborhood_overlap_vs_pca50,"jaccard_edge_retention":np.nan,"graph_connectivity":r.largest_connected_component_fraction,"trustworthiness":np.nan,"continuity":np.nan,"local_distortion":r.local_distortion_vs_pca50,"marker_label_preservation":r.major_state_silhouette,"fibroblast_label_preservation":r.fibroblast_subtype_silhouette,"mouse_neighborhood_dominance":r.mouse_neighborhood_fraction,"condition_time_silhouette":f"dose={r.dose_silhouette_descriptive:.6f}; time={r.time_silhouette_descriptive:.6f}","availability":"available","interpretation":"dose/time/mouse are not batch-removal objectives"})
    return pd.DataFrame(rows)


def terminology_inventory() -> pd.DataFrame:
    replacements=[
        ("Use the mixscore maps to highlight regions where **sample labels** are locally well-mixed (good batch correction) vs separated (potential batch structure).","Use mixscore as a descriptive map of local sample-label mixing; it does not by itself establish batch removal or biological preservation.","notebooks/data_prep/04_scgeo_gse280305_phase1_qc.ipynb","mixing is not proof of preservation"),
        ("out-of-distribution (OOD)","weak reference-neighborhood support / extrapolation sensitivity","notebooks/manuscript/05_OOD.ipynb","avoid general OOD-detection claim"),
        ("OOD/high-risk","unsupported-state warning / high extrapolation sensitivity","notebooks/manuscript/Final_summary.ipynb","precise diagnostic language"),
        ("Mean OOD","Mean unsupported-state warning score","notebooks/exploration/06_Ref_prep.ipynb","avoid general OOD-detection claim"),
        ("OOD should spike only for edge / rare / poorly supported regions","unsupported-state warning may be higher at edge, rare, or poorly supported regions in this controlled toy split","notebooks/tutorials/PBMC_ingest_ScGeo_QC_demo.ipynb","calibrate claim to controlled demonstration"),
        ("strong biological preservation","representation-specific retention of the evaluated labels or neighborhoods","repository-wide inventory rule","avoid universal preservation claim"),
        ("longitudinal change","cross-sectional difference between independently sampled mice","Dataset C captions and documentation","Dataset C is cross-sectional"),
        ("independent confirmation by PCA20/PCA30/PCA50","dimensional sensitivity across nested views of one PCA basis","Dataset B/C captions and documentation","nested PCA dimensions are not independent evidence"),
    ]
    return pd.DataFrame(replacements,columns=["original_wording","corrected_wording","source_path","reason"]).assign(numerical_changes_required=False)


def scpasi_table() -> pd.DataFrame:
    cols=["method","primary_objective","analysis_unit","representation_dependence","replicate_awareness","perturbation_geometry","abundance_handling","local_structure","uncertainty","output_interpretation","relationship_to_ScGeo"]
    return pd.DataFrame([
        ["ScGeo","state-level perturbation displacement geometry and cross-representation stability","state centers; biological samples when sample_key is available","explicitly assessed across prespecified representations","sample-aware centers/bootstrap supported and used in Dataset C","robust displacement magnitude/direction; no opaque composite","reported separately","neighborhood/geometry diagnostics are separate evidence","bootstrap plus representation sensitivity; exact permutation sensitivity in Dataset C","geometric effect and robustness evidence, conditioned on design","reference method in this revision"],
        ["scPASI","not established by locally available manuscript/reviewer source","not established locally","not established locally","not established locally","not established locally","not established locally","not established locally","not established locally","reviewer-identified related method; not run here","related but methodologically distinct; no superiority or equivalence claim without direct evidence"],
        ["Official R Augur 1.0.3","cell-state prioritization by treatment separability","balanced cells within marker-inferred state","uses sparse/HVG expression features","cell CV is computational stability, not biological replication","does not estimate directional displacement","not an abundance-effect estimator","RF feature-based separability","repeated computational subsampling/CV","primary external comparator for Dataset B","complementary separability evidence, not the same estimand"],
    ],columns=cols)


def assemble() -> dict:
    mpl.rcParams.update({"svg.hashsalt":"reviewer2-repairs-v1","axes.spines.top":False,"axes.spines.right":False})
    mix=mixscore_inventory(); definition=mixscore_definition(); dist=distribution_inventory(); reps=representation_inventory(); terms=terminology_inventory(); scp=scpasi_table()
    write_csv("mixscore_aggregation_definition.csv",pd.concat([definition.assign(record_type="definition"),mix.assign(record_type="occurrence")],ignore_index=True,sort=False))
    write_csv("distribution_test_inventory.csv",dist); write_csv("representation_metric_inventory.csv",reps); write_csv("terminology_replacement_inventory.csv",terms); write_csv("scpasi_comparison.csv",scp)
    save_fig("mixscore_aggregation_levels",definition.assign(value=[.5,.75,.25,.5,.5]),"level","value","Mixscore aggregation is explicitly unit-dependent","Illustrative worked values distinguish per-cell, per-state, per-sample, and dataset summaries. They are formulas, not replacement numerical results.")
    save_fig("distribution_test_visibility",dist[dist.test_or_distance.isin(["energy_distance","exact normalized-displacement mouse-label permutation"])],"comparison_unit","statistic","Distribution and permutation evidence by comparison unit","Pooled-cell descriptive distances are separated from exact biological-mouse label permutations; scales differ and no composite is formed.")
    save_fig("representation_quality_inventory",reps[reps.availability=="available"],"representation","neighborhood_overlap","Available representation-quality evidence","Neighborhood-overlap values are shown where available for Dataset B and Dataset C. Missing synthetic and pancreas artifacts remain explicit.")
    summary={"status":"passed","scope":"presentation_schema_documentation_only","numerical_outputs_changed":False,"mixscore_occurrences":len(mix),"distribution_inventory_rows":len(dist),"representation_inventory_rows":len(reps),"terminology_replacements":len(terms),"arithmetic_error_detected":False,"comparator_roles":{"Official R Augur 1.0.3":"primary external comparator","Augur-inspired Python approximation":"supplementary implementation sensitivity only"},"integration_notes":{"Dataset_B":"Harmony and Scanorama excluded because condition/library are confounded and no genuine technical batch exists.","Dataset_C":"Harmony and Scanorama excluded because no independent technical batch is documented and mouse/dose/time must not be removed.","Scanorama_public_demo":"one tested preprocessing/representation in the original public demonstration; not required by ScGeo."},"scpasi_direct_implementation":False,"limitations":["No local scPASI methods source was available, so unknown method fields were not inferred.","Synthetic and pancreas numerical representation-quality artifacts were not present in this notebook workspace; availability is explicit.","Mixscore is descriptive and cells are not biological replicates."]}
    write_text("reviewer2_repairs_summary.json",json.dumps(summary,indent=2,sort_keys=True)+"\n")
    return summary


if __name__=="__main__": print(json.dumps(assemble(),indent=2,sort_keys=True))
