#!/usr/bin/env python3
"""Assemble checksum-pinned Dataset C manuscript artifacts without recomputation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/public_validation/gse211713_dataset_c"
SCGEO = BASE / "scgeo_v2"
OUT = BASE / "manuscript"

PINNED = {
    "scgeo_v2/full_state_evidence.csv": "82f2d51b53aad9d996bbeac9ebdc2125f282ca5809c7a629869d18dae85abea4",
    "scgeo_v2/primary_consensus_state_evidence.csv": "48f9ad636c971cc55690f0edcd9ba5129b77b1fa8aa47dce18008d1797b0f566",
    "scgeo_v2/all_representation_sensitivity.csv": "b3ab9dcec88b4385af28239bb8fcd9f1752332c62d0b908d94c858a0403b83c7",
    "scgeo_v2/negative_neutral_or_unstable_states.csv": "7af27cb2163e5f89bdcdb4a0403d7a538c96c57f1870a78a30f9dd2b68df085c",
    "scgeo_v2/individual_mouse_centers.csv": "aa706e5286716a0bf5a1eb59fbdc88bcfbb5b0b764acab5189ab7c0cdfe9c8f3",
    "scgeo_v2/biological_mouse_bootstrap_intervals.csv": "0eaec3e051a7afd6cf338c3b050872d24b9e81cd38141b4a6b962cf9e5b3fe94",
    "scgeo_v2/biological_mouse_bootstrap_draws.csv": "f9ccc3a94e8c11a677b9ff5b2ccd57f3aa354fbbb503256b72f1ea5a6b522ce1",
    "scgeo_v2/exact_mouse_permutations.csv": "c25d4fef233c88c44cb32cc4acfbb04eaa5aa355851a4cb2b9e3761e1ed87bfd",
    "scgeo_v2/exact_mouse_permutation_values.csv": "0b4f48adc8730a5a80b1b714db0e84a59f1cf02b10369ba416d12a07d9699895",
    "scgeo_v2/exact_mouse_permutation_bh_sensitivity.csv": "7807fe3683b07b337d472c69b2ee2cf42d16b17e15d5383f04f477e06f36c46e",
    "scgeo_v2/mouse_level_abundance.csv": "8f1e64ab0d4bc78b7c2d4da9f92fa3f455dae93b88045c780d11c0e2765bc5ae",
    "scgeo_v2/distribution_metrics.csv": "36c181f70380cb39e2aac4b784009cffaca8f5dcf406da699246b70466adcdb8",
    "scgeo_v2/eligibility_and_coverage_used.csv": "b22c233556dabc31389ef630dd74fcd527bdb34a0941e6c0fef11fa35f1c808a",
    "scgeo_v2/output_schema.json": "e04f4bfff760b83092269ec95764a7105743f47de0e7d7f6aba48ce23827896f",
    "scgeo_v2/final_acceptance_gate.json": "918c0b63fafa073b57acfe59a4d7abf9f7dd427464122e267c1fe87da002e0ac",
    "c6_v2/metadata/representation_metrics.csv": "696f491a49ed9c80e3985272b15e17aa0c43e9ad2d7425dfd44d541ad10a5fef",
    "qc/qc_summary_by_mouse.csv": "2ed7bb4d073914f27058031cd38f70b65b44ae2f240eceba353c6f53092bfd4d",
    "annotation/annotation_confidence_table.csv": "56f7c6ec1b212d8bb2db637be6d9d9bc316f7851c2f32f46472a1c218b88b3d9",
    "annotation/marker_definitions.csv": "1379785c9ccc1b19f36c30812a3975561771173307d5bde2be2644d18332fb7b",
    "coverage/contrast_state_eligibility.csv": "d9099bdb45d393d80aef4cf8fa093db04038d48043cecaf5d6eb8e4c4948432b",
    "audit/study_design.csv": "6f521c29ce8a90944192a2bd3a95372372058613adaf2724490c1308b3709895",
    "invalidated_runs/c6_c7_model_switch_run/invalidation_manifest.json": "4f388bf8dfa44faf4e4b49dcdc66059717ff52a3345aa6fcf7a4f8fa7fe5d633",
}

CONTRAST_LABELS = {
    "control_vs_17gy_early": "Control vs 17 Gy early",
    "control_vs_17gy_late": "Control vs 17 Gy late",
    "17gy_early_vs_late": "17 Gy early vs late",
    "control_vs_all_17gy": "Control vs all 17 Gy",
    "control_vs_all_10gy": "Control vs all 10 Gy",
}
MAJOR = ["Endothelial", "Epithelial", "Fibroblast/stromal", "Lymphoid", "Myeloid", "Proliferating"]
FIBRO = ["Col13a1-like matrix fibroblast", "Col14a1-like matrix fibroblast", "Myofibroblast"]
REPS = ["X_pca30", "X_pca50", "X_scvi", "X_pca20", "X_diffmap"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False)
    tmp.replace(path)


def pin_check() -> dict[str, str]:
    found = {}
    for rel, expected in PINNED.items():
        path = BASE / rel
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"Checksum mismatch for {rel}: {actual} != {expected}")
        found[rel] = actual
    return found


def save_figure(fig: plt.Figure, stem: str, source: pd.DataFrame, alt: str) -> None:
    figdir, srcdir, altdir = OUT / "figures", OUT / "figure_sources", OUT / "alt_text"
    for directory in (figdir, srcdir, altdir):
        directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / f"{stem}.png", dpi=300, bbox_inches="tight", metadata={"Software": "matplotlib"})
    fig.savefig(figdir / f"{stem}.svg", bbox_inches="tight", metadata={"Date": None, "Creator": "matplotlib"})
    plt.close(fig)
    atomic_csv(srcdir / f"{stem}.csv", source)
    atomic_text(altdir / f"{stem}.txt", alt.strip() + "\n")


def effect_panel(ax, consensus: pd.DataFrame, contrast: str, title: str) -> None:
    d = consensus[(consensus.contrast == contrast) & (consensus.hierarchy == "major")].copy()
    d["state"] = pd.Categorical(d.state, MAJOR[::-1], ordered=True)
    d = d.sort_values("state")
    colors = ["#167c80" if x == "representation_stable_effect" else "#c46b45" for x in d.consensus_label]
    ax.barh(d.state.astype(str), d.normalized_magnitude_median, color=colors, alpha=.88)
    ax.errorbar(d.normalized_magnitude_median, np.arange(len(d)),
                xerr=[d.normalized_magnitude_median-d.normalized_magnitude_q25,
                      d.normalized_magnitude_q75-d.normalized_magnitude_median],
                fmt="none", ecolor="#222222", capsize=2, lw=.8)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("Normalized displacement (primary median)")
    ax.grid(axis="x", alpha=.2)


def manuscript_figure(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    cons, full, centers, abundance, perm = (data[k] for k in ("consensus", "full", "centers", "abundance", "perm"))
    fig = plt.figure(figsize=(17, 20), layout="constrained")
    gs = fig.add_gridspec(4, 2, height_ratios=[.75, 1, 1.1, 1])
    ax = fig.add_subplot(gs[0, 0]); ax.axis("off")
    ax.set_title("A  Experimental design", loc="left", fontweight="bold")
    ax.text(.02,.74,"Control\n5 mice",ha="left",va="center",fontsize=14,bbox=dict(boxstyle="round",fc="#e9ecef",ec="#555"))
    ax.text(.40,.74,"17 Gy early\nmonths 1–2, 4 mice",ha="left",va="center",fontsize=14,bbox=dict(boxstyle="round",fc="#f7d6c4",ec="#555"))
    ax.text(.78,.74,"17 Gy late\nmonths 4–5, 4 mice",ha="left",va="center",fontsize=14,bbox=dict(boxstyle="round",fc="#d9c8ea",ec="#555"))
    ax.text(.02,.23,"Independent mouse/GSM libraries; cross-sectional design, not longitudinal.\nPrimary contrasts use mouse-level centers and biological-mouse resampling.",fontsize=12)
    ax = fig.add_subplot(gs[0, 1]); ax.axis("off")
    ax.set_title("F  Biological-replicate evidence", loc="left", fontweight="bold")
    subset = centers[(centers.contrast=="control_vs_17gy_late")&(centers.hierarchy=="major")&(centers.state=="Epithelial")&(centers.representation=="X_pca30")]
    inset=ax.inset_axes([.04,.15,.58,.72])
    for group,color in [("control","#4c78a8"),("17gy_late","#e45756")]:
        z=subset[subset.group==group]; inset.scatter(z.coordinate_1,z.coordinate_2,label=group,s=38,color=color)
        for row in z.itertuples(): inset.annotate(row.mouse_id.replace("GSM649", "…"),(row.coordinate_1,row.coordinate_2),fontsize=6)
    inset.set_xlabel("PCA30 mouse-center coordinate 1"); inset.set_ylabel("coordinate 2"); inset.legend(frameon=False,fontsize=7)
    ax.text(.68,.72,"500\nmouse-bootstrap\ndraws",ha="center",va="center",fontsize=11,bbox=dict(boxstyle="round",fc="#eef5f5"))
    ax.text(.68,.38,"Exact label grids\n5 vs 4: 1/126\n4 vs 4: 1/70\n(complement symmetry)",ha="center",va="center",fontsize=10)
    effect_panel(fig.add_subplot(gs[1,0]),cons,"control_vs_17gy_early","B  Early 17 Gy vs control: no frozen stable consensus")
    effect_panel(fig.add_subplot(gs[1,1]),cons,"control_vs_17gy_late","C  Late 17 Gy vs control")
    effect_panel(fig.add_subplot(gs[2,0]),cons,"17gy_early_vs_late","D  Early vs late 17 Gy (cross-sectional)")
    ax=fig.add_subplot(gs[2,1]); ax.set_title("E  Representation robustness",loc="left",fontweight="bold")
    hm=full[(full.hierarchy=="major")&full.contrast.isin(list(CONTRAST_LABELS)[:3])].pivot_table(index=["contrast","state"],columns="representation",values="normalized_delta_norm").reindex(columns=REPS)
    im=ax.imshow(hm.values,aspect="auto",cmap="viridis"); ax.set_xticks(range(5),["PCA30\nprimary","PCA50\nprimary","scVI\nprimary","PCA20\ndimensional","Diffusion\nexploratory"])
    ax.set_yticks(range(len(hm)),[f"{CONTRAST_LABELS[c]} | {s}" for c,s in hm.index],fontsize=7); fig.colorbar(im,ax=ax,label="Normalized displacement")
    ax=fig.add_subplot(gs[3,0]); ax.set_title("G  Mouse-level abundance (separate evidence)",loc="left",fontweight="bold")
    ab=abundance[(abundance.contrast=="control_vs_17gy_late")&(abundance.hierarchy=="major")]
    xpos={s:i for i,s in enumerate(MAJOR)}
    for group,color,offset in [("control","#4c78a8",-.13),("17gy_late","#e45756",.13)]:
        z=ab[ab.group==group]; ax.scatter([xpos[s]+offset for s in z.state],z.state_fraction,s=20,color=color,alpha=.75,label=group)
    ax.set_xticks(range(6),MAJOR,rotation=30,ha="right"); ax.set_ylabel("State proportion per mouse"); ax.legend(frameon=False); ax.grid(axis="y",alpha=.2)
    ax=fig.add_subplot(gs[3,1]); ax.set_title("H  Fibroblast subtype sensitivity",loc="left",fontweight="bold")
    fb=cons[(cons.hierarchy=="fibroblast")&cons.contrast.isin(list(CONTRAST_LABELS)[:3])].pivot(index="state",columns="contrast",values="normalized_magnitude_median").reindex(index=FIBRO,columns=list(CONTRAST_LABELS)[:3])
    im=ax.imshow(fb.values,aspect="auto",cmap="magma"); ax.set_yticks(range(3),FIBRO); ax.set_xticks(range(3),["Early vs control","Late vs control","Early vs late"],rotation=20,ha="right")
    for i,s in enumerate(FIBRO):
        for j,c in enumerate(list(CONTRAST_LABELS)[:3]):
            lab=cons[(cons.state==s)&(cons.contrast==c)].consensus_label.iloc[0]
            ax.text(j,i,"stable" if lab=="representation_stable_effect" else "unstable",ha="center",va="center",color="white",fontsize=8)
    fig.colorbar(im,ax=ax,label="Normalized displacement")
    fig.suptitle("GSE211713 Dataset C: replicate-aware, representation-robust radiation geometry",fontsize=17,fontweight="bold")
    source=pd.concat([
        cons.assign(source_table="primary_consensus_state_evidence"),
        full.assign(source_table="full_state_evidence"),
        centers.assign(source_table="individual_mouse_centers"),
        abundance.assign(source_table="mouse_level_abundance"),
        perm.assign(source_table="exact_mouse_permutations"),
    ],ignore_index=True,sort=False)
    save_figure(fig,"dataset_c_main_figure",source,"Eight-panel Dataset C figure. The design panel identifies five controls, four early 17 Gy mice, and four late 17 Gy mice as independent cross-sectional libraries. Early displacement is representation-sensitive in all six major states. Late displacement is stable in five states while fibroblast/stromal is unstable. All six major states are stable in the cross-sectional early-versus-late comparison. Remaining panels show representation-specific effects, individual mouse centers, mouse-level abundance, and fibroblast subtype sensitivity. No composite score is used.")
    return source


def simple_supplements(data: dict[str, pd.DataFrame]) -> None:
    specs=[]
    qc=data["qc"].copy(); specs.append(("dataset_c_supp01_qc_and_design",qc,"retained_cells","gsm","irradiation_group","Retained cells per independent library","QC and design summary for 20 independent mouse/GSM libraries."))
    ann=data["annotation"].copy(); specs.append(("dataset_c_supp02_annotation_confidence",ann,"cells","annotation_state","annotation_confidence","Annotation confidence and marker-inferred states","Cell counts by conservative marker-inferred label and confidence."))
    reps=data["reps"].copy(); specs.append(("dataset_c_supp03_representation_quality",reps,"neighborhood_overlap_vs_pca50","representation","role","Representation-quality metrics","Neighborhood overlap with PCA50 for five representations; dose and time silhouettes are descriptive diagnostics, not correction objectives."))
    sens=data["full"][data["full"].representation.isin(["X_pca20","X_diffmap"])].copy(); specs.append(("dataset_c_supp04_pca20_diffusion_sensitivity",sens,"normalized_delta_norm","state","representation","PCA20 and diffusion sensitivity","Normalized displacement in the dimensional and exploratory sensitivity representations."))
    cov=data["coverage"].copy(); specs.append(("dataset_c_supp05_refined_state_coverage",cov,"eligible_mice_group_b","state","contrast","Eligible refined-state coverage","Frozen C1–C5 eligibility counts; this panel does not add numerical ScGeo analyses."))
    ctr=data["centers"].copy(); specs.append(("dataset_c_supp06_mouse_center_heterogeneity",ctr,"coordinate_1","mouse_id","group","Individual mouse-center heterogeneity","First coordinate of each stored geometric-median mouse center; coordinates are representation-specific."))
    boot=data["bootdraws"].copy(); specs.append(("dataset_c_supp07_mouse_bootstrap",boot,"normalized_displacement","state","contrast","Biological-mouse bootstrap distributions","Stored normalized displacement draws from 500 biological-mouse bootstrap iterations per eligible state and representation."))
    per=data["perm_values"].copy(); specs.append(("dataset_c_supp08_exact_permutation",per,"normalized_displacement","state","contrast","Exact mouse-label permutation values","All exact mouse-label assignments; observed assignments are retained, ties count greater than or equal, and BH values are sensitivity only."))
    ab=data["abundance"].copy(); specs.append(("dataset_c_supp09_mouse_abundance",ab,"state_fraction","state","group","Mouse-level abundance","Individual mouse state proportions, maintained separately from displacement."))
    dist=data["distribution"].copy(); specs.append(("dataset_c_supp10_distribution_metrics",dist,"energy_distance","state","comparison_unit","Distribution metrics by comparison unit","Energy distance with pooled-cell and biological-mouse-center comparison units explicitly separated; these are descriptive distances."))
    sec=data["full"][data["full"].contrast.isin(["control_vs_all_17gy","control_vs_all_10gy"])].copy(); specs.append(("dataset_c_supp11_secondary_contrasts",sec,"normalized_delta_norm","state","contrast","Secondary time-heterogeneous contrasts","Secondary all-17 Gy and all-10 Gy contrasts collapse heterogeneous times and are labeled replicate-aware secondary time-heterogeneous."))
    neg=data["negative"].copy(); specs.append(("dataset_c_supp12_neutral_unstable",neg,"normalized_magnitude_median","state","contrast","Neutral, inconclusive, and unstable states","Every canonical neutral or representation-unstable state remains visible."))
    invalid=pd.DataFrame([json.loads((BASE/"invalidated_runs/c6_c7_model_switch_run/invalidation_manifest.json").read_text())]).astype(str)
    specs.append(("dataset_c_supp13_invalidated_run_provenance",invalid,None,None,None,"Invalidated-run provenance","Reproducibility-only summary of the quarantined model-switch run; none of its numerical artifacts are manuscript sources."))
    for stem,frame,value,label,hue,title,alt in specs:
        fig,ax=plt.subplots(figsize=(10,5.5),layout="constrained")
        ax.set_title(title,loc="left",fontweight="bold")
        if value is None or value not in frame or frame.empty:
            ax.axis("off"); ax.text(.02,.92,alt,va="top",wrap=True,fontsize=11)
        else:
            plot=frame.dropna(subset=[value]).copy()
            if len(plot)>2500: plot=plot.iloc[np.linspace(0,len(plot)-1,2500,dtype=int)]
            labs=plot[label].astype(str) if label in plot else pd.Series(np.arange(len(plot)).astype(str))
            top=labs.value_counts().head(16).index; plot=plot[labs.isin(top)].copy(); labs=plot[label].astype(str)
            cats=list(dict.fromkeys(labs)); xmap={v:i for i,v in enumerate(cats)}
            if hue and hue in plot:
                hues=list(dict.fromkeys(plot[hue].astype(str))); cmap=plt.get_cmap("tab10")
                for j,h in enumerate(hues):
                    z=plot[plot[hue].astype(str)==h]; ax.scatter([xmap[x] for x in z[label].astype(str)],z[value],s=13,alpha=.55,label=h,color=cmap(j%10))
                if len(hues)<=10: ax.legend(frameon=False,fontsize=7,ncol=2)
            else: ax.scatter([xmap[x] for x in labs],plot[value],s=13,alpha=.55,color="#32688e")
            ax.set_xticks(range(len(cats)),cats,rotation=35,ha="right",fontsize=7); ax.set_ylabel(value.replace("_"," ")); ax.grid(axis="y",alpha=.2)
        save_figure(fig,stem,frame,alt)


def ledgers(data: dict[str, pd.DataFrame], source_hashes: dict[str,str]) -> None:
    cons=data["consensus"]
    claims=[
        ("Early radiation remodeling was representation-sensitive under the frozen consensus rule.","control_vs_17gy_early","major","all six major states: representation_unstable","Results; main Fig. B","Small mouse groups and representation dependence; association only."),
        ("Late radiation showed representation-stable remodeling in five of six major compartments.","control_vs_17gy_late","major","stable: Epithelial, Endothelial, Proliferating, Myeloid, Lymphoid; Fibroblast/stromal unstable","Results; main Fig. C","Cross-sectional association; no causal claim."),
        ("All six major compartments differed stably between early and late irradiated mice.","17gy_early_vs_late","major","six of six representation_stable_effect","Results; main Fig. D","Different mice were sampled; this is not longitudinal change."),
        ("Fibroblast subtype conclusions were more representation-sensitive.","three primary contrasts","fibroblast","early: all unstable; late: all unstable; early-vs-late: Col13a1-like and myofibroblast stable, Col14a1-like unstable","Results/Supplement","Marker-inferred subtypes and lower 20-cell exploratory coverage threshold."),
        ("Mouse-level replication supports replicate-aware association.","three primary contrasts","major and fibroblast","5 vs 4, 5 vs 4, and 4 vs 4 independent mice; 500 mouse bootstraps","Methods","Small-n uncertainty and discrete exact-permutation grids."),
        ("Exact permutation resolution is limited by design.","three primary contrasts","major","126, 126, and 70 assignments; minimum grids 1/126 and 1/70","Limitations/Supplement","Equal 4-vs-4 assignments have complement symmetry."),
    ]
    rows=[]
    for claim,contrast,hier,value,loc,lim in claims:
        rows.append({"proposed_claim":claim,"supporting_output":"scgeo_v2/primary_consensus_state_evidence.csv; exact_mouse_permutations.csv","exact_value":value,"contrast_scope":contrast,"representation_scope":"PCA30, PCA50, scVI primary; PCA20/diffusion sensitivity","replicate_scope":"mouse_id/GSM biological samples","inferential_status":"replicate-aware association; no causality","limitation":lim,"appropriate_manuscript_location":loc,"reviewer_concern_addressed":"real-data generalization; representation robustness; biological replication"})
    atomic_csv(OUT/"dataset_c_evidence_ledger.csv",pd.DataFrame(rows))
    limitations={"claims":[r[0] for r in claims],"limitations":["Cross-sectional design; no longitudinal reversal, persistence, or causality.","Marker-inferred annotations are uncertain, especially fibroblast subtypes.","Exact permutation resolution is limited by 5-vs-4 and 4-vs-4 designs.","Related contrasts reuse control mice and are not independent.","PCA20/PCA30/PCA50 are nested views, not independent confirmations.","Distribution distances are descriptive and are not substitutes for biological-replicate inference."],"inference_status":{"primary":"replicate_aware_primary","secondary":"replicate_aware_secondary_time_heterogeneous"}}
    atomic_text(OUT/"dataset_c_claims_and_limitations.json",json.dumps(limitations,indent=2)+"\n")
    methods=pd.DataFrame([{"stage":"C0-C5","method":"public-study audit, sparse QC, condition-blind marker annotation","source":"frozen audit/QC/annotation/coverage artifacts","role":"design, cell retention, annotation and prespecified eligibility"},{"stage":"C6","method":"PCA20/30/50, diffusion, scVI; UMAP display only","source":"gse211713_revision_representations_v2.h5ad","role":"representation ensemble; no dose/time/mouse batch removal"},{"stage":"C7","method":"ScGeo mouse centers, 500 mouse bootstraps, exact mouse-label permutations","source":"scgeo_v2 canonical tables","role":"replicate-aware association and sensitivity"},{"stage":"assembly","method":"checksum validation and plotting only","source":"06 manuscript notebook and assembly script","role":"no numerical recomputation"}])
    atomic_csv(OUT/"dataset_c_methods_provenance.csv",methods)
    reviewer=pd.DataFrame([{"reviewer_concern":"single-dataset/generalization","Dataset_C_relevance":"public radiation-response mouse-lung data extend evaluation beyond recovery and Dataset B","boundary":"association, not causality"},{"reviewer_concern":"representation dependence","Dataset_C_relevance":"five representations with prespecified primary and sensitivity roles","boundary":"nested PCA views are not independent evidence"},{"reviewer_concern":"biological replication","Dataset_C_relevance":"mouse/GSM centers and mouse-level resampling","boundary":"small group sizes and reused controls"},{"reviewer_concern":"local geometry/topology","Dataset_C_relevance":"neighborhood overlap, distortion, graph connectivity and label preservation","boundary":"diagnostics do not establish universal preservation"}])
    atomic_csv(OUT/"dataset_c_reviewer_relevance.csv",reviewer)
    caption="""# Dataset C main figure caption\n\n**GSE211713 replicate-aware radiation geometry.** (A) Cross-sectional design with five independent controls, four early 17 Gy mice (months 1–2), and four late 17 Gy mice (months 4–5). (B) No major state met the frozen representation-stable consensus for early 17 Gy versus control. (C) Epithelial, Endothelial, Proliferating, Myeloid, and Lymphoid states met the late-versus-control stable consensus; Fibroblast/stromal remained representation-unstable. (D) All six major states met the stable consensus in the cross-sectional early-versus-late comparison. (E) Primary consensus used PCA30, PCA50, and scVI; PCA20 and diffusion were dimensional and exploratory sensitivities. (F) Individual geometric-median mouse centers, 500 biological-mouse bootstrap draws, and exhaustive mouse-label permutations expose sample-level evidence and its discrete resolution. (G) Mouse-level abundance is shown separately from displacement. (H) Fibroblast subtypes were more representation-sensitive. Results are replicate-aware associations, not longitudinal or causal effects; related contrasts reuse control mice and are not independent.\n"""
    atomic_text(OUT/"dataset_c_main_figure_caption.md",caption)
    captions = ["QC and independent-library design.","Conservative annotation confidence and marker provenance.","Representation-quality diagnostics; condition mixing is not an optimization target.","PCA20 dimensional and diffusion exploratory sensitivity.","Frozen refined-state coverage qualification.","Individual mouse-center heterogeneity.","Biological-mouse bootstrap distributions.","Exact permutation values and within-contrast/representation BH sensitivity.","Individual mouse state proportions.","Descriptive distribution distances with comparison units exposed.","Secondary time-heterogeneous all-17 Gy and all-10 Gy contrasts.","Neutral, inconclusive and representation-unstable findings.","Invalidated-run provenance; no invalidated numerical artifact enters manuscript evidence."]
    supp = "# Dataset C supplementary captions\n\n" + "\n".join(
        f"{i}. **Supplement {i}.** {caption}" for i, caption in enumerate(captions, 1)
    ) + "\n"
    atomic_text(OUT/"dataset_c_supplementary_captions.md",supp)


def assemble() -> dict:
    mpl.rcParams.update({"svg.hashsalt":"gse211713-dataset-c-v1","font.size":9,"axes.spines.top":False,"axes.spines.right":False})
    before=pin_check()
    data={
        "full":pd.read_csv(SCGEO/"full_state_evidence.csv"),
        "consensus":pd.read_csv(SCGEO/"primary_consensus_state_evidence.csv"),
        "negative":pd.read_csv(SCGEO/"negative_neutral_or_unstable_states.csv"),
        "centers":pd.read_csv(SCGEO/"individual_mouse_centers.csv"),
        "boot":pd.read_csv(SCGEO/"biological_mouse_bootstrap_intervals.csv"),
        "bootdraws":pd.read_csv(SCGEO/"biological_mouse_bootstrap_draws.csv"),
        "perm":pd.read_csv(SCGEO/"exact_mouse_permutations.csv"),
        "perm_values":pd.read_csv(SCGEO/"exact_mouse_permutation_values.csv"),
        "perm_bh":pd.read_csv(SCGEO/"exact_mouse_permutation_bh_sensitivity.csv"),
        "abundance":pd.read_csv(SCGEO/"mouse_level_abundance.csv"),
        "distribution":pd.read_csv(SCGEO/"distribution_metrics.csv"),
        "eligibility":pd.read_csv(SCGEO/"eligibility_and_coverage_used.csv"),
        "reps":pd.read_csv(BASE/"c6_v2/metadata/representation_metrics.csv"),
        "qc":pd.read_csv(BASE/"qc/qc_summary_by_mouse.csv"),
        "annotation":pd.read_csv(BASE/"annotation/annotation_confidence_table.csv"),
        "markers":pd.read_csv(BASE/"annotation/marker_definitions.csv"),
        "coverage":pd.read_csv(BASE/"coverage/contrast_state_eligibility.csv"),
        "design":pd.read_csv(BASE/"audit/study_design.csv"),
    }
    manuscript_figure(data)
    simple_supplements(data)
    ledgers(data,before)
    after=pin_check()
    if before != after: raise RuntimeError("Pinned numerical inputs changed during assembly")
    outputs={str(p.relative_to(OUT)):sha256(p) for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="dataset_c_manuscript_assembly_validation.json"}
    validation={"status":"passed","assembly_kind":"presentation_only_no_numerical_recomputation","source_hashes_before":before,"source_hashes_after":after,"numerical_sources_unchanged":before==after,"invalidated_numerical_sources_used":False,"invalidated_manifest_used_for_provenance_only":True,"figures_png":len(list((OUT/"figures").glob("*.png"))),"figures_svg":len(list((OUT/"figures").glob("*.svg"))),"figure_source_csv":len(list((OUT/"figure_sources").glob("*.csv"))),"alt_text_files":len(list((OUT/"alt_text").glob("*.txt"))),"output_checksums":outputs}
    atomic_text(OUT/"dataset_c_manuscript_assembly_validation.json",json.dumps(validation,indent=2,sort_keys=True)+"\n")
    return validation


if __name__ == "__main__":
    print(json.dumps(assemble(),indent=2,sort_keys=True))
