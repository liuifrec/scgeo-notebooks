# ScGeo 0.2.0 — executable tutorials

Learn the published [ScGeo 0.2.0 API](https://github.com/liuifrec/scgeo/tree/v0.2.0)
using small, reproducible examples. This track is separate from the original
paper's notebooks, environments and numerical evidence. It does not run private
data, download manuscript data, or recreate published biological results.

**First, understand the input:** [How the tutorial data are generated](../../../docs/SIMULATION_DESIGN.md) explains the synthetic geometry, gene programs, graph fixtures and sample-count design in pictures and equations. Lessons 00–10 use simulated inputs, not experimental sequencing data. Lesson 11 is optional public-data rendering.

## Install in a separate environment

From the companion repository root, using Python 3.12:

```bash
python -m venv .venv-tutorials
# Linux/macOS:
source .venv-tutorials/bin/activate
# Windows PowerShell: .venv-tutorials\Scripts\Activate.ps1
python -m pip install -r environment/scgeo_v020_tutorials.txt
python -m pip check
python -m pip install "scgeo[notebook]==0.2.0"
python -m ipykernel install --sys-prefix --name scgeo-v020-tutorial --display-name "ScGeo 0.2.0 tutorials"
jupyter lab
```

Select the new tutorial kernel. The notebooks verify the installed package
version and the exact release example sources. They write a fresh ignored
results directory each time the setup cell is run. Do not substitute the new
environment for one used to reproduce the paper.

For headless execution (JupyterLab is not required):

```bash
python scripts/execute_scgeo_v020_tutorials.py --validate-only
python scripts/execute_scgeo_v020_tutorials.py --kernel scgeo-v020-tutorial
# Select a subset if desired:
python scripts/execute_scgeo_v020_tutorials.py --lesson 08 09 10 --kernel scgeo-v020-tutorial
```

For the headless commands, register the named kernel first using the ipykernel
command above, or use `--kernel python3` if that kernel points to this environment.
The runner requires each kernel to use its own Python interpreter, validates
notebook schemas and code, keeps partial executed copies on failure, and checks
all ScGeo PlotResult manifests and their PNG/SVG/PDF and CSV SHA256 hashes.
Run each notebook independently; none depends on output from another lesson.

## Learning sequence

| Notebook | Main topics | Dataset |
| --- | --- | --- |
| [00 Input contract](00_input_contract.ipynb) | Read-only input coverage, sample-aware analysis, six 2-vs-2 allocations | Synthetic |
| [01 Robust shift](01_robust_shift_report.ipynb) | State evidence, sample centers, replicate effects | Synthetic |
| [02 Representation and contrasts](02_representation_and_contrasts.ipynb) | Heatmaps, local distortion, robustness, stored alignment, comparison intervals | Synthetic |
| [03 PAGA maps](03_paga_shift_maps.ipynb) | Graph order, endpoints, stored alignment and focus | Constructed graph |
| [04 Velocity fields](04_velocity_fields.ipynb) | Grids, streams, local support, matched-condition comparisons | Constructed vectors |
| [05 Associations and density](05_associations_density_legacy.ipynb) | Alluvial, chord, Sankey, density, legacy mapping/consensus plots | Synthetic labels |
| [06 Facets and composition](06_facets_and_composition.ipynb) | Shared scales, two-way faceting, external axes, combined figure sources | Synthetic |
| [07 Atlas](07_geometry_dynamics_atlas.ipynb) | Minimal and full integrated geometry–dynamics atlas | Synthetic |
| [08 Shift explanations](08_shift_explanations.ipynb) | Affine signed attribution, genes/programs, local responses, curves and connectors | Synthetic expression |
| [09 Sparse states](09_sparse_states.ipynb) | Zero/singleton control coverage, observed expression, PAGA pies, parent analysis | Synthetic |
| [10 Compass](10_stored_compass.ipynb) | Stored evidence, concordance, sample abundance, coverage, display invariance | Synthetic |
| [11 Optional public data](11_public_data_optional.ipynb) | Checksum-verified PBMC/pancreas rendering, faceting and atlas | Local public files only |

Lessons 00–10 are download-free after the environment is installed. Lesson 11
records `not_run` when neither public file path is supplied; that does not count
as executed public-data validation. To opt in, obtain the exact files listed in
the [release public-data guide](https://github.com/liuifrec/scgeo/blob/v0.2.0/docs/REAL_DATA_PLOTTING.md)
and set `SCGEO_TUTORIAL_PBMC_H5AD` and/or `SCGEO_TUTORIAL_PANCREAS_H5AD` to their
local paths before starting the kernel, or edit the notebook path dictionary.
The preparation helper verifies the file size and SHA256 and never downloads
by default. Only public files matching the pinned checksums are accepted.

## Source and output provenance

The companion uses the PyPI release [0.2.0](https://pypi.org/project/scgeo/0.2.0/),
source commit `46ba9017ceef17255292edd0140602ed13ac5305`. The
[example snapshots](../../../scripts/scgeo_v020_examples/) preserve the exact
13 example scripts from that commit (see the two SOURCE_PART JSON inventories
for Git blob identities). Their original development-era CLI wording is
retained to keep the bytes auditable; the tutorial requirements use the
published wheel. These are examples, not a vendored ScGeo library.

Figures and tables are generated locally under
`results/scgeo_v020_tutorials/`; the runner writes separate executed notebook
copies and an execution summary. Notebook sources are output-free. Each lesson
records package versions, interpreter, companion commit and example provenance.
Source snapshot identity and a version pin do not lock transitive dependencies;
retain the resolved environment for reproducible work. CI additionally verifies
the exact production wheel SHA256 before executing the lessons.

## Interpret the results

Use quantitative representations for estimates and embeddings for display.
Observe biological sampling units, distinguish composition from displacement,
do not infer a missing reference or causal driver, and distinguish a supplied
velocity field from a fitted dynamics model. Synthetic outputs are software
examples, not biological findings. For paper reproduction, use the original
[frozen checkpoints and environments](../../../docs/REPRODUCIBILITY_GUIDE.md).
