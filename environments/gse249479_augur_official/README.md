# Isolated official Augur validation environment

This environment validates the GSE249479 Dataset B Augur-inspired Python
approximation against the official R package without modifying the Python
environment or any Phase 3C result.

- R: 4.1.2
- installation mechanism: `renv`
- official package source: `neurorestore/Augur`
- pinned source commit: `b252b84e4af687d9817813b1db409267eb44ec3f`
- package version at that commit: 1.0.3
- runtime library: ignored validation output subtree

The committed `renv.lock` pins packages installed in the isolated library.
Several already-installed workstation site-library dependencies were reused by
R 4.1.2; their exact versions and library paths are recorded in the generated
`official_augur_dependency_versions.csv` and JSON validation artifacts. Restore
into a separate project library and do not activate it in the repository root.

All validation results remain `descriptive_only`. Augur cross-validation and
subsampling quantify computational cell-level stability, not biological
replication.
