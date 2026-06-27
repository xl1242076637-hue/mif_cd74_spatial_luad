# Early LUAD MIF-CD74 Epithelial-Myeloid Spatial Niche

This project records a reproducible pure-bioinformatics workflow for studying early lung adenocarcinoma (LUAD) progression, with a focus on specificity-audited epithelial-myeloid spatial communication.

## Working Hypothesis

Current analyses nominate `MIF-CD74` as a leading epithelial-myeloid communication candidate during early LUAD progression. `SPP1/TREM2/PLA2G7` is retained as a macrophage-state readout, while `IL1B/TNF/CXCL8` is treated as a benchmark inflammatory axis because broad epithelial-inflammatory niches have already been described in recent early-LUAD spatial-omics studies.

The project uses public datasets, cell-type specificity auditing, spatial neighborhood scoring, multi-cohort evidence ranking, and score-level in-silico perturbation. Virtual perturbation results are interpreted as target prioritization, not causal knockout proof.

## Current Scope

This phase is pure bioinformatics:

1. Collect public early LUAD single-cell, single-nucleus, spatial transcriptomic, and bulk RNA-seq datasets.
2. Build a reproducible data inventory and download workflow.
3. Audit candidate signatures against a large LUAD single-cell reference.
4. Score epithelial-myeloid candidate axes in ordered spatial data.
5. Rank axes across spatial, specificity, bulk, scRNA/snRNA, and perturbation evidence.
6. Use score-level in-silico perturbation to prioritize follow-up targets.

## Current Outputs

- `docs/manuscript_draft.md`: integrated internal manuscript draft.
- `docs/manuscript_next_actions.md`: prioritized manuscript and analysis gaps.
- `docs/results_draft.md`: manuscript-style results draft.
- `docs/methods_draft.md`: methods draft.
- `docs/figure_legends.md`: figure legends draft.
- `docs/reproducibility.md`: command-level reproducibility notes.
- `docs/gse307534_spatial_statistics_summary.md`: patient-aware spatial statistics summary.
- `docs/gse307534_mif_spatial_controls.md`: expression-matched and tissue-density control summary.
- `docs/supplementary_tables_index.md`: packaged supplementary-table index.
- `results/figures/figure1_*` through `results/figures/figure5_*`: current main figure drafts.
- `results/tables/main_axis_evidence_matrix.csv`: integrated axis evidence matrix.
- `results/supplementary_tables/`: reproducibly packaged supplementary CSV files.

## Project Layout

- `config/datasets.yaml`: dataset manifest and source URLs.
- `data/raw/`: downloaded original files.
- `data/interim/`: converted or normalized intermediate files.
- `data/processed/`: final reusable analysis objects.
- `docs/decision_log.md`: human-readable rationale for major choices.
- `docs/analysis_log.md`: chronological command and result log.
- `docs/reproducibility.md`: command-level route for rebuilding major outputs.
- `scripts/`: command-line utilities.
- `src/luad_niche/`: reusable Python helpers.
- `tests/`: tests for reproducible utilities.
- `results/figures/`: output figures.
- `results/tables/`: output tables and inventories.

## Current Environment

The current machine has Python 3.13 with pandas, numpy, scipy, sklearn, matplotlib, seaborn, requests, PyYAML, and pytest. `scanpy`, `anndata`, and `Rscript` are not currently installed, so the first implementation avoids depending on them.
