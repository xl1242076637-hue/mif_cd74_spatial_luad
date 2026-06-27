# Reproducibility Notes

This file records the current command-level route for reproducing the main analysis products. Run commands from the project root:

```powershell
cd '<PROJECT_ROOT>'
$env:PYTHONPATH='<PROJECT_ROOT>\src'
```

## Current main outputs

| Output | Path |
|---|---|
| Dataset inventory | `docs/data_inventory.md` |
| Dataset completeness summary | `docs/dataset_completeness_summary_2026-05-30.md` |
| Novelty audit | `docs/novelty_audit_2026-05-30.md` |
| Literature overlap update | `docs/literature_overlap_update_2026-05-30.md` |
| Reference audit | `docs/reference_audit_2026-05-31.md` |
| Integrated manuscript draft | `docs/manuscript_draft.md` |
| Submission-style manuscript draft | `docs/manuscript_submission_draft.md` |
| Compact submission draft | `docs/manuscript_compact_submission_draft.md` |
| Prioritized next actions | `docs/manuscript_next_actions.md` |
| Communications Biology working draft | `docs/manuscript_communications_biology_draft.md` |
| Communications Biology Chinese review draft | `docs/manuscript_communications_biology_draft_zh.md` |
| Communications Biology Word review draft | `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx` |
| Word export notes | `docs/word_export_notes.md` |
| GRN virtual perturbation method note | `docs/grn_virtual_perturbation_2026-06-03.md` |
| Perturb-seq/GEARS/scGen/CPA feasibility note | `docs/perturbseq_gears_scgen_cpa_feasibility_2026-06-03.md` |
| Results draft | `docs/results_draft.md` |
| Methods draft | `docs/methods_draft.md` |
| Figure legends | `docs/figure_legends.md` |
| Claim-evidence-boundary map | `docs/manuscript_claim_evidence_map.md` |
| Editorial summary and cover-letter seed | `docs/editorial_summary_cover_letter_seed.md` |
| Target-journal formatting route | `docs/target_journal_formatting_route.md` |
| Communications Biology submission-readiness audit | `docs/communications_biology_submission_readiness_audit.md` |
| Patient-aware spatial statistics | `docs/gse307534_spatial_statistics_summary.md` |
| MIF spatial controls | `docs/gse307534_mif_spatial_controls.md` |
| MIF covariate sensitivity | `docs/gse307534_mif_covariate_sensitivity.md` |
| Supplementary-table index | `docs/supplementary_tables_index.md` |
| Main evidence matrix | `results/tables/main_axis_evidence_matrix.csv` |
| Communications Biology reference list | `results/tables/communications_biology_reference_list.csv` |
| Zotero-ready RIS reference export | `results/references/communications_biology_references.ris` |
| GSE308103 GRN-level virtual perturbation tables | `results/tables/gse308103_grn_virtual_perturbation_*.csv` and `results/tables/gse308103_grn_virtual_perturbation_config.json` |
| Main figure source data | `results/tables/figure*_source.csv` |
| Main figures | `results/figures/figure1_*` through `results/figures/figure5_*`; redesigned main figures under `results/figures/nature_redesign/nature_figure1_*` through `nature_figure5_*` |
| Supplementary GRN virtual perturbation figure | `results/figures/supplementary_figure_grn_virtual_perturbation.*` |
| Supplementary GRN virtual perturbation source data | `results/tables/supplementary_figure_grn_virtual_perturbation_*_source.csv` |
| GRN cross-dataset validation tables | `results/tables/grn_cross_dataset_signature_validation*.csv` |
| Supplementary specificity-refinement figure | `results/figures/supplementary_figure_spp1_signature_refinement.*` |
| Supplementary specificity-refinement source data | `results/tables/supplementary_figure_spp1_signature_refinement_*_source.csv` |

## Rebuild selected analyses

These commands assume the raw and interim data have already been downloaded and extracted.

```powershell
python scripts\score_gse308103_snrna_states.py
python scripts\score_gse131907_selected_signatures.py
python scripts\audit_gse189357_signatures_with_gse131907.py
python scripts\score_gse164789_scrna_states.py
python scripts\summarize_gse282617_markers.py
python scripts\summarize_gse308103_candidate_genes.py
python scripts\analyze_gse307534_candidate_axis_spatial.py
python scripts\summarize_gse307534_spatial_statistics.py
python scripts\analyze_gse307534_mif_controls.py
python scripts\analyze_gse307534_mif_covariates.py
python scripts\rank_candidate_mechanisms.py
python scripts\virtual_perturb_gse307534_axes.py
python scripts\continuous_perturb_gse307534_axes.py
python scripts\run_gse308103_grn_virtual_perturbation.py
python scripts\plot_supplementary_grn_virtual_perturbation.py
python scripts\run_gse308103_grn_robustness.py
python scripts\build_grn_cross_dataset_validation.py
python scripts\build_main_evidence_matrix.py
python scripts\plot_main_figures.py
python scripts\plot_main_figures_nature.py
python scripts\plot_supplementary_signature_refinement.py
python scripts\plot_supplementary_orthogonal_validation.py
python scripts\export_supplementary_tables.py
```

## Rebuild main manuscript figures only

```powershell
python scripts\build_main_evidence_matrix.py
python scripts\plot_main_figures.py
python scripts\plot_main_figures_nature.py
```

Figure outputs are written to `results/figures/` as SVG, PDF, TIFF, and PNG. The redesigned main figure set is written to `results/figures/nature_redesign/`. Figure source-data CSV files are written to `results/tables/`.

## Rebuild supplementary validation figures only

```powershell
python scripts\summarize_gse308103_candidate_genes.py
python scripts\plot_supplementary_signature_refinement.py
python scripts\plot_supplementary_orthogonal_validation.py
python scripts\plot_supplementary_grn_virtual_perturbation.py
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' --vanilla scripts\run_original_sctenifoldknk_expanded.R
python scripts\summarize_original_sctenifoldknk_expanded.py
python scripts\export_supplementary_tables.py
```

Specificity-refinement outputs are written as `results/figures/supplementary_figure_spp1_signature_refinement.*` and source-data CSVs under `results/tables/`. Focused orthogonal-validation outputs are written as `results/figures/supplementary_figure_focused_orthogonal_validation.*` and source-data CSVs under `results/tables/`. GRN-level virtual perturbation outputs are written as `results/figures/supplementary_figure_grn_virtual_perturbation.*` and source-data CSVs under `results/tables/`.

The original CRAN `scTenifoldKnk` expanded sensitivity run writes per-target
`diffRegulation` tables under `results/tables/sctenifoldknk_original_expanded/`
and the manuscript-facing boundary summary
`results/tables/sctenifoldknk_original_expanded_interpretation.csv`.

## Rebuild Chinese Word review draft

```powershell
python scripts\export_word_manuscript.py
```

This writes `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx` and `results/references/communications_biology_references.ris`. The RIS export is built from the current reference registry in `results/tables/communications_biology_reference_list.csv`, which now contains 36 references after adding GRN and perturbation-response model method papers.

After Zotero MCP Web API configuration, the same command also writes:

```text
docs/manuscript_communications_biology_draft_zh_zotero_fields.docx
```

The existing Zotero-field draft uses imported Zotero item URIs for references 1-30. References 31-36 were added later for GRN and perturbation-response model context; import the updated RIS into Zotero before regenerating a fully Zotero-linked Word draft. Inspecting `word/document.xml` should show `ADDIN ZOTERO_ITEM CSL_CITATION` field instructions and `http://zotero.org/users/13378991/items/...` URIs for mapped items.

## Run tests

```powershell
$env:PYTHONPATH='src'; pytest tests -q
```

Latest recorded full test result:

```text
126 passed in 13.54s
```

## Current interpretation checkpoint

The current computational conclusion is:

1. `MIF-CD74/CXCR4` is the top-ranked epithelial-myeloid communication candidate.
2. `SPP1/TREM2/PLA2G7` is retained as the macrophage-state readout, not the primary communication novelty.
3. Paired-patient GSE307534 statistics support a source-side `MIF` progression increase.
4. Expression-matched controls retain `MIF` as a prioritized candidate but do not show unique specificity.
5. Covariate sensitivity models retain a positive source-side `MIF` late-lesion effect after adjustment for spot count, permutation-null mean, and radius.
6. Score-level in-silico target prioritization ranks receptor-side `CD74` above `CD44` and `CXCR4`; source-side `MIF` dropout is expected because it collapses a single-gene source score.
7. GSE308103 GRN-level virtual perturbation prioritization stratifies target neighborhoods: `CD74` maps to a C1Q/APOE/ACP5 macrophage network, while `CXCR4` and `PLA2G7` map more strongly to SPP1-like macrophage signatures; this is exploratory network prioritization, not wet-lab perturbation.
8. Original CRAN `scTenifoldKnk` runs locally and has been added as a reduced-panel sensitivity/boundary analysis. It supports epithelial `MIF` network-context sensitivity, but several macrophage-target runs include epithelial markers among top non-target genes, so it should not replace the transparent Python GRN-level implementation.
9. `IL1B/TNF/CXCL8` should be treated as a benchmark inflammatory axis because similar early-LUAD inflammatory niche biology has already been published.

## Environment and Dependency Notes

1. Python Analysis Environment

   All Python scripts in this project are developed and tested under Python 3.10.14.
   All required Python packages with fixed versions are listed in the root  requirements.txt  file.
   Execute the command below in the repository root directory to install all dependencies in one step:

   `bash`

   `pip install -r requirements.txt`


   The installed environment supports all pipelines including GEO data download, spatial transcriptomics h5ad file processing, cell signature scoring, spatial coupling calculation, GRN virtual perturbation, marker gene identification, and manuscript figure generation.

2. R Script Supplementary Description

   Two R scripts are provided for scTenifoldKnk related gene network analysis:
    run_original_sctenifoldknk_expanded.R 
    run_original_sctenifoldknk_smoke.R 
   Required R packages: Seurat, sctenifoldknk, ggplot2, dplyr, tidyr, patchwork
   These R dependencies are not recorded in  requirements.txt ; please install them manually before running R workflows.



