# Public Data Inventory

Updated: 2026-05-30

## First-Wave Datasets

| Accession | Role | Local status | Remote files | Notes |
|---|---|---|---|---|
| `GSE189357` | scRNA-seq reference for AIS/MIA/IAC epithelial and macrophage states | SOFT metadata, `filelist.txt`, full `GSE189357_RAW.tar`, and extracted raw 10x files downloaded | `GSE189357_RAW.tar` approx. 623.90 MB | 9 samples: AIS n=3, MIA n=3, IAC n=3 |
| `GSE189487` | spatial transcriptomics for AIS/MIA/IAC niche mapping | SOFT metadata, `filelist.txt`, full `GSE189487_RAW.tar`, and extracted raw 10x files downloaded | `GSE189487_RAW.tar` approx. 209.54 MB | 6 samples: AIS n=2, MIA n=2, IAC n=2 |
| `GSE282617` | bulk RNA-seq validation cohort | processed CSV and SOFT metadata downloaded | `GSE282617_processed_data.csv.gz` approx. 7.95 MB | 70 samples: Normal n=5, AIS n=20, MIA n=17, IAC n=23, LUSC control n=5 |
| `STDS0000125/GSM5702474` | processed spatial h5ad fallback for one IAC sample from `GSE189487` | processed h5ad and cluster marker TSV downloaded | `GSM5702474_10x_Visium_processed.h5ad` approx. 269.34 MB | Useful to develop h5ad parsing, marker scoring, and spatial adjacency code before full GEO tar download finishes |

## Expansion Datasets for the Comprehensive Study

| Accession | Role | Local status | Remote files | Parsed sample/record count |
|---|---|---|---|---:|
| `GSE308103` | Primary external snRNA-seq validation for epithelial progenitor and myeloid state definitions | SOFT metadata, `filelist.txt`, full `GSE308103_RAW.tar`, and extracted raw count files downloaded | `GSE308103_RAW.tar` 1,561,128,960 bytes | 75 |
| `GSE307534` | Primary external spatial validation for recurrence of the niche across Normal/AAH/AIS/MIA/LUAD | SOFT metadata, `filelist.txt`, full `GSE307534_RAW.tar`, extracted Space Ranger sample directories, and refined spatial signature outputs generated | `GSE307534_RAW.tar` 10,090,321,920 bytes | 56 |
| `GSE164789` | External precursor-LUAD scRNA/scTCR validation for tumor-adjacent and tumor immune states | SOFT metadata, `filelist.txt`, full `GSE164789_RAW.tar`, extracted raw 10x files, and scRNA signature outputs generated | `GSE164789_RAW.tar` 2,117,703,680 bytes | 78 |
| `GSE131907` | Broad LUAD scRNA reference for macrophage/TME specificity | SOFT metadata, feature summary, cell annotation, and raw UMI text matrix downloaded | raw UMI text matrix `408,736,818` bytes; RDS and normalized matrices not downloaded | 58 |

Current harmonized expansion metadata:

| Dataset | Parsed groups |
|---|---|
| `GSE308103` | Normal n=24, AAH n=9, AIS n=14, MIA n=4, LUAD n=24 |
| `GSE307534` | Normal n=1, AAH n=11, AIS n=14, MIA n=4, LUAD n=26 |
| `GSE164789` | Adjacent n=39, Tumor n=39 |
| `GSE131907` | Normal lung n=11, Primary tumor n=15, Metastasis/effusion n=22, Normal lymph node n=10 |

## Important Metadata Harmonization

`GSE282617` uses group abbreviations in GEO:

| GEO group | Interpreted label | Use in LUAD progression analysis |
|---|---|---|
| `ZCF` | Normal lung | Include |
| `YWA` | AIS | Include |
| `WJR` | MIA | Include |
| `JRX` | IAC | Include |
| `FLA` | Lung squamous cell carcinoma control | Exclude from the main LUAD progression trend; retain as a non-LUAD cancer control |

This mapping is implemented in `src/luad_niche/metadata.py` and exported in `results/tables/geo_sample_metadata_annotated.csv`.

## Generated Tables

- `results/tables/geo_file_inventory.csv`: first-wave supplementary file URLs.
- `results/tables/geo_file_inventory_with_sizes.csv`: same inventory with remote file sizes.
- `results/tables/geo_sample_metadata.csv`: parsed GEO SOFT sample metadata.
- `results/tables/geo_sample_metadata_annotated.csv`: harmonized sample stage labels.
- `results/tables/geo_filelists_inventory.csv`: filelist-only download plan for `GSE189357` and `GSE189487`.

## Current Download Decision

The smallest useful validation dataset, `GSE282617_processed_data.csv.gz`, has already been downloaded and used for candidate marker sanity checks.

`GSE189487_RAW.tar` is now fully downloaded and extracted:

- Tar: `data/raw/GSE189487/GSE189487_RAW.tar`
- Local size: `219719680` bytes
- Extracted directory: `data/interim/GSE189487/raw_10x/`
- Extracted contents: 24 files, covering six Visium samples with barcodes, features, MatrixMarket count matrices, and tissue positions.

`GSE189357_RAW.tar` is now fully downloaded and extracted:

- Tar: `data/raw/GSE189357/GSE189357_RAW.tar`
- Local size: `654202880` bytes
- Extracted directory: `data/interim/GSE189357/raw_10x/`
- Extracted contents: 27 files, covering nine scRNA-seq samples with barcodes, features, and MatrixMarket count matrices.

Both first-wave matched datasets (`GSE189357` scRNA-seq and `GSE189487` spatial transcriptomics) are now local. The next data-expansion target is optional rather than blocking: the larger LungPCA datasets (`GSE307534` spatial and `GSE308103` snRNA-seq) can be added after the first scRNA-to-spatial signature mapping works.

The comprehensive-study expansion is now in progress. Because NCBI single-connection downloads are slow in the current environment, large tar downloads use `scripts/download_geo_parallel.py`, which downloads independent HTTP byte ranges into `.parts/` folders before assembling the final tar.

`GSE308103_RAW.tar` is fully downloaded and extracted:

- Tar: `data/raw/GSE308103/GSE308103_RAW.tar`
- Local size: `1561128960` bytes
- Extracted directory: `data/interim/GSE308103/raw_counts/`
- Extracted contents: 75 gzipped raw-count matrix text files.

`GSE307534_RAW.tar` is fully downloaded and extracted:

- Tar: `data/raw/GSE307534/GSE307534_RAW.tar`
- Local size: `10,090,321,920` bytes.
- Extracted outer archives: `data/interim/GSE307534/raw_visium_archives/`
- Extracted Space Ranger directories: `data/interim/GSE307534/raw_visium/`
- Extracted sample directories: 56.
- Temporary `.parts/` directory has been removed after successful assembly and extraction.
- Corrected refined-signature spatial outputs are available under `results/tables/gse307534_refined_signature_*` and `results/figures/gse307534_refined_signature_adjacency_by_stage.png`.

`GSE164789_RAW.tar` is fully downloaded and extracted:

- Tar: `data/raw/GSE164789/GSE164789_RAW.tar`
- Local size: `2,117,703,680` bytes.
- Extracted directory: `data/interim/GSE164789/raw_10x/`
- Tar entries / extracted files: 202.
- Expression matrices discovered: 62 scRNA-seq samples.
- TCR contig files present: 16.
- Temporary `.parts/` directory has been removed after successful assembly and extraction.
- Signature-scoring outputs are available under `results/tables/gse164789_scrna_*` and `results/figures/gse164789_scrna_*`.

`GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` is downloaded:

- Matrix: `data/raw/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- Local size: `408,736,818` bytes.
- Annotation: `data/raw/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz`
- Feature summary: `data/raw/GSE131907/GSE131907_Lung_Cancer_Feature_Summary.xlsx`
- Matrix shape at header inspection: 208,506 cells plus gene index column.
- Temporary `.parts/` directory has been removed after successful assembly.
- Selected-gene specificity outputs are available under `results/tables/gse131907_selected_signature_*` and `results/tables/gse131907_selected_gene_*`.
- GSE131907-audited filtered discovery signatures are available at `results/tables/gse189357_refined_signature_genes_gse131907_specificity_filtered.json`.
- Specificity-filtered spatial reruns are stored in `results/tables/specificity_filtered/` and `results/figures/specificity_filtered/`.

## Current Mechanism-Ranking Outputs

Candidate mechanism configuration:

- `config/candidate_mechanisms.yaml`

Direct spatial-axis evidence:

- `results/tables/gse307534_candidate_axis_spatial_adjacency.csv`
- `results/tables/gse307534_candidate_axis_spatial_by_stage.csv`
- `results/tables/gse307534_candidate_axis_spatial_score_by_stage.csv`
- `results/tables/gse307534_candidate_axis_spatial_genes_used.json`

Mechanism and perturbation ranking:

- `results/tables/candidate_mechanism_axis_ranking.csv`
- `results/tables/candidate_perturbation_gene_ranking.csv`
- `results/tables/candidate_mechanism_ranking_weights.json`
- `results/tables/gse307534_virtual_perturbation_effects.csv`
- `results/tables/gse307534_virtual_perturbation_by_stage.csv`
- `results/tables/gse307534_virtual_perturbation_mia_luad_ranking.csv`
- `results/tables/gse307534_virtual_perturbation_genes_used.json`
- `results/tables/gse307534_continuous_perturbation_effects.csv`
- `results/tables/gse307534_continuous_perturbation_mia_luad_ranking.csv`
- `results/tables/gse307534_continuous_perturbation_genes_used.json`
- `results/tables/main_axis_evidence_matrix.csv`
- `docs/main_axis_evidence_matrix.md`
- `docs/manuscript_results_skeleton.md`

Current top-ranked axes:

| Rank | Axis | Priority score |
|---:|---|---:|
| 1 | `mif_cd74_cxcr4` | 0.789 |
| 2 | `spp1_trem2_macrophage_epithelial` | 0.748 |
| 3 | `cxcl9_cxcl10_cxcr3` | 0.614 |
| 4 | `c1q_apoe_trem2_lgals3` | 0.602 |
| 5 | `inflammatory_il1_tnf_cxcl8` | 0.506 |

## Manuscript Figure Draft Outputs

Figure 1 workflow and dataset composition:

- `results/figures/figure1_workflow_dataset_composition.svg`
- `results/figures/figure1_workflow_dataset_composition.pdf`
- `results/figures/figure1_workflow_dataset_composition.tiff`
- `results/figures/figure1_workflow_dataset_composition.png`
- Source data: `results/tables/figure1_dataset_composition_source.csv`

Figure 2 axis evidence and score-level target prioritization:

- `results/figures/figure2_axis_evidence_perturbation.svg`
- `results/figures/figure2_axis_evidence_perturbation.pdf`
- `results/figures/figure2_axis_evidence_perturbation.tiff`
- `results/figures/figure2_axis_evidence_perturbation.png`
- Source data: `results/tables/figure2_priority_source.csv`
- Source data: `results/tables/figure2_evidence_heatmap_source.csv`
- Source data: `results/tables/figure2_perturbation_source.csv`

Reproduction script:

- `scripts/plot_main_figures.py`

Figure contract:

- `docs/figure_plan.md`

Figure 3 specificity audit:

- `results/figures/figure3_specificity_audit.svg`
- `results/figures/figure3_specificity_audit.pdf`
- `results/figures/figure3_specificity_audit.tiff`
- `results/figures/figure3_specificity_audit.png`
- Source data: `results/tables/figure3_signature_celltype_heatmap_source.csv`
- Source data: `results/tables/figure3_specificity_status_source.csv`
- Source data: `results/tables/figure3_candidate_gene_specificity_source.csv`

Figure 4 spatial axis progression:

- `results/figures/figure4_spatial_axis_progression.svg`
- `results/figures/figure4_spatial_axis_progression.pdf`
- `results/figures/figure4_spatial_axis_progression.tiff`
- `results/figures/figure4_spatial_axis_progression.png`
- Source data: `results/tables/figure4_axis_stage_heatmap_source.csv`
- Source data: `results/tables/figure4_focus_axis_trend_source.csv`
- Source data: `results/tables/figure4_late_axis_summary_source.csv`

Figure 5 score-level target-prioritization priority:

- `results/figures/figure5_virtual_perturbation_priority.svg`
- `results/figures/figure5_virtual_perturbation_priority.pdf`
- `results/figures/figure5_virtual_perturbation_priority.tiff`
- `results/figures/figure5_virtual_perturbation_priority.png`
- Source data: `results/tables/figure5_dose_response_source.csv`
- Source data: `results/tables/figure5_stage_loss_source.csv`
- Source data: `results/tables/figure5_method_concordance_source.csv`
