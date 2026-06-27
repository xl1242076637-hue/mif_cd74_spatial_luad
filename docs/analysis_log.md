# Analysis Log

## 2026-05-29

Initial project setup started.

Context:

- Working directory: `D:\空间转录`
- Python: 3.13.9
- Rscript: not found
- Available Python packages checked: pandas, numpy, scipy, requests, PyYAML, pytest, sklearn, matplotlib, seaborn
- Missing packages checked: scanpy, anndata, GEOparse

Initial dataset priority:

1. `GSE189357`, `GSE189487`, and `GSE282617` for a first runnable early LUAD progression workflow.
2. `GSE307534` and `GSE308103` for later expansion to the larger Cancer Cell LungPCA cohort.

Commands run:

```powershell
python luad_epithelial_macrophage_niche\scripts\download_geo.py --dataset GSE189357 --dataset GSE189487 --dataset GSE282617 --list-only
python luad_epithelial_macrophage_niche\scripts\download_geo.py --dataset GSE282617 --download
python luad_epithelial_macrophage_niche\scripts\download_geo.py --dataset GSE189357 --dataset GSE189487 --include-regex "filelist\.txt$" --download --output luad_epithelial_macrophage_niche\results\tables\geo_filelists_download_plan.jsonl --csv-output luad_epithelial_macrophage_niche\results\tables\geo_filelists_inventory.csv
python luad_epithelial_macrophage_niche\scripts\parse_soft_metadata.py
python luad_epithelial_macrophage_niche\scripts\annotate_sample_metadata.py
```

Files downloaded:

- `data/raw/GSE282617/GSE282617_processed_data.csv.gz`，8,336,900 bytes.
- `data/raw/GSE189357/filelist.txt`，1,894 bytes.
- `data/raw/GSE189487/filelist.txt`，1,766 bytes.
- `data/raw/GSE189357/GSE189357_family.soft.gz`，3,122 bytes.
- `data/raw/GSE189487/GSE189487_family.soft.gz`，3,026 bytes.
- `data/raw/GSE282617/GSE282617_family.soft.gz`，8,134 bytes.

Key observations:

- `GSE189357` provides 9 scRNA-seq samples: AIS n=3, MIA n=3, IAC n=3.
- `GSE189487` provides 6 spatial transcriptomics samples: AIS n=2, MIA n=2, IAC n=2.
- `GSE282617` provides 70 bulk RNA-seq samples. Harmonized labels are Normal n=5, AIS n=20, MIA n=17, IAC n=23, and LUSC control n=5.
- The `FLA` group in `GSE282617` is interpreted as lung squamous cell carcinoma control, not normal lung; it should be excluded from the main LUAD progression trend.
- First runnable computational validation should start with `GSE282617`, while spatial niche discovery requires downloading `GSE189487_RAW.tar` and single-cell annotation requires `GSE189357_RAW.tar`.

Candidate marker sanity check:

```powershell
python luad_epithelial_macrophage_niche\scripts\summarize_gse282617_markers.py
```

Outputs:

- `results/tables/gse282617_candidate_marker_stage_means.csv`
- `results/tables/gse282617_candidate_marker_trends.csv`
- `results/figures/gse282617_candidate_marker_stage_means.png`

All 39 candidate genes in `config/candidate_genes.yaml` were present in `GSE282617`.

Preliminary bulk-level trend observations:

- Epithelial/tumor-associated markers such as `EPCAM`, `KRT8`, `KRT18`, `KRT19`, `MUC1`, `CLDN4`, and `TOP2A` are higher in IAC than normal lung.
- `SPP1` and `APOE` show strong IAC-vs-normal increases among macrophage/niche markers.
- `SFTPC` decreases from normal lung to IAC, consistent with loss of mature AT2-like identity during progression.
- Several macrophage/stromal markers (`C1QA`, `C1QB`, `C1QC`, `CD68`, `VIM`, `S100A4`) are lower in IAC than normal lung in this bulk dataset. This should not be interpreted as contradicting a spatial niche hypothesis; bulk RNA-seq averages cell composition across tissue and cannot resolve local macrophage-epithelial adjacency.

Next practical step:

- Download `GSE189487_RAW.tar` first, because it is the smaller spatial dataset and directly tests whether epithelial and macrophage programs colocalize in AIS/MIA/IAC tissue.
- Download `GSE189357_RAW.tar` afterward to annotate cell states for spatial deconvolution.

## 2026-05-29 Continued

Attempted full `GSE189487_RAW.tar` download from NCBI GEO:

```powershell
python luad_epithelial_macrophage_niche\scripts\download_geo.py --dataset GSE189487 --include-regex "GSE189487_RAW\.tar$" --download --output luad_epithelial_macrophage_niche\results\tables\gse189487_raw_download_plan.jsonl --csv-output luad_epithelial_macrophage_niche\results\tables\gse189487_raw_inventory.csv
```

Result:

- Timed out after 30 minutes.
- Local `data/raw/GSE189487/GSE189487_RAW.tar` was zero bytes.
- NCBI HEAD request confirmed `Content-Length: 219719680` and `Accept-Ranges: bytes`.
- `curl.exe` initially failed because of Windows Schannel certificate revocation checking; adding `--ssl-no-revoke` fixed that issue.
- A one-megabyte range request succeeded, but the transfer speed from NCBI was slow.

Fallback data source:

- Downloaded `GSM5702474_10x_Visium_processed.h5ad` from CNGB/STOMICS.
- Downloaded `GSM5702474_10x_Visium_processed.top100.cluster.markers.tsv` from CNGB/STOMICS.
- These files represent one processed IAC spatial sample from `GSE189487`; they are used for proof-of-code and should not replace the full six-sample GEO analysis.

Commands:

```powershell
curl.exe --ssl-no-revoke -L -C - --connect-timeout 60 --retry 5 --retry-delay 10 -o D:\空间转录\luad_epithelial_macrophage_niche\data\raw\STDS0000125\GSM5702474\GSM5702474_10x_Visium_processed.h5ad https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000125/stomics/GSM5702474/GSM5702474_10x_Visium_processed.h5ad
python luad_epithelial_macrophage_niche\scripts\score_stomics_h5ad.py
python luad_epithelial_macrophage_niche\scripts\test_spatial_adjacency.py --radius-multiplier 1.0 --table-output luad_epithelial_macrophage_niche\results\tables\gsm5702474_epithelial_macrophage_adjacency_radius1nn.csv --figure-output luad_epithelial_macrophage_niche\results\figures\gsm5702474_epithelial_macrophage_adjacency_radius1nn.png
```

New outputs:

- `results/tables/gsm5702474_spatial_panel_scores.csv`
- `results/tables/gsm5702474_panel_genes_used.json`
- `results/figures/gsm5702474_spatial_panel_scores.png`
- `results/tables/gsm5702474_epithelial_macrophage_adjacency.csv`
- `results/tables/gsm5702474_epithelial_macrophage_adjacency_radius1nn.csv`
- `results/figures/gsm5702474_epithelial_macrophage_adjacency.png`
- `results/figures/gsm5702474_epithelial_macrophage_adjacency_radius1nn.png`

Preliminary spatial proof-of-code observations:

- `GSM5702474` contains 4,992 spots in the processed h5ad; one spot has non-finite spatial coordinates and is dropped only for spatial-neighbor calculations.
- Candidate panel scores were computed for epithelial progenitor, proliferation/EMT, macrophage niche, and ligand-receptor axis gene panels.
- In this single IAC sample, macrophage-niche and proliferation/EMT panel scores are moderately correlated across spots, but high epithelial-progenitor and high macrophage-niche spots are not enriched for immediate adjacency under the current top-quartile/radius test.
- This should be treated as a pipeline test, not a biological conclusion, because it is one processed IAC sample and not the full AIS/MIA/IAC spatial cohort.

## 2026-05-29 Full GSE189487 Raw 10x Analysis

Completed the full `GSE189487_RAW.tar` download and verified that the local tar is readable.

Checks:

```powershell
Get-Item luad_epithelial_macrophage_niche\data\raw\GSE189487\GSE189487_RAW.tar
tar -tf luad_epithelial_macrophage_niche\data\raw\GSE189487\GSE189487_RAW.tar
```

Result:

- Local size: `219719680` bytes.
- Tar listed 24 files: six samples, each with `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`, and `tissue_positions_list.csv.gz`.
- Extracted files to `data/interim/GSE189487/raw_10x/`.

Extraction command:

```powershell
New-Item -ItemType Directory -Force -Path luad_epithelial_macrophage_niche\data\interim\GSE189487\raw_10x
tar -xf luad_epithelial_macrophage_niche\data\raw\GSE189487\GSE189487_RAW.tar -C luad_epithelial_macrophage_niche\data\interim\GSE189487\raw_10x
```

Implemented raw 10x support:

- `src/luad_niche/tenx.py`: discovers extracted 10x samples, reads barcodes/features/tissue positions, reads selected marker genes from MatrixMarket matrices, and returns all-feature total counts per spot for correct library-size normalization.
- `src/luad_niche/expression.py`: added `normalize_log1p_counts_by_totals()` so selected marker genes are normalized by full transcriptome spot totals rather than marker-subset totals.
- `scripts/score_gse189487_raw_10x.py`: batch-scores all six spatial samples and runs epithelial-high to macrophage-high spatial adjacency permutation tests.

Main run:

```powershell
python luad_epithelial_macrophage_niche\scripts\score_gse189487_raw_10x.py --radius-multiplier 1.0 --permutations 500
```

Main output tables:

- `results/tables/gse189487_raw_spatial_panel_scores.csv`
- `results/tables/gse189487_raw_panel_genes_used.json`
- `results/tables/gse189487_raw_epithelial_macrophage_adjacency.csv`
- `results/tables/gse189487_raw_epithelial_macrophage_adjacency_by_stage.csv`

Main output figures:

- `results/figures/gse189487_gsm5702473_raw_spatial_panel_scores.png`
- `results/figures/gse189487_gsm5702474_raw_spatial_panel_scores.png`
- `results/figures/gse189487_gsm5702475_raw_spatial_panel_scores.png`
- `results/figures/gse189487_gsm5702476_raw_spatial_panel_scores.png`
- `results/figures/gse189487_gsm5702477_raw_spatial_panel_scores.png`
- `results/figures/gse189487_gsm5702478_raw_spatial_panel_scores.png`
- `results/figures/gse189487_gsm5702473_raw_epithelial_macrophage_adjacency.png`
- `results/figures/gse189487_gsm5702474_raw_epithelial_macrophage_adjacency.png`
- `results/figures/gse189487_gsm5702475_raw_epithelial_macrophage_adjacency.png`
- `results/figures/gse189487_gsm5702476_raw_epithelial_macrophage_adjacency.png`
- `results/figures/gse189487_gsm5702477_raw_epithelial_macrophage_adjacency.png`
- `results/figures/gse189487_gsm5702478_raw_epithelial_macrophage_adjacency.png`
- `results/figures/gse189487_raw_epithelial_macrophage_adjacency_by_stage.png`

Observed per-sample adjacency at `radius_multiplier=1.0`:

| Sample | Stage | Spots used | Observed | Permuted null mean | Delta | p greater | p less |
|---|---:|---:|---:|---:|---:|---:|---:|
| `GSM5702473` | IAC | 4095 | 0.226 | 0.544 | -0.318 | 1.000 | 0.002 |
| `GSM5702474` | IAC | 4416 | 0.450 | 0.596 | -0.145 | 1.000 | 0.002 |
| `GSM5702475` | MIA | 1140 | 0.386 | 0.525 | -0.139 | 1.000 | 0.002 |
| `GSM5702476` | AIS | 1700 | 0.289 | 0.568 | -0.278 | 1.000 | 0.002 |
| `GSM5702477` | MIA | 3911 | 0.402 | 0.580 | -0.179 | 1.000 | 0.002 |
| `GSM5702478` | AIS | 1760 | 0.405 | 0.560 | -0.155 | 1.000 | 0.002 |

Stage-level summary:

| Stage | n samples | Observed mean | Null mean | Mean delta | p less median |
|---|---:|---:|---:|---:|---:|
| AIS | 2 | 0.347 | 0.564 | -0.217 | 0.002 |
| MIA | 2 | 0.394 | 0.553 | -0.159 | 0.002 |
| IAC | 2 | 0.338 | 0.570 | -0.232 | 0.002 |

Interpretation checkpoint:

- The current simple marker-panel definition does not support enriched immediate adjacency between top-quartile epithelial-progenitor spots and top-quartile macrophage-niche spots in raw `GSE189487`; observed fractions are below permutation expectations in all six samples.
- This does not reject the broader project idea. It means the next pure-bioinformatics step should refine the epithelial and macrophage states using `GSE189357` scRNA-seq rather than relying only on broad marker-panel means.
- A first sensitivity check with `radius_multiplier=2.0` made the null adjacency fraction nearly saturated (`~0.97-0.98`) and was therefore less interpretable for a local-niche question. The main analysis uses `radius_multiplier=1.0`, equivalent to approximately one median nearest-neighbor spot spacing.

## 2026-05-29 Full GSE189357 Raw scRNA Download and First Panel Scoring

Downloaded and verified the matched `GSE189357` scRNA-seq raw 10x dataset.

Initial single-connection attempt:

```powershell
curl.exe --ssl-no-revoke -L -C - --connect-timeout 60 --retry 20 --retry-delay 10 --retry-all-errors --speed-time 600 --speed-limit 1024 -o data\raw\GSE189357\GSE189357_RAW.tar https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar
```

Result:

- Too slow from NCBI GEO as a single connection; about 70-80 MB after roughly one hour.
- The GEO `filelist.txt` lists internal tar members, not individually downloadable files; direct attempts to fetch those filenames from the supplement directory returned small HTML error pages and were removed.

Working strategy:

- Used NCBI byte ranges to download the tar in 8 chunks.
- Seven chunks completed directly; the stalled chunk was split into four subchunks and then reassembled.
- Reassembled tar was moved to the standard path: `data/raw/GSE189357/GSE189357_RAW.tar`.

Verification:

- Local size: `654202880` bytes.
- `tar -tf` reports 27 entries.
- Extracted to `data/interim/GSE189357/raw_10x/`.

Extracted scRNA sample structure:

| Sample | Stage | Cells/barcodes | Features | Matrix shape | Nonzero entries |
|---|---:|---:|---:|---:|---:|
| `GSM5699777_TD1` | IAC | 15,216 | 33,538 | 33,538 x 15,216 | 19,291,970 |
| `GSM5699778_TD2` | IAC | 18,064 | 33,538 | 33,538 x 18,064 | 20,653,837 |
| `GSM5699779_TD3` | MIA | 12,658 | 33,538 | 33,538 x 12,658 | 20,832,966 |
| `GSM5699780_TD4` | MIA | 11,756 | 33,538 | 33,538 x 11,756 | 20,893,076 |
| `GSM5699781_TD5` | AIS | 19,079 | 33,538 | 33,538 x 19,079 | 23,384,361 |
| `GSM5699782_TD6` | MIA | 8,999 | 33,538 | 33,538 x 8,999 | 15,080,571 |
| `GSM5699783_TD7` | AIS | 3,507 | 33,538 | 33,538 x 3,507 | 3,907,898 |
| `GSM5699784_TD8` | AIS | 15,521 | 33,538 | 33,538 x 15,521 | 22,755,039 |
| `GSM5699785_TD9` | IAC | 17,573 | 33,538 | 33,538 x 17,573 | 21,486,601 |

Implemented expression-only 10x support:

- `src/luad_niche/tenx.py` now supports both spatial 10x samples with `tissue_positions_list.csv.gz` and scRNA 10x samples without spatial positions.
- Added tests for expression-only 10x discovery and barcode-only obs tables.

First scRNA panel scoring command:

```powershell
python luad_epithelial_macrophage_niche\scripts\score_gse189357_scrna_panels.py
```

Outputs:

- `results/tables/gse189357_scrna_panel_scores.csv`
- `results/tables/gse189357_scrna_panel_score_sample_summary.csv`
- `results/tables/gse189357_scrna_panel_score_stage_summary.csv`
- `results/tables/gse189357_scrna_panel_genes_used.json`
- `results/figures/gse189357_scrna_panel_score_stage_summary.png`

Stage-level first-pass panel means:

| Stage | n samples | total cells | epithelial progenitor | proliferation/EMT | macrophage niche | ligand-receptor axes |
|---|---:|---:|---:|---:|---:|---:|
| AIS | 3 | 38,107 | 0.418 | 1.038 | 0.482 | 0.353 |
| MIA | 3 | 33,413 | 0.227 | 0.929 | 0.408 | 0.391 |
| IAC | 3 | 50,853 | 0.332 | 0.874 | 0.342 | 0.311 |

Interpretation checkpoint:

- This is a panel-score sanity check, not a true scRNA annotation yet.
- Broad candidate panels are heterogeneous across samples. AIS appears high for several broad scores, partly driven by `TD7`; this reinforces that broad panels are not specific enough for the final biological claim.
- Next step should be cell-state refinement: identify epithelial-like and macrophage-like cells from marker expression, derive subtype/state signatures, and then map those refined signatures back to `GSE189487` spatial spots.

## 2026-05-29 scRNA Marker-Score State Refinement and Spatial Remapping

Implemented lightweight marker-score cell-state refinement without requiring Scanpy/Seurat.

New marker config:

- `config/cell_state_markers.yaml`

New code:

- `src/luad_niche/cell_states.py`
- `src/luad_niche/differential.py`
- `scripts/refine_gse189357_cell_states.py`
- `scripts/extract_gse189357_state_markers.py`
- `scripts/score_gse189487_refined_signatures.py`

Commands:

```powershell
python luad_epithelial_macrophage_niche\scripts\refine_gse189357_cell_states.py
python luad_epithelial_macrophage_niche\scripts\extract_gse189357_state_markers.py
python luad_epithelial_macrophage_niche\scripts\score_gse189487_refined_signatures.py
```

Quality-control iteration:

- First macrophage subtype marker extraction showed mast-cell genes (`TPSAB1`, `CPA3`, `TPSB2`, `MS4A2`, `KIT`, `HDC`) contaminating the initial `spp1_macrophage` group.
- Added `mast_cell`, `neutrophil`, and `dendritic` as broad marker-score classes and reran assignments/markers.
- After this correction, `spp1_macrophage` markers became more plausible, with `SPP1`, `CHI3L1`, `CHIT1`, `CTSK`, `PLA2G7`, `CD9`, and `ADAMDEC1` among top genes.

Main scRNA refined signatures:

| Signature | Representative top genes |
|---|---|
| epithelial progenitor-like | `KRT7`, `TACSTD2`, `CLDN3`, `KRT19`, `KRT8`, `CLDN4`, `GDF15`, `EPCAM`, `CD24`, `KRT18` |
| proliferating epithelial | `PCNA`, `TOP2A`, `PCLAF`, `MKI67`, `STMN1`, `CKS1B`, `TK1`, `TYMS`, `UBE2T`, `UBE2C` |
| SPP1 macrophage-like | `SPP1`, `CHI3L1`, `CHIT1`, `CTSK`, `PLA2G7`, `CD9`, `ADAMDEC1` |
| C1Q macrophage-like | `RBP4`, `NUPR1`, `FABP4`, `MARCO`, `C1QB`, `CYP27A1`, `LPL`, `SERPING1`, `APOC1` |
| inflammatory macrophage-like | `IL1B`, `CCL3L1`, `CCL3`, `G0S2`, `GPR183`, `NLRP3`, `PTGS2`, `CXCL8` |

Low-confidence signature note:

- `resident_macrophage` remains low confidence in this lightweight workflow; its extracted genes include neutrophil/monocyte-like markers (`FCN1`, `VCAN`, `AZU1`, `ELANE`). Do not use it as a primary claim without a stronger clustering/decontamination step.

Key outputs:

- `results/tables/gse189357_scrna_cell_state_assignments.csv`
- `results/tables/gse189357_scrna_broad_class_stage_summary.csv`
- `results/tables/gse189357_scrna_epithelial_state_stage_summary.csv`
- `results/tables/gse189357_scrna_macrophage_state_stage_summary.csv`
- `results/tables/gse189357_refined_state_markers.csv`
- `results/tables/gse189357_refined_state_top_markers.csv`
- `results/tables/gse189357_refined_state_signature_genes.json`
- `results/tables/gse189487_refined_signature_spatial_scores.csv`
- `results/tables/gse189487_refined_signature_adjacency.csv`
- `results/tables/gse189487_refined_signature_adjacency_by_stage.csv`

Refined spatial remapping checkpoint:

| Target macrophage signature | AIS mean delta | MIA mean delta | IAC mean delta | Most notable result |
|---|---:|---:|---:|---|
| `spp1_macrophage` | -0.144 | +0.058 | -0.064 | MIA positive adjacency, median `p_greater` about 0.006 |
| `c1q_macrophage` | -0.180 | -0.130 | -0.207 | Negative adjacency across stages |
| `inflammatory_macrophage` | -0.074 | +0.041 | -0.025 | MIA weak positive trend |
| `resident_macrophage` | -0.088 | -0.022 | -0.065 | Low-confidence signature, not a primary result |

Working interpretation:

- The broad epithelial-macrophage panel was negative, but refined scRNA-derived signatures reveal a more specific possible niche: epithelial progenitor-like spots are locally enriched near `SPP1 macrophage-like` signature spots in MIA samples.
- This is a better candidate story than a generic epithelial-macrophage niche.
- Next computational step should focus on robustness: top-N signature sensitivity, radius/quantile sensitivity, and ligand-receptor axes around the MIA epithelial progenitor-SPP1 macrophage niche.

## 2026-05-29 SPP1 Niche Robustness Grid

Implemented a parameter-sensitivity analysis for the prioritized MIA epithelial progenitor-like to `SPP1 macrophage-like` spatial niche.

New script:

- `scripts/robustness_gse189487_spp1_niche.py`

Command:

```powershell
python scripts\robustness_gse189487_spp1_niche.py
```

Grid:

- epithelial and macrophage signatures rebuilt from `GSE189357` ranked markers at top-N = 10, 20, 30 genes.
- high-score spot cutoffs = top 30%, 25%, and 20% (`quantile=0.70,0.75,0.80`).
- spatial radius = 0.75, 1.00, and 1.25 times the median nearest-neighbor distance.
- 200 spatial-label permutations per sample/parameter set.

Key outputs:

- `results/tables/gse189487_spp1_niche_robustness.csv`
- `results/tables/gse189487_spp1_niche_robustness_by_stage.csv`
- `results/tables/gse189487_spp1_niche_robustness_mia_summary.csv`
- `results/figures/gse189487_spp1_niche_robustness_mia.png`

Stage-level result:

| Stage | Positive parameter sets | Positive and `p_greater < 0.05` sets | Mean delta across grid |
|---|---:|---:|---:|
| AIS | 0/27 | 0/27 | -0.158 |
| MIA | 13/27 | 6/27 | +0.002 |
| IAC | 0/27 | 0/27 | -0.098 |

MIA signal details:

- The statistically strongest MIA results are all top30 signature settings.
- Significant MIA parameter sets are `top_n=30` with `radius_multiplier=0.75` or `1.0` across all tested score quantiles.
- The best mean enrichment delta is `+0.096` at `top_n=30`, `quantile=0.75`, `radius_multiplier=0.75`, with median `p_greater=0.00498`.
- In the significant top30 settings, both MIA samples (`GSM5702475` and `GSM5702477`) have positive enrichment, so the result is not driven by only one MIA sample.

Interpretation checkpoint:

- The candidate niche is stage-specific in the current six-sample spatial cohort: AIS and IAC do not show positive enrichment across the tested parameter grid.
- The signal is real enough to keep as the main computational story, but it is parameter-sensitive. It depends on a broader top30 `SPP1 macrophage-like` signature and a local one-hop or sub-one-hop spatial radius.
- The cautious claim should be: MIA samples show a local spatial enrichment between epithelial progenitor-like and `SPP1 macrophage-like` signatures under refined scRNA-derived signatures, especially with top30 marker panels. Do not claim a universally robust epithelial-macrophage niche.

Verification:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
pytest tests -q
```

Result: `39 passed`.

## 2026-05-29 Expansion Cohort Setup for the Comprehensive Study

Expanded the project from a minimal discovery workflow to a multi-cohort design.

New or updated code:

- `config/datasets.yaml`: promoted `GSE308103` and `GSE307534` to primary expansion cohorts and added `GSE164789` and `GSE131907`.
- `src/luad_niche/geo.py`: added GEO family SOFT URL construction.
- `scripts/download_geo_soft.py`: downloads GEO family SOFT metadata for manifest datasets.
- `src/luad_niche/download.py`: added local/remote download status columns, resumable range downloads, and parallel HTTP byte-range downloading.
- `scripts/download_geo_parallel.py`: downloads large GEO supplementary files as independent range chunks and assembles the final tar.
- `src/luad_niche/metadata.py`: added stage/condition harmonization for LungPCA (`GSE307534`, `GSE308103`), `GSE164789`, and `GSE131907`.

Commands:

```powershell
python scripts\download_geo.py --dataset GSE308103 --dataset GSE307534 --dataset GSE164789 --dataset GSE131907 --fetch-sizes --output results\tables\geo_expansion_download_plan.jsonl --csv-output results\tables\geo_expansion_file_inventory.csv
python scripts\download_geo_soft.py --dataset GSE308103 --dataset GSE307534 --dataset GSE164789 --dataset GSE131907 --download --output results\tables\geo_expansion_soft_download_plan.jsonl --csv-output results\tables\geo_expansion_soft_file_inventory.csv
python scripts\download_geo.py --dataset GSE308103 --dataset GSE307534 --dataset GSE164789 --dataset GSE131907 --include-regex 'filelist\.txt$|cell_annotation|Feature_Summary' --fetch-sizes --download --output results\tables\geo_expansion_small_download_plan.jsonl --csv-output results\tables\geo_expansion_small_file_inventory.csv
python scripts\parse_soft_metadata.py
python scripts\annotate_sample_metadata.py
```

Expansion file inventory:

| Dataset | Main remote file | Remote size |
|---|---|---:|
| `GSE308103` | `GSE308103_RAW.tar` | 1,561,128,960 bytes |
| `GSE307534` | `GSE307534_RAW.tar` | 10,090,321,920 bytes |
| `GSE164789` | `GSE164789_RAW.tar` | 2,117,703,680 bytes |
| `GSE131907` | raw/normalized matrix files plus annotation | up to 2.9 GB per matrix file |

Parsed and harmonized metadata now cover 352 GEO sample records:

| Dataset | Harmonized groups |
|---|---|
| `GSE308103` | Normal n=24, AAH n=9, AIS n=14, MIA n=4, LUAD n=24 |
| `GSE307534` | Normal n=1, AAH n=11, AIS n=14, MIA n=4, LUAD n=26 |
| `GSE164789` | Adjacent n=39, Tumor n=39 |
| `GSE131907` | Normal lung n=11, Primary tumor n=15, Metastasis/effusion n=22, Normal lymph node n=10 |

Download issue and resolution:

- Direct NCBI single-connection download of `GSE308103_RAW.tar` was too slow and unstable in the current environment.
- The first attempt stopped after `37,748,736` bytes with an `IncompleteRead`.
- A resumed single-connection attempt reached `89,128,960` bytes but remained too slow for practical use.
- Implemented `scripts/download_geo_parallel.py` to download many HTTP byte ranges in parallel.
- Started background parallel download for `GSE308103_RAW.tar` with 16 workers and 2 MB chunks.

Current background command:

```powershell
python scripts\download_geo_parallel.py --dataset GSE308103 --include-regex 'GSE308103_RAW\.tar$' --workers 16 --range-chunk-mb 2 --max-attempts 100 --timeout 360 --output results\tables\gse308103_parallel_download_plan.jsonl --csv-output results\tables\gse308103_parallel_file_inventory.csv
```

Logs:

- `results/logs/gse308103_parallel_download.out.log`
- `results/logs/gse308103_parallel_download.err.log`
- `results/logs/gse308103_parallel_download.pid`

Verification:

```powershell
pytest tests -q
```

Result: `46 passed`.

## 2026-05-30 GSE308103 Download Recovery and Extraction

User noticed that the `GSE308103_RAW.tar` download appeared stopped. Checked the background parallel download process and confirmed it had stalled at 90.06%:

- Old PID: `18776`
- Complete parts: 671
- Part bytes: `1,405,939,712`
- No byte growth over a 60-second observation window.

Action:

- Stopped the stalled process.
- Restarted `scripts/download_geo_parallel.py` with 24 workers, 2 MB chunks, and 200 max attempts.
- The resumed download advanced from 90.06% to 98.93% within one minute and then completed.

Verification and extraction:

```powershell
tar -tf data\raw\GSE308103\GSE308103_RAW.tar
tar -xf data\raw\GSE308103\GSE308103_RAW.tar -C data\interim\GSE308103\raw_counts
```

Result:

- Final tar size: `1,561,128,960` bytes.
- Tar entries: 75.
- Extracted files: 75.
- Extracted directory: `data/interim/GSE308103/raw_counts/`.
- Cleaned temporary `data/raw/GSE308103/GSE308103_RAW.tar.parts/` after successful verification/extraction.

Started the next expansion download:

```powershell
python scripts\download_geo_parallel.py --dataset GSE307534 --include-regex 'GSE307534_RAW\.tar$' --workers 32 --range-chunk-mb 4 --max-attempts 200 --timeout 300 --output results\tables\gse307534_parallel_download_plan.jsonl --csv-output results\tables\gse307534_parallel_file_inventory.csv
```

Initial status after launch:

- PID file: `results/logs/gse307534_parallel_download.pid`
- PID: `7632`
- Complete parts after about 1 minute: 56
- Downloaded bytes: `234,881,024`
- Progress: 2.33%

## 2026-05-30 GSE307534 Spatial Validation Rerun

Completed the LungPCA Visium expansion cohort and reran refined signature spatial mapping.

Data status:

- `GSE307534_RAW.tar` fully downloaded: `10,090,321,920` bytes.
- Tar verified and extracted into 56 inner sample archives.
- Space Ranger directories extracted to `data/interim/GSE307534/raw_visium/`.
- Parsed cohort: Normal n=1, AAH n=11, AIS n=14, MIA n=4, LUAD n=26.

Code correction before rerun:

- Fixed duplicate-lesion sample mapping from GEO filelist names to extracted Space Ranger directories.
- Examples: `P4_AAH`/`P4_AAH-1` now map to `P4_AAH1`/`P4_AAH2`; `P21_AIS`/`P21_AIS-1` now map to `P21_AIS1`/`P21_AIS2`.
- Added explicit `status` handling for adjacency tests where source-high or target-high spots are absent.
- Invalid adjacency rows now report `status=insufficient_high_spots` and are excluded from stage-level mean enrichment.

Verification:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
pytest tests -q
python scripts\score_gse307534_refined_signatures.py --permutations 100
```

Results:

- Tests: `54 passed`.
- Rerun completed across all 56 spatial samples.
- Main outputs:
  - `results/tables/gse307534_refined_signature_spatial_scores.csv`
  - `results/tables/gse307534_refined_signature_adjacency.csv`
  - `results/tables/gse307534_refined_signature_adjacency_by_stage.csv`
  - `results/tables/gse307534_refined_signature_score_by_stage.csv`
  - `results/figures/gse307534_refined_signature_adjacency_by_stage.png`

Invalid high-spot samples:

- `GSM9226189` / `P10_MIA`
- `GSM9226172` / `P3_AIS`
- `GSM9226176` / `P4_AAH2`
- `GSM9226181` / `P6_LUAD`

These samples had zero high-scoring epithelial progenitor-like and macrophage-signature spots under the current top-quantile rule, so they are not used for adjacency effect-size means.

Corrected stage-level epithelial progenitor-like to macrophage-signature adjacency deltas:

| Target macrophage signature | Normal | AAH | AIS | MIA | LUAD |
|---|---:|---:|---:|---:|---:|
| `spp1_macrophage` | 0.205 | 0.155 | 0.230 | 0.316 | 0.292 |
| `c1q_macrophage` | 0.070 | 0.080 | 0.143 | 0.265 | 0.203 |
| `inflammatory_macrophage` | 0.200 | 0.166 | 0.206 | 0.239 | 0.204 |
| `resident_macrophage` | 0.164 | 0.124 | 0.183 | 0.230 | 0.201 |

Interpretation:

- The external GSE307534 Visium cohort supports recurrence of local epithelial progenitor-like to macrophage-like spatial coupling.
- The strongest corrected signal is still around MIA/LUAD, especially for `spp1_macrophage`, but it is not MIA-exclusive.
- AIS already shows a positive signal, and Normal/AAH are not completely negative, so the next analysis must test specificity against generic spatial gradients.
- Immediate next controls should include broad epithelial/macrophage scores, random gene-set controls, library-size or tissue-density covariates, and patient/lesion-paired summaries.

## 2026-05-30 GSE164789 Download and scRNA Validation

Downloaded and extracted the precursor-LUAD scRNA/scTCR expansion cohort.

Data status:

- `GSE164789_RAW.tar` fully downloaded: `2,117,703,680` bytes.
- Tar entries: 202.
- Extracted directory: `data/interim/GSE164789/raw_10x/`.
- Expression matrices discovered: 62 scRNA-seq 10x MatrixMarket samples.
- Remaining GEO records are TCR contig annotation files, not expression matrices.
- Temporary `.parts/` directory removed after successful tar assembly and extraction.

Code updates:

- Extended `src/luad_niche/tenx.py` to support dot-named 10x files:
  - `sample.barcodes.tsv.gz`
  - `sample.genes.tsv.gz`
  - `sample.matrix.mtx.gz`
- Added tests for dot-named two-column `genes.tsv.gz` files.
- Added `scripts/score_gse164789_scrna_states.py` for Adjacent/Tumor scRNA signature scoring.

Verification:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\score_gse164789_scrna_states.py --limit-samples 2
python scripts\score_gse164789_scrna_states.py
pytest tests -q
```

Results:

- Tests: `56 passed`.
- Full scRNA scoring completed for 62 expression samples:
  - Adjacent: 31 samples, 180,251 cells.
  - Tumor: 31 samples, 176,340 cells.
- Main outputs:
  - `results/tables/gse164789_scrna_cell_state_assignments.csv`
  - `results/tables/gse164789_scrna_broad_class_stage_summary.csv`
  - `results/tables/gse164789_scrna_epithelial_state_stage_summary.csv`
  - `results/tables/gse164789_scrna_macrophage_state_stage_summary.csv`
  - `results/tables/gse164789_scrna_refined_signature_stage_summary.csv`
  - `results/figures/gse164789_scrna_*`

Key Tumor vs Adjacent observations:

| Readout | Adjacent | Tumor |
|---|---:|---:|
| Epithelial broad-class mean fraction | 0.108 | 0.210 |
| Macrophage broad-class mean fraction | 0.239 | 0.236 |
| Epithelial `progenitor_like` mean fraction among epithelial cells | 0.155 | 0.126 |
| Epithelial `proliferating_epithelial` mean fraction among epithelial cells | 0.189 | 0.222 |
| Macrophage `spp1_macrophage` mean fraction among macrophages | 0.035 | 0.059 |
| Macrophage `inflammatory_macrophage` mean fraction among macrophages | 0.224 | 0.248 |
| Refined SPP1 macrophage score in macrophages | 0.120 | 0.143 |
| Refined inflammatory macrophage score in macrophages | 0.479 | 0.580 |

Interpretation:

- `GSE164789` supports tumor-side strengthening of SPP1/inflammatory macrophage programs.
- It does not simply reproduce an epithelial progenitor fraction increase; tumor samples have more epithelial cells overall and a higher proliferating epithelial fraction.
- This cohort is best used as immune-state and tumor/adjacent validation, while the stage-specific early-progression claim should continue to rely on `GSE189357/GSE189487`, `GSE308103`, `GSE307534`, and `GSE282617`.

## 2026-05-30 GSE131907 Matrix Download

Downloaded the broad LUAD single-cell reference matrix for later cell-type specificity analysis.

Data status:

- Existing local files:
  - `data/raw/GSE131907/GSE131907_family.soft.gz`
  - `data/raw/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz`
  - `data/raw/GSE131907/GSE131907_Lung_Cancer_Feature_Summary.xlsx`
- Newly downloaded:
  - `data/raw/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
  - Local size: `408,736,818` bytes.
- Temporary `.parts/` directory removed after successful assembly.

File inspection:

- Matrix is gene-by-cell, tab-delimited, gzip-compressed text.
- Header has 208,506 cell columns plus the `Index` gene column.
- Annotation columns include `Barcode`, `Sample`, `Sample_Origin`, `Cell_type`, `Cell_type.refined`, and `Cell_subtype`.

Decision for next analysis:

- Use `GSE131907` for specificity, not stage discovery.
- Next code step should implement a streaming selected-gene reader for the wide gene-by-cell matrix, then aggregate candidate signature scores by `Sample_Origin`, `Cell_type.refined`, and `Cell_subtype`.

## 2026-05-30 GSE131907 Specificity Scoring and Signature Audit

Implemented wide-matrix selected-gene scoring for the GSE131907 lung cancer single-cell atlas.

Code updates:

- Added `src/luad_niche/wide_matrix.py` for streaming selected genes from gene-by-cell gzip text matrices.
- Added `scripts/score_gse131907_selected_signatures.py`.
- Added `src/luad_niche/signature_specificity.py`.
- Added `scripts/audit_gse189357_signatures_with_gse131907.py`.

Verification:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
pytest tests -q
python scripts\score_gse131907_selected_signatures.py
python scripts\audit_gse189357_signatures_with_gse131907.py
```

Results:

- GSE131907 scoring completed for 208,506 cells and 198 selected genes.
- Outputs:
  - `results/tables/gse131907_selected_signature_cell_scores.csv`
  - `results/tables/gse131907_selected_signature_celltype_summary.csv`
  - `results/tables/gse131907_selected_signature_origin_summary.csv`
  - `results/tables/gse131907_selected_signature_subtype_summary.csv`
  - `results/tables/gse131907_selected_gene_celltype_summary.csv`
  - `results/tables/gse131907_selected_gene_top_celltype.csv`
  - `results/figures/gse131907_selected_signature_celltype_heatmap.png`

Key specificity checks:

- Broad `spp1_macrophage` score is highest in `Myeloid cells` (mean 1.366), supporting the broad marker panel as a myeloid-relevant readout.
- `refined_inflammatory_macrophage` is also highest in `Myeloid cells` (mean 0.496), with top subtypes including `CD163+CD14+ DCs`, monocytes, and mo-Mac.
- `refined_epithelial_progenitor_like` is highest in `Epithelial cells` (mean 1.075), supporting epithelial specificity.
- However, `refined_spp1_macrophage` is highest in `Epithelial cells` (mean 0.276) rather than `Myeloid cells` (mean 0.187), indicating epithelial/tumor-program contamination.

SPP1 macrophage refined-signature audit:

- 14/30 genes have GSE131907 top expression in `Myeloid cells`.
- 15/30 genes are off-target.
- 1/30 genes is missing from the GSE131907 selected matrix.
- Off-target high contributors include epithelial-enriched `WFDC2`, `OCIAD2`, `CST6`, `CD9`, `RAMP1`, and fibroblast-enriched `CTSK`, `RARRES1`, `SDC2`, `TNS1`.

Generated specificity-filtered signatures:

- `results/tables/gse189357_refined_signature_gse131907_specificity_audit.csv`
- `results/tables/gse189357_refined_signature_gse131907_specificity_summary.csv`
- `results/tables/gse189357_refined_signature_genes_gse131907_specificity_filtered.json`

Filtered `spp1_macrophage_vs_other_macrophage` genes:

`SPP1`, `CHIT1`, `OTOA`, `HAMP`, `HBB`, `PLA2G7`, `ADAMDEC1`, `SDS`, `FABP3`, `TREM2`, `RAB42`, `ATP6V0D2`, `CAMK1`, `EEPD1`.

## 2026-05-30 Specificity-Filtered Spatial Rerun

Reran spatial adjacency using the GSE131907 specificity-filtered signatures.

Commands:

```powershell
python scripts\score_gse189487_refined_signatures.py --signatures results\tables\gse189357_refined_signature_genes_gse131907_specificity_filtered.json --table-dir results\tables\specificity_filtered --figure-dir results\figures\specificity_filtered
python scripts\score_gse307534_refined_signatures.py --signatures results\tables\gse189357_refined_signature_genes_gse131907_specificity_filtered.json --table-dir results\tables\specificity_filtered --figure-dir results\figures\specificity_filtered --permutations 100
```

Key result:

- In the small discovery spatial cohort `GSE189487`, specificity-filtered macrophage panels no longer show positive epithelial-progenitor adjacency. All filtered target panels are negative in AIS/MIA/IAC.
- For filtered `spp1_macrophage` in `GSE189487`:
  - AIS delta -0.280
  - MIA delta -0.149
  - IAC delta -0.152
- In the larger external spatial cohort `GSE307534`, specificity-filtered `spp1_macrophage` adjacency remains positive but much weaker:
  - Normal delta 0.009
  - AAH delta 0.014
  - AIS delta 0.030
  - MIA delta 0.085
  - LUAD delta 0.088

Interpretation:

- The original refined SPP1 macrophage spatial signal was partly driven by epithelial/stromal genes in the refined marker list.
- The stronger, safer conclusion is not "a robust MIA-specific SPP1 macrophage niche" from the first discovery pair.
- The revised direction should be: epithelial plasticity and myeloid inflammatory/SPP1 programs are coupled during progression, but macrophage-specific spatial claims require specificity-filtered signatures and larger-cohort support.

## 2026-05-30 Candidate Mechanism Ranking and Direct Spatial-Axis Evidence

Updated the mechanism-prioritization layer after the GSE131907 specificity audit.

Code updates:

- `scripts/score_gse131907_selected_signatures.py` now includes genes from `config/candidate_mechanisms.yaml` when building the selected-gene list.
- `scripts/score_gse308103_snrna_states.py` and `scripts/score_gse164789_scrna_states.py` now retain `c1q_macrophage_score` and `refined_c1q_macrophage_score` in the compact cell and stage summaries.
- Added `src/luad_niche/spatial_axis.py`.
- Added `scripts/analyze_gse307534_candidate_axis_spatial.py` to test whether each candidate mechanism-axis source/target panel is spatially adjacent to epithelial progenitor-like spots in GSE307534.
- `scripts/rank_candidate_mechanisms.py` now includes optional direct candidate-axis spatial evidence from `results/tables/gse307534_candidate_axis_spatial_by_stage.csv`.

Commands:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\score_gse131907_selected_signatures.py
python scripts\audit_gse189357_signatures_with_gse131907.py
python scripts\score_gse308103_snrna_states.py
python scripts\score_gse164789_scrna_states.py
python scripts\analyze_gse307534_candidate_axis_spatial.py
python scripts\rank_candidate_mechanisms.py
```

GSE131907 rerun:

- Scored 208,506 cells with 211 selected genes.
- Candidate mechanism genes such as `MIF`, `CD74`, `CXCR4`, `CD44`, `ITGAV`, `ITGB1`, `ITGA5`, `IL1R1`, `IL1RAP`, `LRP1`, `CXCL9`, `CXCL10`, `CXCR3`, and `MARCO` are now represented in the gene top-celltype table.
- Remaining missing selected genes are `SELENOP`, `PCLAF`, `CENPX`, `JPT1`, `GSDME`, and `CCL3L1`; these do not block the current mechanism axes.

Direct GSE307534 candidate-axis spatial evidence:

Mean MIA/LUAD enrichment delta, averaged across valid samples:

| Axis/evidence | MIA/LUAD mean delta | Valid tests |
|---|---:|---:|
| `mif_cd74_cxcr4`, source near epithelial progenitor | 0.204 | 30 |
| `spp1_trem2_macrophage_epithelial`, target near epithelial progenitor | 0.195 | 30 |
| `mif_cd74_cxcr4`, target near epithelial progenitor | 0.168 | 30 |
| `c1q_apoe_trem2_lgals3`, target near epithelial progenitor | 0.161 | 30 |
| `inflammatory_il1_tnf_cxcl8`, target near epithelial progenitor | 0.143 | 30 |
| `spp1_trem2_macrophage_epithelial`, source near epithelial progenitor | 0.113 | 30 |
| `c1q_apoe_trem2_lgals3`, source near epithelial progenitor | 0.071 | 30 |
| `cxcl9_cxcl10_cxcr3`, source near epithelial progenitor | 0.055 | 30 |
| `inflammatory_il1_tnf_cxcl8`, source near epithelial progenitor | 0.036 | 30 |
| `cxcl9_cxcl10_cxcr3`, target near epithelial progenitor | 0.025 | 30 |

Updated multi-cohort mechanism ranking:

| Rank | Axis | Priority score | Notes |
|---:|---|---:|---|
| 1 | `mif_cd74_cxcr4` | 0.789 | Strong epithelial-source/myeloid-target specificity and strongest direct GSE307534 source-near-epithelial spatial evidence. |
| 2 | `spp1_trem2_macrophage_epithelial` | 0.748 | Best macrophage-state perturbation anchor; strong bulk/snRNA support and positive filtered spatial recurrence, but receptor-side epithelial specificity is modest. |
| 3 | `cxcl9_cxcl10_cxcr3` | 0.614 | Cell-type specificity is good, but this is more immune-recruitment than epithelial-macrophage niche biology. |
| 4 | `c1q_apoe_trem2_lgals3` | 0.602 | Myeloid-specific and spatially recurrent, but snRNA/scRNA c1q macrophage scores do not show a strong late/tumor increase. |
| 5 | `inflammatory_il1_tnf_cxcl8` | 0.506 | Useful inflammatory readout; weaker bulk and receptor-side specificity under the current scoring. |

Top perturbation candidates from the updated ranking:

- Mechanistic epithelial-to-myeloid axis: `MIF`, `CD74`, `CXCR4`.
- Macrophage-state axis: `SPP1`, `TREM2`, `PLA2G7`.
- Secondary immune axes: `CXCL9`, `CXCL10`, `CXCR3`, `APOE`, `LGALS3`, `IL1B`, `TNF`, `NLRP3`, `CXCL8`.

Interpretation:

- The computational story should now separate two linked layers:
  - `MIF-CD74/CXCR4` as the leading ligand-receptor mechanism hypothesis for epithelial progenitor-like/stress epithelial spots communicating with myeloid compartments.
  - `SPP1/TREM2/PLA2G7` as the leading macrophage-state and later perturbation-readout axis.
- This is a better framing than forcing the whole story to be a pure `SPP1 macrophage niche`, because the earlier unfiltered refined SPP1 panel had specificity issues.

## 2026-05-30 GSE307534 Score-Level Target Prioritization

Implemented a first score-level in-silico target-prioritization layer for the GSE307534 spatial cohort.

Code updates:

- Added `src/luad_niche/perturbation.py`.
- Added `scripts/virtual_perturb_gse307534_axes.py`.
- Added `tests/test_perturbation.py`.

Command:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\virtual_perturb_gse307534_axes.py
```

Outputs:

- `results/tables/gse307534_virtual_perturbation_effects.csv`
- `results/tables/gse307534_virtual_perturbation_by_stage.csv`
- `results/tables/gse307534_virtual_perturbation_mia_luad_ranking.csv`
- `results/tables/gse307534_virtual_perturbation_genes_used.json`

Scope:

- 56 GSE307534 Visium samples.
- 1,232 sample-level perturbation effect rows.
- Focus genes: `MIF`, `CD74`, `CXCR4`, `SPP1`, `TREM2`, `PLA2G7`.
- Perturbation factors: `0` score dropout and `0.5` score down-weighting.

Key MIA/LUAD ranking results:

| Perturbed genes | Axis | Evidence | Factor | Panel relative delta | Observed adjacency delta | Interpretation |
|---|---|---|---:|---:|---:|---|
| `MIF,CD74,CXCR4` | `mif_cd74_cxcr4` | source near epithelial progenitor | 0 | -1.000 | -0.565 | Full source-side score dropout collapses the MIF source score and the source-near-epithelial spatial signal. |
| `MIF` | `mif_cd74_cxcr4` | source near epithelial progenitor | 0 | -1.000 | -0.565 | Same as above because the source panel is currently single-gene `MIF`. |
| `CD74` | `mif_cd74_cxcr4` | target near epithelial progenitor | 0 | -0.714 | -0.004 | `CD74` is the dominant contributor to the receptor-side panel score; spatial top-rank adjacency changes little. |
| `CXCR4` | `mif_cd74_cxcr4` | target near epithelial progenitor | 0 | -0.059 | +0.005 | `CXCR4` contributes much less than `CD74` in the current GSE307534 panel. |
| `SPP1,TREM2,PLA2G7` | `spp1_trem2_macrophage_epithelial` | source near epithelial progenitor | 0 | -0.127 | -0.008 | SPP1/TREM2/PLA2G7 perturbation weakens the macrophage-state panel, but less dramatically than MIF-axis source perturbation. |
| `SPP1` | `spp1_trem2_macrophage_epithelial` | source near epithelial progenitor | 0 | -0.062 | -0.010 | `SPP1` is the strongest single-gene contributor among the current SPP1/TREM2/PLA2G7 perturbation set. |
| `TREM2` | `spp1_trem2_macrophage_epithelial` | source near epithelial progenitor | 0 | -0.039 | +0.001 | Smaller source-panel effect than `SPP1`. |
| `PLA2G7` | `spp1_trem2_macrophage_epithelial` | source near epithelial progenitor | 0 | -0.026 | +0.001 | Smallest effect among the three in the current spatial scoring setup. |

Important caveat:

- This is not causal validation. It is a score-level perturbation of public spatial transcriptomic expression matrices.
- Multiplicative knockdown, especially `0.5`, often changes panel magnitude but not top-quantile spot ranking, so spatial adjacency may remain unchanged.
- For single-gene panels, score dropout can collapse the high-spot set. This is useful for sensitivity analysis but should be described conservatively as score-level in-silico target prioritization.

Interpretation:

- The perturbation layer supports `MIF` as the most sensitive source-side computational target in the `MIF-CD74/CXCR4` axis.
- On the receptor side, `CD74` is much stronger than `CXCR4` in GSE307534.
- `SPP1` remains the main macrophage-state readout gene among `SPP1/TREM2/PLA2G7`.

## 2026-05-30 Continuous Spatial-Coupling Perturbation

Added a continuous spatial-coupling perturbation analysis to complement the top-quantile adjacency perturbation above.

Reason:

- The first score-level target-prioritization analysis uses high-spot/top-quantile adjacency.
- This is useful for discrete spatial-niche enrichment, but `0.5` knockdown can leave the top-quantile spot ranking unchanged.
- A continuous coupling score is better for dose-response-like perturbation because it uses score magnitude directly.

Code updates:

- Added `src/luad_niche/spatial_coupling.py`.
- Added `scripts/continuous_perturb_gse307534_axes.py`.
- Added `tests/test_spatial_coupling.py`.

Command:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\continuous_perturb_gse307534_axes.py
```

Outputs:

- `results/tables/gse307534_continuous_perturbation_effects.csv`
- `results/tables/gse307534_continuous_perturbation_mia_luad_ranking.csv`
- `results/tables/gse307534_continuous_perturbation_genes_used.json`

Method:

- For each spot, compute the mean candidate target-panel score in local neighboring spots within one median-nearest-neighbor radius.
- Compute a continuous coupling score as `epithelial_progenitor_like_score * local_neighbor_target_score`.
- Rerun this score after in-silico gene/axis perturbation.
- Summarize MIA/LUAD relative coupling changes.

Key MIA/LUAD continuous perturbation results:

| Perturbed genes | Axis | Evidence | Factor | Coupling relative delta | Interpretation |
|---|---|---|---:|---:|---|
| `MIF,CD74,CXCR4` | `mif_cd74_cxcr4` | source near epithelial progenitor | 0 | -1.000 | Full MIF-axis source dropout collapses continuous coupling. |
| `MIF` | `mif_cd74_cxcr4` | source near epithelial progenitor | 0 | -1.000 | Same as full source dropout because the current source panel is single-gene `MIF`. |
| `MIF,CD74,CXCR4` | `mif_cd74_cxcr4` | target near epithelial progenitor | 0 | -0.800 | Full axis dropout strongly reduces receptor-side coupling. |
| `CD74` | `mif_cd74_cxcr4` | target near epithelial progenitor | 0 | -0.740 | `CD74` explains most receptor-side coupling. |
| `MIF` | `mif_cd74_cxcr4` | source near epithelial progenitor | 0.5 | -0.500 | Continuous score captures the expected 50% dose response. |
| `CD74` | `mif_cd74_cxcr4` | target near epithelial progenitor | 0.5 | -0.370 | Receptor-side knockdown response remains strong. |
| `SPP1,TREM2,PLA2G7` | `spp1_trem2_macrophage_epithelial` | source near epithelial progenitor | 0 | -0.173 | Macrophage-state axis coupling decreases, but less than MIF-CD74. |
| `SPP1` | `spp1_trem2_macrophage_epithelial` | source near epithelial progenitor | 0 | -0.121 | `SPP1` is the strongest single-gene contributor within this macrophage-state panel. |
| `CXCR4` | `mif_cd74_cxcr4` | target near epithelial progenitor | 0 | -0.060 | `CXCR4` is much weaker than `CD74` in the current spatial data. |

Interpretation:

- The continuous analysis confirms the top-quantile perturbation conclusion.
- `MIF` and `CD74` remain the strongest computational perturbation candidates.
- `SPP1` remains the leading macrophage-state readout candidate, but its perturbation effect is smaller than the MIF-CD74 axis.
- This gives us a stronger and cleaner virtual-perturbation section because it includes both thresholded spatial niche loss and continuous dose-response-like coupling loss.

## 2026-05-30 Main Evidence Matrix and Manuscript Skeleton

Built a compact manuscript-facing evidence matrix to unify the current results across spatial, specificity, bulk, snRNA/scRNA, and score-level target-prioritization evidence.

Code updates:

- Added `src/luad_niche/evidence_matrix.py`.
- Added `scripts/build_main_evidence_matrix.py`.
- Added `tests/test_evidence_matrix.py`.

Command:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\build_main_evidence_matrix.py
```

Outputs:

- `results/tables/main_axis_evidence_matrix.csv`
- `docs/main_axis_evidence_matrix.md`
- `docs/manuscript_results_skeleton.md`

Current manuscript-facing axis interpretation:

| Rank | Axis | Grade | Interpretation |
|---:|---|---|---|
| 1 | `mif_cd74_cxcr4` | lead | Lead epithelial-to-myeloid communication hypothesis; prioritize `MIF` and `CD74`. |
| 2 | `spp1_trem2_macrophage_epithelial` | strong | Macrophage-state readout axis; prioritize `SPP1`, with `TREM2/PLA2G7` secondary. |
| 3 | `cxcl9_cxcl10_cxcr3` | supporting | Immune-recruitment supporting axis, less central to epithelial-myeloid mechanism. |
| 4 | `c1q_apoe_trem2_lgals3` | supporting | Myeloid immunoregulatory supporting axis with mixed late/tumor support. |
| 5 | `inflammatory_il1_tnf_cxcl8` | benchmark/secondary | Published inflammatory-niche benchmark/positive control, not primary novelty. |

The manuscript skeleton now organizes the results into five sections:

1. Public multi-cohort framework.
2. Specificity audit reframing the original SPP1-niche hypothesis.
3. Spatial candidate-axis analysis nominating MIF-CD74/CXCR4.
4. Multi-cohort separation of communication axis versus macrophage-state readout.
5. In-silico perturbation prioritizing MIF and CD74.

## 2026-05-30 Manuscript Figure Drafts

Generated the first manuscript-facing figure drafts from existing result tables.

Code updates:

- Added `src/luad_niche/figure_data.py`.
- Added `scripts/plot_main_figures.py`.
- Added `tests/test_figure_data.py`.
- Added `docs/figure_plan.md`.

Command:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\plot_main_figures.py
```

Figure exports:

- `results/figures/figure1_workflow_dataset_composition.svg`
- `results/figures/figure1_workflow_dataset_composition.pdf`
- `results/figures/figure1_workflow_dataset_composition.tiff`
- `results/figures/figure1_workflow_dataset_composition.png`
- `results/figures/figure2_axis_evidence_perturbation.svg`
- `results/figures/figure2_axis_evidence_perturbation.pdf`
- `results/figures/figure2_axis_evidence_perturbation.tiff`
- `results/figures/figure2_axis_evidence_perturbation.png`

Source-data exports:

- `results/tables/figure1_dataset_composition_source.csv`
- `results/tables/figure2_priority_source.csv`
- `results/tables/figure2_evidence_heatmap_source.csv`
- `results/tables/figure2_perturbation_source.csv`

Figure interpretation:

- Figure 1 frames the analysis as a seven-cohort public-data workflow with ordered spatial and single-nucleus progression cohorts as the central evidence.
- Figure 2 shows that `MIF-CD74/CXCR4` is the top integrated axis and that continuous score-level perturbation prioritizes `MIF` and `CD74`.
- `SPP1/TREM2/PLA2G7` is retained as the macrophage-state readout rather than the primary communication novelty.

QA notes:

- SVG export keeps text editable through `svg.fonttype = none`.
- TIFF and PNG previews were generated from the same Python/matplotlib backend.
- The perturbation panel is explicitly labeled as target-prioritization evidence, not causal proof.

## 2026-05-30 Specificity Audit and Spatial Axis Figure Drafts

Extended the manuscript figure script to generate Figure 3 and Figure 4.

Code updates:

- Extended `src/luad_niche/figure_data.py` with specificity-audit, candidate-gene specificity, and spatial-axis stage helpers.
- Extended `tests/test_figure_data.py` for the new figure-data helpers.
- Extended `scripts/plot_main_figures.py` to export Figure 3 and Figure 4.
- Updated `docs/figure_plan.md`.

Command:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\plot_main_figures.py
```

New figure exports:

- `results/figures/figure3_specificity_audit.svg`
- `results/figures/figure3_specificity_audit.pdf`
- `results/figures/figure3_specificity_audit.tiff`
- `results/figures/figure3_specificity_audit.png`
- `results/figures/figure4_spatial_axis_progression.svg`
- `results/figures/figure4_spatial_axis_progression.pdf`
- `results/figures/figure4_spatial_axis_progression.tiff`
- `results/figures/figure4_spatial_axis_progression.png`

New source-data exports:

- `results/tables/figure3_signature_celltype_heatmap_source.csv`
- `results/tables/figure3_specificity_status_source.csv`
- `results/tables/figure3_candidate_gene_specificity_source.csv`
- `results/tables/figure4_axis_stage_heatmap_source.csv`
- `results/tables/figure4_focus_axis_trend_source.csv`
- `results/tables/figure4_late_axis_summary_source.csv`

Key figure findings:

- GSE131907 specificity audit shows broad `spp1_macrophage` is myeloid-enriched, but the refined SPP1 macrophage marker set has high epithelial/fibroblast relative signal.
- The original SPP1 macrophage marker audit keeps 14/30 expected myeloid genes, flags 15/30 as off-target, and has 1/30 missing gene.
- Candidate-gene top-cell-type calls support epithelial `MIF` and myeloid/immune `CD74`, `SPP1`, `TREM2`, `PLA2G7`, `IL1B`, `TNF`, and `CXCL8`.
- GSE307534 spatial-axis analysis shows `MIF-CD74/CXCR4` source-side enrichment increases across progression and has the top MIA/LUAD mean enrichment delta among displayed axis rows.
- `SPP1/TREM2/PLA2G7` receptor-side enrichment remains strong, supporting its interpretation as a macrophage-state readout rather than the main communication-axis novelty.

## 2026-05-30 Figure 5 Score-Level Target Prioritization Priority

Generated a dedicated Figure 5 for score-level in-silico target prioritization.

Code updates:

- Extended `src/luad_niche/figure_data.py` with perturbation dose-response, stage-loss, and method-concordance helpers.
- Extended `tests/test_figure_data.py` for the new perturbation helpers.
- Extended `scripts/plot_main_figures.py` to export Figure 5.
- Updated `docs/figure_plan.md`.

Command:

```powershell
$env:PYTHONPATH='D:\空间转录\luad_epithelial_macrophage_niche\src'
python scripts\plot_main_figures.py
```

New figure exports:

- `results/figures/figure5_virtual_perturbation_priority.svg`
- `results/figures/figure5_virtual_perturbation_priority.pdf`
- `results/figures/figure5_virtual_perturbation_priority.tiff`
- `results/figures/figure5_virtual_perturbation_priority.png`

New source-data exports:

- `results/tables/figure5_dose_response_source.csv`
- `results/tables/figure5_stage_loss_source.csv`
- `results/tables/figure5_method_concordance_source.csv`

Key figure findings:

- MIF score-level full dropout collapses source-side continuous coupling (`coupling_remaining = 0`).
- CD74 receptor-side score-level full dropout retains only about 26% of baseline coupling in the MIA/LUAD summary.
- CXCR4, SPP1, TREM2, and PLA2G7 have smaller coupling-loss effects than MIF/CD74.
- Continuous-coupling and top-quantile dropout perturbation approaches both prioritize MIF and CD74.

Interpretation boundary:

- This figure supports score-level in-silico target prioritization, not causal proof.
- Recommended downstream validation priority remains MIF/CD74 perturbation, with SPP1/TREM2/PLA2G7 as macrophage-state readouts.

## 2026-05-30 Literature Recheck, Automation, and Results Draft

Rechecked whether the revised `MIF-CD74` early-LUAD niche framing overlaps with existing work.

Key literature conclusion:

- The broad epithelial/alveolar progenitor plus proinflammatory macrophage niche story is highly overlapped with recent Cancer Cell work and should not be claimed as our main novelty.
- `MIF-CD74` has prior LUAD/NSCLC biological support, including inflammation-driven LUAD progression and broader TME communication studies.
- The defensible angle is narrower: public-data, specificity-audited spatial-axis prioritization that separates a `MIF-CD74` communication hypothesis from an `SPP1/TREM2/PLA2G7` macrophage-state readout and uses score-level in-silico target prioritization to prioritize `MIF` and `CD74`.

New documentation:

- `docs/literature_overlap_update_2026-05-30.md`
- `docs/results_draft.md`

Automation:

- Created heartbeat automation `luad-30`.
- Schedule: every 30 minutes.
- Purpose: inspect project status and continue the next analysis/manuscript step in the current thread, reporting only substantive progress.

Interpretation boundary:

- Keep perturbation language as score-level in-silico target prioritization, not causal proof.

## 2026-05-30 Manuscript Methods, Figure Legends, and Reproducibility Notes

Converted the current analysis state into manuscript-support documents.

New documentation:

- `docs/figure_legends.md`
- `docs/methods_draft.md`
- `docs/reproducibility.md`

Updated documentation:

- `README.md`

Key content:

- Figure legends now cover Figure 1-5 and the planned supplementary tables.
- The methods draft records the current public-data design, harmonized metadata, expression scoring, GSE131907 specificity audit, GSE307534 spatial candidate-axis analysis, evidence matrix, score-level in-silico target prioritization, figure generation, software, testing, and interpretation boundaries.
- The reproducibility note records command-level rebuild routes for core analyses, main figures, and tests.
- README now reflects the revised `MIF-CD74` lead story rather than the earlier broad epithelial progenitor-macrophage niche framing.

No code changes were made in this step.

## 2026-05-30 Integrated Manuscript Draft and Next-Actions Triage

Assembled the current manuscript components into a full internal working draft.

New documentation:

- `docs/manuscript_draft.md`
- `docs/manuscript_next_actions.md`

Updated documentation:

- `README.md`
- `docs/reproducibility.md`

Draft sections now present:

- title and short title;
- abstract;
- introduction;
- results;
- discussion;
- limitations;
- condensed methods;
- data and code availability;
- working reference list;
- interpretation guardrails.

Next analysis priority:

- Add sample-level and patient-aware statistical support for the main `GSE307534` spatial-axis trends before expanding to additional large datasets.

No code changes were made in this step.

## 2026-05-30 Patient-Aware Spatial Statistics, CD44 Perturbation, and Supplementary Tables

Added patient-aware statistics for GSE307534 and expanded receptor-side perturbation to include `CD44`.

Code updates:

- Added `src/luad_niche/spatial_statistics.py`.
- Added `scripts/summarize_gse307534_spatial_statistics.py`.
- Added `tests/test_spatial_statistics.py`.
- Added `scripts/export_supplementary_tables.py`.
- Added `tests/test_export_supplementary_tables.py`.
- Updated `config/candidate_mechanisms.yaml`.
- Updated `scripts/virtual_perturb_gse307534_axes.py`.
- Updated `src/luad_niche/figure_data.py`.
- Updated `scripts/plot_main_figures.py`.

Commands:

```powershell
$env:PYTHONPATH='<PROJECT_ROOT>\src'
python scripts\summarize_gse307534_spatial_statistics.py
python scripts\rank_candidate_mechanisms.py
python scripts\virtual_perturb_gse307534_axes.py
python scripts\continuous_perturb_gse307534_axes.py
python scripts\build_main_evidence_matrix.py
python scripts\plot_main_figures.py
python scripts\export_supplementary_tables.py
pytest tests -q
```

New outputs:

- `docs/gse307534_spatial_statistics_summary.md`
- `docs/supplementary_tables_index.md`
- `results/tables/gse307534_candidate_axis_spatial_with_patient.csv`
- `results/tables/gse307534_candidate_axis_late_vs_precursor_sample_stats.csv`
- `results/tables/gse307534_candidate_axis_late_vs_precursor_patient_stats.csv`
- `results/tables/gse307534_candidate_axis_paired_patient_differences.csv`
- `results/tables/gse307534_candidate_axis_paired_patient_stats.csv`
- `results/supplementary_tables/`

Key spatial-statistics result:

- Paired precursor-versus-late analysis included 20 patients.
- Source-side `MIF` mean paired enrichment delta: 0.141.
- 95% bootstrap CI: 0.094 to 0.188.
- Positive paired-patient fraction: 0.950.
- BH-adjusted Wilcoxon q-value: 0.00019.
- Combined receptor-side paired difference: 0.014, 95% CI -0.055 to 0.088; do not overstate it as a robust stage increase.

Updated receptor-side continuous perturbation ranking:

- `CD74` full dropout coupling loss: 74.0%.
- `CD44` full dropout coupling loss: 20.0%.
- `CXCR4` full dropout coupling loss: 6.0%.

Interpretation:

- The strongest patient-aware spatial progression evidence is epithelial/source-side `MIF`.
- `CD74` remains the top receptor-side follow-up target because the receptor score depends on it most strongly.
- `CD44` is an intermediate receptor-side candidate.
- `CXCR4` is supporting rather than central in the current model.
- Perturbation results remain score-level in-silico target prioritization, not causal proof.

QA:

- Figure 1-5 exports rebuilt.
- Figure 5 SVG contains editable `CD44` text.
- Supplementary-table exporter packaged 12 stable CSV tables.
- Full test suite: `91 passed in 4.74s`.

## 2026-05-30 MIF Expression-Matched and Tissue-Density Controls

Added source-side MIF spatial controls to test whether the paired progression signal behaves like generic high-expression or tissue-area effects.

Code updates:

- Added `src/luad_niche/spatial_controls.py`.
- Added `scripts/analyze_gse307534_mif_controls.py`.
- Added `tests/test_spatial_controls.py`.
- Extended `scripts/export_supplementary_tables.py` with control tables.

Command:

```powershell
$env:PYTHONPATH='<PROJECT_ROOT>\src'
python scripts\analyze_gse307534_mif_controls.py
python scripts\export_supplementary_tables.py
pytest tests -q
```

Design:

- Used one reference Visium section per ordered stage to choose 20 expression-matched single-gene controls.
- Excluded candidate-axis genes, discovery-signature genes, mitochondrial genes, and ribosomal genes.
- Recomputed epithelial-progenitor-neighborhood adjacency for `MIF` and matched genes across all 56 GSE307534 sections.
- Reused broad macrophage-signature adjacency and tissue-geometry metrics as additional controls.

New outputs:

- `docs/gse307534_mif_spatial_controls.md`
- `results/tables/gse307534_mif_expression_matched_controls.csv`
- `results/tables/gse307534_mif_random_control_adjacency.csv`
- `results/tables/gse307534_mif_random_control_paired_differences.csv`
- `results/tables/gse307534_mif_random_control_paired_stats.csv`
- `results/tables/gse307534_mif_random_control_summary.csv`
- `results/tables/gse307534_mif_density_control_summary.csv`
- `results/tables/gse307534_mif_broad_signature_control_paired_stats.csv`
- `results/figures/supplementary_figure_mif_spatial_controls.*`

Key results:

- `MIF` paired late-minus-precursor enrichment delta: 0.141.
- `MIF` percentile versus 20 expression-matched single-gene controls: 0.950.
- Empirical upper-tail p-value: 0.095.
- `PRDX1` is slightly above MIF; `BSG` is also near the top.
- MIF enrichment versus in-tissue spot count: Spearman r=0.122, p=0.369.
- Paired MIF enrichment change versus spot-count change: r=-0.310, p=0.184.
- MIF enrichment versus permutation-null mean: r=0.352, p=0.0077.

Interpretation update:

- Keep MIF as a prioritized epithelial/source-side candidate.
- Do not claim unique MIF specificity.
- Treat full MIF dropout cautiously: it is expected because the source panel contains only MIF.
- The receptor-side multi-gene ranking (`CD74 > CD44 > CXCR4`) is more discriminative.

QA:

- Supplementary-table exporter now packages 17 stable CSV tables.
- Supplementary control figure exports: SVG, PDF, PNG, and TIFF.
- Supplementary control SVG retains editable text labels.
- Full test suite: `94 passed in 4.80s`.

## 2026-05-30 Focused GSE308103 snRNA and GSE282617 Bulk Orthogonal Validation

Added focused non-spatial candidate-gene summaries and a manuscript-facing orthogonal-validation supplementary figure.

Code updates:

- Added `src/luad_niche/orthogonal_validation.py`.
- Added `scripts/summarize_gse308103_candidate_genes.py`.
- Added `scripts/plot_supplementary_orthogonal_validation.py`.
- Added `tests/test_orthogonal_validation.py`.
- Added a targeted `read_tabular_selected_genes()` reader to `src/luad_niche/sn_matrix.py`.
- Extended `scripts/export_supplementary_tables.py` to package `ST18-ST21`.

Performance note:

- The first targeted snRNA pass reused the older reader that recalculated per-cell totals from every raw matrix row.
- The optimized pass reuses existing assignment-level `total_counts` and parses numeric values only for the 10 requested genes.
- All 75 GSE308103 samples completed successfully after this optimization.

Commands:

```powershell
$env:PYTHONPATH='<PROJECT_ROOT>\src'
python scripts\summarize_gse308103_candidate_genes.py
python scripts\plot_supplementary_orthogonal_validation.py
python scripts\export_supplementary_tables.py
pytest tests -q
```

New outputs:

- `results/tables/gse308103_snrna_candidate_gene_sample_summary.csv`
- `results/tables/gse308103_snrna_candidate_gene_stage_summary.csv`
- `results/tables/gse308103_snrna_candidate_gene_genes_used.json`
- `results/tables/supplementary_figure_focused_orthogonal_validation_snrna_source.csv`
- `results/tables/supplementary_figure_focused_orthogonal_validation_bulk_source.csv`
- `results/figures/supplementary_figure_focused_orthogonal_validation.*`

Key snRNA results:

- Epithelial `MIF`: Normal 0.227 to LUAD 0.577.
- Macrophage `CD74`: Normal 4.495 to LUAD 4.766.
- Macrophage `CXCR4`: Normal 0.580 to LUAD 0.897.
- Macrophage `SPP1`: Normal 0.069 to LUAD 0.602.
- Macrophage `TREM2`: Normal 0.510 to LUAD 0.755.
- Macrophage `PLA2G7`: Normal 0.305 to LUAD 0.507.
- `CD44`, `IL1B`, `TNF`, and `CXCL8` did not show a uniform Normal-to-LUAD increase in the focused macrophage summaries.

Key bulk results:

- `MIF`: Normal 67.17 to IAC 123.60.
- `CXCR4`: Normal 27.36 to IAC 105.67.
- `SPP1`: Normal 3.53 to IAC 56.62.
- `CD74`: Normal 870.28 to IAC 893.11.
- `CD44`: Normal 22.63 to IAC 11.79.

Interpretation update:

- Non-spatial cohorts support source-side `MIF` progression and macrophage-state remodeling.
- Receptor-side genes do not form a uniform monotonic progression program.
- `CD74` remains the lead receptor-side follow-up target because of score-level spatial dependency ranking, not a strong monotonic expression trend.
- The focused panels are orthogonal trend support and boundary-setting, not evidence of cell-cell contact or causal perturbation.

QA:

- Supplementary Figure 2 exports were generated as SVG, PDF, PNG, and TIFF.
- SVG text labels remain editable.
- Supplementary-table exporter now packages 21 stable CSV tables.
- Full test suite: `101 passed in 4.56s`.

## 2026-05-31 Reference Audit and Manuscript Claim Tightening

Added a structured literature audit to support the current narrow manuscript frame.

New outputs:

- `docs/reference_audit_2026-05-31.md`
- `results/tables/reference_audit_2026-05-31.csv`

Updated documents:

- `docs/manuscript_draft.md`
- `docs/results_draft.md`
- `docs/decision_log.md`

Key literature conclusions:

- The closest overlap remains Peng et al. *Cancer Cell* 2026 on alveolar progenitor/proinflammatory niches in lung precursor lesions.
- MIF-CD74 has prior LUAD biological support, especially TNF-dependent inflammation upregulating MIF-CD74, so novelty must be framed as specificity-audited public-data prioritization rather than de novo pathway discovery.
- CD44 and CXCR4 are mechanistically defensible receptor-side comparators through MIF receptor-complex biology.
- CD74 behavior in early LUAD spatial literature is complex, reinforcing the current boundary that receptor-side genes are target-prioritization candidates rather than a proven monotonic progression program.
- `SPP1+` and `TREM2+` macrophage literature supports retaining `SPP1/TREM2/PLA2G7` as a macrophage-state readout.

Interpretation update:

- No change to the main computational conclusion.
- The reference audit strengthens the manuscript guardrails and makes the novelty claim more defensible.
- Score-level perturbation remains described as score-level in-silico target prioritization, not wet-lab perturbation evidence.

QA:

- Reference-audit CSV parsed successfully with 9 unique entries.
- Full test suite after document/table updates: `101 passed in 4.64s`.

## 2026-05-31 MIF Spatial Covariate Sensitivity

Added covariate-adjusted sensitivity models for source-side `MIF` spatial enrichment in GSE307534.

Code updates:

- Extended `src/luad_niche/spatial_controls.py` with `fit_ols_sensitivity()` and `build_paired_change_table()`.
- Added `scripts/analyze_gse307534_mif_covariates.py`.
- Extended `tests/test_spatial_controls.py`.
- Extended `scripts/export_supplementary_tables.py` with `ST22-ST23`.

Commands:

```powershell
$env:PYTHONPATH='<PROJECT_ROOT>\src'
python scripts\analyze_gse307534_mif_covariates.py
python scripts\export_supplementary_tables.py
pytest tests -q
```

New outputs:

- `docs/gse307534_mif_covariate_sensitivity.md`
- `results/tables/gse307534_mif_covariate_sensitivity_models.csv`
- `results/tables/gse307534_mif_covariate_paired_changes.csv`

Key results:

- Sample-level late-versus-precursor coefficient: 0.133 unadjusted.
- Sample-level coefficient after adjusting for in-tissue spot count, permutation-null mean, and neighborhood radius: 0.124.
- Paired-patient late-minus-precursor effect: 0.141.
- Paired-patient effect remained 0.141 after adjustment for changes in spot count, null mean, and radius.

Interpretation update:

- Source-side MIF spatial progression remains directionally stable after basic geometry and null-background adjustment.
- This strengthens robustness, but remains exploratory and non-causal.
- Score-level perturbation remains score-level in-silico target prioritization.

QA:

- Supplementary-table exporter now packages 23 stable CSV tables.
- Full test suite: `104 passed in 4.72s`.

## 2026-05-31 Submission-Style Manuscript Draft

Created a cleaner reader-facing manuscript draft while preserving the internal working draft.

New output:

- `docs/manuscript_submission_draft.md`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Drafting changes:

- Reframed the abstract around the precise gap: public-data prioritization after specificity auditing, rather than generic inflammatory niche discovery.
- Organized Results as an evidence ladder: multi-cohort design, specificity audit, patient-aware spatial statistics, score-level in-silico target prioritization, controls, and orthogonal validation.
- Strengthened Discussion guardrails for `MIF-CD74` novelty, receptor-side interpretation, and the limits of score-level gene dropout.
- Preserved the required boundary that perturbation outputs are score-level in-silico target prioritization, not wet-lab perturbation or causal validation.

Next recommended implementation step:

- Generate a supplementary figure comparing the original contaminated `SPP1` macrophage signature with the specificity-filtered version, with source data.

QA:

- Checked wording around perturbation claims in the new draft and updated documents.
- Full test suite after document updates: `104 passed in 4.75s`.

## 2026-05-31 SPP1 Signature Specificity-Refinement Figure

Generated a supplementary figure and source-data package documenting why the original refined `SPP1` macrophage signature was demoted from a primary communication mechanism to a macrophage-state readout.

Code updates:

- Added `prepare_signature_refinement_source()` and `prepare_signature_refinement_status_summary()` to `src/luad_niche/figure_data.py`.
- Added `scripts/plot_supplementary_signature_refinement.py`.
- Extended `tests/test_figure_data.py`.
- Extended `scripts/export_supplementary_tables.py` and `tests/test_export_supplementary_tables.py` with `ST24-ST26`.

New outputs:

- `results/figures/supplementary_figure_spp1_signature_refinement.svg`
- `results/figures/supplementary_figure_spp1_signature_refinement.pdf`
- `results/figures/supplementary_figure_spp1_signature_refinement.tiff`
- `results/figures/supplementary_figure_spp1_signature_refinement.png`
- `results/tables/supplementary_figure_spp1_signature_refinement_gene_source.csv`
- `results/tables/supplementary_figure_spp1_signature_refinement_status_source.csv`
- `results/tables/supplementary_figure_spp1_signature_refinement_celltype_source.csv`

Updated documents:

- `docs/figure_legends.md`
- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`
- `docs/supplementary_tables_index.md`

Key result:

- The original 30-gene refined `SPP1` macrophage signature contains 14 expected myeloid genes, 15 off-target genes, and 1 missing gene; specificity filtering retains 14 expected myeloid genes.

Interpretation update:

- This strengthens the manuscript's hypothesis-refinement logic and supports treating `SPP1/TREM2/PLA2G7` as a macrophage-state readout rather than the primary communication novelty.

QA:

- Re-rendered main figures after terminology cleanup.
- Checked docs, scripts, source, and tests for stale perturbation-overclaim wording; no positive causal-claim wording remained.
- Full test suite after figure/code updates: `107 passed in 4.83s`.

## 2026-05-31 Manuscript Claim-Evidence Map

Added a reviewer-facing claim-evidence-boundary map to make the manuscript easier to judge and safer to revise.

New output:

- `docs/manuscript_claim_evidence_map.md`

Updated documents:

- `docs/manuscript_submission_draft.md`
- `docs/manuscript_next_actions.md`
- `docs/decision_log.md`
- `docs/reproducibility.md`

Key content:

- Mapped each major claim to the relevant figures, supplementary figures, supplementary tables, and interpretation boundary.
- Added explicit callouts to Supplementary Fig. 1 for MIF controls and Supplementary Fig. 2 for non-spatial orthogonal validation in the submission-style Results.
- Updated next actions to prioritize journal-format tightening rather than additional broad data downloads.

QA:

- Full test suite after documentation updates: `107 passed in 4.74s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements, not positive causal claims.

## 2026-05-31 Editorial Summary and Cover-Letter Seed

Converted the claim-evidence map into a journal-neutral editorial package for later submission adaptation.

New output:

- `docs/editorial_summary_cover_letter_seed.md`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Key content:

- Editor-facing summary of the manuscript's bounded contribution.
- Short significance statement.
- Cover-letter seed.
- Graphical abstract logic.
- Target-journal framing options and do-not-overstate checklist.

Interpretation guardrail:

- The document keeps `MIF-CD74` as a follow-up candidate and describes perturbation results only as score-level in-silico target prioritization.

QA:

- Full test suite after editorial-package documentation updates: `107 passed in 4.55s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements, not positive causal claims.

## 2026-05-31 Compact Submission Draft

Created a shorter journal-neutral manuscript draft to bridge the full reader-facing draft and a future target-journal submission.

New output:

- `docs/manuscript_compact_submission_draft.md`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Drafting changes:

- Compressed the Abstract, Introduction, Results, Discussion, and Methods Summary while retaining the key quantitative evidence.
- Preserved the core interpretation: `MIF-CD74` is a follow-up candidate, `SPP1/TREM2/PLA2G7` is a macrophage-state readout, and score-level in-silico target prioritization is model-dependency ranking rather than causal validation.
- Updated next actions so the next manuscript step is exact target-journal formatting rather than more broad expansion analysis.

QA:

- Full test suite after compact-draft documentation updates: `107 passed in 4.56s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements, not positive causal claims.

## 2026-05-31 Target-Journal Formatting Route

Added a target-journal decision and formatting route for the current pure-bioinformatics manuscript package.

New output:

- `docs/target_journal_formatting_route.md`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Key content:

- Recommended Communications Biology as the first target-specific formatting route.
- Kept npj Precision Oncology and iScience as higher-risk alternatives, Scientific Reports as a fallback, and Cancer Research Communications as a later validation-oriented option.
- Marked BMC Cancer as not recommended for the current pure-bioinformatics manuscript because its scope guidance explicitly discourages public-database-only prediction papers without validation.
- Preserved cautious wording: `MIF-CD74` is a follow-up candidate, and the gene-dependency analysis remains score-level in-silico target prioritization.

Next recommended implementation step:

- Create `docs/manuscript_communications_biology_draft.md` from `docs/manuscript_compact_submission_draft.md`.

QA:

- Full test suite after target-route documentation updates: `107 passed in 4.78s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements, not positive causal claims.

## 2026-05-31 Communications Biology Working Draft

Converted the compact journal-neutral draft into a Communications Biology-oriented working manuscript.

New output:

- `docs/manuscript_communications_biology_draft.md`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Key content:

- Added a target-specific title, abstract, Introduction, Results, Discussion, Methods, Data Availability, Code Availability, current figure set, and author notes.
- Kept the central conclusion bounded to public-data prioritization of a follow-up `MIF-CD74` candidate.
- Preserved the distinction between source-side `MIF` spatial progression and receptor-side `CD74` score dependence.
- Described the gene-dependency analysis as score-level in-silico target prioritization only.

Next recommended implementation step:

- Add final formatted references, then audit figure callouts, source-data labels, data availability, and code availability for Communications Biology submission readiness.

QA:

- Communications Biology working draft size: 3,001 words and 23,148 characters.
- Full test suite after draft-generation updates: `107 passed in 4.66s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements, not positive causal claims.

## 2026-05-31 Communications Biology Reference Pass

Added the first target-specific numbered citation pass to the Communications Biology working draft.

New output:

- `results/tables/communications_biology_reference_list.csv`

Updated documents:

- `docs/manuscript_communications_biology_draft.md`
- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`
- `docs/decision_log.md`

Key content:

- Inserted bracketed numbered citations into the Abstract, Introduction, Results, and Discussion.
- Added eight core references covering closest-overlap early-LUAD spatial biology, discovery scRNA/spatial context, MIF-CD74/CD44/CXCR4 receptor biology, CD74 spatial-boundary interpretation, and TREM2/SPP1 macrophage-state context.
- Exported a machine-readable reference-role table for later citation-manager or submission-system cleanup.

Next recommended implementation step:

- Audit figure callouts, source-data labels, data availability, and code availability for Communications Biology submission readiness.

QA:

- Communications Biology working draft after references: 3,390 words and 25,640 characters.
- Citation-marker audit found references 1-8 all used in the manuscript.
- Full test suite after reference-pass updates: `107 passed in 4.83s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements, not positive causal claims.

## 2026-05-31 Communications Biology Submission-Readiness Audit

Audited the current target-specific manuscript package for submission readiness.

New output:

- `docs/communications_biology_submission_readiness_audit.md`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Key content:

- Confirmed that Figures 1-5 and Supplementary Figures 1-3 are present in SVG, PDF, TIFF, and PNG.
- Confirmed that all main-figure source-data CSVs are present.
- Confirmed that Supplementary Figures 2-3 have uniformly named source-data CSVs.
- Noted that Supplementary Figure 1 is supported by existing packaged MIF-control tables but could benefit from uniformly named panel-level source-data files.
- Identified public code availability as the main remaining submission blocker because the draft currently points to local `scripts/` and `src/luad_niche/` paths.

Next recommended implementation step:

- Build a clean public-code/package plan for Communications Biology submission.

QA:

- Full test suite after submission-readiness audit updates: `107 passed in 4.77s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements, not positive causal claims.

## 2026-05-31 Communications Biology Chinese Review Draft

Translated the Communications Biology working draft into a synchronized Chinese review version for internal reading and discussion.

New output:

- `docs/manuscript_communications_biology_draft_zh.md`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Key content:

- Preserved the same section structure as the English working draft.
- Retained gene symbols, accession IDs, figure labels, statistics, and numbered references.
- Kept the interpretation boundary: `MIF-CD74` is a follow-up candidate, and score-level in-silico target prioritization is not causal or wet-lab perturbation evidence.

Next recommended implementation step:

- Use the Chinese review draft to review scientific logic with collaborators, then update the English Communications Biology draft where wording or framing needs tightening.

QA:

- Chinese review draft size: 14,674 characters across 105 lines.
- Chinese guardrail scan found no positive claims of virtual knockout, causal proof, or completed experimental validation.
- Full test suite after Chinese-draft documentation updates: `107 passed in 4.65s`.

## 2026-05-31 Chinese Word Review Export

Exported the synchronized Chinese Communications Biology review draft to Word with manuscript text, figures, tables, and Zotero-ready references.

New outputs:

- `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`
- `docs/word_export_notes.md`
- `results/references/communications_biology_references.ris`
- `scripts/export_word_manuscript.py`

Updated documents:

- `docs/manuscript_next_actions.md`
- `docs/reproducibility.md`

Key content:

- The Word file contains the Chinese manuscript text, numbered in-text citations, a reference list, 8 figure images, and 5 tables.
- The tables include dataset composition, integrated candidate-axis ranking, GSE307534 paired-patient spatial statistics, score-level in-silico target-prioritization ranking, and supplementary table manifest.
- The RIS file is ready for Zotero import.
- Local Zotero API was not enabled, so citations are exported as readable numbered citations rather than live Zotero Word-plugin field codes.

QA:

- Word document inspection: 133 paragraphs, 5 tables, and 8 inline images.
- Zotero-ready RIS file generated from the DOI/reference list.
- Full test suite after Word-export script and documentation updates: `107 passed in 4.74s`.
- Guardrail-wording scan found only explicit "avoid/do not" boundary statements and QA notes, not positive causal claims.

Next recommended implementation step:

- If live Zotero citations are required, import the RIS into Zotero and use the Word Zotero plugin to replace numbered citations with live fields.

## 2026-05-31 Zotero MCP Import and Word Field Export

After Codex restart, Zotero MCP write access was available through the Zotero Web API credentials. The 8 manuscript references were imported from DOI metadata into the configured Zotero user library and mapped into the Word export script.

Updated outputs:

- `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`
- `docs/manuscript_communications_biology_draft_zh_zotero_fields.docx`
- `results/references/communications_biology_references.ris`
- `scripts/export_word_manuscript.py`
- `docs/word_export_notes.md`
- `docs/reproducibility.md`
- `docs/reference_expansion_2026-05-31.md`

Key content:

- The Zotero-field Word draft contains Zotero CSL field instructions for the manuscript's numbered citations.
- The field URIs now use the configured Zotero user library item URIs rather than local-only `users/0` URIs.
- The stable Word draft remains available as a readable numbered-citation version.
- The Codex Zotero MCP config was changed to Web API mode to avoid local `userID=0` and Web user ID conflicts.

QA:

- Zotero Web library DOI audit found exactly one matching item for each of the 8 references.
- Zotero-field Word document inspection: 5 tables, 8 inline images, 8 media files, 15 Zotero citation fields, 23 configured-user Zotero URI hits, and 0 local `users/0` URI hits.
- Full test suite with `PYTHONPATH=src`: `107 passed in 5.01s`.
- Guardrail-wording scan found only explicit avoid/do-not boundary statements and QA notes, not positive causal claims.

## 2026-05-31 Reference Expansion to 30 Papers

Expanded the manuscript reference set from 8 to 30 DOI-verified papers after literature review across early LUAD progression, single-cell/spatial LUAD ecology, MIF-CD74 biology, macrophage-state context, spatial transcriptomics and cell-cell communication methodology.

Updated outputs:

- `docs/manuscript_communications_biology_draft.md`
- `docs/manuscript_communications_biology_draft_zh.md`
- `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`
- `docs/manuscript_communications_biology_draft_zh_zotero_fields.docx`
- `results/tables/communications_biology_reference_list.csv`
- `results/references/communications_biology_references.ris`
- `scripts/export_word_manuscript.py`
- `docs/word_export_notes.md`
- `docs/reproducibility.md`

Key content:

- Added 22 references covering `GSE131907` context, early-stage LUAD single-cell ecology, preneoplasia-to-invasive immune evolution, subsolid/ground-glass nodule biology, MIF/CD74 lung-cancer evidence, TREM2 tumor-myeloid biology and spatial/cell-cell communication methods.
- Added title and author metadata columns to the reference registry so Word/RIS export no longer depends on hard-coded fallback metadata for new references.
- Imported all 30 DOI references into the configured Zotero Web library and mapped item keys for Zotero field-code export.
- Updated English and Chinese drafts with expanded citation markers rather than only appending uncited references.

QA:

- Zotero Web library DOI audit found exactly one matching item for each of the 30 references.
- Reference registry contains 30 rows; RIS export contains 30 `TY  - JOUR` entries.
- Zotero-field Word document inspection: 5 tables, 8 inline images, 8 media files, 16 Zotero citation fields, 51 configured-user Zotero URI hits, and 0 local `users/0` URI hits.
- Full test suite with `PYTHONPATH=src`: `107 passed in 4.68s`.
- Guardrail-wording scan found only explicit avoid/do-not boundary statements and QA notes, not positive causal claims.

## 2026-06-03 GRN-Level Virtual Perturbation Prioritization

Reviewed scTenifoldKnk and CellOracle-style GRN perturbation options and selected a scTenifoldKnk-inspired first implementation for this project. The reason is practical and biological: CellOracle is strongest for TF-driven cell-identity perturbation with base-GRN/motif or chromatin-accessibility priors, while the current target list is dominated by non-TF ligand, receptor, and macrophage-state genes.

Implemented and ran a transparent GSE308103 workflow:

- Code: `src/luad_niche/grn_perturbation.py`
- Tests: `tests/test_grn_perturbation.py`
- Runner: `scripts/run_gse308103_grn_virtual_perturbation.py`
- Method note: `docs/grn_virtual_perturbation_2026-06-03.md`
- Outputs: `results/tables/gse308103_grn_virtual_perturbation_*.csv` and `results/tables/gse308103_grn_virtual_perturbation_config.json`

Run configuration:

- Stages: `MIA,LUAD`
- Classes: `epithelial,macrophage`
- Max cells per class: 6,000
- Network: positive Pearson coexpression GRN, correlation threshold 0.05, row-normalized edges
- Propagation: 3 steps, restart 0.15, decay 0.5

Network summary:

- Epithelial: 6,000 cells, 116 input genes, 111 network genes, 3,268 edges; tested target `MIF`.
- Macrophage: 6,000 cells, 116 input genes, 115 network genes, 3,940 edges; tested targets `CD74`, `CD44`, `CXCR4`, `SPP1`, `TREM2`, `PLA2G7`.

Top GRN-level target neighborhoods:

- `CD44`: C1Q macrophage signature, mean impact 0.01774; top impacted genes `MARCO`, `MCEMP1`, `INHBA`, `LTA4H`, `PPARG`.
- `TREM2`: C1Q macrophage signature, mean impact 0.01773; top impacted genes `APOE`, `APOC1`, `ACP5`, `CYP27A1`, `LGALS3`.
- `PLA2G7`: SPP1 macrophage signature, mean impact 0.01650; top impacted genes `APOE`, `ACP5`, `APOC1`, `CHIT1`, `CYP27A1`.
- `CD74`: C1Q macrophage signature, mean impact 0.01508; top impacted genes `C1QA`, `C1QB`, `C1QC`, `APOE`, `ACP5`.
- `MIF`: proliferating epithelial signature, mean impact 0.01247; top impacted genes `MDK`, `EPCAM`, `KRT8`, `CD24`, `PERP`.
- `CXCR4`: SPP1 macrophage signature, mean impact 0.01210; top impacted genes `PRDM1`, `GPR183`, `ADAMDEC1`, `PLA2G7`, `APOE`.
- `SPP1`: SPP1 macrophage signature, mean impact 0.01022; top impacted genes `MIF`, `ANXA2`, `TREM2`, `G0S2`, `CAMK1`.

Interpretation update:

- The GRN layer does not simply repeat the spatial score-dropout ranking.
- It supports keeping `CD74` as a receptor-side candidate linked to a C1Q/APOE/ACP5 macrophage network.
- It suggests `CXCR4` and `PLA2G7` are better interpreted as SPP1-like macrophage-state network modulators/readouts.
- This is GRN-level virtual perturbation prioritization only; it does not establish causal knockout biology.

QA:

- Re-ran `python scripts\run_gse308103_grn_virtual_perturbation.py`: wrote 794 gene-effect rows, 26 signature-impact rows, and 7 target-ranking rows.
- Full test suite with `PYTHONPATH=src`: `112 passed in 4.76s`.
- Guardrail scan found only implementation names or explicit boundary statements, not positive claims that computational perturbation proves causality.

## 2026-06-03 Supplementary GRN Figure and Table

Generated a manuscript-facing supplementary figure and packaged target-ranking table for the GSE308103 GRN-level virtual perturbation layer.

Outputs:

- Figure script: `scripts/plot_supplementary_grn_virtual_perturbation.py`
- Figure exports: `results/figures/supplementary_figure_grn_virtual_perturbation.svg`, `.pdf`, `.tiff`, and `.png`
- Figure source data: `results/tables/supplementary_figure_grn_virtual_perturbation_target_source.csv`
- Figure source data: `results/tables/supplementary_figure_grn_virtual_perturbation_top_gene_source.csv`
- Packaged supplementary table: `results/supplementary_tables/ST27_gse308103_grn_virtual_perturbation_target_ranking.csv`

Figure logic:

- Panel a ranks targets by top-signature mean propagated impact score.
- Panel b displays the top five propagated gene neighborhoods per target.
- Panel c records the GSE308103 epithelial and macrophage GRN sizes.

Interpretation:

- `CD74` remains linked to a C1Q/APOE/ACP5 macrophage network.
- `CXCR4` and `PLA2G7` are more connected to SPP1-like macrophage-state neighborhoods.
- The figure caption explicitly frames the result as GRN-level virtual perturbation prioritization, not wet-lab perturbation or causal validation.

QA:

- Visual check opened `results/figures/supplementary_figure_grn_virtual_perturbation.png`; the figure is non-empty and text is readable.
- Re-ran `python scripts\export_supplementary_tables.py`: packaged 27 supplementary tables including `ST27`.
- Full test suite with `PYTHONPATH=src`: `115 passed in 4.87s`.
- Guardrail scan found only explicit boundary statements, not positive causal claims.

## 2026-06-03 GRN Robustness and Cross-Dataset Validation

Addressed the dataset-scope question for the GRN-level virtual perturbation layer. The primary GRN perturbation analysis remains anchored in `GSE308103` because it has ordered early-LUAD snRNA raw-count matrices and enough epithelial/macrophage cells. Other datasets were assigned external validation roles rather than treated as parallel GRN-discovery cohorts.

Robustness outputs:

- `scripts/run_gse308103_grn_robustness.py`
- `results/tables/gse308103_grn_virtual_perturbation_robustness_detail.csv`
- `results/tables/gse308103_grn_virtual_perturbation_robustness_summary.csv`
- `results/tables/gse308103_grn_virtual_perturbation_robustness_network_summary.csv`
- Packaged table: `results/supplementary_tables/ST28_gse308103_grn_virtual_perturbation_robustness_summary.csv`

Robustness configuration:

- Thresholds: 0.03, 0.05, 0.08
- Seeds: 7, 11, 23
- Cells: up to 6,000 per class
- Runs: 9 target-ranking runs

Key robustness results:

- `CD74`: median rank 4, rank range 3-5, top signature always `c1q_macrophage_vs_other_macrophage`.
- `CXCR4`: median rank 6, rank range 5-7, top signature always `spp1_macrophage_vs_other_macrophage`.
- `PLA2G7`: median rank 3, rank range 1-4, top signature always `spp1_macrophage_vs_other_macrophage`.
- `MIF`: median rank 5, rank range 5-6, top signature always `proliferating_epithelial_vs_other_epithelial`.

Cross-dataset validation outputs:

- `scripts/build_grn_cross_dataset_validation.py`
- `results/tables/grn_cross_dataset_signature_validation.csv`
- `results/tables/grn_cross_dataset_signature_validation_summary.csv`
- Packaged table: `results/supplementary_tables/ST29_grn_cross_dataset_signature_validation.csv`
- Packaged table: `results/supplementary_tables/ST30_grn_cross_dataset_signature_validation_summary.csv`

Cross-dataset interpretation:

- `MIF`/proliferating epithelial state is higher in comparison groups across GSE307534, GSE189357, and GSE164789, and is epithelial-specific in GSE131907.
- `SPP1`-like macrophage state linked to `CXCR4`, `PLA2G7`, and `SPP1` is higher in comparison groups in the two scRNA external datasets, but not in the GSE307534 late-minus-precursor spatial score summary; this should be written as supportive but not uniformly monotonic.
- `C1Q` macrophage state linked to `CD74`, `CD44`, and `TREM2` is myeloid-specific in GSE131907 but trends lower in the comparison groups across the tested state-fraction/score summaries; this is a boundary-setting result and argues against saying the C1Q state itself expands monotonically.

Packaging:

- `python scripts\export_supplementary_tables.py` now packages 30 supplementary tables through `ST30`.

QA:

- Full test suite with `PYTHONPATH=src`: `121 passed in 4.89s`.
- Guardrail scan found only explicit boundary statements, not positive causal claims.

## 2026-06-03 Manuscript Integration and Perturb-seq Feasibility Decision

Integrated the GRN-level virtual perturbation and robustness results into the Communications Biology working draft and recorded the perturbation-response model decision.

Manuscript updates:

- Added a Results subsection, `GRN-level virtual perturbation stratifies target-linked expression neighborhoods`, to `docs/manuscript_communications_biology_draft.md`.
- Added a Methods subsection for GSE308103 GRN-level virtual perturbation prioritization, robustness settings, and cross-dataset validation.
- Added a Discussion boundary paragraph explaining that GEARS/scGen/CPA are not used as current primary evidence because a matched public Perturb-seq dataset for the early-LUAD epithelial-myeloid MIF-CD74 context was not identified.
- Added Supplementary Figure 4 to the current figure set and updated figure/reproducibility wording.

Perturb-seq feasibility output:

- New note: `docs/perturbseq_gears_scgen_cpa_feasibility_2026-06-03.md`.
- Current decision: do not add GEARS/scGen/CPA unless future public or in-house perturbation data contain relevant targets and cellular context.
- Future trigger targets: `MIF`, `CD74`, `CD44`, `CXCR4`, `SPP1`, `TREM2`, `PLA2G7`, or close pathway genes in lung epithelial, macrophage, co-culture, organoid, or tumor-microenvironment context.

Reference updates:

- Expanded `results/tables/communications_biology_reference_list.csv` from 30 to 36 references.
- Updated `results/references/communications_biology_references.ris` with scTenifoldKnk, CellOracle, scPerturb, GEARS, scGen, and CPA references.

Interpretation:

- The current strongest computational conclusion remains `MIF-CD74` candidate prioritization, not causal validation.
- The GRN layer adds network context: `CD74` maps to a C1Q/APOE/ACP5 macrophage neighborhood, while `CXCR4` and `PLA2G7` map more strongly to SPP1-like macrophage-state neighborhoods.
- Public Perturb-seq model use should be conditional, not automatic.

## 2026-06-03 Manuscript and Figure Storyline Recheck

Re-summarized the current manuscript and figure architecture after adding the GRN-level virtual perturbation and Perturb-seq feasibility decisions.

Output:

- `docs/manuscript_figure_storyline_2026-06-03.md`

Key recommendation:

- Keep five main figures and four supplementary figures.
- The strongest Results order is: public framework -> specificity audit -> integrated ranking -> paired-patient spatial `MIF` result -> `MIF` controls/sensitivity -> score-level target prioritization -> GRN-level virtual perturbation prioritization -> orthogonal cohort support.
- Figure 1 should be regenerated or lightly revised before final submission because the current figure files were produced before the GRN-level prioritization layer was added, while the current legend now mentions it.

## 2026-06-03 Nature/Cell-Style Figure Revision

Revised the manuscript figure style and regenerated publication-format exports.

Code/style changes:

- Updated `src/luad_niche/nature_figure_style.py` with a more unified low-saturation Nature/Cell-style palette.
- Updated `scripts/plot_main_figures_nature.py` so Figure 1 now includes both score-level and GRN-level target-prioritization layers.
- Updated `scripts/plot_main_figures_nature.py` so Figure 5 uses the visible label `Score reduction`.
- Updated `scripts/plot_supplementary_grn_virtual_perturbation.py` to use the shared Nature-style palette/export helpers.
- Updated `tests/test_figure_data.py` so the visual contract tests semantic color mapping and valid hex colors rather than locking the older palette hex values.

Regenerated exports:

- `results/figures/nature_redesign/nature_figure1_workflow_dataset_composition.*`
- `results/figures/nature_redesign/nature_figure2_axis_evidence_perturbation.*`
- `results/figures/nature_redesign/nature_figure3_specificity_audit.*`
- `results/figures/nature_redesign/nature_figure4_spatial_axis_progression.*`
- `results/figures/nature_redesign/nature_figure5_virtual_perturbation_priority.*`
- `results/figures/supplementary_figure_grn_virtual_perturbation.*`

Visual QA:

- Previewed Figure 1, Figure 2, Figure 5 and Supplementary Figure 4 PNG outputs.
- Figure 1 now clearly shows `GRN` as a separate GSE308103 target-neighborhood prioritization layer.
- The main red/teal/blue/green/gold/violet palette is more consistent across main and supplementary figures.
- Figure 5 visible wording now aligns with score-level in-silico target-prioritization boundaries.

QA:

- Full test suite with `PYTHONPATH=src`: `121 passed in 4.74s`.
- Guardrail scan found boundary statements only; no positive claim that computational perturbation proves causality.

## 2026-06-03 Results Order and Legend Synchronization

Updated the Communications Biology manuscript narrative after the Nature/Cell-style figure refresh.

Manuscript updates:

- Reordered `docs/manuscript_communications_biology_draft.md` so the `MIF` control and sensitivity section now appears immediately after the paired-patient spatial `MIF` result.
- Preserved the evidence order recommended in `docs/manuscript_figure_storyline_2026-06-03.md`: public framework -> specificity audit -> integrated ranking -> paired spatial `MIF` -> `MIF` controls -> score-level target prioritization -> GRN-level virtual perturbation prioritization -> orthogonal cohorts.

Figure legend updates:

- Revised Figure 1 to describe the six-step workflow now shown in the regenerated figure: cohort harmonization, specificity audit, spatial candidate-axis scoring, score-level target prioritization, GRN-level prioritization and bounded `MIF-CD74` candidate nomination.
- Revised Figure 2 and Figure 5 wording from broad perturbation language toward score-reduction, dropout and down-weighting language.
- Kept all score-level and GRN-level analyses framed as computational target prioritization rather than causal validation.

## 2026-06-03 Submission Consistency Audit

Performed a manuscript-facing consistency pass across the Communications Biology working draft, figure legends, supplementary-table index and reproducibility notes.

New audit note:

- `docs/submission_consistency_audit_2026-06-03.md`

Edits made:

- Changed the figure-legend working title from stronger axis phrasing to a bounded communication candidate claim.
- Revised Figure 4 wording so the main spatial progression claim is source-side `MIF` enrichment, while receptor-side and macrophage-state signals are framed as support or boundary-setting readouts.
- Revised Figure 5 wording so score-level target prioritization emphasizes receptor-side `CD74`; source-side `MIF` dropout is treated as the expected collapse of a single-gene score.
- Reordered the Methods narrative so `Control and sensitivity analyses` follows `Patient-aware progression statistics`, matching the Results evidence ladder.
- Updated supplementary-table and reproducibility wording from generic perturbation ranking toward score-reduction/dropout language.
- Updated `scripts/plot_main_figures.py` and `scripts/plot_main_figures_nature.py` so the visible Figure 5 title and panel annotation match the receptor-side `CD74` interpretation.
- Regenerated the standard main figures and Nature-redesigned main figures after the Figure 5 title synchronization.

Current remaining submission risks:

- Final journal citation formatting is still needed.
- Zotero-linked Word export should be regenerated after importing references 31-36 if a fully field-linked draft is required.
- Data Availability still needs a public repository or Supplementary Data destination before submission.

## 2026-06-03 Chinese Word Export Synchronization

Synchronized the Chinese review manuscript and regenerated Word review outputs after the GRN and Figure 5 interpretation updates.

Files updated:

- Rebuilt `docs/manuscript_communications_biology_draft_zh.md` as a 2026-06-03 synchronized Chinese review draft.
- Updated `scripts/export_word_manuscript.py` so Figure 5 caption emphasizes receptor-side `CD74` and Supplementary Figure 4 is inserted.
- Regenerated `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`.
- Regenerated `docs/manuscript_communications_biology_draft_zh_zotero_fields.docx`.
- Regenerated `results/references/communications_biology_references.ris`.
- Updated `docs/word_export_notes.md`.

Export QA:

- Stable Word draft contains the updated wording `受体侧 CD74`, `补充图 4`, and `GRN-level virtual perturbation`.
- Stable Word draft no longer contains the earlier overbroad Figure 5 wording.
- Zotero-field Word draft contains `ADDIN ZOTERO_ITEM CSL_CITATION`, Zotero user-item URIs for mapped references 1-30, CSL itemData, and the GRN-related reference content.
- RIS export contains 36 journal entries.

## 2026-06-04 Bilingual Manuscript Review Export

Created a formatted English-Chinese side-by-side review version of the Communications Biology working draft, with figures and key tables included for manuscript review.

Files updated or added:

- Updated the manuscript dates in `docs/manuscript_communications_biology_draft.md` and `docs/manuscript_communications_biology_draft_zh.md` to 2026-06-04.
- Added `scripts/export_bilingual_word_manuscript.py`.
- Added `tests/test_export_bilingual_word_manuscript.py`.
- Generated `docs/manuscript_communications_biology_bilingual_with_figures_tables.docx`.
- Updated `docs/word_export_notes.md`.

Export contents:

- Landscape Word layout with block-level English/Chinese manuscript text.
- Main Figures 1-5 from `results/figures/nature_redesign/`.
- Supplementary Figures 1-4, including the GRN-level virtual perturbation prioritization figure.
- Five key manuscript-facing tables: dataset composition, integrated candidate-axis ranking, paired spatial statistics, score-level in-silico target prioritization, and supplementary table manifest.

QA:

- The bilingual Word export contains `English`, `中文`, `score-level in-silico target prioritization`, `受体侧 CD74`, `Supplementary Figure 4`, `补充图 4`, and `GRN-level virtual perturbation`.
- The bilingual Word export contains 9 embedded media files.
- Guardrail scan confirmed that the old broad wording `MIF 和 CD74 排名靠前`, `因果证明`, and `实验验证了` are absent from the exported Word XML.
- Targeted exporter tests: `3 passed`.
- Full test suite after export updates: `124 passed in 14.14s`.

## 2026-06-04 Separate Chinese and English Word Exports

Generated separate single-language Word manuscripts in addition to the bilingual side-by-side review file.

Files updated or added:

- Added `scripts/export_english_word_manuscript.py`.
- Added `tests/test_export_english_word_manuscript.py`.
- Regenerated `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`.
- Regenerated `docs/manuscript_communications_biology_draft_zh_zotero_fields.docx`.
- Generated `docs/manuscript_communications_biology_draft_en_with_figures_tables.docx`.
- Regenerated `results/references/communications_biology_references.ris`.
- Updated `docs/word_export_notes.md`.

Export contents:

- Chinese Word export: Chinese manuscript text, numbered citations, 36-reference list, Main Figures 1-5, Supplementary Figures 1-4, and five key tables.
- English Word export: English manuscript text, numbered citations, 36-reference list, Main Figures 1-5, Supplementary Figures 1-4, and five key tables.

QA:

- Both single-language Word files contain 9 embedded media files.
- Both exports contain `score-level in-silico target prioritization` and `GRN-level virtual perturbation`.
- The Chinese export contains `补充图 4` and `受体侧 CD74`; the English export contains `Supplementary Figure 4` and `receptor-side CD74`.
- Guardrail scan confirmed that `ranks MIF and CD74`, `MIF 和 CD74 排名靠前`, `因果证明`, and `实验验证了` are absent from the single-language Word XML outputs.

## 2026-06-10 Original scTenifoldKnk Smoke Test

Tested whether the original CRAN `scTenifoldKnk` package can run on the
project's reduced GSE308103 matrices.

Environment:

- Rscript: `C:\Program Files\R\R-4.6.0\bin\Rscript.exe`
- R: `R version 4.6.0 (2026-04-24 ucrt)`
- Package: `scTenifoldKnk 1.0.3`

Code added:

- `scripts/run_original_sctenifoldknk_smoke.R`

Command:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' --vanilla scripts\run_original_sctenifoldknk_smoke.R
```

Outputs added:

- `results/tables/sctenifoldknk_smoke/gse308103_epithelial_mif_sctenifoldknk_diffregulation.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_epithelial_mif_sctenifoldknk_top_non_target_genes.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_macrophage_cd74_sctenifoldknk_diffregulation.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_macrophage_cd74_sctenifoldknk_top_non_target_genes.csv`
- `results/tables/sctenifoldknk_smoke/sctenifoldknk_original_smoke_summary.csv`
- `docs/sctenifoldknk_original_smoke_2026-06-10.md`

Result summary:

- `MIF` epithelial smoke test: 90 cells, 116 genes, 116 diffRegulation rows.
  Excluding `MIF` itself, the top network-affected genes included `KRT8`,
  `ANXA2`, `PERP`, `CD24`, `LGALS3`, `CLDN4`, `CLDN3`, and `EPCAM`.
- `CD74` macrophage smoke test: 90 cells, 116 genes, 116 diffRegulation rows.
  Excluding `CD74` itself, the top network-affected genes included `CYP27A1`,
  `TREM2`, `APOC1`, `PLA2G7`, `LPL`, `LRP1`, `LGALS3`, and `APOE`.

Interpretation boundary:

- The original package is locally usable and can be incorporated as an optional
  sensitivity analysis.
- Because this was a reduced smoke test and the best non-target adjusted
  P values were 0.305553 for `MIF` and 0.191556 for `CD74`, this result should
  remain GRN-level virtual perturbation prioritization / network-context
  evidence, not causal knockout validation.

## 2026-06-16 Original scTenifoldKnk Expanded Sensitivity And Manuscript Refresh

Extended the original CRAN `scTenifoldKnk` sensitivity analysis from the two
smoke-test targets to all seven candidate targets used in the manuscript:
`MIF`, `CD74`, `CD44`, `CXCR4`, `SPP1`, `TREM2`, and `PLA2G7`.

Code and commands:

- Added `scripts/run_original_sctenifoldknk_expanded.R`.
- Added `scripts/summarize_original_sctenifoldknk_expanded.py`.
- Ran:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' --vanilla scripts\run_original_sctenifoldknk_expanded.R
python scripts\summarize_original_sctenifoldknk_expanded.py
python scripts\export_supplementary_tables.py
python scripts\export_word_manuscript.py
python scripts\export_english_word_manuscript.py
```

New or refreshed outputs:

- `results/tables/sctenifoldknk_original_expanded/sctenifoldknk_original_expanded_summary.csv`
- `results/tables/sctenifoldknk_original_expanded/*_diffregulation.csv`
- `results/tables/sctenifoldknk_original_expanded/*_top_non_target_genes.csv`
- `results/tables/sctenifoldknk_original_expanded_interpretation.csv`
- `results/supplementary_tables/ST31_sctenifoldknk_original_expanded_interpretation.csv`
- `results/supplementary_tables/supplementary_table_manifest.csv`
- `docs/supplementary_tables_index.md`
- `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`
- `docs/manuscript_communications_biology_draft_zh_zotero_fields.docx`
- `docs/manuscript_communications_biology_draft_en_with_figures_tables.docx`

Result boundary:

- The original package was locally feasible for all seven reduced-panel runs.
- `MIF` showed an epithelial network-context sensitivity signal.
- Several macrophage target runs returned top non-target genes containing
  epithelial markers, consistent with a reduced-panel boundary.
- The original-package outputs are therefore included as supplementary
  sensitivity/boundary analysis, not causal knockout validation and not a
  replacement for the transparent Python GRN-level implementation.
