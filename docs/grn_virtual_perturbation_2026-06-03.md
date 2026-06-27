# GRN-Level Virtual Perturbation Note

Date: 2026-06-03

## Question

Can we add a GRN-level virtual perturbation layer to the early-LUAD epithelial-myeloid niche project, beyond the existing GSE307534 score-level in-silico target prioritization?

## Method comparison

| Method | Strength | Main fit to this project | Limitation for this project |
|---|---|---|---|
| scTenifoldKnk | Uses scRNA/snRNA expression to infer a single-cell GRN, removes a target gene's outgoing network edges, and compares WT versus virtual-KO networks. Good for expression-only target prioritization, including non-TF genes. | Best first fit. Our local `GSE308103` snRNA raw-count matrices are available, and our candidates include non-TFs such as `MIF`, `CD74`, `CD44`, `CXCR4`, `SPP1`, `TREM2`, and `PLA2G7`. | The full R/Python package is heavier than needed for a quick reproducible project layer; local `scTenifoldpy` smoke testing was slow. |
| CellOracle | Strong for TF-centered cell-identity perturbation. It builds context-dependent GRNs and simulates TF perturbation effects as cell-state transition vectors. | Good future option if we add TF regulators or matched/public chromatin priors. | Less direct for ligand/receptor genes and macrophage-state markers. It generally benefits from a base GRN built from scATAC-seq, bulk ATAC-seq, promoter databases, or a trusted TF-target list. |

Decision: use a transparent scTenifoldKnk-inspired route now, and reserve CellOracle for a later TF-focused extension.

## Implemented workflow

Files:

- `src/luad_niche/grn_perturbation.py`
- `tests/test_grn_perturbation.py`
- `scripts/run_gse308103_grn_virtual_perturbation.py`

Input:

- `data/interim/GSE308103/raw_counts/*.raw_counts.mtx.txt.gz`
- `results/tables/gse308103_snrna_cell_state_assignments.csv`
- `config/candidate_mechanisms.yaml`
- `results/tables/gse189357_refined_signature_genes_gse131907_specificity_filtered.json`

Command:

```powershell
python scripts\run_gse308103_grn_virtual_perturbation.py
```

Core parameters:

- Stages: `MIA,LUAD`
- Classes: `epithelial,macrophage`
- Max cells per class: 6,000
- Network: positive Pearson correlation GRN, correlation threshold 0.05
- Propagation: 3 steps, restart 0.15, decay 0.5

Output:

- `results/tables/gse308103_grn_virtual_perturbation_gene_effects.csv`
- `results/tables/gse308103_grn_virtual_perturbation_signature_impacts.csv`
- `results/tables/gse308103_grn_virtual_perturbation_target_ranking.csv`
- `results/tables/gse308103_grn_virtual_perturbation_network_summary.csv`
- `results/tables/gse308103_grn_virtual_perturbation_config.json`

## Results snapshot

| Target | Class | Top impacted signature | Mean impact | Top impacted genes |
|---|---|---|---:|---|
| `CD44` | macrophage | `c1q_macrophage_vs_other_macrophage` | 0.01774 | `MARCO`, `MCEMP1`, `INHBA`, `LTA4H`, `PPARG` |
| `TREM2` | macrophage | `c1q_macrophage_vs_other_macrophage` | 0.01773 | `APOE`, `APOC1`, `ACP5`, `CYP27A1`, `LGALS3` |
| `PLA2G7` | macrophage | `spp1_macrophage_vs_other_macrophage` | 0.01650 | `APOE`, `ACP5`, `APOC1`, `CHIT1`, `CYP27A1` |
| `CD74` | macrophage | `c1q_macrophage_vs_other_macrophage` | 0.01508 | `C1QA`, `C1QB`, `C1QC`, `APOE`, `ACP5` |
| `MIF` | epithelial | `proliferating_epithelial_vs_other_epithelial` | 0.01247 | `MDK`, `EPCAM`, `KRT8`, `CD24`, `PERP` |
| `CXCR4` | macrophage | `spp1_macrophage_vs_other_macrophage` | 0.01210 | `PRDM1`, `GPR183`, `ADAMDEC1`, `PLA2G7`, `APOE` |
| `SPP1` | macrophage | `spp1_macrophage_vs_other_macrophage` | 0.01022 | `MIF`, `ANXA2`, `TREM2`, `G0S2`, `CAMK1` |

Interpretation:

- The GRN result is not identical to the spatial score-dropout result, which is useful because it contributes an orthogonal network view.
- `CD74` remains a credible receptor-side candidate, now linked to a C1Q/APOE/ACP5 macrophage network.
- `CXCR4` and `PLA2G7` appear more connected to SPP1-like macrophage-state signatures.
- `MIF` connects to an epithelial proliferative/progenitor-like expression neighborhood, but this remains a network-prioritization signal rather than causal evidence.

## Boundary language

Use:

- GRN-level virtual perturbation prioritization
- scTenifoldKnk-inspired outgoing-edge perturbation
- candidate network-impact ranking

Avoid:

- causal knockout validation
- virtual knockout proves
- experimentally validated
- CellOracle/scTenifoldKnk full workflow, unless the full package workflow is actually run

## Sources consulted

- scTenifoldKnk paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9058914/
- scTenifoldKnk code/package page: https://github.com/cailab-tamu/scTenifoldKnk
- scTenifoldKnk CRAN/METACRAN package summary: https://www.r-pkg.org/pkg/scTenifoldKnk
- CellOracle paper: https://www.nature.com/articles/s41586-022-05688-9
- CellOracle documentation, base-GRN options: https://morris-lab.github.io/CellOracle.documentation/tutorials/base_grn.html
- CellOracle documentation, perturbation simulation: https://morris-lab.github.io/CellOracle.documentation/tutorials/simulation.html

## 2026-06-10 original scTenifoldKnk smoke test

The original CRAN `scTenifoldKnk` package was installed and tested locally:

- Rscript: `C:\Program Files\R\R-4.6.0\bin\Rscript.exe`
- R: `R version 4.6.0 (2026-04-24 ucrt)`
- Package: `scTenifoldKnk 1.0.3`
- Reproducible script: `scripts/run_original_sctenifoldknk_smoke.R`

Command:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' --vanilla scripts\run_original_sctenifoldknk_smoke.R
```

Inputs were reduced GSE308103 matrices for epithelial `MIF` and macrophage
`CD74`, each with 90 cells and 116 genes.

Outputs:

- `results/tables/sctenifoldknk_smoke/gse308103_epithelial_mif_sctenifoldknk_diffregulation.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_epithelial_mif_sctenifoldknk_top_non_target_genes.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_macrophage_cd74_sctenifoldknk_diffregulation.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_macrophage_cd74_sctenifoldknk_top_non_target_genes.csv`
- `results/tables/sctenifoldknk_smoke/sctenifoldknk_original_smoke_summary.csv`

Result snapshot after excluding the perturbed target itself:

- `MIF` epithelial: top non-target genes included `KRT8`, `ANXA2`, `PERP`,
  `CD24`, `LGALS3`, `CLDN4`, `CLDN3`, and `EPCAM`; best non-target adjusted
  P value was 0.305553.
- `CD74` macrophage: top non-target genes included `CYP27A1`, `TREM2`,
  `APOC1`, `PLA2G7`, `LPL`, `LRP1`, `LGALS3`, and `APOE`; best non-target
  adjusted P value was 0.191556.

Interpretation: the original package is usable in the local project
environment and gives network-context signals that are directionally compatible
with the current epithelial progenitor-like and APOE/TREM2/PLA2G7 macrophage
axes. Because this was a reduced smoke test and non-target FDR values were not
significant, it should be presented only as optional sensitivity analysis /
GRN-level virtual perturbation prioritization, not causal knockout validation.
