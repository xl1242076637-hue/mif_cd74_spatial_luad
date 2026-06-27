# Original scTenifoldKnk smoke test, 2026-06-10

## Purpose

This note records a feasibility test of the original CRAN `scTenifoldKnk`
implementation on the LUAD epithelial progenitor-macrophage niche project.
The result should be interpreted as GRN-level virtual perturbation
prioritization and network-context evidence, not causal knockout validation.

## Environment

- Rscript: `C:\Program Files\R\R-4.6.0\bin\Rscript.exe`
- R: `R version 4.6.0 (2026-04-24 ucrt)`
- Package: `scTenifoldKnk 1.0.3`

`Rscript` is available via the absolute path above, but it is not currently
assumed to be on the system `PATH`.

## Reproducible command

Run from the project root:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' --vanilla scripts\run_original_sctenifoldknk_smoke.R
```

## Inputs

The smoke test uses reduced gene-by-cell matrices derived from GSE308103:

| Class | Target | Cells | Genes | Input |
|---|---:|---:|---:|---|
| epithelial | `MIF` | 90 | 116 | `results/tables/sctenifoldknk_smoke/gse308103_epithelial_mif_counts_for_sctenifoldknk.csv` |
| macrophage | `CD74` | 90 | 116 | `results/tables/sctenifoldknk_smoke/gse308103_macrophage_cd74_counts_for_sctenifoldknk.csv` |

## Parameters

- `qc = FALSE`
- `gKO = target_gene`
- `nc_nNet = 2`
- `nc_nCells = 60`
- `nc_nComp = 2`
- `td_K = 2`
- `td_maxIter = 80`
- `ma_nDim = 2`
- `nCores = 1`
- `set.seed(7)`

These are intentionally lightweight smoke-test settings. A manuscript-level
sensitivity run should increase network repetitions and assess target stability.

## Outputs

- `results/tables/sctenifoldknk_smoke/gse308103_epithelial_mif_sctenifoldknk_diffregulation.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_epithelial_mif_sctenifoldknk_top_non_target_genes.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_macrophage_cd74_sctenifoldknk_diffregulation.csv`
- `results/tables/sctenifoldknk_smoke/gse308103_macrophage_cd74_sctenifoldknk_top_non_target_genes.csv`
- `results/tables/sctenifoldknk_smoke/sctenifoldknk_original_smoke_summary.csv`
- `scripts/run_original_sctenifoldknk_smoke.R`

## Result snapshot

After excluding the perturbed target itself from the ranking:

| Target | Class | Top non-target genes by adjusted P | Top non-target adjusted P |
|---|---|---|---:|
| `MIF` | epithelial | `KRT8`, `ANXA2`, `PERP`, `IFT57`, `CYP4B1`, `CD24`, `LGALS3`, `CLDN4`, `CLDN3`, `EPCAM` | 0.305553 |
| `CD74` | macrophage | `CYP27A1`, `TREM2`, `APOC1`, `PLA2G7`, `LPL`, `LRP1`, `LGALS3`, `APOE`, `FABP3`, `RAB42` | 0.191556 |

## Interpretation

The original package runs successfully on the project-formatted GSE308103
matrices and returns the expected `diffRegulation` fields: `gene`, `distance`,
`Z`, `FC`, `p.value`, and `p.adj`.

The smoke-test signal is directionally compatible with the current project
axis. `MIF` perturbation highlights epithelial structural/progenitor-associated
genes such as `KRT8`, `CD24`, `CLDN3`, `CLDN4`, and `EPCAM`. `CD74`
perturbation highlights macrophage lipid/phagocytic and disease-associated
genes including `TREM2`, `APOC1`, `PLA2G7`, `LRP1`, `LGALS3`, and `APOE`.

The adjusted P values for non-target genes are not significant in this reduced
smoke test. Therefore, these outputs should be used as feasibility and
network-context support only. They should not be described as causal knockout
validation or experimental perturbation evidence.

## Recommendation

Keep the existing transparent Python GRN-level perturbation as the main
reproducible analysis for now. Use the original `scTenifoldKnk` package as an
optional sensitivity analysis after expanding settings and targets. Recommended
next targets are `MIF`, `CD74`, `CD44`, `CXCR4`, `SPP1`, `TREM2`, and `PLA2G7`.
