# Figure Legends Draft

Working title:

**Specificity-audited spatial transcriptomic mining nominates a MIF-CD74 epithelial-myeloid communication candidate during early LUAD progression**

## Figure 1. Public multi-cohort framework for early LUAD epithelial-myeloid niche mining

**a**, Study workflow. Public early lung adenocarcinoma (LUAD) datasets were organized into discovery, ordered spatial progression, ordered single-nucleus progression, bulk trend, specificity-audit, and tumor-adjacent support layers. The workflow proceeds through six analytical steps: cohort harmonization, specificity audit, spatial candidate-axis scoring, score-level in-silico target prioritization, GRN-level virtual perturbation prioritization, and bounded `MIF-CD74` candidate nomination. **b**, Stage composition of the seven public cohorts used in the current analysis. `GSE307534` and `GSE308103` provide the main ordered Normal-AAH-AIS-MIA-LUAD progression evidence, while `GSE131907` provides a large cell-type specificity reference. **c**, Evidence role assigned to each dataset. Sample counts reflect harmonized metadata; for `GSE164789`, 62 expression matrices were analyzed from the tumor-adjacent expression files.

## Figure 2. Integrated evidence prioritizes MIF-CD74/CXCR4 over broader macrophage-state axes

**a**, Integrated priority ranking of candidate epithelial-myeloid axes. The `MIF-CD74/CXCR4` axis ranked first, followed by the `SPP1/TREM2/PLA2G7` macrophage-state readout axis. **b**, Normalized evidence-component heatmap summarizing spatial, specificity, bulk, single-nucleus/single-cell, and score-level target-prioritization support for each candidate axis. **c**, Continuous score-reduction analysis in MIA/LUAD spatial samples from `GSE307534` (`n=30` samples). Full score dropout of `MIF` produced the strongest source-side coupling loss, and full score dropout of `CD74` produced the strongest receptor-side coupling loss. These effects are interpreted as score-level in-silico target prioritization rather than causal evidence.

## Figure 3. Specificity auditing reframes the original SPP1 macrophage niche hypothesis

**a**, Relative cell-type scores of selected discovery and candidate signatures in the large `GSE131907` LUAD single-cell reference. The broad macrophage programs are myeloid-enriched, but the original refined `SPP1` macrophage-like signature shows substantial epithelial and stromal off-target signal. **b**, Specificity status of the original refined `SPP1` macrophage marker set. Of 30 audited genes, 14 matched the expected myeloid direction, 15 were off-target, and 1 was absent from the reference matrix. **c**, Top-expressing reference cell type for selected candidate genes. `MIF` is epithelial-enriched, whereas `CD74`, `SPP1`, `TREM2`, `PLA2G7`, `IL1B`, `TNF`, and `CXCL8` are assigned to myeloid or immune compartments. This audit supports separating the `MIF-CD74` communication hypothesis from the `SPP1/TREM2/PLA2G7` macrophage-state readout.

## Figure 4. Spatial progression analysis supports source-side MIF epithelial-neighborhood enrichment

**a**, Stage-wise observed-minus-null spatial enrichment of candidate source and receptor programs near epithelial progenitor-like spots in `GSE307534`. **b**, Focused progression trends for the `MIF-CD74/CXCR4` and `SPP1/TREM2/PLA2G7` axes. `MIF-CD74/CXCR4` source-side enrichment increases into MIA/LUAD, while receptor-side and macrophage-state signals are interpreted as supporting or boundary-setting readouts rather than monotonic spatial progression markers. **c**, Mean MIA/LUAD enrichment ranking across the strongest spatial-axis rows. In a paired sensitivity analysis of 20 patients with precursor and late lesions, the `MIF` source-side late-minus-precursor enrichment delta was 0.141 (95% bootstrap CI 0.094 to 0.188; BH-adjusted Wilcoxon q=0.00019). Spatial analyses are Visium spot-level neighborhood enrichments and should not be interpreted as single-cell physical contact.

## Figure 5. Score-level in-silico target prioritization ranks receptor-side CD74

**a**, Schematic of the score-level target-prioritization strategy. Individual candidate genes were down-weighted or dropped from source or receptor axis scores, and epithelial-neighborhood coupling was recalculated. **b**, Dose-response-like coupling retention after 50% down-weighting or full score dropout. `MIF` dropout collapsed its single-gene source score, while receptor-side loss ranked `CD74` (74.0%) above `CD44` (20.0%) and `CXCR4` (6.0%). **c**, Stage-wise full-dropout coupling loss. **d**, Concordance between continuous score-reduction and top-quantile dropout analyses. Receptor-side comparisons prioritize `CD74`, while `CD44`, `SPP1`, `TREM2`, `CXCR4`, and `PLA2G7` show weaker effects. MIF dropout is expected for a single-gene source panel. These analyses nominate follow-up validation targets but do not establish causal biology.

## Supplementary Figure 1. Expression-matched and tissue-density controls for source-side MIF enrichment

**a**, Paired late-minus-precursor enrichment deltas for 20 expression-matched single-gene controls. The horizontal line marks the `MIF` paired delta. `MIF` lies at the 95th percentile of matched controls, but `PRDX1` shows a slightly larger effect and the empirical upper-tail p-value is 0.095. **b**, Paired changes for broad macrophage-signature controls. Positive shifts in broader signatures place MIF within a wider spatial-remodeling context. **c**, Source-side MIF enrichment versus in-tissue spot count across GSE307534 sections. The weak association (Spearman r=0.122, p=0.369) argues against a simple tissue-area explanation. These controls support prioritization of MIF but do not establish unique specificity or causal signaling.

## Supplementary Figure 2. Focused orthogonal expression support in snRNA and bulk progression cohorts

**a**, Compartment-matched candidate-gene expression trends in `GSE308103` snRNA-seq. `MIF` is summarized in epithelial cells, while receptor and macrophage-state genes are summarized in macrophages; values are row z-scores of stage-level mean normalized expression. Epithelial `MIF` increased from 0.227 in Normal samples to 0.577 in LUAD. **b**, Focused bulk RNA-seq expression trends in `GSE282617` across Normal, AIS, MIA, and IAC, shown as row z-scores. **c**, Normal-to-IAC bulk expression deltas for the focused genes. Bulk `MIF`, `CXCR4`, and `SPP1` increased, while `CD44` decreased. **d**, Selected snRNA trends for `MIF`, `CD74`, `SPP1`, and benchmark `IL1B`. These panels provide non-spatial orthogonal trend support and boundary-setting; they do not demonstrate ligand-receptor contact or causal perturbation.

## Supplementary Figure 3. Specificity refinement of the original SPP1 macrophage signature

**a**, Specificity-status counts before and after auditing the original refined `SPP1` macrophage signature in the large `GSE131907` LUAD single-cell reference. The original 30-gene signature contained 14 expected myeloid genes, 15 off-target genes, and 1 missing gene; after filtering, 14 expected myeloid genes were retained. **b**, Gene-level audit across the original marker ranks. Retained genes are emphasized with black outlines, while removed genes show epithelial, fibroblast, endothelial, or missing-reference behavior. **c**, Top reference cell type of retained and removed genes. Removed genes expose mixed-lineage signal, supporting the decision to treat `SPP1/TREM2/PLA2G7` as a macrophage-state readout rather than the primary communication novelty.

## Supplementary Figure 4. GRN-level virtual perturbation stratifies target-linked expression neighborhoods

**a**, GSE308103 MIA/LUAD GRN-level virtual perturbation prioritization of candidate targets. Bars show the mean propagated impact score for the top affected epithelial or macrophage signature after outgoing-edge perturbation of each target in a positive-correlation GRN. **b**, Top five propagated gene neighborhoods for each target. `CD74` is linked to a C1Q/APOE/ACP5 macrophage neighborhood, whereas `CXCR4` and `PLA2G7` are linked to SPP1-like macrophage-state neighborhoods. **c**, Network summary for the epithelial and macrophage GRNs. This analysis is computational target prioritization and does not represent wet-lab perturbation or causal validation.

## Supplementary Table Drafts

| Table | Content | Source files |
|---|---|---|
| Supplementary Table 1 | Dataset inventory and local completeness | `docs/data_inventory.md`; `docs/dataset_completeness_summary_2026-05-30.md` |
| Supplementary Table 2 | Specificity-audited signature genes | `results/tables/gse189357_refined_signature_gse131907_specificity_audit.csv` |
| Supplementary Table 3 | Full candidate-axis spatial statistics | `results/tables/gse307534_candidate_axis_spatial_adjacency.csv`; `results/tables/gse307534_candidate_axis_spatial_by_stage.csv` |
| Supplementary Table 4 | Sample-level, patient-aggregated, and paired-patient spatial contrasts | `results/tables/gse307534_candidate_axis_late_vs_precursor_*stats.csv` |
| Supplementary Table 5 | Integrated axis evidence matrix | `results/tables/main_axis_evidence_matrix.csv` |
| Supplementary Table 6 | Score-level in-silico target-prioritization outputs | `results/tables/gse307534_virtual_perturbation_effects.csv`; `results/tables/gse307534_continuous_perturbation_effects.csv` |
| Supplementary Table 7 | SPP1 signature specificity-refinement source data | `results/tables/supplementary_figure_spp1_signature_refinement_*_source.csv` |
| Supplementary Table 8 | GSE308103 GRN-level virtual perturbation target ranking | `results/supplementary_tables/ST27_gse308103_grn_virtual_perturbation_target_ranking.csv` |
| Supplementary Table 9 | GSE308103 GRN-level virtual perturbation robustness summary | `results/supplementary_tables/ST28_gse308103_grn_virtual_perturbation_robustness_summary.csv` |
| Supplementary Table 10 | Cross-dataset validation rows for GRN-prioritized signatures | `results/supplementary_tables/ST29_grn_cross_dataset_signature_validation.csv` |
| Supplementary Table 11 | Cross-dataset validation summary for GRN-prioritized signatures | `results/supplementary_tables/ST30_grn_cross_dataset_signature_validation_summary.csv` |

Packaged manuscript-facing CSV files are indexed in `docs/supplementary_tables_index.md`.
