# Manuscript Figure Plan

This file records the current figure logic so the manuscript draft remains reproducible across long analysis sessions.

## Figure 1: Workflow and Dataset Composition

Core conclusion:

The study is a multi-cohort, specificity-audited public-data framework for early LUAD epithelial-myeloid niche mining.

Figure archetype:

Schematic-led composite.

Target journal/output:

Nature-style manuscript draft, double-column width. Outputs are SVG, PDF, TIFF, and PNG.

Backend:

Python/matplotlib only.

Final size:

7.2 x 5.1 inches before tight export.

Panel map:

- a: Analysis workflow from public cohorts to specificity audit, spatial scoring, evidence ranking, and score-level perturbation.
- b: Dataset stage-composition stacked bar chart.
- c: Evidence-role table for the seven public cohorts.

Evidence hierarchy:

- Hero evidence: the workflow and cohort composition establish that the project is comprehensive rather than a single-dataset observation.
- Validation evidence: sample composition shows the ordered progression cohorts used for main spatial and snRNA support.
- Controls/robustness: specificity reference, bulk trend, and tumor-adjacent cohorts are shown as support datasets rather than primary spatial evidence.

Statistics needed:

No inferential statistics in this figure. Panel b shows metadata sample counts by interpreted stage.

Source data needed:

- `results/tables/geo_sample_metadata_annotated.csv`
- `results/tables/figure1_dataset_composition_source.csv`

Image-integrity notes:

This is vector schematic and tabular plotting only. No microscopy or spatial tissue image processing is used in this figure.

Reviewer risk:

`GSE164789` has 78 tumor/adjacent metadata rows but 62 expression matrices used for scoring. The figure labels the bar axis as metadata samples and the endpoint label as 62 expression matrices to avoid overclaiming expression-level completeness.

## Figure 2: Axis Evidence and Score-Level Target Prioritization

Core conclusion:

Integrated multi-cohort evidence and score-level in-silico target prioritization prioritize the MIF-CD74 axis over the SPP1 macrophage-state readout for follow-up validation.

Figure archetype:

Quantitative grid.

Target journal/output:

Nature-style manuscript draft, double-column width. Outputs are SVG, PDF, TIFF, and PNG.

Backend:

Python/matplotlib only.

Final size:

7.2 x 5.6 inches before tight export.

Panel map:

- a: Integrated priority ranking of candidate epithelial-myeloid axes.
- b: Normalized evidence-component heatmap.
- c: Continuous score-level target-prioritization effect sizes in MIA/LUAD spatial samples.

Evidence hierarchy:

- Hero evidence: MIF-CD74/CXCR4 has the highest integrated priority score.
- Validation evidence: component heatmap separates spatial, specificity, bulk, snRNA/scRNA, and perturbation evidence.
- Controls/robustness: SPP1/TREM2/PLA2G7 is retained as a macrophage-state readout; IL1B/TNF/CXCL8 remains a benchmark inflammatory axis.

Statistics needed:

Panel c reports `n=30` MIA/LUAD spatial samples from GSE307534 for the displayed perturbation summaries. Current figure draft does not display p-values because the perturbation is a score-level target-prioritization analysis rather than causal inference.

Source data needed:

- `results/tables/main_axis_evidence_matrix.csv`
- `results/tables/gse307534_continuous_perturbation_mia_luad_ranking.csv`
- `results/tables/figure2_priority_source.csv`
- `results/tables/figure2_evidence_heatmap_source.csv`
- `results/tables/figure2_perturbation_source.csv`

Display-score normalization:

- Spatial niche: mean of source/target spatial deltas divided by 0.20, clipped to 0-1.
- Specificity audit: mean of source/target specificity fractions.
- Bulk trend: positive bulk delta divided by 40, clipped to 0-1.
- snRNA program: positive snRNA delta divided by 0.30, clipped to 0-1.
- Tumor-adjacent scRNA: positive tumor-adjacent scRNA delta divided by 0.10, clipped to 0-1.
- Target prioritization: top continuous perturbation priority score.

Image-integrity notes:

This is generated from quantitative tables only. No image manipulation is involved.

Reviewer risk:

The perturbation panel must be described as score-level in-silico target prioritization, not as causal validation. MIF and CD74 are prioritized for downstream experimental validation; SPP1 remains the macrophage-state readout.

## Figure 3: Specificity Audit

Core conclusion:

GSE131907 specificity auditing shows that the original SPP1 macrophage niche hypothesis needs reframing: broad SPP1 macrophage signal is myeloid-enriched, but the refined SPP1 marker set contains substantial epithelial/stromal off-target signal, whereas the MIF-CD74 candidate genes show cleaner epithelial-myeloid directionality.

Figure archetype:

Quantitative grid.

Target journal/output:

Nature-style manuscript draft, double-column width. Outputs are SVG, PDF, TIFF, and PNG.

Backend:

Python/matplotlib only.

Final size:

7.2 x 5.9 inches before tight export.

Panel map:

- a: GSE131907 relative cell-type score heatmap for selected discovery/refined signatures.
- b: Specificity audit fractions for refined GSE189357 marker signatures.
- c: Top-expressing GSE131907 cell type for candidate genes.

Evidence hierarchy:

- Hero evidence: refined SPP1 macrophage has high epithelial/fibroblast relative signal in the large scRNA reference.
- Validation evidence: only 14/30 original SPP1 macrophage markers pass the expected myeloid top-cell-type audit; 15/30 are off-target and 1/30 is missing.
- Controls/robustness: MIF is top-expressed in epithelial cells, while CD74/SPP1/TREM2/PLA2G7/IL1B/TNF/CXCL8 are top-expressed in myeloid or immune cells.

Statistics needed:

No inferential statistics in this figure. Panels show reference-cell mean signature scores and marker top-cell-type audit counts.

Source data needed:

- `results/tables/gse131907_selected_signature_celltype_summary.csv`
- `results/tables/gse189357_refined_signature_gse131907_specificity_summary.csv`
- `results/tables/gse131907_selected_gene_top_celltype.csv`
- `results/tables/figure3_signature_celltype_heatmap_source.csv`
- `results/tables/figure3_specificity_status_source.csv`
- `results/tables/figure3_candidate_gene_specificity_source.csv`

Image-integrity notes:

This is generated from quantitative tables only. No image manipulation is involved.

Reviewer risk:

This figure should be described as a specificity audit and hypothesis refinement, not as proof that SPP1 biology is irrelevant. SPP1 remains the macrophage-state readout.

## Figure 4: Spatial Axis Progression

Core conclusion:

In GSE307534 spatial transcriptomics, MIF-CD74/CXCR4 shows increasing epithelial-neighborhood spatial enrichment across early LUAD progression and ranks among the strongest MIA/LUAD spatial axes.

Figure archetype:

Quantitative grid.

Target journal/output:

Nature-style manuscript draft, double-column width. Outputs are SVG, PDF, TIFF, and PNG.

Backend:

Python/matplotlib only.

Final size:

7.2 x 6.0 inches before tight export.

Panel map:

- a: Stage-wise heatmap of observed-minus-null adjacency enrichment for candidate axes.
- b: Focused stage trend for MIF-CD74 and SPP1/TREM2/PLA2G7 source/receptor sides.
- c: Mean MIA/LUAD enrichment ranking for the strongest spatial-axis rows.

Evidence hierarchy:

- Hero evidence: MIF-CD74 source-side enrichment rises from normal/AAH/AIS into MIA/LUAD and is the top MIA/LUAD mean axis row.
- Validation evidence: patient-aware sensitivity analysis supports the source-side MIF progression signal; receptor-side MIF-CD74/CXCR4 remains a target-prioritization signal rather than a robust paired stage-increase claim.
- Controls/robustness: SPP1/TREM2/PLA2G7 receptor-side enrichment is strong, supporting its use as a macrophage-state readout.

Statistics needed:

Panel c uses MIA/LUAD stage summaries across 30 spatial samples from GSE307534. A separate paired sensitivity analysis includes 20 patients with precursor and late lesions: MIF source-side mean late-minus-precursor delta 0.141, 95% bootstrap CI 0.094 to 0.188, BH-adjusted Wilcoxon q=0.00019. Empirical permutation p-value medians are available in the source spatial table but are not yet displayed on the figure.

Source data needed:

- `results/tables/gse307534_candidate_axis_spatial_by_stage.csv`
- `results/tables/figure4_axis_stage_heatmap_source.csv`
- `results/tables/figure4_focus_axis_trend_source.csv`
- `results/tables/figure4_late_axis_summary_source.csv`

Image-integrity notes:

This is generated from quantitative spatial-neighborhood summary tables only. No histology image processing is used in this figure.

Reviewer risk:

The spatial signal is Visium spot-level and should be described as spatial coupling/enrichment, not direct physical cell-cell contact at single-cell resolution.

## Figure 5: Score-Level In-Silico Target Prioritization

Core conclusion:

Continuous score-level in-silico target prioritization ranks MIF and CD74 as the strongest follow-up targets within the MIF-CD74 epithelial-myeloid axis. Receptor-side dropout effects rank CD74 above CD44 and CXCR4, while SPP1, TREM2, and PLA2G7 show weaker macrophage-state coupling-loss effects.

Figure archetype:

Asymmetric mixed-modality figure: schematic plus quantitative perturbation panels.

Target journal/output:

Nature-style manuscript draft, double-column width. Outputs are SVG, PDF, TIFF, and PNG.

Backend:

Python/matplotlib only.

Final size:

7.2 x 5.8 inches before tight export.

Panel map:

- a: Score-level in-silico target-prioritization model.
- b: Dose-response-like coupling retention for x0.5 down-weighting and x0 full score dropout.
- c: Full-dropout relative coupling loss by stage.
- d: Concordance between continuous coupling perturbation and top-quantile dropout perturbation.

Evidence hierarchy:

- Hero evidence: MIF full score dropout collapses source-side continuous coupling, and CD74 receptor-side dropout produces the next strongest effect.
- Validation evidence: MIF and CD74 remain strong across stages and across two perturbation-scoring approaches.
- Controls/robustness: CD44, CXCR4, SPP1, TREM2, and PLA2G7 are shown as weaker comparators.

Statistics needed:

Panels b-d use GSE307534 spatial perturbation summaries. Panel c includes the stage-specific sample counts inherited from the source table; MIA/LUAD summaries include 30 samples. No p-values are shown because this figure is a target-prioritization analysis rather than a causal experiment.

Source data needed:

- `results/tables/gse307534_continuous_perturbation_mia_luad_ranking.csv`
- `results/tables/gse307534_continuous_perturbation_effects.csv`
- `results/tables/gse307534_virtual_perturbation_mia_luad_ranking.csv`
- `results/tables/figure5_dose_response_source.csv`
- `results/tables/figure5_stage_loss_source.csv`
- `results/tables/figure5_method_concordance_source.csv`

Image-integrity notes:

This is generated from quantitative spatial-perturbation tables only. No microscopy or histology image processing is used.

Reviewer risk:

Avoid causal wording. Use "score-level dropout", "score-level perturbation", or "score-level in-silico target prioritization". This figure nominates MIF/CD74 for follow-up wet-lab validation but does not establish causal sufficiency.

Additional interpretation note:

- Full MIF dropout is expected because the source-side selected panel contains only `MIF`.
- Treat MIF dropout as model dependency, not independent causal support.
- The receptor-side `CD74 > CD44 > CXCR4` comparison is more discriminative.

## Supplementary Figure 1: MIF Spatial Controls

Core conclusion:

Source-side MIF enrichment is not explained by section spot count alone and lies near the top of an expression-matched single-gene control distribution, but it is not uniquely specific.

Panel map:

- a: Paired late-minus-precursor deltas for MIF and 20 expression-matched single-gene controls.
- b: Broad macrophage-signature paired controls.
- c: MIF enrichment versus in-tissue spot count.

Reviewer risk:

The matched-control empirical p-value is 0.095 and `PRDX1` is slightly above MIF. State this transparently. The control supports prioritization, not proof that MIF is uniquely exceptional.

## Supplementary Figure 2: Focused Orthogonal Validation

Core conclusion:

Focused non-spatial cohorts provide orthogonal expression-trend support for selected MIF-axis and macrophage-state genes, but these panels do not replace the spatial-coupling evidence.

Figure archetype:

Quantitative grid.

Target journal/output:

Nature-style manuscript draft, double-column width. Outputs are SVG, PDF, TIFF, and PNG.

Backend:

Python/matplotlib only.

Final size:

7.2 x 6.7 inches before tight export.

Panel map:

- a: GSE308103 snRNA compartment-matched candidate-gene expression heatmap. `MIF` is summarized in epithelial cells; receptor and macrophage-state genes are summarized in macrophages.
- b: GSE282617 bulk candidate-gene expression heatmap across Normal, AIS, MIA, and IAC.
- c: GSE282617 Normal-to-IAC bulk deltas for the focused genes.
- d: Selected GSE308103 snRNA trend lines for `MIF`, `CD74`, `SPP1`, and benchmark `IL1B`.

Evidence hierarchy:

- Hero evidence: non-spatial cohorts are used as orthogonal trend checks for the genes prioritized by spatial analysis.
- Validation evidence: bulk trends show whether candidate genes rise, stay flat, or decrease with progression.
- Controls/robustness: inflammatory benchmark genes remain visible so that MIF-CD74 is not conflated with a generic IL1B/TNF/CXCL8 niche.

Statistics needed:

Panel a summarizes 75 GSE308103 snRNA samples at the sample-averaged stage level. Panel b-c summarize 70 GSE282617 bulk samples at the stage-mean level. No p-values are displayed because these panels are trend-support visualizations.

Source data needed:

- `results/tables/gse308103_snrna_candidate_gene_sample_summary.csv`
- `results/tables/gse308103_snrna_candidate_gene_stage_summary.csv`
- `results/tables/gse282617_candidate_marker_stage_means.csv`
- `results/tables/supplementary_figure_focused_orthogonal_validation_snrna_source.csv`
- `results/tables/supplementary_figure_focused_orthogonal_validation_bulk_source.csv`

Image-integrity notes:

This is generated from quantitative expression-summary tables only. No image manipulation is involved.

Reviewer risk:

These are non-spatial expression trends. They should be described as orthogonal support and boundary-setting, not as evidence of cell-cell communication or spatial niche formation.
