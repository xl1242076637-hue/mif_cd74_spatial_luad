# Decision Log

## 2026-05-29: Split the Work Into Pure Bioinformatics First, Wet Validation Later

The first project will be pure bioinformatics. Wet-lab validation will be treated as a second project after the computational analysis nominates a small set of candidate niches and perturbation targets.

Rationale:

- The public early LUAD datasets are sufficient to build a complete computational story.
- Tissue spatial validation and cell/animal functional validation require different resources, timelines, and collaborators.
- Separating the projects prevents the first manuscript from becoming too broad.

## 2026-05-29: Prioritize Early LUAD Progression Over Immunotherapy Response

The first computational project will focus on early LUAD progression and epithelial progenitor-macrophage spatial niches.

Rationale:

- The early LUAD progression direction has public single-cell, spatial, and bulk validation datasets.
- The disease stages are biologically coherent: normal, AAH, AIS, MIA, IAC/LUAD.
- The story naturally connects to later wet validation through epithelial-macrophage colocalization and ligand-receptor axes.
- Immunotherapy-response datasets are valuable but often have fewer spatial samples and more raw-data processing burden.

## 2026-05-29: Use Python-First Reproducibility

The current environment has Python but no `Rscript`, `scanpy`, or `anndata`. The project will start with Python utilities for provenance, downloading, inventory, and tabular processing. Later single-cell analysis can add `scanpy`/`anndata` or an R/Seurat environment when needed.

## 2026-05-29: Use Raw GSE189487 as the First Spatial Cohort

The full `GSE189487_RAW.tar` is now available locally and contains all six AIS/MIA/IAC Visium samples. The first spatial pipeline will therefore use these raw 10x MatrixMarket files instead of relying on the one-sample processed STOMICS h5ad fallback.

Rationale:

- Raw matrices let the same scoring and adjacency logic run across all six samples.
- This avoids drawing biological conclusions from the single processed IAC proof-of-code sample.
- The current environment lacks R/Seurat and Scanpy/AnnData, so a lightweight Python reader based on `gzip`, `pandas`, and `scipy.io.mmread` is the most reproducible path for now.

## 2026-05-29: Normalize Raw Spatial Scores by Full Spot Library Size

Candidate marker counts from raw 10x data will be normalized by all-feature per-spot total counts before `log1p` transformation and panel scoring.

Rationale:

- Normalizing by only the selected candidate marker subset would make scores depend on marker-panel composition.
- Full-matrix spot totals better represent standard spatial transcriptomics library-size normalization.

## 2026-05-29: Use One-Hop Spatial Radius for Main Adjacency Test

The main epithelial-high to macrophage-high adjacency test uses `radius_multiplier=1.0`, equivalent to approximately one median nearest-neighbor spot spacing.

Rationale:

- A sensitivity run with `radius_multiplier=2.0` made the permuted null nearly saturated (`~0.97-0.98`), which is too broad for an immediate spatial-niche question.
- A one-hop radius better matches the idea of local epithelial-macrophage proximity.
- The broad marker panels did not show positive adjacency enrichment in the six raw `GSE189487` samples. The next computational step should refine epithelial progenitor and macrophage subtype signatures using the matched `GSE189357` scRNA-seq dataset.

## 2026-05-29: Treat GSE189357 Panel Scoring as a Sanity Check, Not Annotation

The first scRNA step scores the same broad candidate panels in raw `GSE189357` 10x matrices, but this is only a sanity check. It should not be treated as cell-state annotation.

Rationale:

- The broad panels mix cell identity, tumor epithelial state, proliferation, macrophage programs, and ligand-receptor axes.
- Sample-level heterogeneity is strong, so stage means alone are not enough to support a mechanistic claim.
- The next biologically meaningful step is to define epithelial-like and macrophage-like cells, derive more specific state signatures, and then test those refined signatures in spatial data.

Decision:

- Do not force the project story around the current broad panel results.
- Use `GSE189357` to build refined epithelial progenitor and macrophage subtype signatures before repeating spatial colocalization and ligand-receptor analyses.

## 2026-05-29: Add Mast/Neutrophil/Dendritic Broad Classes Before Trusting Macrophage Subtypes

Initial marker extraction from the first macrophage subtype labels showed mast-cell genes among the top `spp1_macrophage` markers. This indicated marker-score contamination rather than a valid macrophage subtype result.

Decision:

- Add `mast_cell`, `neutrophil`, and `dendritic` marker-score classes to the broad classifier.
- Rerun scRNA state assignment and marker extraction after adding these contaminant-removal classes.
- Treat the corrected `SPP1 macrophage-like`, `C1Q macrophage-like`, and `inflammatory macrophage-like` signatures as usable for spatial remapping.
- Treat the current `resident_macrophage` signature as low confidence because its extracted top genes still show neutrophil/monocyte-like contamination.

## 2026-05-29: Prioritize the MIA Epithelial Progenitor-SPP1 Macrophage Niche

After remapping refined scRNA-derived signatures to `GSE189487` spatial samples, the strongest positive signal is in MIA: epithelial progenitor-like spots show positive adjacency to `SPP1 macrophage-like` signature spots.

Decision:

- Prioritize `epithelial_progenitor_like` to `spp1_macrophage` spatial adjacency as the main candidate niche for the next phase.
- Keep `inflammatory_macrophage` as a secondary exploratory axis.
- Do not lead with generic macrophage or C1Q/resident macrophage adjacency, because those are not positively enriched in the current spatial tests.

## 2026-05-29: Keep the SPP1 Niche as Main Candidate, but Treat It as Parameter-Sensitive

The robustness grid supports the MIA epithelial progenitor-like to `SPP1 macrophage-like` niche as the main candidate, but not as a fully parameter-invariant result.

Evidence:

- AIS and IAC show no positive parameter sets across the tested grid.
- MIA has 6/27 parameter sets with positive enrichment and median `p_greater < 0.05`.
- The positive/significant MIA signal concentrates in top30 signatures and local radii (`0.75` or `1.0` times median nearest-neighbor distance).
- Both MIA spatial samples are positive in the significant top30 settings.

Decision:

- Continue with this niche as the primary computational story.
- Use top30 refined signatures for downstream niche-focused analyses, while reporting top-N sensitivity transparently.
- Frame the result as a refined, MIA-enriched local spatial niche rather than a generic or universally robust epithelial-macrophage interaction.
- Next priority: nominate ligand-receptor axes and perturbation candidates specifically within the MIA epithelial progenitor-like and `SPP1 macrophage-like` context.

## 2026-05-29: Expand From a Minimal Feasible Project to a Comprehensive Multi-Cohort Study

The project should not be optimized for the smallest possible download. The stronger design is a comprehensive, multi-modal public-data study of early LUAD precursor progression.

Decision:

- Keep `GSE189357`/`GSE189487` as the current discovery pair because they already nominate the MIA epithelial progenitor-like to `SPP1 macrophage-like` niche.
- Add `GSE308103` snRNA-seq as the primary external single-nucleus validation cohort for epithelial progenitor and myeloid state definitions.
- Add `GSE307534` Visium as the primary external spatial validation cohort for the full Normal/AAH/AIS/MIA/LUAD sequence.
- Add `GSE164789` as a precursor-LUAD scRNA/scTCR immune-state validation cohort.
- Add `GSE131907` as a broad LUAD single-cell reference to test macrophage/TME specificity beyond precursor-only data.
- Use `GSE282617` bulk RNA-seq as a compact expression-trend validation cohort.

Scientific framing:

- The main story should be broader than one `SPP1` signal: alveolar/epithelial progenitor plasticity and proinflammatory or immunoregulatory myeloid niches may co-evolve during AAH/AIS/MIA/LUAD progression.
- `SPP1`, `IL1B`, `TREM2`, `CXCL9/CXCL10`, and `MIF-CD74/CXCR4` axes should be treated as candidate mechanisms to rank, not assumptions.
- Score-level in-silico target prioritization should be used after cross-cohort and spatial evidence to prioritize mechanisms for later wet-lab validation, not as standalone proof.

## 2026-05-30: Reframe GSE307534 as External Spatial Recurrence, Not MIA-Specific Proof

After correcting duplicate-lesion accession mapping and excluding invalid zero-high-spot adjacency tests, `GSE307534` shows recurrent epithelial progenitor-like to macrophage-like spatial coupling across the larger Normal/AAH/AIS/MIA/LUAD Visium cohort.

Decision:

- Keep `epithelial_progenitor_like` to `spp1_macrophage` as the lead candidate niche because it remains strongest in MIA/LUAD and was nominated independently by `GSE189487`.
- Do not claim that the niche is MIA-exclusive based on `GSE307534`.
- Frame the external result as progression-associated spatial recurrence or strengthening, with AIS already showing signal and LUAD remaining strong.
- Add specificity controls before mechanism claims: broad epithelial/macrophage controls, random gene-set controls, total-count/tissue-density controls, and patient-paired lesion summaries where possible.

## 2026-05-30: Use GSE164789 as Tumor/Adjacent Immune-State Validation

`GSE164789` has 62 expression matrices and 16 TCR contig files. The first expression scoring shows tumor-side enrichment of SPP1 and inflammatory macrophage programs, while epithelial progenitor-like fraction is not higher in tumor than adjacent under the current marker-score classification.

Decision:

- Use `GSE164789` mainly to validate macrophage-state direction and tumor/adjacent immune remodeling.
- Do not use it as primary evidence for stage-specific AAH/AIS/MIA progression because its harmonized labels are Adjacent/Tumor rather than ordered precursor stages.
- Treat its epithelial result as a caution: the project should separate epithelial progenitor/plasticity, proliferation, and tumor-cell abundance instead of merging them into one undifferentiated "progenitor" signal.

## 2026-05-30: Use GSE131907 for Specificity Rather Than Discovery

`GSE131907` is a large LUAD single-cell reference with broad lung cancer, normal lung, metastasis/effusion, and lymph-node contexts. The raw UMI text matrix and annotation are now local.

Decision:

- Use `GSE131907` to ask whether candidate macrophage programs are macrophage/TME-specific and whether epithelial progenitor-like signatures concentrate in epithelial/tumor compartments.
- Do not treat it as an ordered early-progression validation cohort.
- Analyze it with aggregated selected-gene signature scoring by annotated cell type and sample origin, rather than trying to rerun the full per-cell state pipeline immediately.

## 2026-05-30: Downgrade the Original Refined SPP1 Niche Claim After Specificity Audit

GSE131907 specificity scoring showed that the original `spp1_macrophage_vs_other_macrophage` refined signature from GSE189357 contains substantial epithelial/stromal contamination. The broad SPP1 macrophage marker panel is myeloid-enriched, but the refined 30-gene SPP1 macrophage signature is highest in epithelial cells in GSE131907.

Decision:

- Do not use the unfiltered refined SPP1 macrophage panel as a primary macrophage-specific spatial niche marker.
- Keep SPP1/TREM2/PLA2G7/CHIT1-type myeloid genes as a candidate macrophage program, but separate them from epithelial genes such as `WFDC2`, `OCIAD2`, `CST6`, and `RAMP1`.
- Treat the original GSE189487 MIA SPP1 result as hypothesis-generating only.
- Lead with the larger and more conservative pattern: epithelial plasticity is coupled to myeloid inflammatory/SPP1 programs across progression, with MIA/LUAD enrichment strongest in GSE307534 after specificity filtering.
- Before score-level target prioritization, build mechanism candidates from specificity-filtered myeloid and epithelial programs rather than from the original contaminated refined panel.

## 2026-05-30: Split the Lead Story Into Mechanism Axis and Macrophage-State Axis

After adding direct GSE307534 candidate-axis spatial adjacency and rerunning the multi-cohort mechanism ranking, `mif_cd74_cxcr4` ranks slightly above `spp1_trem2_macrophage_epithelial`.

Decision:

- Treat `MIF-CD74/CXCR4` as the top ligand-receptor mechanism hypothesis.
- Treat `SPP1/TREM2/PLA2G7` as the top macrophage-state and perturbation-readout hypothesis.
- Keep the project title/framing centered on early LUAD epithelial progenitor-like and myeloid spatial niche formation, rather than on a single `SPP1` marker.

Rationale:

- `MIF` is epithelial-enriched in GSE131907, while `CD74`/`CD44` are myeloid-enriched and the MIF-axis source/target panels show strong GSE307534 spatial adjacency to epithelial progenitor-like spots in MIA/LUAD.
- `SPP1/TREM2/PLA2G7` retains strong macrophage specificity and cross-cohort support, but its receptor-side epithelial specificity is weaker under GSE131907 top-celltype auditing.
- The two axes are complementary: `MIF-CD74/CXCR4` is a plausible upstream communication axis, whereas `SPP1/TREM2/PLA2G7` captures the macrophage phenotype that may be remodeled.

Next priority:

- Build a score-level in-silico target-prioritization layer that reports how much the epithelial-progenitor/myeloid spatial niche score changes after computationally dropping or down-weighting `MIF`, `CD74`, `CXCR4`, `SPP1`, `TREM2`, and `PLA2G7`.
- Keep perturbation language conservative: these analyses nominate validation targets; they do not prove causal biology without later cell or animal experiments.

## 2026-05-30: Novelty Audit Finds High Overlap With Cancer Cell 2026

A targeted literature check found a close, high-impact overlap:

Peng F, Sinjab A, Dai Y, et al. **Multimodal spatial-omics reveal co-evolution of alveolar progenitors and proinflammatory niches in progression of lung precursor lesions.** *Cancer Cell*. 2026;44(2):321-339.e13. DOI: `10.1016/j.ccell.2025.10.004`.

Decision:

- Do not pursue the project as a generic "alveolar/epithelial progenitor-proinflammatory macrophage spatial niche in early LUAD" story.
- Treat `IL1B-high macrophage/proinflammatory niche` biology as an already-published benchmark rather than our main novelty.
- Pivot to a narrower, non-duplicate computational question: specificity-audited public-data mining and score-level in-silico target prioritization of `MIF-CD74/CXCR4`, with `SPP1/TREM2/PLA2G7` as a macrophage-state readout.

Revised working title:

**Specificity-audited spatial transcriptomic mining nominates a MIF-CD74/CXCR4 epithelial-myeloid communication axis during early LUAD progression**

Full audit note:

- `docs/novelty_audit_2026-05-30.md`

## 2026-05-30: Score-Level In-Silico Target Prioritization Ranks MIF and CD74 for Follow-Up

The first GSE307534 score-level in-silico target-prioritization analysis supports the revised mechanism split.

Decision:

- Prioritize `MIF` as the epithelial/source-side computational perturbation candidate.
- Prioritize `CD74` over `CXCR4` as the receptor-side myeloid/target candidate in the current spatial data.
- Keep `SPP1` as the macrophage-state readout candidate, with `TREM2` and `PLA2G7` as secondary macrophage-state genes.

Rationale:

- `MIF` score dropout collapses the `mif_cd74_cxcr4` source score and strongly reduces source-near-epithelial observed adjacency in MIA/LUAD.
- `CD74` explains most of the receptor-side `MIF-CD74/CXCR4` panel score; `CXCR4` has a much smaller contribution under the current selected-gene spatial scoring.
- `SPP1` has the largest single-gene effect within the `SPP1/TREM2/PLA2G7` macrophage-state perturbation set.

Caveat:

- Call this "score-level in-silico target prioritization", not causal validation.
- Wet-lab follow-up would still require actual MIF knockdown/blockade or CD74 perturbation in epithelial-myeloid co-culture or animal models.

## 2026-05-30: Continuous Coupling Confirms the MIF-CD74 Priority

The continuous spatial-coupling perturbation analysis confirms the top-quantile adjacency perturbation result and is better suited for knockdown-like dose response.

Decision:

- Use both perturbation readouts in the final computational story:
  - top-quantile adjacency perturbation for intuitive spatial niche loss;
  - continuous coupling perturbation for dose-response-like sensitivity.
- In figures/tables, lead with `MIF` and `CD74`.
- Treat `CXCR4` as a supporting receptor gene rather than the main receptor-side perturbation target.

Rationale:

- `MIF` dropout gives a -1.000 relative continuous coupling change for the source side, and `MIF` 50% down-weighting gives -0.500.
- `CD74` dropout gives a -0.740 relative continuous coupling change on the receptor side, much stronger than `CXCR4` (-0.060).
- `SPP1` remains the strongest macrophage-state readout gene within the `SPP1/TREM2/PLA2G7` panel, but its continuous coupling effect is smaller than the MIF-CD74 axis.

## 2026-05-30: Proceed With a Narrow MIF-CD74 Draft After Literature Recheck

A follow-up literature search found related work on early-LUAD inflammatory spatial niches, `MIF-CD74` in LUAD biology, and `SPP1+` macrophage spatial programs.

Decision:

- Continue writing the initial manuscript draft, but keep the novelty claim narrow.
- Do not frame the study as discovering epithelial-inflammatory niches in early LUAD.
- Do not frame the study as discovering `MIF-CD74` de novo in lung cancer.
- Frame the contribution as a reproducible, public-data, specificity-audited prioritization workflow that nominates `MIF-CD74` as a leading epithelial-myeloid communication candidate in ordered early-LUAD spatial data.
- Use `SPP1/TREM2/PLA2G7` as macrophage-state readouts and `IL1B/TNF/CXCL8` as a benchmark inflammatory axis.

Rationale:

- The broad progenitor-inflammatory niche is already strongly covered by recent multimodal spatial-omics work.
- Generic `SPP1+` macrophage spatial crosstalk is crowded in LUAD/NSCLC.
- The current analysis still has a defendable difference through cross-cohort specificity auditing, candidate-axis ranking, and score-level in-silico target prioritization of `MIF` and `CD74`.

## 2026-05-30: Patient-Aware Statistics Strengthen the Source-Side MIF Result

Added sample-level, patient-aggregated, and paired-patient sensitivity analyses for GSE307534 spatial-axis enrichment.

Decision:

- Lead the spatial progression result with source-side `MIF` enrichment.
- Describe receptor-side `MIF-CD74/CXCR4` as a target-prioritization signal, not as a proven paired stage-increasing program.
- Keep lesion-aware and paired-patient summaries in the supplementary tables.

Rationale:

- Among 20 patients with both precursor (`AAH/AIS`) and late (`MIA/LUAD`) lesions, source-side `MIF` late-minus-precursor enrichment increased by 0.141 on average.
- The 95% bootstrap confidence interval was 0.094 to 0.188.
- Nineteen of 20 paired patients showed a positive change.
- The BH-adjusted Wilcoxon q-value was 0.00019.
- The combined receptor-side panel showed a smaller and non-robust paired progression difference: 0.014, 95% CI -0.055 to 0.088.

## 2026-05-30: Receptor-Side Perturbation Ranking Is CD74 Above CD44 Above CXCR4

Extended score-level perturbation to include `CD44`, which was already part of the MIF-axis receptor panel but was not included in the first focused perturbation run.

Decision:

- Prioritize `CD74` as the main receptor-side follow-up target.
- Retain `CD44` as an intermediate receptor-side candidate.
- Treat `CXCR4` as a weaker supporting co-receptor candidate in the current spatial-score model.

Rationale:

- Full score-dropout continuous-coupling loss:
  - `CD74`: 74.0%.
  - `CD44`: 20.0%.
  - `CXCR4`: 6.0%.
- The result is a dependency ranking within the current score model, not causal wet-lab validation.

## 2026-05-30: Matched-Gene Controls Retain MIF as a Candidate but Temper Specificity Claims

Added expression-matched single-gene, broad-signature, and tissue-density controls for source-side MIF spatial enrichment.

Decision:

- Keep `MIF` as the lead source-side follow-up candidate.
- Do not claim that the source-side MIF trend is uniquely specific.
- Treat the full MIF score-dropout result as expected model dependency because the selected source panel contains only `MIF`.
- Give more interpretive weight to the receptor-side `CD74 > CD44 > CXCR4` comparison.

Rationale:

- In paired-patient analysis, `MIF` remains at the 95th percentile of 20 expression-matched single-gene controls.
- The matched-control empirical upper-tail p-value is 0.095, so the result supports prioritization but not uniqueness.
- `PRDX1` has a slightly larger matched-control paired effect, and `BSG` also ranks near MIF.
- MIF enrichment is not strongly correlated with in-tissue spot count across sections: Spearman r=0.122, p=0.369.
- MIF enrichment is associated with permutation-null mean: r=0.352, p=0.0077, so broader background spatial structure should remain explicit in interpretation.

## 2026-05-30: Focused Non-Spatial Cohorts Support MIF but Constrain Receptor-Side Interpretation

Added candidate-gene summaries for the ordered `GSE308103` snRNA cohort and focused validation panels combining `GSE308103` with `GSE282617` bulk trends.

Decision:

- Keep source-side `MIF` as the lead progression-associated follow-up candidate.
- Keep `CD74` as the main receptor-side experimental follow-up target because of its receptor-score dependency, not because it shows a large monotonic expression increase.
- Treat `CXCR4` as supporting and `CD44` as an intermediate receptor-score contributor with non-monotonic expression behavior.
- Use `SPP1/TREM2/PLA2G7` as macrophage-state readouts.
- Describe snRNA and bulk panels as orthogonal trend support and boundary-setting, not spatial communication evidence.

Rationale:

- In `GSE308103`, epithelial `MIF` mean normalized expression increased from 0.227 in Normal samples to 0.577 in LUAD.
- Macrophage `CD74` increased modestly from 4.495 to 4.766, while `CXCR4` increased from 0.580 to 0.897.
- `SPP1`, `TREM2`, and `PLA2G7` increased in macrophages, whereas `CD44` and the inflammatory benchmark genes did not show a uniform LUAD increase.
- In `GSE282617`, bulk `MIF`, `CXCR4`, and `SPP1` increased from Normal to IAC; `CD74` increased modestly and `CD44` decreased.

## 2026-05-31: Reference Audit Narrows the Manuscript Claim

Added a structured reference audit to separate novelty-boundary, mechanism-support, macrophage-state, dataset, and background roles.

Decision:

- Keep the headline focused on specificity-audited public-data prioritization of `MIF-CD74`, not discovery of early-LUAD inflammatory niches.
- Cite the Peng *Cancer Cell* 2026 paper as the closest overlap and benchmark for alveolar progenitor/proinflammatory niche biology.
- Cite MIF-CD74 lung-cancer and receptor-complex papers as biological plausibility, while explicitly stating the axis is not de novo.
- Cite CD74 part-solid LUAD spatial work and our own non-spatial validation as reasons not to claim a simple receptor-side monotonic progression program.
- Cite `TREM2+` and `SPP1+` macrophage papers as support for treating `SPP1/TREM2/PLA2G7` as macrophage-state readouts.

Output:

- `docs/reference_audit_2026-05-31.md`
- `results/tables/reference_audit_2026-05-31.csv`

## 2026-05-31: Covariate Sensitivity Supports Source-Side MIF Robustness

Added linear sensitivity models for source-side `MIF` spatial enrichment.

Decision:

- Keep the source-side MIF spatial progression signal as robust to basic spatial covariates.
- Continue treating this as robustness support rather than causal evidence.
- Keep receptor-side `CD74` prioritization tied to score-level spatial dependency rather than monotonic expression progression.

Rationale:

- Sample-level late-versus-precursor coefficient was 0.133 unadjusted and 0.124 after adjustment for in-tissue spot count, permutation-null mean, and neighborhood radius.
- Paired-patient late-minus-precursor effect remained 0.141 after adjustment for changes in those covariates.
- These models address tissue-geometry and null-background concerns without changing the main interpretation.

## 2026-05-31: Communications Biology Reference Strategy

Added the first target-specific numbered reference pass to the Communications Biology working draft.

Decision:

- Use eight core references in the first target-specific draft rather than a broad review-style bibliography.
- Keep Peng et al. as the closest-overlap and dataset-context reference.
- Keep Zhu et al. as discovery scRNA/spatial context.
- Use Cao et al., Shi et al., and Schwartz et al. for MIF-CD74/CD44/CXCR4 biological plausibility.
- Use Zhang et al. iScience as the receptor-side CD74 boundary reference.
- Use TREM2 and SPP1 macrophage references only to support macrophage-state readout wording.

Rationale:

- The manuscript's novelty is a specificity-audited public-data prioritization workflow, not pathway discovery.
- A short, role-mapped reference list reduces the risk of diluting the main claim.
- The reference pass preserves cautious wording around score-level in-silico target prioritization.

Output:

- `docs/manuscript_communications_biology_draft.md`
- `results/tables/communications_biology_reference_list.csv`

## 2026-06-03: GRN-Level Virtual Perturbation Method Choice

Decision:

- Use a scTenifoldKnk-inspired GRN-level virtual perturbation prioritization layer on `GSE308103` rather than CellOracle for the first implementation.
- Keep the existing GSE307534 score-dropout/down-weighting analysis separate and label it as score-level in-silico target prioritization.
- Label the new GSE308103 analysis as GRN-level virtual perturbation prioritization, with an explicit boundary that it is not wet-lab knockout validation or causal proof.

Rationale:

- CellOracle is better suited to TF-driven cell-identity perturbation, particularly when a base GRN can be constructed from scATAC-seq, bulk ATAC-seq, promoter databases, or a trusted TF-target prior.
- The current candidate list contains mostly non-TF ligand, receptor, and macrophage-state genes: `MIF`, `CD74`, `CD44`, `CXCR4`, `SPP1`, `TREM2`, and `PLA2G7`.
- `GSE308103` provides local snRNA raw-count matrices and broad-class assignments, allowing a transparent expression-derived GRN layer to be added without requiring matched chromatin data.
- A local `scTenifoldpy` smoke test was too slow for immediate project use, so the implemented workflow keeps the core scTenifoldKnk idea of removing target outgoing edges but uses a compact positive-correlation network and signal-propagation summary.

Output:

- `src/luad_niche/grn_perturbation.py`
- `tests/test_grn_perturbation.py`
- `scripts/run_gse308103_grn_virtual_perturbation.py`
- `docs/grn_virtual_perturbation_2026-06-03.md`
- `results/tables/gse308103_grn_virtual_perturbation_gene_effects.csv`
- `results/tables/gse308103_grn_virtual_perturbation_signature_impacts.csv`
- `results/tables/gse308103_grn_virtual_perturbation_target_ranking.csv`
- `results/tables/gse308103_grn_virtual_perturbation_network_summary.csv`
- `results/tables/gse308103_grn_virtual_perturbation_config.json`

## 2026-06-03: GRN Dataset Scope and External Validation

Decision:

- Keep `GSE308103` as the primary GRN-level virtual perturbation dataset.
- Use other public cohorts as external validation or boundary-setting layers, not as parallel GRN-discovery cohorts.

Rationale:

- GRN-level virtual perturbation needs single-cell or single-nucleus expression matrices, broad cell-class assignments, and sufficient cells in the compartments being modeled.
- `GSE308103` is the best local dataset for this purpose because it is an ordered early-LUAD single-nucleus progression cohort with local raw-count matrices and enough epithelial/macrophage cells.
- `GSE307534` is better used for spatial support, `GSE131907` for specificity, `GSE189357` for early-LUAD scRNA state context, and `GSE164789` for tumor-adjacent scRNA state context.
- Running every dataset through the same GRN perturbation workflow would look more comprehensive but would mix different biological scopes and input quality; for the current manuscript, cross-dataset validation of the GRN-prioritized states is more defensible.

Output:

- `scripts/run_gse308103_grn_robustness.py`
- `scripts/build_grn_cross_dataset_validation.py`
- `results/tables/gse308103_grn_virtual_perturbation_robustness_summary.csv`
- `results/tables/grn_cross_dataset_signature_validation.csv`
- `results/tables/grn_cross_dataset_signature_validation_summary.csv`
- `results/supplementary_tables/ST28_gse308103_grn_virtual_perturbation_robustness_summary.csv`
- `results/supplementary_tables/ST29_grn_cross_dataset_signature_validation.csv`
- `results/supplementary_tables/ST30_grn_cross_dataset_signature_validation_summary.csv`

## 2026-06-03: Perturb-seq / GEARS / scGen / CPA Feasibility Decision

Decision:

- Do not add GEARS, scGen, or CPA as primary manuscript evidence in the current version.
- Keep the GSE307534 score-dropout/down-weighting analysis as score-level in-silico target prioritization.
- Keep the GSE308103 network analysis as GRN-level virtual perturbation prioritization.
- Treat Perturb-seq response-prediction modeling as a conditional future extension.

Rationale:

- GEARS, scGen, and CPA are strongest when trained or evaluated against matched perturbation-response data.
- The current project uses observational public early-LUAD spatial, snRNA, scRNA, and bulk datasets.
- A feasibility screen found broad public perturbation resources and some lung-cancer-adjacent perturbation datasets, but not a matched early-LUAD epithelial-myeloid dataset perturbing `MIF`, `CD74`, `CD44`, `CXCR4`, `SPP1`, `TREM2`, or `PLA2G7`.
- Adding a mismatched perturbation-response model would weaken the causal boundary of the manuscript.

Output:

- `docs/perturbseq_gears_scgen_cpa_feasibility_2026-06-03.md`
- `docs/manuscript_communications_biology_draft.md`
- `results/tables/communications_biology_reference_list.csv`
- `results/references/communications_biology_references.ris`
