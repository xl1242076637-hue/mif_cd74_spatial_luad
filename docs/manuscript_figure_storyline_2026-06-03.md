# Manuscript and Figure Storyline

Date: 2026-06-03

## One-Sentence Argument

In early LUAD progression, a specificity-audited public-data workflow prioritizes epithelial `MIF` and myeloid `CD74` as a bounded epithelial-myeloid communication candidate, while treating `SPP1/TREM2/PLA2G7` as a macrophage-state readout and computational perturbation outputs as target-prioritization evidence rather than causal proof.

## Recommended Results Logic

1. Public multi-cohort framework.
   - Purpose: establish that each dataset has a predefined role.
   - Evidence: seven public datasets, with `GSE307534` as spatial backbone, `GSE308103` as snRNA progression backbone, and `GSE131907` as specificity reference.
   - Figure: Figure 1.

2. Specificity audit reframes the original hypothesis.
   - Purpose: show that the original `SPP1` macrophage niche idea was too broad as a primary novelty claim.
   - Evidence: the 30-gene refined `SPP1` macrophage signature contains 14 expected myeloid genes, 15 off-target genes and 1 missing gene in `GSE131907`.
   - Figure: Figure 3 and Supplementary Figure 3.

3. Integrated evidence ranks `MIF-CD74/CXCR4` first.
   - Purpose: move from broad macrophage niche to a sharper candidate axis.
   - Evidence: `MIF-CD74/CXCR4` priority score 0.789; `SPP1/TREM2/PLA2G7` second at 0.748.
   - Figure: Figure 2.

4. Patient-aware spatial progression identifies source-side `MIF`.
   - Purpose: present the strongest primary result.
   - Evidence: 20 paired patients; source-side `MIF` late-minus-precursor delta 0.141; 95% CI 0.094-0.188; q=0.00019; 19/20 positive.
   - Figure: Figure 4.

5. Controls and sensitivity analyses temper the `MIF` claim.
   - Purpose: show that `MIF` is prioritized but not uniquely specific.
   - Evidence: `MIF` lies at the 95th percentile of 20 expression-matched controls, empirical upper-tail p=0.095; `PRDX1` slightly exceeds MIF; spatial geometry does not remove the positive effect.
   - Figure: Supplementary Figure 1.
   - Recommendation: this section should appear immediately after the spatial progression section, before perturbation-prioritization results.

6. Score-level in-silico target prioritization ranks receptor-side `CD74`.
   - Purpose: prioritize follow-up targets within the current spatial score model.
   - Evidence: receptor-side full-dropout loss ranks `CD74` 74.0%, `CD44` 20.0%, `CXCR4` 6.0%; `MIF` dropout collapses its single-gene source score and should be interpreted conservatively.
   - Figure: Figure 5.

7. GRN-level virtual perturbation adds network context.
   - Purpose: separate target-linked expression neighborhoods in `GSE308103`.
   - Evidence: `CD74` maps to a C1Q/APOE/ACP5 macrophage neighborhood; `CXCR4` and `PLA2G7` map more strongly to SPP1-like macrophage neighborhoods; `MIF` maps to a proliferating epithelial neighborhood; top signature assignments are stable across nine threshold/seed settings.
   - Figure: Supplementary Figure 4.

8. Orthogonal non-spatial cohorts support and constrain interpretation.
   - Purpose: show cross-dataset expression support without overclaiming spatial communication.
   - Evidence: `GSE308103` epithelial `MIF` increases from Normal to LUAD; bulk `MIF`, `CXCR4` and `SPP1` increase in `GSE282617`; receptor-side `CD74` is not a simple monotonic spatial progression marker.
   - Figure: Supplementary Figure 2.

## Main Figures

| Figure | Current file | Main message | Manuscript role |
|---|---|---|---|
| Figure 1 | `results/figures/nature_redesign/nature_figure1_workflow_dataset_composition.*` | Seven public cohorts were assigned predefined evidence roles | Orientation and reproducibility |
| Figure 2 | `results/figures/nature_redesign/nature_figure2_axis_evidence_perturbation.*` | Integrated evidence ranks `MIF-CD74/CXCR4` first | Candidate-axis selection |
| Figure 3 | `results/figures/nature_redesign/nature_figure3_specificity_audit.*` | Specificity audit separates `MIF-CD74` from broad `SPP1` macrophage signal | Novelty boundary and signature quality |
| Figure 4 | `results/figures/nature_redesign/nature_figure4_spatial_axis_progression.*` | Paired-patient spatial analysis supports source-side `MIF` progression | Primary biological signal |
| Figure 5 | `results/figures/nature_redesign/nature_figure5_virtual_perturbation_priority.*` | Score-level target prioritization ranks `CD74` above `CD44` and `CXCR4` on receptor side | Follow-up target prioritization |

## Supplementary Figures

| Figure | Current file | Main message | Manuscript role |
|---|---|---|---|
| Supplementary Figure 1 | `results/figures/supplementary_figure_mif_spatial_controls.*` | `MIF` is high-ranking among matched controls but not uniquely specific | Robustness and claim calibration |
| Supplementary Figure 2 | `results/figures/supplementary_figure_focused_orthogonal_validation.*` | snRNA and bulk cohorts support source-side `MIF` and macrophage-state trends while constraining receptor interpretation | Orthogonal support |
| Supplementary Figure 3 | `results/figures/supplementary_figure_spp1_signature_refinement.*` | The original refined `SPP1` signature contains substantial off-target signal | Signature-audit transparency |
| Supplementary Figure 4 | `results/figures/supplementary_figure_grn_virtual_perturbation.*` | GRN-level virtual perturbation prioritization stratifies target-linked epithelial and macrophage neighborhoods | Network-context prioritization |

## Figure-Level Issues to Fix Before Final Submission

1. Figure 1 has been regenerated with the GRN-level prioritization layer included in the workflow panel.
2. Figure 2 currently carries integrated evidence and score-level target prioritization. If the main manuscript keeps GRN as Supplementary Figure 4, Figure 2 does not need to absorb GRN, but its caption should not imply that score-level dropout is the only target-prioritization layer.
3. Figure 5 title/file name uses `virtual_perturbation_priority`; in text, call this `score-level in-silico target prioritization` to avoid overclaiming.
4. Supplementary Figure 4 is currently the only figure showing GRN-level results. This is appropriate unless the final manuscript wants to elevate GRN as a main analytical layer.
5. The textual order has been updated to spatial result -> MIF controls -> score-level target prioritization -> GRN-level target prioritization -> orthogonal cohort support, which should be more reviewer-friendly because the primary `MIF` claim is calibrated before computational target-prioritization layers are introduced.

## Recommended Final Figure Architecture

Keep five main figures:

- Figure 1: study design and dataset evidence roles.
- Figure 2: integrated evidence ranking.
- Figure 3: specificity audit.
- Figure 4: paired-patient spatial progression.
- Figure 5: score-level target prioritization.

Keep four supplementary figures:

- Supplementary Figure 1: MIF controls and tissue-geometry sensitivity.
- Supplementary Figure 2: orthogonal snRNA/bulk expression support.
- Supplementary Figure 3: detailed SPP1 signature refinement.
- Supplementary Figure 4: GRN-level virtual perturbation prioritization and network context.

Current conclusion: the paper is strongest as a public-data prioritization manuscript, not a perturbation-response modeling manuscript. The figures should therefore emphasize audit, patient-aware spatial signal, controlled target prioritization and boundary-setting.

## 2026-06-03 Figure Revision Status

Completed after the storyline recheck:

- Regenerated Figure 1 with a six-step workflow: cohorts -> audit -> spatial -> score-level priority -> GRN-level priority -> MIF-CD74 target.
- Updated the shared Nature-style palette so `MIF`/epithelial signals, macrophage/GRN signals, support axes and neutral context use a consistent low-saturation manuscript palette.
- Regenerated all Nature-redesigned main figures and the GRN supplementary figure as SVG, PDF, TIFF and PNG.
- Changed Figure 5 visible wording to `Score reduction` to keep the figure aligned with score-level target-prioritization language.
- Reordered the Communications Biology Results narrative so `MIF` control and sensitivity analyses appear immediately after the paired-patient spatial result.
- Updated Figure 1, Figure 2 and Figure 5 legends to align with the revised six-step workflow and score-level dropout/down-weighting language.
