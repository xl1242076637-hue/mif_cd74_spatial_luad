# Word Export Notes

Date: 2026-06-04

## Generated files

| File | Purpose |
|---|---|
| `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx` | Chinese Word review draft with synchronized manuscript text, expanded numbered citations, 36 references, 9 figures, and 5 tables |
| `docs/manuscript_communications_biology_draft_zh_zotero_fields.docx` | Experimental Zotero-field Word draft using Zotero CSL field instructions and imported Zotero item URIs |
| `docs/manuscript_communications_biology_draft_en_with_figures_tables.docx` | English Word review draft with manuscript text, expanded numbered citations, 36 references, 9 figures, and 5 tables |
| `docs/manuscript_communications_biology_bilingual_with_figures_tables.docx` | English-Chinese side-by-side review draft with manuscript text, 9 figures, and 5 key tables |
| `results/references/communications_biology_references.ris` | Zotero-ready RIS reference file generated from the project DOI/reference list |
| `scripts/export_word_manuscript.py` | Reproducible Word/RIS export script |
| `scripts/export_english_word_manuscript.py` | Reproducible English Word export script |
| `scripts/export_bilingual_word_manuscript.py` | Reproducible bilingual Word review export script |

## Contents of the Word file

- Chinese synchronized manuscript text updated through the GRN-level prioritization layer.
- Numbered in-text citations and 36-reference list.
- Main Figures 1-5.
- Supplementary Figures 1-4.
- Table 1: dataset composition and evidence roles.
- Table 2: integrated candidate-axis priority ranking.
- Table 3: GSE307534 paired-patient spatial progression statistics.
- Table 4: score-level in-silico target-prioritization continuous-coupling ranking.
- Table 5: supplementary table manifest.

## Bilingual review file

The bilingual review file is intended for internal editing rather than journal submission. It uses a landscape Word page and a two-column body layout:

- Left column: current English Communications Biology working draft.
- Right column: synchronized Chinese review draft.
- Figure section: bilingual captions for Main Figures 1-5 and Supplementary Figures 1-4.
- Table section: bilingual titles and column labels for the five key manuscript-facing tables.

The bilingual file was checked at the Word XML level after export. It contains the expected bilingual headers, `score-level in-silico target prioritization`, `受体侧 CD74`, `Supplementary Figure 4`, `补充图 4`, and `GRN-level virtual perturbation`; it embeds 9 media files; and it does not contain the older overbroad wording `MIF 和 CD74 排名靠前`, `因果证明`, or `实验验证了`.

## Separate single-language files

The separate Chinese and English Word files are intended for cleaner reading and editing when side-by-side layout is not needed:

- `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`: Chinese-only review manuscript.
- `docs/manuscript_communications_biology_draft_en_with_figures_tables.docx`: English-only review manuscript.

Both single-language files include 9 embedded figure images and the same five manuscript-facing tables. Word XML QA confirmed that both contain the expected score-level and GRN-level boundary wording, while the older overbroad phrasing is absent.

## Zotero status

After restarting Codex, Zotero MCP write access was available through the configured Zotero Web API credentials. The manuscript reference set was expanded from 8 to 30 DOI-verified references and imported into the user's Zotero Web library. References 31-36 were added later for GRN and perturbation-response model context and are included in the updated CSV/RIS registry.

The full reference registry is stored in:

```text
results/tables/communications_biology_reference_list.csv
```

It includes reference number, role, short citation, year, journal, DOI, URL, title and authors. Zotero item keys for the originally imported references 1-30 are mapped in `scripts/export_word_manuscript.py`; references 31-36 are included in the generated RIS and embedded as CSL itemData in the Zotero-field Word draft, but they should still be imported into Zotero before a fully URI-linked Zotero Word regeneration.

The Codex Zotero MCP config was changed to Web API mode (`ZOTERO_LOCAL=false`) so future Zotero MCP searches and writes target the configured user library instead of the local-only `users/0` API.

Two Word outputs are kept:

- `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx` is the stable readable review file with numbered citations and 36 references.
- `docs/manuscript_communications_biology_draft_zh_zotero_fields.docx` contains Zotero CSL field instructions. References 1-30 are linked to imported Zotero item URIs; references 31-36 currently carry CSL itemData and need Zotero import for full URI linkage. This is a best-effort automated field-code export; Word/Zotero may still require opening the document and clicking Zotero refresh to fully activate or normalize the fields.

Fallback route if Zotero does not recognize the automated fields:

1. Import `results/references/communications_biology_references.ris` into Zotero if needed.
2. Open the stable Word draft.
3. Use the Zotero Word plugin to replace bracketed numbered citations with live Zotero citations.
4. Insert a Zotero bibliography in the References section and choose the required journal citation style.

## Rebuild command

Run from the project root:

```powershell
python scripts\export_word_manuscript.py
python scripts\export_english_word_manuscript.py
python scripts\export_bilingual_word_manuscript.py
```

## Interpretation guardrail

The Word export preserves the manuscript wording that score-level in-silico target prioritization is a score-dependency analysis and that GRN-level virtual perturbation prioritization is a network-context target-ranking layer, not a wet-lab perturbation experiment or causal validation.
