# Manuscript tables — source-file guide

This document is a navigation guide for the manuscript's retained analysis files. It does not assert that every manuscript table has been independently regenerated or checked against the latest Word manuscript.

## Concordant Hallmark pathways (13)

Relevant committed scripts:

- `scripts/22_enriquecimiento_Hallmark_universo_comun_GSE108022.R`
- `scripts/23_enriquecimiento_Hallmark_GSE28511_universo_comun_9428.R`
- `scripts/24_comparar_Hallmark_GSE108022_GSE28511.R`
- `scripts/25_leading_edge_Hallmarks_replicados.R`
- `scripts/28_mapear_genes_leading_edges_replicados.R`

The repository's reproducibility inventory records 13 concordant Hallmarks. Consult the committed source results and the manuscript's Supplementary Table S1 for the final published presentation; this guide does not reproduce or change the numerical table.

## Shared metabolic core (16 genes)

The reproducibility inventory records a 16-gene intersection between shared leading edges of `HALLMARK_OXIDATIVE_PHOSPHORYLATION` and `HALLMARK_FATTY_ACID_METABOLISM`, corresponding to manuscript Table 2.

Relevant committed scripts:

- `scripts/25_leading_edge_Hallmarks_replicados.R`
- `scripts/28_mapear_genes_leading_edges_replicados.R`
- `scripts/43_verificar_seleccion_16_genes_metabolicos.R`

The retained shared leading-edge CSV is stored under `data/` (see `REPRODUCIBILITY_AUDIT.md` for its provenance). The gene-selection script documents the intersection; no additional selection or recalculation is performed by this guide.

## CKM and myogenic-marker association results

- `results/resultados_originales_Figuras_4_5.txt` preserves the original six PC1–myogenic-marker Spearman associations and 16 individual metabolic-gene–CKM associations.
- `results/RMS_Figuras_4_5_objetos_originales.RData` preserves compact original R objects.
- `scripts/45_exportar_resultados_originales_Figuras_4_5.R` is available to export already computed objects, when execution is later requested.

These results must not be mistaken for the differentiation-score correlations shown in the three-panel manuscript Figure 4.

## Adjusted regression

- The original `summary(modelo_final)` is retained in `results/resultados_originales_Figuras_4_5.txt`.
- `scripts/47_generar_Figura_5_desde_modelo_original.R` uses the archived fitted model for a coefficient plot, without refitting it.

## Important boundaries

- The table numbering and exact table-to-file mapping must follow the current final manuscript, not assumptions made in this guide.
- This document is not a substitute for the source CSV files or the manuscript's supplementary tables.
- No analyses, validations, numerical changes, or manuscript edits were performed in preparing this guide.
