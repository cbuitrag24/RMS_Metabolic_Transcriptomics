# Reproducibility audit — RMS manuscript

This is a **working inventory**, not a claim that the complete manuscript has been reproduced. It records the files present in the repository as checked on 2026-11-01. The original analyses and scripts must not be silently altered.

## Scope

Manuscript analyses: GSE108022 (64 RMS tumours and 5 normal skeletal-muscle samples on GPL16791) and GSE28511 (18 RMS and 3 normal samples); differential expression, Hallmark enrichment, 16-gene metabolic selection, within-tumour PCA, myogenic-marker correlations, an external differentiation-related score, adjusted regression and sensitivity analyses.

## Already present — source code

**GSE108022 RNA-seq**
- `scripts/14_preparar_diseno_DE.R`
- `15_preparar_counts_edgeR.R` (repository root)
- `scripts/22_preparar_analisis_primario_GPL16791.R`
- `23_dispersion_QL_GPL16791.R` (repository root)
- `scripts/24_DE_RMS_vs_Normal_GPL16791.R`
- `scripts/26_corregir_tabla_DE_genes.R`

**GSE28511 microarray**
- `scripts/03_preparar_log2_GSE28511.R`
- `scripts/07_preparar_diseno_GSE28511.R`
- `scripts/09_preparar_limma_primario_GSE28511.R`
- `scripts/10_ajuste_limma_primario_GSE28511.R`
- `scripts/11_integrar_anotacion_limma_GSE28511.R`

**Cross-dataset enrichment and leading edges**
- `scripts/22_enriquecimiento_Hallmark_universo_comun_GSE108022.R`
- `scripts/23_enriquecimiento_Hallmark_GSE28511_universo_comun_9428.R`
- `scripts/24_comparar_Hallmark_GSE108022_GSE28511.R`
- `scripts/25_leading_edge_Hallmarks_replicados.R`
- `scripts/28_mapear_genes_leading_edges_replicados.R`
- `scripts/29_auditar_solapamiento_Hallmarks_entre_ejes.R`

Presence is verified; end-to-end execution and scientific outputs are **not yet verified**.

## Still to locate and verify against the manuscript

- [ ] Exact original code and intermediate inputs producing the final 16-gene list and Figure 2 (check whether existing leading-edge scripts suffice).
- [x] Figure 3: documented PCA reconstruction added as `scripts/42_reconstruccion_PCA_metabolica_RMS_PENDIENTE_VALIDACION.R`. Executed locally by the investigator: 16 genes, 64 RMS tumours, PC1 36.481659%, PC2 13.131902%, all PC1 loadings negative. Reconstructed PDF visually compared with original; this is not the historical script. End-to-end reproduction from raw public data remains unverified.
- [ ] Spearman analyses for the six myogenic markers and the PC1–CKM association.
- [ ] Construction of the external differentiation-related score from the cited published gene lists.
- [ ] Correlations of CKM, metabolic PC1 and the differentiation score (Figure 4).
- [ ] Multivariable regression, fusion-category coding and sensitivity analyses (Figure 5).
- [ ] Scripts that **generate** the manuscript figures and tables from analysis outputs.
- [ ] Required data download instructions, source versions, dependency versions and complete run order.
- [ ] Reproduction test from original public data to final numerical results.

## Deliberately excluded from the primary analysis

The earlier `16_filtrar_genes_edgeR.R`, `17_normalizar_TMM_edgeR.R`, and `20_estimar_dispersion_edgeR.R` scripts refer to the broader, mixed-platform GSE108022 workflow, not the manuscript's GPL16791-only primary comparison. The exploratory Python scripts `10_pca_solo_RMS.py` and `11_grafico_pca_solo_RMS.py` do not implement the manuscript's 16-gene PCA. Do not substitute these for the manuscript analyses.

## Data and outputs policy

Prefer accession-based downloads of public raw data and scripts that regenerate intermediate objects and figures. Do not upload bulky raw datasets or final figures merely as substitutes for reproducible code. Include small indispensable metadata or gene lists only after confirming provenance and licensing.

## Next action

Figure 3 reconstruction is now in `scripts/42_reconstruccion_PCA_metabolica_RMS_PENDIENTE_VALIDACION.R` and was validated locally against the published figure. Continue locating/verifying original scripts for remaining manuscript analyses; do not claim end-to-end reproduction from raw data.
