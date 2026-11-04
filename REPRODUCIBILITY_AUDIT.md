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

Presence is verified. **Local investigator-reported verification (2026-11-01):** `scripts/25_leading_edge_Hallmarks_replicados.R` was executed successfully using the original `Hallmark_fgsea_GSE108022_universo_comun_9428.rds`, `Hallmark_fgsea_GSE28511_universo_comun_9428.rds`, and `Hallmark_replicados_GSE108022_GSE28511.csv` in `validacion_GSE28511/`. The console reported 13 concordant Hallmarks, 495 shared Hallmark–gene associations, and 374 unique genes. The script generated `Hallmark_leadingEdge_concordancia_resumen.csv`, `Hallmark_leadingEdge_genes_detalle.csv`, `Hallmark_leadingEdge_genes_compartidos.csv`, `Hallmark_leadingEdge_genes_frecuencia.csv`, and `Hallmark_leadingEdge_concordancia.rds`. The original `Hallmark_leadingEdge_genes_compartidos.csv` has since been inspected and committed to `data/`. The intersection of the shared leading edges for `HALLMARK_OXIDATIVE_PHOSPHORYLATION` (89 genes) and `HALLMARK_FATTY_ACID_METABOLISM` (47 genes) contains exactly the 16 genes listed in manuscript Table 2. `scripts/43_verificar_seleccion_16_genes_metabolicos.R` documents this comparison; its full execution in a clean environment has not been independently demonstrated. Other generated intermediate files remain local. End-to-end execution from raw data remains unverified.

## Still to locate and verify against the manuscript

- [x] Provenance of the final 16-gene selection: original shared leading-edge CSV committed under `data/`, with direct intersection matching Table 2 exactly; documented verification script 43 added. This does not yet establish the provenance of Figure 2 or an end-to-end rerun.
- [x] Figure 3: documented PCA reconstruction added as `scripts/42_reconstruccion_PCA_metabolica_RMS_PENDIENTE_VALIDACION.R`. Executed locally by the investigator: 16 genes, 64 RMS tumours, PC1 36.481659%, PC2 13.131902%, all PC1 loadings negative. Reconstructed PDF visually compared with original; this is not the historical script. End-to-end reproduction from raw public data remains unverified.
- [x] Original PC1-versus-six-myogenic-marker Spearman results recovered from `resultado_PC1_miogenesis` in the investigator's retained R session. CKM: rho = -0.41341575, p = 0.0006852737, FDR = 0.004111642; MYH3: rho = -0.34166667, p = 0.0057247432, FDR = 0.017174229. MYOD1, MYOG, MYH1 and MYH2 have FDR >= 0.05. This is recovery of existing output, not a new calculation; original analysis script and a standalone committed output file remain to be documented. **Partially recovered:** the original `tabla_cor_CKM` object was printed from the investigator's retained R session, giving 16 gene-wise Spearman correlations with CKM and Benjamini–Hochberg FDR. Six genes meet FDR < 0.05: ECH1 (rho 0.53946886, FDR 6.767733e-05), PDHA1 (0.52284799, 7.462989e-05), ACO2 (0.50805861, 9.718995e-05), PDHB (0.39418498, 0.004107318), SDHD (0.39381868, 0.004107318), and MDH1 (0.34285714, 0.01479306). These are individual gene–CKM associations, **not** the PC1–CKM or CKM–differentiation-score correlations. The original object is retained locally; no code or output file has yet been committed for these correlations.
- [ ] Construction of the external differentiation-related score from the cited published gene lists.
- [ ] Correlations of CKM, metabolic PC1 and the differentiation score (Figure 4).
- [x] Original Figure 5 multivariable model object located in the investigator's retained R session (`sesion_RMS_manuscrito_2026-10-17.RData`), and its `summary(modelo_final)` inspected: `lm(PC1 ~ CKM + diferenciacion + fusion)`, n = 64 (58 residual df, six fitted coefficients); CKM beta = -0.38654, SE = 0.10067, p = 0.000307; differentiation beta = 0.38814, p = 0.415550; R-squared = 0.2436, adjusted R-squared = 0.1784, overall F(5,58) = 3.735, p = 0.005332. Fusion factor levels shown include PAX3-FOXO1, PAX3-NCOA1 and PAX7-FOXO1 relative to an unspecified reference category. The original model object is retained locally, not committed; source code, upstream data preparation and sensitivity analyses still require provenance documentation. This is recovery of an existing model, not a new computation.
- [ ] Scripts that **generate** the manuscript figures and tables from analysis outputs.
- [ ] Required data download instructions, source versions, dependency versions and complete run order.
- [ ] Reproduction test from original public data to final numerical results.

## Deliberately excluded from the primary analysis

The earlier `16_filtrar_genes_edgeR.R`, `17_normalizar_TMM_edgeR.R`, and `20_estimar_dispersion_edgeR.R` scripts refer to the broader, mixed-platform GSE108022 workflow, not the manuscript's GPL16791-only primary comparison. The exploratory Python scripts `10_pca_solo_RMS.py` and `11_grafico_pca_solo_RMS.py` do not implement the manuscript's 16-gene PCA. Do not substitute these for the manuscript analyses.

## Data and outputs policy

Prefer accession-based downloads of public raw data and scripts that regenerate intermediate objects and figures. Do not upload bulky raw datasets or final figures merely as substitutes for reproducible code. Include small indispensable metadata or gene lists only after confirming provenance and licensing.

## Next action

Figure 3 reconstruction is now in `scripts/42_reconstruccion_PCA_metabolica_RMS_PENDIENTE_VALIDACION.R` and was validated locally against the published figure. The original leading-edge concordance script has now run successfully on retained enrichment results. The investigator's original R session also retains `modelo_final`, `modelo_sensibilidad`, `rho_LOO_CKM`, `correlaciones_CKM`, `tabla_cor_CKM`, `score_diferenciacion`, and `pca_metabolico`. The Figure 5 model summary has been recovered; the other objects have not yet been inspected. The 16-gene metabolic-core selection is now traceable to the committed shared leading-edge CSV and documented in script 43. Next priority: locate or document the existing source code and inputs for manuscript Figures 2, 4 and 5, correlations and regression; do not claim end-to-end reproduction from raw data.

## Archived original Figure 4 and Figure 5 outputs (2026-11-04)

`results/resultados_originales_Figuras_4_5.txt` is a text export from the investigator's retained R session (`sesion_RMS_manuscrito_2026-10-17.RData`). It preserves the six PC1–myogenic-marker Spearman correlations with p values and FDR, 16 gene-wise CKM correlations with p values and FDR, and `summary(modelo_final)` for the adjusted Figure 5 model. The output was inspected against the user-provided text before inclusion. The binary RData backup was deliberately **not** uploaded. This export does **not** include the CKM–differentiation-score or PC1–differentiation-score correlations, nor the original code to regenerate the correlations, differentiation score, regression inputs or figures. It therefore supports result provenance but does not demonstrate full computational reproducibility.
