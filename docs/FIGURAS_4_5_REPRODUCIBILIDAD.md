# Manuscript Figures 4 and 5 — archived data and plotting scripts

This note identifies the committed files supporting the manuscript figures. It does **not** claim that the figure scripts have been executed or validated.

## Original archived outputs

- `results/resultados_originales_Figuras_4_5.txt`: exported original Spearman results for PC1 versus six myogenic markers, gene-wise CKM associations, and the summary of the adjusted model.
- `results/RMS_Figuras_4_5_objetos_originales.RData`: compact R workspace exported from the investigator's retained session. Intended to retain original PCA, CKM expression, differentiation score, correlations, fitted models and sensitivity results. The binary file is kept as an archival input, not a replacement for documented raw-data processing.

The original full R workspace (`sesion_RMS_manuscrito_2026-10-17.RData`) is **not** required to run the scripts below.

## Source scripts

1. `scripts/44_cargar_objetos_originales_Figuras_4_5.R` — loads and lists the archived objects; does not perform statistical calculations.
2. `scripts/45_exportar_resultados_originales_Figuras_4_5.R` — exports archived result tables, model summaries and selected R objects; does not refit models.
3. `scripts/46_generar_Figura_4_desde_objetos_originales.R` — draws three scatterplot panels from archived observations (CKM vs metabolic PC1; CKM vs differentiation score; differentiation score vs metabolic PC1). Correlation annotations are transcribed from the investigator-supplied original Figure 4 PDF; they are **not** recalculated by this script.
4. `scripts/47_generar_Figura_5_desde_modelo_original.R` — draws the original model coefficients with approximate normal-theory 95% confidence intervals using the archived `lm` object. The model is **not** refitted.

## Intended run order

From the repository root in R:

```r
source("scripts/44_cargar_objetos_originales_Figuras_4_5.R")
source("scripts/45_exportar_resultados_originales_Figuras_4_5.R")
source("scripts/46_generar_Figura_4_desde_objetos_originales.R")
source("scripts/47_generar_Figura_5_desde_modelo_original.R")
```

The scripts were committed as source files. Their execution and output comparison with the original manuscript PDFs remain **pending**; do not describe the generated figures as validated reproductions.

## Expected generated files (not yet confirmed)

- `results/exportacion_objetos_originales/` — exports from script 45
- `figures/Figura_4_Metabolismo_Diferenciacion_RMS_recreada.pdf`
- `figures/Figura_5_Modelo_Ajustado_RMS_recreada.pdf`

## Interpretation and limitations

- Figure 4 is a layout reconstruction using archived observation vectors; it is not necessarily identical to the historical plotting code or PDF.
- Figure 5 uses `confint.default()`, which produces normal-theory intervals. The original interval method has not been established from the original plotting script.
- The full upstream workflow from public GEO data to the archived objects is not established by these scripts alone.
- Do not replace the original manuscript figures or assert full computational reproducibility without a separate, deliberate check.
