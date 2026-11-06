# Execution and dependency guide — RMS metabolic transcriptomics

This guide is based on the scripts and documentation currently retained in this repository. It is **not** evidence of a clean end-to-end run. The repository contains historical analysis scripts, documented reconstructions and archived outputs; these must remain distinguishable.

## Study scope

The manuscript compares rhabdomyosarcoma and normal skeletal muscle using GEO series GSE108022 and GSE28511. The analyses include differential expression, concordant Hallmark pathways, a 16-gene metabolic core, within-tumour PCA, myogenic-marker associations, differentiation-related scores and adjusted regression.

## Data acquisition

GEO series:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108022
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28511

The RNA-seq workflow described in the root README requires `GSE108022_counts_matrix.txt.gz`, `GSE108022_metadata.txt.gz` and `GSE108022_family.xml`, plus intermediate files from preceding steps. Consult individual scripts for exact local paths and preprocessing assumptions; this document does not assert that all upstream input files are committed.

## Documented RNA-seq workflow (GSE108022, GPL16791)

The existing README describes the following logical order. Some older scripts are in the repository root and others in `scripts/`; follow the actual paths when running:

1. `14_preparar_diseno_DE.R` — sample metadata and design
2. `15_preparar_counts_edgeR.R` — count-matrix DGEList
3. `22_preparar_analisis_primario_GPL16791.R` — GPL16791 selection, filtering, TMM
4. `23_dispersion_QL_GPL16791.R` — robust quasi-likelihood fit
5. `24_DE_RMS_vs_Normal_GPL16791.R` — differential-expression comparison
6. `26_corregir_tabla_DE_genes.R` — corrected definitive DE tables (`_CORREGIDO`)

The preliminary exports from step 5 should not be substituted for the corrected DE tables from step 6.

## GSE28511 microarray workflow

The retained script sequence is:

1. `scripts/03_preparar_log2_GSE28511.R`
2. `scripts/07_preparar_diseno_GSE28511.R`
3. `scripts/09_preparar_limma_primario_GSE28511.R`
4. `scripts/10_ajuste_limma_primario_GSE28511.R`
5. `scripts/11_integrar_anotacion_limma_GSE28511.R`

This is a navigation sequence, not a claim that all dependencies and paths have been tested in a clean environment.

## Cross-dataset Hallmark analysis

Relevant retained scripts:

- `scripts/22_enriquecimiento_Hallmark_universo_comun_GSE108022.R`
- `scripts/23_enriquecimiento_Hallmark_GSE28511_universo_comun_9428.R`
- `scripts/24_comparar_Hallmark_GSE108022_GSE28511.R`
- `scripts/25_leading_edge_Hallmarks_replicados.R`
- `scripts/28_mapear_genes_leading_edges_replicados.R`

The 16-gene core can be derived from the retained `data/Hallmark_leadingEdge_genes_compartidos.csv` using `scripts/43_verificar_seleccion_16_genes_metabolicos.R`. The root README provides its invocation and limitations.

## Figures

- Figure 3: `scripts/42_reconstruccion_PCA_metabolica_RMS_PENDIENTE_VALIDACION.R` — documented reconstruction, not the lost original script. The historical filename is retained.
- Figures 4–5: archived source objects and plotting scripts 44–47. See `docs/FIGURAS_4_5_REPRODUCIBILIDAD.md`.

The scripts 46 and 47 are plotting reconstructions, not verified exact reproductions of the supplied original PDFs.

## Software and packages

The workflows use R. The RNA-seq scripts refer to edgeR and the microarray workflow to limma; pathway-enrichment work uses fgsea outputs. Exact R version, complete package versions and all additional dependencies have **not** been pinned in a tested lockfile. Install and document them from the individual scripts before any future independent rerun. Do not infer package versions from the dates of the files.

## Reproducibility boundaries

- The committed source files and compact archived R objects support transparency and future reproduction work.
- No complete raw-data-to-manuscript reproduction has been demonstrated by this guide.
- Figure-generation scripts 46–47 have not been run or compared with the original PDFs as part of this documentation step.
- No new statistical analyses, figure validations or manuscript edits were performed to prepare this file.

For the detailed existing inventory, consult `REPRODUCIBILITY_AUDIT.md`.
