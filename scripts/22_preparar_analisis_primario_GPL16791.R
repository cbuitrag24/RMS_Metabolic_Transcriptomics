# ============================================================
# 22_preparar_analisis_primario_GPL16791.R
# Preparación del análisis primario:
# RMS vs músculo normal dentro de GPL16791
# GSE108022
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar datos sin filtrar y diseño de muestras
# ------------------------------------------------------------

dge_all <- readRDS("GSE108022_edgeR_sin_filtrar.rds")

design_data <- read.csv(
  "DE_diseno_muestras.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!identical(colnames(dge_all$counts), design_data$samplename)) {
  stop("ERROR: el orden de las muestras no coincide.")
}

cat("Dataset completo:\n")
cat("Genes:", nrow(dge_all), "\n")
cat("Muestras:", ncol(dge_all), "\n")


# ------------------------------------------------------------
# 2. Seleccionar exclusivamente GPL16791
# ------------------------------------------------------------

keep_samples <- design_data$platform == "GPL16791"

dge <- dge_all[, keep_samples, keep.lib.sizes = FALSE]
design_sub <- design_data[keep_samples, , drop = FALSE]

if (!identical(colnames(dge$counts), design_sub$samplename)) {
  stop("ERROR: el orden de las muestras GPL16791 no coincide.")
}

cat("\nSubconjunto GPL16791:\n")
cat("Muestras:", ncol(dge), "\n")

cat("\nMuestras por condición:\n")
print(table(design_sub$condition))


# ------------------------------------------------------------
# 3. Definir condición
# ------------------------------------------------------------

design_sub$condition <- factor(
  design_sub$condition,
  levels = c("Normal", "RMS")
)

design <- model.matrix(
  ~ condition,
  data = design_sub
)

cat("\nColumnas del diseño:\n")
print(colnames(design))

cat("\nDimensiones del diseño:\n")
print(dim(design))

cat("\nRango del diseño:", qr(design)$rank, "\n")

if (qr(design)$rank != ncol(design)) {
  stop("ERROR: la matriz de diseño no es de rango completo.")
}


# ------------------------------------------------------------
# 4. Filtrado específico para este análisis
# ------------------------------------------------------------

genes_antes <- nrow(dge)

keep_genes <- filterByExpr(
  dge,
  design = design
)

dge <- dge[
  keep_genes,
  ,
  keep.lib.sizes = FALSE
]

cat("\nFiltrado con filterByExpr:\n")
cat("Genes antes:", genes_antes, "\n")
cat("Genes conservados:", nrow(dge), "\n")
cat("Genes eliminados:", genes_antes - nrow(dge), "\n")
cat(
  "Porcentaje conservado:",
  round(100 * nrow(dge) / genes_antes, 2),
  "%\n"
)


# ------------------------------------------------------------
# 5. Normalización TMM dentro de GPL16791
# ------------------------------------------------------------

dge <- normLibSizes(
  dge,
  method = "TMM"
)

cat("\nResumen de factores TMM:\n")
print(summary(dge$samples$norm.factors))

if (anyNA(dge$samples$norm.factors)) {
  stop("ERROR: existen factores TMM NA.")
}

if (any(dge$samples$norm.factors <= 0)) {
  stop("ERROR: existen factores TMM no positivos.")
}


# ------------------------------------------------------------
# 6. Mostrar factores TMM por condición
# ------------------------------------------------------------

control_TMM <- data.frame(
  samplename = design_sub$samplename,
  condition = design_sub$condition,
  lib_size = dge$samples$lib.size,
  norm_factor = dge$samples$norm.factors,
  effective_lib_size =
    dge$samples$lib.size * dge$samples$norm.factors,
  stringsAsFactors = FALSE
)

cat("\nFactores TMM por condición:\n")

print(
  aggregate(
    norm_factor ~ condition,
    data = control_TMM,
    FUN = function(x) {
      c(
        n = length(x),
        min = min(x),
        median = median(x),
        mean = mean(x),
        max = max(x)
      )
    }
  )
)


# ------------------------------------------------------------
# 7. Guardar objetos
# ------------------------------------------------------------

saveRDS(
  dge,
  "GSE108022_GPL16791_filtrado_TMM.rds"
)

saveRDS(
  design,
  "GSE108022_GPL16791_design_matrix.rds"
)

write.csv(
  design_sub,
  "GPL16791_diseno_muestras.csv",
  row.names = FALSE
)

write.csv(
  control_TMM,
  "GPL16791_control_TMM.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("- GSE108022_GPL16791_filtrado_TMM.rds\n")
cat("- GSE108022_GPL16791_design_matrix.rds\n")
cat("- GPL16791_diseno_muestras.csv\n")
cat("- GPL16791_control_TMM.csv\n")