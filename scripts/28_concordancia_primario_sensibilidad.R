# ============================================================
# 28_concordancia_primario_sensibilidad.R
# Concordancia entre:
#   1) análisis primario GPL16791
#   2) sensibilidad 106 muestras ajustada por plataforma
# GSE108022
# ============================================================

# ------------------------------------------------------------
# 1. Cargar tablas completas
# ------------------------------------------------------------

primario <- read.csv(
  "GPL16791_DE_RMS_vs_Normal_todos_los_genes_CORREGIDO.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sensibilidad <- read.csv(
  "ALL106_DE_RMS_vs_Normal_ajustado_plataforma_todos.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Tablas cargadas.\n")

cat("\nGenes análisis primario:\n")
print(nrow(primario))

cat("\nGenes análisis de sensibilidad:\n")
print(nrow(sensibilidad))


# ------------------------------------------------------------
# 2. Controles de identidad
# ------------------------------------------------------------

columnas_necesarias <- c("gene", "logFC", "PValue", "FDR")

if (!all(columnas_necesarias %in% colnames(primario))) {
  stop("ERROR: faltan columnas necesarias en la tabla primaria.")
}

if (!all(columnas_necesarias %in% colnames(sensibilidad))) {
  stop("ERROR: faltan columnas necesarias en la tabla de sensibilidad.")
}

if (anyNA(primario$gene) || any(primario$gene == "")) {
  stop("ERROR: hay genes vacíos o NA en el análisis primario.")
}

if (anyNA(sensibilidad$gene) || any(sensibilidad$gene == "")) {
  stop("ERROR: hay genes vacíos o NA en el análisis de sensibilidad.")
}

if (anyDuplicated(primario$gene) > 0) {
  stop("ERROR: hay genes duplicados en el análisis primario.")
}

if (anyDuplicated(sensibilidad$gene) > 0) {
  stop("ERROR: hay genes duplicados en el análisis de sensibilidad.")
}

cat("\nVERIFICACIÓN: identificadores génicos válidos y no duplicados.\n")


# ------------------------------------------------------------
# 3. Emparejar exclusivamente por símbolo génico
# ------------------------------------------------------------

comparacion <- merge(
  primario[, c("gene", "logFC", "PValue", "FDR")],
  sensibilidad[, c("gene", "logFC", "PValue", "FDR")],
  by = "gene",
  suffixes = c("_primario", "_sensibilidad"),
  all = FALSE,
  sort = FALSE
)

cat("\nGenes comunes entre ambos análisis:\n")
print(nrow(comparacion))

cat("\nGenes exclusivos del análisis primario:\n")
print(sum(!primario$gene %in% comparacion$gene))

cat("\nGenes exclusivos del análisis de sensibilidad:\n")
print(sum(!sensibilidad$gene %in% comparacion$gene))

if (nrow(comparacion) == 0) {
  stop("ERROR: no existen genes comunes.")
}

if (anyNA(comparacion$logFC_primario) ||
    anyNA(comparacion$logFC_sensibilidad)) {
  stop("ERROR: hay NA en los logFC después del emparejamiento.")
}


# ------------------------------------------------------------
# 4. Correlación de tamaños de efecto
# ------------------------------------------------------------

pearson_logFC <- cor(
  comparacion$logFC_primario,
  comparacion$logFC_sensibilidad,
  method = "pearson"
)

spearman_logFC <- cor(
  comparacion$logFC_primario,
  comparacion$logFC_sensibilidad,
  method = "spearman"
)

cat("\nCorrelación Pearson de logFC:\n")
print(pearson_logFC)

cat("\nCorrelación Spearman de logFC:\n")
print(spearman_logFC)


# ------------------------------------------------------------
# 5. Diferencia absoluta de logFC
# ------------------------------------------------------------

comparacion$diferencia_logFC <- abs(
  comparacion$logFC_primario -
  comparacion$logFC_sensibilidad
)

cat("\nResumen de |diferencia de logFC|:\n")
print(summary(comparacion$diferencia_logFC))


# ------------------------------------------------------------
# 6. Concordancia de dirección
#
# Se excluyen solamente logFC exactamente iguales a cero
# para este cálculo.
# ------------------------------------------------------------

no_cero <- comparacion[
  comparacion$logFC_primario != 0 &
  comparacion$logFC_sensibilidad != 0,
  ,
  drop = FALSE
]

misma_direccion <- sign(no_cero$logFC_primario) ==
                   sign(no_cero$logFC_sensibilidad)

cat("\nGenes evaluables para concordancia de dirección:\n")
print(nrow(no_cero))

cat("\nGenes con la misma dirección de efecto:\n")
print(sum(misma_direccion))

cat("\nProporción con la misma dirección de efecto:\n")
print(mean(misma_direccion))


# ------------------------------------------------------------
# 7. Estabilidad de genes con
# FDR < 0.05 y |logFC| >= 1
# ------------------------------------------------------------

comparacion$criterio_primario <-
  comparacion$FDR_primario < 0.05 &
  abs(comparacion$logFC_primario) >= 1

comparacion$criterio_sensibilidad <-
  comparacion$FDR_sensibilidad < 0.05 &
  abs(comparacion$logFC_sensibilidad) >= 1

ambos <- comparacion$criterio_primario &
         comparacion$criterio_sensibilidad

solo_primario <- comparacion$criterio_primario &
                 !comparacion$criterio_sensibilidad

solo_sensibilidad <- !comparacion$criterio_primario &
                     comparacion$criterio_sensibilidad

ninguno <- !comparacion$criterio_primario &
           !comparacion$criterio_sensibilidad

cat("\nGenes comunes que cumplen criterio en primario:\n")
print(sum(comparacion$criterio_primario))

cat("\nGenes comunes que cumplen criterio en sensibilidad:\n")
print(sum(comparacion$criterio_sensibilidad))

cat("\nCumplen criterio en ambos:\n")
print(sum(ambos))

cat("\nCumplen criterio sólo en primario:\n")
print(sum(solo_primario))

cat("\nCumplen criterio sólo en sensibilidad:\n")
print(sum(solo_sensibilidad))

cat("\nNo cumplen criterio en ninguno:\n")
print(sum(ninguno))


# ------------------------------------------------------------
# 8. Concordancia direccional dentro de los genes
# que cumplen el criterio en ambos análisis
# ------------------------------------------------------------

comparacion_ambos <- comparacion[
  ambos,
  ,
  drop = FALSE
]

if (nrow(comparacion_ambos) > 0) {

  direccion_ambos <-
    sign(comparacion_ambos$logFC_primario) ==
    sign(comparacion_ambos$logFC_sensibilidad)

  cat("\nGenes que cumplen criterio en ambos y conservan dirección:\n")
  print(sum(direccion_ambos))

  cat("\nProporción direccional entre los que cumplen en ambos:\n")
  print(mean(direccion_ambos))

} else {

  cat("\nNo hay genes que cumplan el criterio en ambos análisis.\n")
}


# ------------------------------------------------------------
# 9. Índice de Jaccard para los conjuntos que cumplen criterio
# ------------------------------------------------------------

union_criterio <- comparacion$criterio_primario |
                  comparacion$criterio_sensibilidad

if (sum(union_criterio) > 0) {

  jaccard <- sum(ambos) / sum(union_criterio)

  cat("\nÍndice de Jaccard de los conjuntos con criterio:\n")
  print(jaccard)

}


# ------------------------------------------------------------
# 10. Guardar tabla completa de comparación
# ------------------------------------------------------------

comparacion <- comparacion[
  order(comparacion$PValue_primario),
  ,
  drop = FALSE
]

rownames(comparacion) <- NULL

write.csv(
  comparacion,
  "concordancia_primario_vs_sensibilidad.csv",
  row.names = FALSE
)

cat("\nArchivo generado:\n")
cat("- concordancia_primario_vs_sensibilidad.csv\n")

cat("\nRESULTADO: comparación cuantitativa completada.\n")