# ============================================================
# 30_preparar_rankings_FN_FOXO1FP_vs_Normal.R
#
# VERSIÓN CORREGIDA
#
# OBJETIVO:
# Construir rankings prerankeados para:
#
#   1) FN vs Normal
#   2) FOXO1-FP vs Normal
#
# utilizando las tablas DE previamente auditadas, que
# contienen los símbolos génicos correctos.
#
# Estadístico:
#   sign(logFC) * sqrt(F)
#
# ranking > 0 : mayor expresión en RMS
# ranking < 0 : mayor expresión en músculo normal
# ============================================================


# ------------------------------------------------------------
# 1. Archivos de entrada
# ------------------------------------------------------------

archivo_FN <-
  "GPL16791_DE_FN_vs_Normal_todos.csv"

archivo_FP <-
  "GPL16791_DE_FOXO1FP_vs_Normal_todos.csv"


archivos <- c(
  archivo_FN,
  archivo_FP
)


faltantes <- archivos[
  !file.exists(archivos)
]


if (length(faltantes) > 0) {

  cat("\nArchivos faltantes:\n")
  print(faltantes)

  stop("ERROR: faltan archivos necesarios.")
}


# ------------------------------------------------------------
# 2. Cargar tablas DE auditadas
# ------------------------------------------------------------

tabla_FN <- read.csv(
  archivo_FN,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


tabla_FP <- read.csv(
  archivo_FP,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\nDimensiones tabla FN vs Normal:\n")
print(dim(tabla_FN))

cat("\nDimensiones tabla FOXO1-FP vs Normal:\n")
print(dim(tabla_FP))


cat("\nColumnas FN:\n")
print(colnames(tabla_FN))

cat("\nColumnas FOXO1-FP:\n")
print(colnames(tabla_FP))


# ------------------------------------------------------------
# 3. Verificar columnas necesarias
# ------------------------------------------------------------

columnas_necesarias <- c(
  "gene",
  "logFC",
  "F"
)


if (!all(
  columnas_necesarias %in% colnames(tabla_FN)
)) {
  stop(
    "ERROR: FN no contiene gene, logFC y F."
  )
}


if (!all(
  columnas_necesarias %in% colnames(tabla_FP)
)) {
  stop(
    "ERROR: FOXO1-FP no contiene gene, logFC y F."
  )
}


# ------------------------------------------------------------
# 4. Limpiar símbolos génicos
# ------------------------------------------------------------

tabla_FN$gene <- trimws(
  as.character(tabla_FN$gene)
)


tabla_FP$gene <- trimws(
  as.character(tabla_FP$gene)
)


# ------------------------------------------------------------
# 5. Auditoría de identidad
# ------------------------------------------------------------

if (
  anyNA(tabla_FN$gene) ||
  any(tabla_FN$gene == "")
) {
  stop(
    "ERROR: existen genes vacíos en FN."
  )
}


if (
  anyNA(tabla_FP$gene) ||
  any(tabla_FP$gene == "")
) {
  stop(
    "ERROR: existen genes vacíos en FOXO1-FP."
  )
}


if (anyDuplicated(tabla_FN$gene)) {
  stop(
    "ERROR: existen genes duplicados en FN."
  )
}


if (anyDuplicated(tabla_FP$gene)) {
  stop(
    "ERROR: existen genes duplicados en FOXO1-FP."
  )
}


cat("\nPrimeros símbolos FN:\n")
print(
  head(
    tabla_FN$gene,
    10
  )
)


cat("\nPrimeros símbolos FOXO1-FP:\n")
print(
  head(
    tabla_FP$gene,
    10
  )
)


# ------------------------------------------------------------
# 6. Verificar estadísticas
# ------------------------------------------------------------

if (
  anyNA(tabla_FN$logFC) ||
  anyNA(tabla_FN$F) ||
  anyNA(tabla_FP$logFC) ||
  anyNA(tabla_FP$F)
) {
  stop(
    "ERROR: existen NA en logFC o F."
  )
}


if (
  any(!is.finite(tabla_FN$logFC)) ||
  any(!is.finite(tabla_FN$F)) ||
  any(!is.finite(tabla_FP$logFC)) ||
  any(!is.finite(tabla_FP$F))
) {
  stop(
    "ERROR: existen estadísticas no finitas."
  )
}


if (
  any(tabla_FN$F < 0) ||
  any(tabla_FP$F < 0)
) {
  stop(
    "ERROR: existen valores F negativos."
  )
}


# ------------------------------------------------------------
# 7. Construir rankings
# ------------------------------------------------------------

ranking_FN <- data.frame(

  gene =
    tabla_FN$gene,

  logFC =
    tabla_FN$logFC,

  F =
    tabla_FN$F,

  ranking_stat =
    sign(tabla_FN$logFC) *
    sqrt(tabla_FN$F),

  stringsAsFactors = FALSE
)


ranking_FP <- data.frame(

  gene =
    tabla_FP$gene,

  logFC =
    tabla_FP$logFC,

  F =
    tabla_FP$F,

  ranking_stat =
    sign(tabla_FP$logFC) *
    sqrt(tabla_FP$F),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. Verificar signo
# ------------------------------------------------------------

if (any(
  sign(ranking_FN$ranking_stat) !=
  sign(ranking_FN$logFC)
)) {
  stop(
    "ERROR: inconsistencia de signo en FN."
  )
}


if (any(
  sign(ranking_FP$ranking_stat) !=
  sign(ranking_FP$logFC)
)) {
  stop(
    "ERROR: inconsistencia de signo en FOXO1-FP."
  )
}


# ------------------------------------------------------------
# 9. Verificar fórmula del ranking
# ------------------------------------------------------------

error_FN <- max(
  abs(
    ranking_FN$ranking_stat -
    sign(ranking_FN$logFC) *
    sqrt(ranking_FN$F)
  )
)


error_FP <- max(
  abs(
    ranking_FP$ranking_stat -
    sign(ranking_FP$logFC) *
    sqrt(ranking_FP$F)
  )
)


cat("\nError máximo fórmula ranking FN:\n")
cat(error_FN, "\n")


cat("\nError máximo fórmula ranking FOXO1-FP:\n")
cat(error_FP, "\n")


if (
  error_FN != 0 ||
  error_FP != 0
) {
  stop(
    "ERROR: inconsistencia en el cálculo del ranking."
  )
}


# ------------------------------------------------------------
# 10. Ordenar rankings
# ------------------------------------------------------------

ranking_FN <- ranking_FN[
  order(
    ranking_FN$ranking_stat,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


ranking_FP <- ranking_FP[
  order(
    ranking_FP$ranking_stat,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


rownames(ranking_FN) <- NULL
rownames(ranking_FP) <- NULL


# ------------------------------------------------------------
# 11. Comparabilidad de universos
# ------------------------------------------------------------

comunes <- intersect(
  ranking_FN$gene,
  ranking_FP$gene
)


solo_FN <- setdiff(
  ranking_FN$gene,
  ranking_FP$gene
)


solo_FP <- setdiff(
  ranking_FP$gene,
  ranking_FN$gene
)


cat("\n====================================================\n")
cat("COMPARABILIDAD DE LOS RANKINGS CORREGIDOS\n")
cat("====================================================\n")


cat("\nGenes FN:\n")
cat(nrow(ranking_FN), "\n")


cat("\nGenes FOXO1-FP:\n")
cat(nrow(ranking_FP), "\n")


cat("\nGenes comunes:\n")
cat(length(comunes), "\n")


cat("\nGenes sólo FN:\n")
cat(length(solo_FN), "\n")


cat("\nGenes sólo FOXO1-FP:\n")
cat(length(solo_FP), "\n")


# ------------------------------------------------------------
# 12. Resumen FN
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RANKING FN vs NORMAL CORREGIDO\n")
cat("====================================================\n")


cat("\nResumen ranking_stat:\n")
print(
  summary(
    ranking_FN$ranking_stat
  )
)


cat("\nPositivos:\n")
cat(
  sum(ranking_FN$ranking_stat > 0),
  "\n"
)


cat("\nNegativos:\n")
cat(
  sum(ranking_FN$ranking_stat < 0),
  "\n"
)


cat("\nCeros:\n")
cat(
  sum(ranking_FN$ranking_stat == 0),
  "\n"
)


# ------------------------------------------------------------
# 13. Resumen FOXO1-FP
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RANKING FOXO1-FP vs NORMAL CORREGIDO\n")
cat("====================================================\n")


cat("\nResumen ranking_stat:\n")
print(
  summary(
    ranking_FP$ranking_stat
  )
)


cat("\nPositivos:\n")
cat(
  sum(ranking_FP$ranking_stat > 0),
  "\n"
)


cat("\nNegativos:\n")
cat(
  sum(ranking_FP$ranking_stat < 0),
  "\n"
)


cat("\nCeros:\n")
cat(
  sum(ranking_FP$ranking_stat == 0),
  "\n"
)


# ------------------------------------------------------------
# 14. Concordancia descriptiva
# ------------------------------------------------------------

comparacion <- merge(

  ranking_FN[
    ,
    c(
      "gene",
      "logFC",
      "ranking_stat"
    )
  ],

  ranking_FP[
    ,
    c(
      "gene",
      "logFC",
      "ranking_stat"
    )
  ],

  by = "gene",

  suffixes = c(
    "_FN",
    "_FOXO1FP"
  )
)


pearson_logFC <- cor(
  comparacion$logFC_FN,
  comparacion$logFC_FOXO1FP,
  method = "pearson"
)


spearman_logFC <- cor(
  comparacion$logFC_FN,
  comparacion$logFC_FOXO1FP,
  method = "spearman"
)


misma_direccion <- mean(
  sign(comparacion$logFC_FN) ==
  sign(comparacion$logFC_FOXO1FP)
)


cat("\n====================================================\n")
cat("CONCORDANCIA DESCRIPTIVA CORREGIDA\n")
cat("====================================================\n")


cat(
  "\nPearson logFC:",
  pearson_logFC,
  "\n"
)


cat(
  "Spearman logFC:",
  spearman_logFC,
  "\n"
)


cat(
  "Misma dirección:",
  misma_direccion,
  "\n"
)


# ------------------------------------------------------------
# 15. Guardar
#
# Sobrescribe los rankings incorrectos generados anteriormente.
# ------------------------------------------------------------

write.csv(
  ranking_FN,
  "ranking_GSE108022_FN_vs_Normal.csv",
  row.names = FALSE
)


saveRDS(
  ranking_FN,
  "ranking_GSE108022_FN_vs_Normal.rds"
)


write.csv(
  ranking_FP,
  "ranking_GSE108022_FOXO1FP_vs_Normal.csv",
  row.names = FALSE
)


saveRDS(
  ranking_FP,
  "ranking_GSE108022_FOXO1FP_vs_Normal.rds"
)


write.csv(
  comparacion,
  "concordancia_rankings_FN_FOXO1FP_vs_Normal.csv",
  row.names = FALSE
)


cat("\nArchivos corregidos generados:\n")

cat(
  "- ranking_GSE108022_FN_vs_Normal.csv/.rds\n"
)

cat(
  "- ranking_GSE108022_FOXO1FP_vs_Normal.csv/.rds\n"
)

cat(
  "- concordancia_rankings_FN_FOXO1FP_vs_Normal.csv\n"
)


cat("\n====================================================\n")
cat("RANKINGS CORREGIDOS COMPLETADOS\n")
cat("====================================================\n")