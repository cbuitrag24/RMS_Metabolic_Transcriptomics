# ============================================================
# 31_restringir_rankings_FN_FOXO1FP_universo_9428.R
#
# VERSIÓN CORREGIDA
#
# OBJETIVO:
# Construir un universo estrictamente comparable para:
#
#   FN vs Normal
#   FOXO1-FP vs Normal
#
# partiendo del universo común de 9.428 genes utilizado para
# la replicación GSE108022 vs GSE28511 y conservando únicamente
# los genes efectivamente disponibles en AMBOS análisis de
# subgrupos.
#
# IMPORTANTE:
# - El tamaño final NO se impone.
# - Se calcula por intersección.
# - No se imputan genes ausentes.
# - No se reutilizan estadísticas de otros contrastes.
# - No se selecciona por logFC, FDR, F ni dirección.
# - ranking_stat no se modifica.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_FN <-
  "ranking_GSE108022_FN_vs_Normal.rds"

archivo_FP <-
  "ranking_GSE108022_FOXO1FP_vs_Normal.rds"

archivo_universo_original <-
  file.path(
    "validacion_GSE28511",
    "ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.rds"
  )


archivos <- c(
  archivo_FN,
  archivo_FP,
  archivo_universo_original
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
# 2. Cargar
# ------------------------------------------------------------

ranking_FN <- readRDS(
  archivo_FN
)

ranking_FP <- readRDS(
  archivo_FP
)

universo_original <- readRDS(
  archivo_universo_original
)


cat("\nDimensiones ranking FN:\n")
print(dim(ranking_FN))

cat("\nDimensiones ranking FOXO1-FP:\n")
print(dim(ranking_FP))

cat("\nDimensiones universo original:\n")
print(dim(universo_original))


# ------------------------------------------------------------
# 3. Verificar estructura
# ------------------------------------------------------------

columnas_necesarias <- c(
  "gene",
  "logFC",
  "F",
  "ranking_stat"
)


if (!all(
  columnas_necesarias %in% colnames(ranking_FN)
)) {
  stop(
    "ERROR: estructura inesperada del ranking FN."
  )
}


if (!all(
  columnas_necesarias %in% colnames(ranking_FP)
)) {
  stop(
    "ERROR: estructura inesperada del ranking FOXO1-FP."
  )
}


if (!("gene" %in% colnames(universo_original))) {
  stop(
    "ERROR: el universo original no contiene gene."
  )
}


# ------------------------------------------------------------
# 4. Normalizar identificadores
# ------------------------------------------------------------

ranking_FN$gene <- trimws(
  as.character(ranking_FN$gene)
)

ranking_FP$gene <- trimws(
  as.character(ranking_FP$gene)
)

universo_original$gene <- trimws(
  as.character(universo_original$gene)
)


# ------------------------------------------------------------
# 5. Verificar identidad y unicidad
# ------------------------------------------------------------

if (
  anyNA(ranking_FN$gene) ||
  any(ranking_FN$gene == "")
) {
  stop(
    "ERROR: existen genes vacíos en FN."
  )
}


if (
  anyNA(ranking_FP$gene) ||
  any(ranking_FP$gene == "")
) {
  stop(
    "ERROR: existen genes vacíos en FOXO1-FP."
  )
}


if (
  anyNA(universo_original$gene) ||
  any(universo_original$gene == "")
) {
  stop(
    "ERROR: existen genes vacíos en el universo original."
  )
}


if (anyDuplicated(ranking_FN$gene)) {
  stop(
    "ERROR: existen genes duplicados en FN."
  )
}


if (anyDuplicated(ranking_FP$gene)) {
  stop(
    "ERROR: existen genes duplicados en FOXO1-FP."
  )
}


if (anyDuplicated(universo_original$gene)) {
  stop(
    "ERROR: existen genes duplicados en el universo original."
  )
}


if (nrow(universo_original) != 9428) {
  stop(
    paste(
      "ERROR: el universo original debería contener 9428 genes;",
      "se encontraron",
      nrow(universo_original)
    )
  )
}


# ------------------------------------------------------------
# 6. Auditar cobertura por separado
# ------------------------------------------------------------

genes_originales <- universo_original$gene


ausentes_FN <- setdiff(
  genes_originales,
  ranking_FN$gene
)

ausentes_FP <- setdiff(
  genes_originales,
  ranking_FP$gene
)


cat("\n====================================================\n")
cat("AUDITORÍA DEL UNIVERSO ORIGINAL\n")
cat("====================================================\n")


cat("\nGenes del universo original:\n")
cat(
  length(genes_originales),
  "\n"
)


cat("\nAusentes en FN:\n")
cat(
  length(ausentes_FN),
  "\n"
)


if (length(ausentes_FN) > 0) {
  print(ausentes_FN)
}


cat("\nAusentes en FOXO1-FP:\n")
cat(
  length(ausentes_FP),
  "\n"
)


if (length(ausentes_FP) > 0) {
  print(ausentes_FP)
}


# ------------------------------------------------------------
# 7. Construir intersección REAL
#
# Un gen se conserva solamente si:
#
# 1. pertenece al universo original de replicación
# 2. está disponible en FN vs Normal
# 3. está disponible en FOXO1-FP vs Normal
#
# Se conserva el orden del universo original.
# ------------------------------------------------------------

disponible_ambos <- (
  genes_originales %in% ranking_FN$gene &
  genes_originales %in% ranking_FP$gene
)


universo_subgrupos <- genes_originales[
  disponible_ambos
]


excluidos_universo <- genes_originales[
  !disponible_ambos
]


cat("\n====================================================\n")
cat("UNIVERSO COMPARABLE DE SUBGRUPOS\n")
cat("====================================================\n")


cat("\nGenes conservados:\n")
cat(
  length(universo_subgrupos),
  "\n"
)


cat("\nGenes excluidos respecto del universo original:\n")
cat(
  length(excluidos_universo),
  "\n"
)


if (length(excluidos_universo) > 0) {

  cat("\nIdentidad de genes excluidos:\n")

  print(
    excluidos_universo
  )
}


# ------------------------------------------------------------
# 8. Verificación lógica de la intersección
# ------------------------------------------------------------

if (!all(
  universo_subgrupos %in% ranking_FN$gene
)) {
  stop(
    "ERROR: el universo final contiene genes ausentes en FN."
  )
}


if (!all(
  universo_subgrupos %in% ranking_FP$gene
)) {
  stop(
    "ERROR: el universo final contiene genes ausentes en FOXO1-FP."
  )
}


if (anyDuplicated(universo_subgrupos)) {
  stop(
    "ERROR: el universo final contiene genes duplicados."
  )
}


# ------------------------------------------------------------
# 9. Extraer ambos rankings en exactamente el mismo orden
# ------------------------------------------------------------

idx_FN <- match(
  universo_subgrupos,
  ranking_FN$gene
)

idx_FP <- match(
  universo_subgrupos,
  ranking_FP$gene
)


if (
  anyNA(idx_FN) ||
  anyNA(idx_FP)
) {
  stop(
    "ERROR: match produjo NA inesperados."
  )
}


FN_final <- ranking_FN[
  idx_FN,
  ,
  drop = FALSE
]


FP_final <- ranking_FP[
  idx_FP,
  ,
  drop = FALSE
]


rownames(FN_final) <- NULL
rownames(FP_final) <- NULL


# ------------------------------------------------------------
# 10. Verificar identidad exacta y orden
# ------------------------------------------------------------

if (!identical(
  FN_final$gene,
  universo_subgrupos
)) {
  stop(
    "ERROR: FN no conserva exactamente el universo final."
  )
}


if (!identical(
  FP_final$gene,
  universo_subgrupos
)) {
  stop(
    "ERROR: FOXO1-FP no conserva exactamente el universo final."
  )
}


if (!identical(
  FN_final$gene,
  FP_final$gene
)) {
  stop(
    "ERROR: los rankings finales no tienen idénticos genes y orden."
  )
}


# ------------------------------------------------------------
# 11. Verificar que ranking_stat no fue alterado
# ------------------------------------------------------------

original_FN_recuperado <- ranking_FN$ranking_stat[
  match(
    FN_final$gene,
    ranking_FN$gene
  )
]


original_FP_recuperado <- ranking_FP$ranking_stat[
  match(
    FP_final$gene,
    ranking_FP$gene
  )
]


error_FN <- max(
  abs(
    FN_final$ranking_stat -
    original_FN_recuperado
  )
)


error_FP <- max(
  abs(
    FP_final$ranking_stat -
    original_FP_recuperado
  )
)


cat("\nError máximo ranking_stat FN:\n")
cat(
  error_FN,
  "\n"
)


cat("\nError máximo ranking_stat FOXO1-FP:\n")
cat(
  error_FP,
  "\n"
)


if (
  error_FN != 0 ||
  error_FP != 0
) {
  stop(
    "ERROR: se alteró algún ranking_stat."
  )
}


# ------------------------------------------------------------
# 12. Resumen FN
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("FN vs NORMAL - UNIVERSO FINAL\n")
cat("====================================================\n")


cat("\nGenes:\n")
cat(
  nrow(FN_final),
  "\n"
)


cat("\nResumen ranking_stat:\n")
print(
  summary(
    FN_final$ranking_stat
  )
)


cat("\nPositivos:\n")
cat(
  sum(FN_final$ranking_stat > 0),
  "\n"
)


cat("\nNegativos:\n")
cat(
  sum(FN_final$ranking_stat < 0),
  "\n"
)


cat("\nCeros:\n")
cat(
  sum(FN_final$ranking_stat == 0),
  "\n"
)


# ------------------------------------------------------------
# 13. Resumen FOXO1-FP
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("FOXO1-FP vs NORMAL - UNIVERSO FINAL\n")
cat("====================================================\n")


cat("\nGenes:\n")
cat(
  nrow(FP_final),
  "\n"
)


cat("\nResumen ranking_stat:\n")
print(
  summary(
    FP_final$ranking_stat
  )
)


cat("\nPositivos:\n")
cat(
  sum(FP_final$ranking_stat > 0),
  "\n"
)


cat("\nNegativos:\n")
cat(
  sum(FP_final$ranking_stat < 0),
  "\n"
)


cat("\nCeros:\n")
cat(
  sum(FP_final$ranking_stat == 0),
  "\n"
)


# ------------------------------------------------------------
# 14. Concordancia descriptiva en el universo FINAL
# ------------------------------------------------------------

pearson_logFC <- cor(
  FN_final$logFC,
  FP_final$logFC,
  method = "pearson"
)


spearman_logFC <- cor(
  FN_final$logFC,
  FP_final$logFC,
  method = "spearman"
)


misma_direccion <- mean(
  sign(FN_final$logFC) ==
  sign(FP_final$logFC)
)


pearson_ranking <- cor(
  FN_final$ranking_stat,
  FP_final$ranking_stat,
  method = "pearson"
)


spearman_ranking <- cor(
  FN_final$ranking_stat,
  FP_final$ranking_stat,
  method = "spearman"
)


cat("\n====================================================\n")
cat("CONCORDANCIA EN EL UNIVERSO FINAL\n")
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
  "Misma dirección logFC:",
  misma_direccion,
  "\n"
)


cat(
  "Pearson ranking_stat:",
  pearson_ranking,
  "\n"
)


cat(
  "Spearman ranking_stat:",
  spearman_ranking,
  "\n"
)


# ------------------------------------------------------------
# 15. Tabla explícita de genes excluidos
# ------------------------------------------------------------

tabla_excluidos <- data.frame(

  gene =
    excluidos_universo,

  presente_universo_original =
    TRUE,

  presente_FN =
    excluidos_universo %in%
    ranking_FN$gene,

  presente_FOXO1FP =
    excluidos_universo %in%
    ranking_FP$gene,

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 16. Auditoría
# ------------------------------------------------------------

auditoria <- data.frame(

  medida = c(
    "n_universo_original",
    "n_universo_final_subgrupos",
    "n_excluidos",
    "n_ausentes_FN",
    "n_ausentes_FOXO1FP",
    "error_max_ranking_FN",
    "error_max_ranking_FOXO1FP",
    "pearson_logFC",
    "spearman_logFC",
    "misma_direccion_logFC",
    "pearson_ranking_stat",
    "spearman_ranking_stat"
  ),

  valor = c(
    length(genes_originales),
    length(universo_subgrupos),
    length(excluidos_universo),
    length(ausentes_FN),
    length(ausentes_FP),
    error_FN,
    error_FP,
    pearson_logFC,
    spearman_logFC,
    misma_direccion,
    pearson_ranking,
    spearman_ranking
  ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 17. Guardar
#
# El nombre usa "universo_comparable" porque el tamaño se
# determina empíricamente y no se fuerza a priori.
# ------------------------------------------------------------

write.csv(
  FN_final,
  "ranking_GSE108022_FN_vs_Normal_universo_comparable.csv",
  row.names = FALSE
)


saveRDS(
  FN_final,
  "ranking_GSE108022_FN_vs_Normal_universo_comparable.rds"
)


write.csv(
  FP_final,
  "ranking_GSE108022_FOXO1FP_vs_Normal_universo_comparable.csv",
  row.names = FALSE
)


saveRDS(
  FP_final,
  "ranking_GSE108022_FOXO1FP_vs_Normal_universo_comparable.rds"
)


write.csv(
  tabla_excluidos,
  "genes_excluidos_universo_subgrupos.csv",
  row.names = FALSE
)


write.csv(
  auditoria,
  "auditoria_universo_comparable_FN_FOXO1FP.csv",
  row.names = FALSE
)


saveRDS(
  list(
    universo_original = genes_originales,
    universo_final = universo_subgrupos,
    excluidos = tabla_excluidos,
    FN = FN_final,
    FOXO1FP = FP_final,
    auditoria = auditoria
  ),
  "universo_comparable_FN_FOXO1FP.rds"
)


cat("\nArchivos generados:\n")

cat(
  "- ranking_GSE108022_FN_vs_Normal_universo_comparable.csv/.rds\n"
)

cat(
  "- ranking_GSE108022_FOXO1FP_vs_Normal_universo_comparable.csv/.rds\n"
)

cat(
  "- genes_excluidos_universo_subgrupos.csv\n"
)

cat(
  "- auditoria_universo_comparable_FN_FOXO1FP.csv\n"
)

cat(
  "- universo_comparable_FN_FOXO1FP.rds\n"
)


cat("\n====================================================\n")
cat("UNIVERSO COMPARABLE CONSTRUIDO Y AUDITADO\n")
cat("====================================================\n")