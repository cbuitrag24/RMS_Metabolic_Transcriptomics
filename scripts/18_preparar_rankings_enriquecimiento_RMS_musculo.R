# ============================================================
# 18_preparar_rankings_enriquecimiento_RMS_musculo.R
#
# OBJETIVO:
# Preparar rankings independientes para análisis funcional
# prerankeado de RMS vs músculo.
#
# GSE108022 (RNA-seq / edgeR):
#   ranking = sign(logFC) * sqrt(F)
#
# GSE28511 (microarray / limma):
#   ranking = estadístico t moderado
#
# DIRECCIÓN:
#   valor positivo = mayor en RMS
#   valor negativo = mayor en músculo
#
# IMPORTANTE:
# - NO se aplica ningún corte de FDR.
# - NO se seleccionan genes según concordancia entre estudios.
# - NO se realiza todavía enriquecimiento funcional.
# - GSE28511 principal: genes con UNA sola sonda.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_rna <- file.path(
  "..",
  "GPL16791_DE_RMS_vs_Normal_todos_los_genes_CORREGIDO.csv"
)

archivo_micro <-
  "GSE28511_primario_limma_todas_las_sondas.rds"

archivo_universo <-
  "universo_comun_GSE108022_GSE28511.csv"


archivos <- c(
  archivo_rna,
  archivo_micro,
  archivo_universo
)

faltantes <- archivos[!file.exists(archivos)]

if (length(faltantes) > 0) {
  stop(
    paste(
      "ERROR: faltan archivos:",
      paste(faltantes, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

rna <- read.csv(
  archivo_rna,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

micro <- readRDS(
  archivo_micro
)

universo <- read.csv(
  archivo_universo,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("Dimensiones GSE108022:\n")
print(dim(rna))

cat("\nDimensiones GSE28511:\n")
print(dim(micro))

cat("\nDimensiones universo común:\n")
print(dim(universo))


# ------------------------------------------------------------
# 3. Verificar columnas
# ------------------------------------------------------------

necesarias_rna <- c(
  "gene",
  "logFC",
  "F",
  "PValue",
  "FDR"
)

necesarias_micro <- c(
  "Probe_ID",
  "Symbol",
  "logFC",
  "t",
  "PValue",
  "FDR"
)

necesarias_universo <- c(
  "gene",
  "n_sondas_GSE28511"
)


if (!all(necesarias_rna %in% colnames(rna))) {
  stop(
    paste(
      "ERROR: faltan columnas en GSE108022:",
      paste(
        setdiff(necesarias_rna, colnames(rna)),
        collapse = ", "
      )
    )
  )
}

if (!all(necesarias_micro %in% colnames(micro))) {
  stop(
    paste(
      "ERROR: faltan columnas en GSE28511:",
      paste(
        setdiff(necesarias_micro, colnames(micro)),
        collapse = ", "
      )
    )
  )
}

if (!all(necesarias_universo %in% colnames(universo))) {
  stop(
    "ERROR: faltan columnas necesarias en el universo común."
  )
}


# ------------------------------------------------------------
# 4. Limpiar identificadores
# ------------------------------------------------------------

rna$gene <- trimws(
  as.character(rna$gene)
)

micro$Symbol <- trimws(
  as.character(micro$Symbol)
)

micro$Symbol[is.na(micro$Symbol)] <- ""

universo$gene <- trimws(
  as.character(universo$gene)
)


# ------------------------------------------------------------
# 5. Controles de GSE108022
# ------------------------------------------------------------

if (anyDuplicated(rna$gene)) {
  stop(
    "ERROR: existen genes duplicados en GSE108022."
  )
}

if (
  anyNA(rna$logFC) ||
  anyNA(rna$F) ||
  any(!is.finite(rna$logFC)) ||
  any(!is.finite(rna$F))
) {
  stop(
    "ERROR: logFC o F inválidos en GSE108022."
  )
}

if (any(rna$F < 0)) {
  stop(
    "ERROR: existen estadísticos F negativos."
  )
}


# ------------------------------------------------------------
# 6. Ranking independiente GSE108022
#
# Para un contraste de 1 grado de libertad:
# sqrt(F) aporta la magnitud del estadístico.
# El signo se recupera del logFC.
# ------------------------------------------------------------

rna$ranking_stat <- (
  sign(rna$logFC) *
  sqrt(rna$F)
)


if (
  anyNA(rna$ranking_stat) ||
  any(!is.finite(rna$ranking_stat))
) {
  stop(
    "ERROR: ranking inválido en GSE108022."
  )
}


ranking_rna <- data.frame(
  gene = rna$gene,
  logFC = as.numeric(rna$logFC),
  F = as.numeric(rna$F),
  ranking_stat = as.numeric(rna$ranking_stat),
  stringsAsFactors = FALSE
)

ranking_rna <- ranking_rna[
  order(
    ranking_rna$ranking_stat,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

rownames(ranking_rna) <- NULL


# ------------------------------------------------------------
# 7. Universo principal GSE28511:
# genes comunes con una única sonda
# ------------------------------------------------------------

genes_sonda_unica <- universo$gene[
  universo$n_sondas_GSE28511 == 1
]


if (length(genes_sonda_unica) != 9428) {
  stop(
    paste(
      "ERROR: se esperaban 9428 genes de sonda única y se encontraron",
      length(genes_sonda_unica)
    )
  )
}


idx_micro <- match(
  genes_sonda_unica,
  micro$Symbol
)

if (anyNA(idx_micro)) {
  stop(
    "ERROR: faltan genes de sonda única en GSE28511."
  )
}


micro_u <- micro[
  idx_micro,
  ,
  drop = FALSE
]


if (!identical(
  micro_u$Symbol,
  genes_sonda_unica
)) {
  stop(
    "ERROR: identidad gen-sonda incorrecta en GSE28511."
  )
}

if (anyDuplicated(micro_u$Symbol)) {
  stop(
    "ERROR: existen Symbols duplicados entre genes de sonda única."
  )
}


# ------------------------------------------------------------
# 8. Controles del estadístico t
# ------------------------------------------------------------

if (
  anyNA(micro_u$t) ||
  any(!is.finite(micro_u$t))
) {
  stop(
    "ERROR: estadísticos t inválidos en GSE28511."
  )
}


# Verificar que el signo de t coincide con el logFC.
# Para este contraste deberían tener exactamente
# la misma dirección.

signos_validos <- (
  sign(micro_u$t) ==
  sign(micro_u$logFC)
)

if (!all(signos_validos)) {
  stop(
    paste(
      "ERROR:",
      sum(!signos_validos),
      "sondas tienen signo de t distinto del logFC."
    )
  )
}


cat(
  "\nVERIFICACIÓN: el signo del estadístico t",
  "coincide con el logFC en los 9428 genes.\n"
)


# ------------------------------------------------------------
# 9. Ranking independiente GSE28511
# ------------------------------------------------------------

ranking_micro <- data.frame(
  gene =
    micro_u$Symbol,

  Probe_ID =
    as.character(micro_u$Probe_ID),

  logFC =
    as.numeric(micro_u$logFC),

  t =
    as.numeric(micro_u$t),

  ranking_stat =
    as.numeric(micro_u$t),

  stringsAsFactors = FALSE
)


ranking_micro <- ranking_micro[
  order(
    ranking_micro$ranking_stat,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

rownames(ranking_micro) <- NULL


# ------------------------------------------------------------
# 10. Resumen de rankings
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RANKING GSE108022\n")
cat("====================================================\n")

cat("\nGenes:\n")
cat(nrow(ranking_rna), "\n")

cat("\nResumen del estadístico firmado:\n")
print(
  summary(
    ranking_rna$ranking_stat
  )
)

cat("\nGenes con ranking > 0:\n")
cat(
  sum(ranking_rna$ranking_stat > 0),
  "\n"
)

cat("\nGenes con ranking < 0:\n")
cat(
  sum(ranking_rna$ranking_stat < 0),
  "\n"
)

cat("\nGenes con ranking = 0:\n")
cat(
  sum(ranking_rna$ranking_stat == 0),
  "\n"
)


cat("\n====================================================\n")
cat("RANKING GSE28511 - SONDA ÚNICA\n")
cat("====================================================\n")

cat("\nGenes:\n")
cat(nrow(ranking_micro), "\n")

cat("\nResumen del estadístico t:\n")
print(
  summary(
    ranking_micro$ranking_stat
  )
)

cat("\nGenes con ranking > 0:\n")
cat(
  sum(ranking_micro$ranking_stat > 0),
  "\n"
)

cat("\nGenes con ranking < 0:\n")
cat(
  sum(ranking_micro$ranking_stat < 0),
  "\n"
)

cat("\nGenes con ranking = 0:\n")
cat(
  sum(ranking_micro$ranking_stat == 0),
  "\n"
)


# ------------------------------------------------------------
# 11. Verificar orientación biológica
#
# Positivo debe significar mayor en RMS en ambos estudios.
# Negativo debe significar mayor en músculo.
# ------------------------------------------------------------

if (!all(
  sign(ranking_rna$ranking_stat) ==
  sign(ranking_rna$logFC)
)) {
  stop(
    "ERROR: orientación incorrecta del ranking GSE108022."
  )
}

if (!all(
  sign(ranking_micro$ranking_stat) ==
  sign(ranking_micro$logFC)
)) {
  stop(
    "ERROR: orientación incorrecta del ranking GSE28511."
  )
}


cat("\n====================================================\n")
cat("ORIENTACIÓN\n")
cat("====================================================\n")

cat(
  "\nRanking positivo = mayor expresión en RMS.\n"
)

cat(
  "Ranking negativo = mayor expresión en músculo.\n"
)

cat(
  "\nVERIFICACIÓN: orientación consistente",
  "en ambos estudios.\n"
)


# ------------------------------------------------------------
# 12. Guardar rankings
# ------------------------------------------------------------

write.csv(
  ranking_rna,
  "ranking_GSE108022_RMS_vs_musculo.csv",
  row.names = FALSE
)

saveRDS(
  ranking_rna,
  "ranking_GSE108022_RMS_vs_musculo.rds"
)


write.csv(
  ranking_micro,
  "ranking_GSE28511_RMS_vs_musculo_sonda_unica.csv",
  row.names = FALSE
)

saveRDS(
  ranking_micro,
  "ranking_GSE28511_RMS_vs_musculo_sonda_unica.rds"
)


cat("\nArchivos generados:\n")

cat(
  "- ranking_GSE108022_RMS_vs_musculo.csv\n"
)

cat(
  "- ranking_GSE108022_RMS_vs_musculo.rds\n"
)

cat(
  "- ranking_GSE28511_RMS_vs_musculo_sonda_unica.csv\n"
)

cat(
  "- ranking_GSE28511_RMS_vs_musculo_sonda_unica.rds\n"
)


# ------------------------------------------------------------
# 13. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("PREPARACIÓN DE RANKINGS COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nTodavía NO se realizó ningún análisis",
  "de enriquecimiento funcional.\n"
)

cat(
  "\nLos rankings fueron construidos",
  "independientemente para cada estudio.\n"
)