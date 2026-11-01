# ============================================================
# 15_concordancia_GSE108022_GSE28511_sonda_unica.R
#
# OBJETIVO:
# Evaluar concordancia de efectos RMS vs músculo entre:
#
# GSE108022:
# RNA-seq, análisis primario GPL16791
# 64 RMS vs 5 Normal
#
# GSE28511:
# microarray GPL6947
# 18 RMS vs 3 Skeletal_muscle
#
# ANÁLISIS PRIMARIO DE CONCORDANCIA:
# solamente genes del universo común representados por
# UNA ÚNICA SONDA en GSE28511.
#
# Se evalúa:
# - correlación Pearson de logFC
# - correlación Spearman de logFC
# - concordancia direccional
# - concordancia en genes con evidencia fuerte en GSE108022
#
# IMPORTANTE:
# - NO se seleccionan genes usando GSE28511.
# - NO se define todavía ningún gen como "validado".
# - NO se utilizan genes con múltiples sondas en este análisis.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_rnaseq <- file.path(
  "..",
  "GPL16791_DE_RMS_vs_Normal_todos_los_genes_CORREGIDO.csv"
)

archivo_micro <-
  "GSE28511_primario_limma_todas_las_sondas.rds"

archivo_universo <-
  "universo_comun_GSE108022_GSE28511.csv"


archivos <- c(
  archivo_rnaseq,
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

rnaseq <- read.csv(
  archivo_rnaseq,
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
print(dim(rnaseq))

cat("\nDimensiones GSE28511:\n")
print(dim(micro))

cat("\nDimensiones universo común:\n")
print(dim(universo))


# ------------------------------------------------------------
# 3. Verificar columnas necesarias
# ------------------------------------------------------------

necesarias_rnaseq <- c(
  "gene",
  "logFC",
  "PValue",
  "FDR"
)

necesarias_micro <- c(
  "Probe_ID",
  "Symbol",
  "logFC",
  "PValue",
  "FDR"
)

necesarias_universo <- c(
  "gene",
  "n_sondas_GSE28511"
)


if (!all(necesarias_rnaseq %in% colnames(rnaseq))) {
  stop(
    "ERROR: faltan columnas necesarias en GSE108022."
  )
}

if (!all(necesarias_micro %in% colnames(micro))) {
  stop(
    "ERROR: faltan columnas necesarias en GSE28511."
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

rnaseq$gene <- trimws(
  as.character(rnaseq$gene)
)

micro$Symbol <- trimws(
  as.character(micro$Symbol)
)

micro$Symbol[is.na(micro$Symbol)] <- ""

universo$gene <- trimws(
  as.character(universo$gene)
)


# ------------------------------------------------------------
# 5. Seleccionar genes con UNA sola sonda
#
# Esta selección depende exclusivamente de la arquitectura
# de la plataforma, NO de los resultados diferenciales.
# ------------------------------------------------------------

genes_unicos <- universo$gene[
  universo$n_sondas_GSE28511 == 1
]


cat("\n====================================================\n")
cat("UNIVERSO DE SONDA ÚNICA\n")
cat("====================================================\n")

cat("\nGenes comunes totales:\n")
cat(nrow(universo), "\n")

cat("\nGenes con una única sonda:\n")
cat(length(genes_unicos), "\n")


if (length(genes_unicos) != 9428) {
  stop(
    paste(
      "ERROR: se esperaban 9428 genes de sonda única y se encontraron",
      length(genes_unicos)
    )
  )
}


# ------------------------------------------------------------
# 6. Extraer RNA-seq
# ------------------------------------------------------------

idx_rna <- match(
  genes_unicos,
  rnaseq$gene
)

if (anyNA(idx_rna)) {
  stop(
    "ERROR: faltan genes de sonda única en GSE108022."
  )
}

rna_u <- rnaseq[
  idx_rna,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Extraer microarray
# ------------------------------------------------------------

idx_micro <- match(
  genes_unicos,
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


# ------------------------------------------------------------
# 8. Verificar identidad exacta
# ------------------------------------------------------------

if (!identical(
  rna_u$gene,
  genes_unicos
)) {
  stop(
    "ERROR: identidad incorrecta en RNA-seq."
  )
}

if (!identical(
  micro_u$Symbol,
  genes_unicos
)) {
  stop(
    "ERROR: identidad incorrecta en microarray."
  )
}


if (anyDuplicated(rna_u$gene)) {
  stop(
    "ERROR: genes duplicados en RNA-seq."
  )
}

if (anyDuplicated(micro_u$Symbol)) {
  stop(
    "ERROR: genes duplicados en microarray de sonda única."
  )
}


cat(
  "\nVERIFICACIÓN: los 9428 genes están",
  "perfectamente alineados entre ambos estudios.\n"
)


# ------------------------------------------------------------
# 9. Construir tabla integrada
# ------------------------------------------------------------

tabla <- data.frame(
  gene = genes_unicos,

  logFC_GSE108022 =
    as.numeric(rna_u$logFC),

  PValue_GSE108022 =
    as.numeric(rna_u$PValue),

  FDR_GSE108022 =
    as.numeric(rna_u$FDR),

  Probe_ID_GSE28511 =
    as.character(micro_u$Probe_ID),

  logFC_GSE28511 =
    as.numeric(micro_u$logFC),

  PValue_GSE28511 =
    as.numeric(micro_u$PValue),

  FDR_GSE28511 =
    as.numeric(micro_u$FDR),

  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 10. Control de valores
# ------------------------------------------------------------

columnas_numericas <- c(
  "logFC_GSE108022",
  "PValue_GSE108022",
  "FDR_GSE108022",
  "logFC_GSE28511",
  "PValue_GSE28511",
  "FDR_GSE28511"
)

for (col in columnas_numericas) {

  x <- tabla[[col]]

  if (anyNA(x)) {
    stop(
      paste(
        "ERROR: existen NA en",
        col
      )
    )
  }

  if (any(!is.finite(x))) {
    stop(
      paste(
        "ERROR: existen valores no finitos en",
        col
      )
    )
  }
}


# ------------------------------------------------------------
# 11. Correlación global de efectos
# ------------------------------------------------------------

pearson <- cor(
  tabla$logFC_GSE108022,
  tabla$logFC_GSE28511,
  method = "pearson"
)

spearman <- cor(
  tabla$logFC_GSE108022,
  tabla$logFC_GSE28511,
  method = "spearman"
)


cat("\n====================================================\n")
cat("CONCORDANCIA GLOBAL DE logFC\n")
cat("====================================================\n")

cat("\nPearson:\n")
cat(pearson, "\n")

cat("\nSpearman:\n")
cat(spearman, "\n")


# ------------------------------------------------------------
# 12. Concordancia direccional global
#
# Solo se evalúan genes cuyo logFC no sea exactamente cero
# en ninguno de los estudios.
# ------------------------------------------------------------

evaluables <- (
  tabla$logFC_GSE108022 != 0 &
  tabla$logFC_GSE28511 != 0
)

misma_direccion <- (
  sign(tabla$logFC_GSE108022[evaluables]) ==
  sign(tabla$logFC_GSE28511[evaluables])
)


cat("\n====================================================\n")
cat("CONCORDANCIA DIRECCIONAL GLOBAL\n")
cat("====================================================\n")

cat("\nGenes evaluables:\n")
cat(sum(evaluables), "\n")

cat("\nGenes con misma dirección:\n")
cat(sum(misma_direccion), "\n")

cat("\nProporción con misma dirección:\n")
cat(mean(misma_direccion), "\n")


# ------------------------------------------------------------
# 13. Conjunto FUERTE definido EXCLUSIVAMENTE en descubrimiento
#
# Criterio ya utilizado en GSE108022:
# FDR < 0.05 y |logFC| >= 1
#
# GSE28511 NO participa en la selección.
# Esto evita seleccionar los genes por el resultado
# de la cohorte de validación.
# ------------------------------------------------------------

fuerte_descubrimiento <- (
  tabla$FDR_GSE108022 < 0.05 &
  abs(tabla$logFC_GSE108022) >= 1
)


cat("\n====================================================\n")
cat("GENES FUERTES EN GSE108022\n")
cat("FDR < 0.05 Y |logFC| >= 1\n")
cat("====================================================\n")

cat("\nGenes de sonda única que cumplen el criterio:\n")
cat(sum(fuerte_descubrimiento), "\n")


# ------------------------------------------------------------
# 14. Dirección en validación de genes fuertes
#
# Aquí NO exigimos significación en GSE28511.
# Primero preguntamos únicamente si el efecto tiene
# la misma dirección en la cohorte independiente.
# ------------------------------------------------------------

fuertes <- tabla[
  fuerte_descubrimiento,
  ,
  drop = FALSE
]

evaluables_fuertes <- (
  fuertes$logFC_GSE108022 != 0 &
  fuertes$logFC_GSE28511 != 0
)

misma_direccion_fuertes <- (
  sign(
    fuertes$logFC_GSE108022[evaluables_fuertes]
  ) ==
  sign(
    fuertes$logFC_GSE28511[evaluables_fuertes]
  )
)


cat("\nGenes fuertes evaluables direccionalmente:\n")
cat(sum(evaluables_fuertes), "\n")

cat("\nGenes fuertes con misma dirección en GSE28511:\n")
cat(sum(misma_direccion_fuertes), "\n")

cat("\nProporción de genes fuertes con misma dirección:\n")
cat(mean(misma_direccion_fuertes), "\n")


# ------------------------------------------------------------
# 15. Correlación de efectos en genes fuertes
#
# Se informa descriptivamente.
# La selección se realizó únicamente en descubrimiento.
# ------------------------------------------------------------

if (nrow(fuertes) >= 3) {

  pearson_fuertes <- cor(
    fuertes$logFC_GSE108022,
    fuertes$logFC_GSE28511,
    method = "pearson"
  )

  spearman_fuertes <- cor(
    fuertes$logFC_GSE108022,
    fuertes$logFC_GSE28511,
    method = "spearman"
  )

  cat("\nPearson en genes fuertes:\n")
  cat(pearson_fuertes, "\n")

  cat("\nSpearman en genes fuertes:\n")
  cat(spearman_fuertes, "\n")
}


# ------------------------------------------------------------
# 16. Añadir variables descriptivas a tabla
# ------------------------------------------------------------

tabla$misma_direccion <- (
  sign(tabla$logFC_GSE108022) ==
  sign(tabla$logFC_GSE28511)
)

tabla$fuerte_GSE108022 <- fuerte_descubrimiento


# ------------------------------------------------------------
# 17. Controles finales
# ------------------------------------------------------------

if (nrow(tabla) != 9428) {
  stop(
    "ERROR: número inesperado de genes en tabla integrada."
  )
}

if (anyDuplicated(tabla$gene)) {
  stop(
    "ERROR: existen genes duplicados."
  )
}


cat("\n====================================================\n")
cat("VERIFICACIÓN FINAL\n")
cat("====================================================\n")

cat("\nGenes integrados:\n")
cat(nrow(tabla), "\n")

cat("\nGenes duplicados:\n")
cat(sum(duplicated(tabla$gene)), "\n")

cat(
  "\nVERIFICACIÓN: tabla integrada construida",
  "sin pérdida ni duplicación de genes.\n"
)


# ------------------------------------------------------------
# 18. Guardar resultados
# ------------------------------------------------------------

write.csv(
  tabla,
  "concordancia_GSE108022_GSE28511_sonda_unica.csv",
  row.names = FALSE
)

saveRDS(
  tabla,
  "concordancia_GSE108022_GSE28511_sonda_unica.rds"
)


cat("\nArchivos generados:\n")
cat(
  "- concordancia_GSE108022_GSE28511_sonda_unica.csv\n"
)
cat(
  "- concordancia_GSE108022_GSE28511_sonda_unica.rds\n"
)


# ------------------------------------------------------------
# 19. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("ANÁLISIS DE CONCORDANCIA COMPLETADO\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: estos resultados evalúan",
  "replicación de dirección y magnitud relativa",
  "entre estudios.\n"
)

cat(
  "\nTodavía NO se ha definido ningún gen",
  "como validado individualmente.\n"
)