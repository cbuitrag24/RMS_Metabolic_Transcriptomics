# ============================================================
# 13_auditar_solapamiento_GSE108022_GSE28511.R
#
# Objetivo:
# - definir el universo génico común entre:
#
#   GSE108022:
#   análisis primario RNA-seq GPL16791
#   64 RMS vs 5 Normal
#
#   GSE28511:
#   microarray GPL6947
#   18 RMS vs 3 Skeletal_muscle
#
# - usar exclusivamente Symbol explícito para el emparejamiento
# - cuantificar cuántos genes de GSE108022 están representados
#   por 1 o múltiples sondas en GSE28511
#
# IMPORTANTE:
# - NO compara logFC
# - NO compara P-values ni FDR
# - NO define genes "validados"
# - NO selecciona una sonda representativa
# ============================================================


# ------------------------------------------------------------
# 1. Archivos necesarios
# ------------------------------------------------------------

archivo_rnaseq <- file.path(
  "..",
  "GPL16791_DE_RMS_vs_Normal_todos_los_genes_CORREGIDO.csv"
)

archivo_microarray <-
  "GSE28511_primario_limma_todas_las_sondas.rds"


archivos <- c(
  archivo_rnaseq,
  archivo_microarray
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
  archivo_microarray
)


cat("Dimensiones GSE108022:\n")
print(dim(rnaseq))

cat("\nDimensiones GSE28511 a nivel de sonda:\n")
print(dim(micro))


# ------------------------------------------------------------
# 3. Inspeccionar nombres de columnas
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("COLUMNAS DISPONIBLES\n")
cat("====================================================\n")

cat("\nGSE108022:\n")
print(colnames(rnaseq))

cat("\nGSE28511:\n")
print(colnames(micro))


# ------------------------------------------------------------
# 4. Identificar columna génica de GSE108022
#
# En la tabla corregida esperamos una columna "gene".
# No intentaremos inferir identidades desde otras columnas.
# ------------------------------------------------------------

if (!("gene" %in% colnames(rnaseq))) {

  stop(
    paste(
      "ERROR: la tabla corregida de GSE108022",
      "no contiene la columna 'gene'."
    )
  )
}


genes_rnaseq <- trimws(
  as.character(rnaseq$gene)
)

genes_rnaseq[is.na(genes_rnaseq)] <- ""


# ------------------------------------------------------------
# 5. Verificar identidad génica de GSE108022
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA GSE108022\n")
cat("====================================================\n")

cat("\nFilas totales:\n")
cat(length(genes_rnaseq), "\n")

cat("\nGenes con nombre no vacío:\n")
cat(sum(genes_rnaseq != ""), "\n")

cat("\nGenes con nombre vacío:\n")
cat(sum(genes_rnaseq == ""), "\n")

cat("\nGenes duplicados no vacíos:\n")

duplicados_rnaseq <- sum(
  duplicated(
    genes_rnaseq[genes_rnaseq != ""]
  )
)

cat(duplicados_rnaseq, "\n")


if (any(genes_rnaseq == "")) {

  stop(
    "ERROR: existen genes sin nombre en GSE108022."
  )
}


if (duplicados_rnaseq != 0) {

  stop(
    "ERROR: existen genes duplicados en GSE108022."
  )
}


# ------------------------------------------------------------
# 6. Verificar Symbol de GSE28511
# ------------------------------------------------------------

if (!("Probe_ID" %in% colnames(micro))) {

  stop(
    "ERROR: GSE28511 no contiene Probe_ID."
  )
}

if (!("Symbol" %in% colnames(micro))) {

  stop(
    "ERROR: GSE28511 no contiene Symbol."
  )
}


symbol_micro <- trimws(
  as.character(micro$Symbol)
)

symbol_micro[is.na(symbol_micro)] <- ""


cat("\n====================================================\n")
cat("AUDITORÍA GSE28511\n")
cat("====================================================\n")

cat("\nSondas totales:\n")
cat(nrow(micro), "\n")

cat("\nSondas con Symbol:\n")
cat(sum(symbol_micro != ""), "\n")

cat("\nSondas sin Symbol:\n")
cat(sum(symbol_micro == ""), "\n")

genes_micro <- unique(
  symbol_micro[symbol_micro != ""]
)

cat("\nSymbols distintos representados:\n")
cat(length(genes_micro), "\n")


# ------------------------------------------------------------
# 7. Construir universo común por coincidencia EXACTA de Symbol
#
# No se actualizan aliases.
# No se convierten nombres históricos.
# No se utilizan Entrez para rescatar coincidencias.
# ------------------------------------------------------------

genes_comunes <- intersect(
  genes_rnaseq,
  genes_micro
)

genes_solo_rnaseq <- setdiff(
  genes_rnaseq,
  genes_micro
)

genes_solo_micro <- setdiff(
  genes_micro,
  genes_rnaseq
)


cat("\n====================================================\n")
cat("UNIVERSO COMÚN POR SYMBOL EXACTO\n")
cat("====================================================\n")

cat("\nGenes GSE108022:\n")
cat(length(genes_rnaseq), "\n")

cat("\nSymbols GSE28511:\n")
cat(length(genes_micro), "\n")

cat("\nGenes comunes:\n")
cat(length(genes_comunes), "\n")

cat("\nGenes de GSE108022 no representados por Symbol exacto:\n")
cat(length(genes_solo_rnaseq), "\n")

cat("\nSymbols de GSE28511 no presentes en GSE108022:\n")
cat(length(genes_solo_micro), "\n")


prop_rnaseq_representada <-
  length(genes_comunes) /
  length(genes_rnaseq)


cat("\nProporción de genes GSE108022 representados en GSE28511:\n")
cat(
  round(
    prop_rnaseq_representada,
    6
  ),
  "\n"
)


# ------------------------------------------------------------
# 8. Multiplicidad de sondas SOLO dentro del universo común
# ------------------------------------------------------------

micro_con_symbol <- micro[
  symbol_micro != "",
  ,
  drop = FALSE
]

micro_con_symbol$Symbol <- symbol_micro[
  symbol_micro != ""
]


micro_comun <- micro_con_symbol[
  micro_con_symbol$Symbol %in% genes_comunes,
  ,
  drop = FALSE
]


multiplicidad_comun <- table(
  micro_comun$Symbol
)


cat("\n====================================================\n")
cat("MULTIPLICIDAD DE SONDAS EN EL UNIVERSO COMÚN\n")
cat("====================================================\n")

cat("\nGenes comunes evaluados:\n")
cat(length(multiplicidad_comun), "\n")

cat("\nGenes comunes con exactamente 1 sonda:\n")
cat(
  sum(multiplicidad_comun == 1),
  "\n"
)

cat("\nGenes comunes con más de 1 sonda:\n")
cat(
  sum(multiplicidad_comun > 1),
  "\n"
)

cat("\nMáximo número de sondas por gen común:\n")
cat(
  max(multiplicidad_comun),
  "\n"
)


# ------------------------------------------------------------
# 9. Distribución completa de multiplicidad
# ------------------------------------------------------------

cat("\nDistribución del número de sondas por gen común:\n")

print(
  table(
    as.integer(multiplicidad_comun)
  )
)


# ------------------------------------------------------------
# 10. Construir tabla del universo común
#
# Una fila por gen.
# No incorpora todavía resultados diferenciales del microarray.
# ------------------------------------------------------------

tabla_universo <- data.frame(
  gene = genes_rnaseq[
    genes_rnaseq %in% genes_comunes
  ],
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Mantener el orden de GSE108022.

tabla_universo$n_sondas_GSE28511 <-
  as.integer(
    multiplicidad_comun[
      tabla_universo$gene
    ]
  )


if (anyNA(tabla_universo$n_sondas_GSE28511)) {

  stop(
    paste(
      "ERROR: algún gen común quedó",
      "sin multiplicidad de sondas."
    )
  )
}


# ------------------------------------------------------------
# 11. Verificaciones
# ------------------------------------------------------------

if (
  nrow(tabla_universo) !=
    length(genes_comunes)
) {

  stop(
    "ERROR: número inesperado de genes en el universo común."
  )
}


if (anyDuplicated(tabla_universo$gene)) {

  stop(
    "ERROR: existen genes duplicados en el universo común."
  )
}


if (
  any(
    tabla_universo$n_sondas_GSE28511 < 1
  )
) {

  stop(
    "ERROR: existe un gen común sin sondas."
  )
}


cat("\n====================================================\n")
cat("VERIFICACIÓN FINAL\n")
cat("====================================================\n")

cat("\nFilas de la tabla del universo común:\n")
cat(nrow(tabla_universo), "\n")

cat("\nGenes duplicados:\n")
cat(
  sum(duplicated(tabla_universo$gene)),
  "\n"
)

cat(
  "\nVERIFICACIÓN: universo común construido",
  "sin duplicados y manteniendo el orden de GSE108022.\n"
)


# ------------------------------------------------------------
# 12. Guardar resultados de auditoría
# ------------------------------------------------------------

write.csv(
  tabla_universo,
  "universo_comun_GSE108022_GSE28511.csv",
  row.names = FALSE
)


write.csv(
  data.frame(
    gene = genes_solo_rnaseq,
    stringsAsFactors = FALSE
  ),
  "genes_GSE108022_sin_Symbol_exacto_GSE28511.csv",
  row.names = FALSE
)


cat("\nArchivos generados:\n")
cat(
  "- universo_comun_GSE108022_GSE28511.csv\n"
)
cat(
  "- genes_GSE108022_sin_Symbol_exacto_GSE28511.csv\n"
)


# ------------------------------------------------------------
# 13. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA DE SOLAPAMIENTO COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: todavía NO se compararon",
  "logFC, P-values ni FDR entre los estudios.\n"
)

cat(
  "Todavía NO se definió ningún gen como validado.\n"
)