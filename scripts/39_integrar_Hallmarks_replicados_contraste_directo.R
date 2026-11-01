
# ============================================================
# SCRIPT 39
# Integracion de 13 Hallmarks replicados
# con el contraste directo FOXO1-FP vs FN
#
# No recalcula fgsea ni modifica resultados originales.
# ============================================================

archivo_replicados <-
  "Hallmark_13_replicados_consistencia_FN_FOXO1FP.csv"

archivo_directo <-
  "Hallmark_fgsea_GSE108022_FOXO1FP_vs_FN_directo_9427.csv"

if (!file.exists(archivo_replicados) ||
    !file.exists(archivo_directo)) {
  stop("ERROR: falta alguno de los archivos de entrada.")
}

# 1. Leer resultados originales

replicados <- read.csv(
  archivo_replicados,
  stringsAsFactors = FALSE
)

directo <- read.csv(
  archivo_directo,
  stringsAsFactors = FALSE
)

# 2. Verificar columnas

columnas_replicados <- c(
  "pathway",
  "NES_replicacion_1",
  "NES_replicacion_2",
  "NES_FN",
  "FDR_FN",
  "NES_FOXO1FP",
  "FDR_FOXO1FP",
  "consistente_ambos",
  "eje_original"
)

columnas_directo <- c(
  "pathway",
  "NES",
  "pval",
  "padj",
  "size"
)

if (!all(columnas_replicados %in% names(replicados))) {
  stop("ERROR: columnas faltantes en replicados.")
}

if (!all(columnas_directo %in% names(directo))) {
  stop("ERROR: columnas faltantes en contraste directo.")
}

# 3. Auditar identidad de los programas

if (nrow(replicados) != 13) {
  stop("ERROR: se esperaban 13 Hallmarks replicados.")
}

if (nrow(directo) != 50) {
  stop("ERROR: se esperaban 50 Hallmarks directos.")
}

if (anyNA(replicados$pathway) ||
    anyDuplicated(replicados$pathway) ||
    anyNA(directo$pathway) ||
    anyDuplicated(directo$pathway)) {
  stop("ERROR: identificadores de Hallmark invalidos.")
}

if (!all(replicados$consistente_ambos %in% TRUE)) {
  stop("ERROR: algun Hallmark no es consistente en ambos.")
}

indices <- match(replicados$pathway, directo$pathway)

if (anyNA(indices)) {
  stop("ERROR: faltan Hallmarks en el contraste directo.")
}

# 4. Integrar conservando el orden original

integrado <- replicados

integrado$NES_directo <- directo$NES[indices]
integrado$PValue_directo <- directo$pval[indices]
integrado$FDR_directo <- directo$padj[indices]
integrado$size_directo <- directo$size[indices]

if (any(!is.finite(integrado$NES_directo)) ||
    anyNA(integrado$FDR_directo)) {
  stop("ERROR: estadisticos directos invalidos.")
}

# 5. Clasificar resultados directos

integrado$significativo_directo <-
  integrado$FDR_directo < 0.05

integrado$direccion_directa <- ifelse(
  !integrado$significativo_directo,
  "No significativo",
  ifelse(
    integrado$NES_directo > 0,
    "FOXO1-FP",
    "FN"
  )
)

# 6. Resumen

cat("\nHallmarks replicados:", nrow(integrado), "\n")

cat(
  "Con diferencia directa significativa:",
  sum(integrado$significativo_directo),
  "\n"
)

cat(
  "Enriquecidos hacia FOXO1-FP:",
  sum(integrado$direccion_directa == "FOXO1-FP"),
  "\n"
)

cat(
  "Enriquecidos hacia FN:",
  sum(integrado$direccion_directa == "FN"),
  "\n"
)

cat("\nTabla integrada:\n")

print(
  integrado[
    ,
    c(
      "pathway",
      "eje_original",
      "NES_FN",
      "NES_FOXO1FP",
      "NES_directo",
      "FDR_directo",
      "direccion_directa"
    )
  ],
  row.names = FALSE
)

# 7. Guardar

write.csv(
  integrado,
  "Hallmark_13_replicados_integracion_directa.csv",
  row.names = FALSE
)

saveRDS(
  integrado,
  "Hallmark_13_replicados_integracion_directa.rds"
)

cat("\nIntegracion completada y guardada.\n")