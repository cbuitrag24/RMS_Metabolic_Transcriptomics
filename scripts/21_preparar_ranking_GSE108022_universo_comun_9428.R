# ============================================================
# 21_preparar_ranking_GSE108022_universo_comun_9428.R
#
# OBJETIVO:
# Restringir el ranking de GSE108022 al mismo universo génico
# utilizado en la validación primaria de GSE28511:
# genes comunes representados por exactamente una sonda.
#
# IMPORTANTE:
# - No se recalcula ninguna estadística.
# - No se seleccionan genes según resultados de GSE28511.
# - No se realiza todavía enriquecimiento funcional.
# - Se conserva el estadístico prerankeado original de GSE108022.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_ranking <-
  "ranking_GSE108022_RMS_vs_musculo.rds"

archivo_validacion <-
  "concordancia_GSE108022_GSE28511_sonda_unica.rds"


if (!file.exists(archivo_ranking)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_ranking
    )
  )
}

if (!file.exists(archivo_validacion)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_validacion
    )
  )
}


# ------------------------------------------------------------
# 2. Cargar objetos
# ------------------------------------------------------------

ranking <- readRDS(
  archivo_ranking
)

validacion <- readRDS(
  archivo_validacion
)


cat("\nDimensiones ranking completo GSE108022:\n")
print(
  dim(ranking)
)

cat("\nDimensiones tabla de universo sonda única:\n")
print(
  dim(validacion)
)

cat("\nColumnas ranking:\n")
print(
  colnames(ranking)
)

cat("\nColumnas tabla sonda única:\n")
print(
  colnames(validacion)
)


# ------------------------------------------------------------
# 3. Verificar estructura del ranking
# ------------------------------------------------------------

if (!all(
  c("gene", "ranking_stat") %in% colnames(ranking)
)) {
  stop(
    "ERROR: el ranking no contiene gene y ranking_stat."
  )
}


ranking$gene <- trimws(
  as.character(ranking$gene)
)

ranking$ranking_stat <- as.numeric(
  ranking$ranking_stat
)


if (anyDuplicated(ranking$gene)) {
  stop(
    "ERROR: existen genes duplicados en el ranking."
  )
}

if (
  anyNA(ranking$gene) ||
  any(ranking$gene == "")
) {
  stop(
    "ERROR: existen genes vacíos en el ranking."
  )
}

if (
  anyNA(ranking$ranking_stat) ||
  any(!is.finite(ranking$ranking_stat))
) {
  stop(
    "ERROR: existen estadísticas inválidas en el ranking."
  )
}


# ------------------------------------------------------------
# 4. Identificar de forma inequívoca la columna del gen
#    de GSE108022 en la tabla de concordancia
#
# Primero se buscan nombres explícitos esperables.
# No se adivina silenciosamente una columna.
# ------------------------------------------------------------

candidatas <- c(
  "gene",
  "gene_GSE108022",
  "Gene_GSE108022",
  "Symbol",
  "symbol"
)

col_gen <- candidatas[
  candidatas %in% colnames(validacion)
]


if (length(col_gen) != 1) {

  cat(
    "\nNo pudo identificarse de forma inequívoca",
    "la columna del gen de GSE108022.\n"
  )

  cat(
    "Columnas disponibles:\n"
  )

  print(
    colnames(validacion)
  )

  stop(
    "ERROR: revisar la columna identificadora del gen."
  )
}


col_gen <- col_gen[1]


cat("\nColumna utilizada como identificador génico:\n")
cat(
  col_gen,
  "\n"
)


genes_universo <- trimws(
  as.character(
    validacion[[col_gen]]
  )
)


# ------------------------------------------------------------
# 5. Auditar universo de validación
# ------------------------------------------------------------

if (
  anyNA(genes_universo) ||
  any(genes_universo == "")
) {
  stop(
    "ERROR: existen identificadores génicos vacíos en el universo."
  )
}


if (anyDuplicated(genes_universo)) {
  stop(
    "ERROR: existen genes duplicados en el universo sonda única."
  )
}


cat("\nNúmero de genes del universo sonda única:\n")
cat(
  length(genes_universo),
  "\n"
)


if (length(genes_universo) != 9428) {
  stop(
    paste(
      "ERROR: se esperaban 9428 genes y se encontraron",
      length(genes_universo)
    )
  )
}


# ------------------------------------------------------------
# 6. Verificar presencia exacta en GSE108022
# ------------------------------------------------------------

faltantes <- setdiff(
  genes_universo,
  ranking$gene
)


cat("\nGenes del universo ausentes del ranking GSE108022:\n")
cat(
  length(faltantes),
  "\n"
)


if (length(faltantes) > 0) {

  cat(
    "\nPrimeros genes faltantes:\n"
  )

  print(
    head(
      faltantes,
      20
    )
  )

  stop(
    "ERROR: el universo sonda única no está completamente representado."
  )
}


# ------------------------------------------------------------
# 7. Restringir el ranking
#
# match() conserva exactamente el conjunto definido por
# genes_universo. Después se vuelve a ordenar por ranking_stat.
# ------------------------------------------------------------

idx <- match(
  genes_universo,
  ranking$gene
)


if (anyNA(idx)) {
  stop(
    "ERROR: match produjo posiciones NA."
  )
}


ranking_comun <- ranking[
  idx,
  ,
  drop = FALSE
]


# Verificación exacta antes de ordenar

if (!identical(
  ranking_comun$gene,
  genes_universo
)) {
  stop(
    "ERROR: no se preservó correctamente el universo génico."
  )
}


# ------------------------------------------------------------
# 8. Ordenar por estadístico decreciente
# ------------------------------------------------------------

ranking_comun <- ranking_comun[
  order(
    ranking_comun$ranking_stat,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

rownames(ranking_comun) <- NULL


# ------------------------------------------------------------
# 9. Auditorías finales
# ------------------------------------------------------------

if (nrow(ranking_comun) != 9428) {
  stop(
    "ERROR: el ranking restringido no contiene 9428 genes."
  )
}


if (anyDuplicated(ranking_comun$gene)) {
  stop(
    "ERROR: hay genes duplicados en el ranking restringido."
  )
}


if (!setequal(
  ranking_comun$gene,
  genes_universo
)) {
  stop(
    "ERROR: el conjunto génico final no coincide con el universo esperado."
  )
}


if (is.unsorted(
  -ranking_comun$ranking_stat
)) {
  stop(
    "ERROR: el ranking restringido no quedó ordenado correctamente."
  )
}


# Comprobar que no se modificó ninguna estadística

idx_original <- match(
  ranking_comun$gene,
  ranking$gene
)

error_stat <- max(
  abs(
    ranking_comun$ranking_stat -
    ranking$ranking_stat[idx_original]
  )
)


cat("\n====================================================\n")
cat("AUDITORÍA RANKING COMÚN\n")
cat("====================================================\n")


cat("\nGenes ranking completo GSE108022:\n")
cat(
  nrow(ranking),
  "\n"
)

cat("\nGenes universo común sonda única:\n")
cat(
  nrow(ranking_comun),
  "\n"
)

cat("\nGenes eliminados únicamente por restricción de universo:\n")
cat(
  nrow(ranking) - nrow(ranking_comun),
  "\n"
)

cat("\nError máximo al comparar ranking_stat con el original:\n")
cat(
  format(
    error_stat,
    scientific = TRUE
  ),
  "\n"
)


if (error_stat != 0) {
  stop(
    "ERROR: alguna estadística fue modificada."
  )
}


cat("\nResumen del estadístico restringido:\n")
print(
  summary(
    ranking_comun$ranking_stat
  )
)


cat("\nEstadísticos positivos:\n")
cat(
  sum(ranking_comun$ranking_stat > 0),
  "\n"
)

cat("\nEstadísticos negativos:\n")
cat(
  sum(ranking_comun$ranking_stat < 0),
  "\n"
)

cat("\nEstadísticos iguales a cero:\n")
cat(
  sum(ranking_comun$ranking_stat == 0),
  "\n"
)


# ------------------------------------------------------------
# 10. Guardar
# ------------------------------------------------------------

write.csv(
  ranking_comun,
  "ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.csv",
  row.names = FALSE
)

saveRDS(
  ranking_comun,
  "ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.rds"
)


cat("\nArchivos generados:\n")

cat(
  "- ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.csv\n"
)

cat(
  "- ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.rds\n"
)


cat("\n====================================================\n")
cat("PREPARACIÓN DEL RANKING COMÚN COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nTodavía NO se realizó un nuevo análisis de enriquecimiento.\n"
)