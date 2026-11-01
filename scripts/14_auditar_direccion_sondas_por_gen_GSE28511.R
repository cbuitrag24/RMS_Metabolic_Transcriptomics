# ============================================================
# 14_auditar_direccion_sondas_por_gen_GSE28511.R
#
# OBJETIVO:
# Evaluar la concordancia direccional entre sondas
# asignadas al mismo Symbol en GSE28511.
#
# Contraste previamente ajustado:
# RMS - Skeletal_muscle
#
# IMPORTANTE:
# - No compara resultados con GSE108022.
# - No selecciona sondas por significación.
# - No define genes diferencialmente expresados.
# - No considera concordancia direccional como
#   evidencia suficiente de validación.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_micro <-
  "GSE28511_primario_limma_todas_las_sondas.rds"

archivo_universo <-
  "universo_comun_GSE108022_GSE28511.csv"

if (!file.exists(archivo_micro)) {
  stop("ERROR: falta la tabla de sondas GSE28511.")
}

if (!file.exists(archivo_universo)) {
  stop("ERROR: falta el universo común.")
}


# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

micro <- readRDS(archivo_micro)

universo <- read.csv(
  archivo_universo,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensiones microarray:\n")
print(dim(micro))

cat("\nDimensiones universo común:\n")
print(dim(universo))


# ------------------------------------------------------------
# 3. Verificar columnas e identidades
# ------------------------------------------------------------

necesarias_micro <- c(
  "Probe_ID",
  "Symbol",
  "logFC"
)

necesarias_universo <- c(
  "gene",
  "n_sondas_GSE28511"
)

if (!all(necesarias_micro %in% colnames(micro))) {
  stop("ERROR: faltan columnas en microarray.")
}

if (!all(necesarias_universo %in% colnames(universo))) {
  stop("ERROR: faltan columnas en universo común.")
}

if (anyDuplicated(micro$Probe_ID)) {
  stop("ERROR: hay sondas duplicadas.")
}

if (anyDuplicated(universo$gene)) {
  stop("ERROR: hay genes duplicados en universo.")
}

if (anyNA(micro$logFC) ||
    any(!is.finite(micro$logFC))) {
  stop("ERROR: logFC contiene valores inválidos.")
}

if (anyNA(universo$gene) ||
    any(trimws(universo$gene) == "")) {
  stop("ERROR: existen genes sin identidad.")
}


# ------------------------------------------------------------
# 4. Seleccionar universo común
# ------------------------------------------------------------

micro$Symbol <- trimws(
  as.character(micro$Symbol)
)

micro$Symbol[is.na(micro$Symbol)] <- ""

micro_comun <- micro[
  micro$Symbol %in% universo$gene,
  ,
  drop = FALSE
]

cat("\n====================================================\n")
cat("UNIVERSO COMÚN\n")
cat("====================================================\n")

cat("\nGenes comunes esperados:\n")
cat(nrow(universo), "\n")

cat("\nGenes comunes encontrados:\n")
cat(length(unique(micro_comun$Symbol)), "\n")

cat("\nSondas correspondientes:\n")
cat(nrow(micro_comun), "\n")

if (
  length(unique(micro_comun$Symbol)) !=
  nrow(universo)
) {
  stop("ERROR: faltan genes del universo común.")
}


# ------------------------------------------------------------
# 5. Agrupar sondas por Symbol
# ------------------------------------------------------------

indices <- split(
  seq_len(nrow(micro_comun)),
  micro_comun$Symbol
)

# Se conservan todas las sondas de cada gen.
# No se elige una sonda representativa.


# ------------------------------------------------------------
# 6. Resumen direccional por gen
# ------------------------------------------------------------

resumen <- vector(
  "list",
  length(indices)
)

names(resumen) <- names(indices)

for (s in names(indices)) {

  idx <- indices[[s]]

  valores <- micro_comun$logFC[idx]

  n_pos <- sum(valores > 0)
  n_neg <- sum(valores < 0)
  n_cero <- sum(valores == 0)

  n <- length(valores)

  # Clasificación puramente descriptiva.
  # No utiliza P-values ni FDR.

  if (n_cero == n) {

    clasificacion <- "todos_cero"

  } else if (n_pos == n) {

    clasificacion <- "todas_positivas"

  } else if (n_neg == n) {

    clasificacion <- "todas_negativas"

  } else if (n_cero > 0 &&
             n_pos > 0 &&
             n_neg == 0) {

    clasificacion <- "positivas_y_cero"

  } else if (n_cero > 0 &&
             n_neg > 0 &&
             n_pos == 0) {

    clasificacion <- "negativas_y_cero"

  } else {

    clasificacion <- "direcciones_mixtas"
  }

  resumen[[s]] <- data.frame(
    gene = s,
    n_sondas = n,
    n_positivas = n_pos,
    n_negativas = n_neg,
    n_cero = n_cero,
    fraccion_positivas = n_pos / n,
    fraccion_negativas = n_neg / n,
    logFC_min = min(valores),
    logFC_mediana = median(valores),
    logFC_max = max(valores),
    clasificacion = clasificacion,
    stringsAsFactors = FALSE
  )
}

tabla <- do.call(
  rbind,
  resumen
)

rownames(tabla) <- NULL


# ------------------------------------------------------------
# 7. Recuperar orden original del universo
# ------------------------------------------------------------

orden <- match(
  universo$gene,
  tabla$gene
)

if (anyNA(orden)) {
  stop("ERROR: no se pudo alinear el universo.")
}

tabla <- tabla[
  orden,
  ,
  drop = FALSE
]

rownames(tabla) <- NULL


# ------------------------------------------------------------
# 8. Verificar multiplicidad
# ------------------------------------------------------------

if (!identical(
  as.integer(tabla$n_sondas),
  as.integer(universo$n_sondas_GSE28511)
)) {
  stop(
    "ERROR: multiplicidad inconsistente con script 13."
  )
}

cat(
  "\nVERIFICACIÓN: multiplicidad de sondas",
  "coincide exactamente con el universo común.\n"
)


# ------------------------------------------------------------
# 9. Resultados descriptivos
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA DIRECCIONAL GLOBAL\n")
cat("====================================================\n")

cat("\nGenes evaluados:\n")
cat(nrow(tabla), "\n")

cat("\nClasificaciones:\n")
print(table(tabla$clasificacion))


# ------------------------------------------------------------
# 10. Genes con una sola sonda
# ------------------------------------------------------------

unica <- tabla[
  tabla$n_sondas == 1,
  ,
  drop = FALSE
]

cat("\n====================================================\n")
cat("GENES CON UNA SOLA SONDA\n")
cat("====================================================\n")

cat("\nNúmero de genes:\n")
cat(nrow(unica), "\n")

cat("\nClasificación:\n")
print(table(unica$clasificacion))


# ------------------------------------------------------------
# 11. Genes con múltiples sondas
# ------------------------------------------------------------

multiple <- tabla[
  tabla$n_sondas > 1,
  ,
  drop = FALSE
]

cat("\n====================================================\n")
cat("GENES CON MÚLTIPLES SONDAS\n")
cat("====================================================\n")

cat("\nNúmero de genes:\n")
cat(nrow(multiple), "\n")

cat("\nClasificación:\n")
print(table(multiple$clasificacion))

cat("\nSondas todas positivas:\n")
cat(sum(
  multiple$clasificacion == "todas_positivas"
), "\n")

cat("\nSondas todas negativas:\n")
cat(sum(
  multiple$clasificacion == "todas_negativas"
), "\n")

cat("\nSondas con direcciones mixtas:\n")
cat(sum(
  multiple$clasificacion == "direcciones_mixtas"
), "\n")


# ------------------------------------------------------------
# 12. Control de consistencia
# ------------------------------------------------------------

if (any(
  tabla$n_positivas +
  tabla$n_negativas +
  tabla$n_cero !=
  tabla$n_sondas
)) {
  stop(
    "ERROR: los recuentos direccionales no coinciden."
  )
}

if (anyDuplicated(tabla$gene)) {
  stop(
    "ERROR: genes duplicados en tabla final."
  )
}

if (nrow(tabla) != 14553) {
  stop(
    "ERROR: número inesperado de genes comunes."
  )
}

cat(
  "\nVERIFICACIÓN: recuentos, identidades",
  "y número de genes correctos.\n"
)


# ------------------------------------------------------------
# 13. Guardar auditoría
# ------------------------------------------------------------

write.csv(
  tabla,
  "GSE28511_auditoria_direccion_sondas_por_gen.csv",
  row.names = FALSE
)

saveRDS(
  tabla,
  "GSE28511_auditoria_direccion_sondas_por_gen.rds"
)

cat("\nArchivos generados:\n")

cat(
  "- GSE28511_auditoria_direccion_sondas_por_gen.csv\n"
)

cat(
  "- GSE28511_auditoria_direccion_sondas_por_gen.rds\n"
)


# ------------------------------------------------------------
# 14. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA DIRECCIONAL COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: la concordancia entre sondas",
  "es descriptiva y NO demuestra significación",
  "estadística ni validación independiente.\n"
)

cat(
  "\nTodavía NO se compararon los resultados",
  "con GSE108022.\n"
)