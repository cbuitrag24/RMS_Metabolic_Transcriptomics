# ============================================================
# 33_concordancia_FN_FOXO1FP_vs_Normal.R
# Concordancia entre:
#   FN vs Normal
#   FOXO1-FP vs Normal
# GPL16791 - GSE108022
# ============================================================

# ------------------------------------------------------------
# 1. Cargar tablas completas
# ------------------------------------------------------------

fn <- read.csv(
  "GPL16791_DE_FN_vs_Normal_todos.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fp <- read.csv(
  "GPL16791_DE_FOXO1FP_vs_Normal_todos.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Tablas cargadas.\n")

cat("\nGenes FN vs Normal:\n")
print(nrow(fn))

cat("\nGenes FOXO1-FP vs Normal:\n")
print(nrow(fp))


# ------------------------------------------------------------
# 2. Verificaciones
# ------------------------------------------------------------

columnas_necesarias <- c(
  "gene",
  "logFC",
  "PValue",
  "FDR"
)

if (!all(columnas_necesarias %in% colnames(fn))) {
  stop("ERROR: faltan columnas necesarias en FN vs Normal.")
}

if (!all(columnas_necesarias %in% colnames(fp))) {
  stop("ERROR: faltan columnas necesarias en FOXO1-FP vs Normal.")
}

if (anyNA(fn$gene) || any(fn$gene == "")) {
  stop("ERROR: genes inválidos en FN vs Normal.")
}

if (anyNA(fp$gene) || any(fp$gene == "")) {
  stop("ERROR: genes inválidos en FOXO1-FP vs Normal.")
}

if (anyDuplicated(fn$gene) > 0) {
  stop("ERROR: genes duplicados en FN vs Normal.")
}

if (anyDuplicated(fp$gene) > 0) {
  stop("ERROR: genes duplicados en FOXO1-FP vs Normal.")
}

cat(
  "\nVERIFICACIÓN: identificadores génicos válidos y no duplicados.\n"
)


# ------------------------------------------------------------
# 3. Emparejar por símbolo génico
# ------------------------------------------------------------

comparacion <- merge(
  fn[, c("gene", "logFC", "PValue", "FDR")],
  fp[, c("gene", "logFC", "PValue", "FDR")],
  by = "gene",
  suffixes = c("_FN", "_FOXO1FP"),
  all = FALSE,
  sort = FALSE
)

cat("\nGenes comunes:\n")
print(nrow(comparacion))

cat("\nGenes exclusivos de FN vs Normal:\n")
print(sum(!fn$gene %in% comparacion$gene))

cat("\nGenes exclusivos de FOXO1-FP vs Normal:\n")
print(sum(!fp$gene %in% comparacion$gene))

if (nrow(comparacion) == 0) {
  stop("ERROR: no existen genes comunes.")
}

if (anyNA(comparacion$logFC_FN) ||
    anyNA(comparacion$logFC_FOXO1FP)) {
  stop("ERROR: existen NA en los logFC.")
}


# ------------------------------------------------------------
# 4. Correlación de tamaños de efecto
# ------------------------------------------------------------

pearson <- cor(
  comparacion$logFC_FN,
  comparacion$logFC_FOXO1FP,
  method = "pearson"
)

spearman <- cor(
  comparacion$logFC_FN,
  comparacion$logFC_FOXO1FP,
  method = "spearman"
)

cat("\nCorrelación Pearson de logFC:\n")
print(pearson)

cat("\nCorrelación Spearman de logFC:\n")
print(spearman)


# ------------------------------------------------------------
# 5. Diferencia absoluta entre logFC
# ------------------------------------------------------------

comparacion$diferencia_abs_logFC <- abs(
  comparacion$logFC_FN -
  comparacion$logFC_FOXO1FP
)

cat("\nResumen de |logFC_FN - logFC_FOXO1FP|:\n")
print(summary(comparacion$diferencia_abs_logFC))


# ------------------------------------------------------------
# 6. Concordancia de dirección global
# ------------------------------------------------------------

evaluables <- comparacion[
  comparacion$logFC_FN != 0 &
  comparacion$logFC_FOXO1FP != 0,
  ,
  drop = FALSE
]

misma_direccion <- (
  sign(evaluables$logFC_FN) ==
  sign(evaluables$logFC_FOXO1FP)
)

cat("\nGenes evaluables para dirección:\n")
print(nrow(evaluables))

cat("\nGenes con la misma dirección:\n")
print(sum(misma_direccion))

cat("\nProporción con la misma dirección:\n")
print(mean(misma_direccion))


# ------------------------------------------------------------
# 7. Aplicar criterio previamente utilizado:
# FDR < 0.05 y |logFC| >= 1
# ------------------------------------------------------------

comparacion$criterio_FN <- (
  comparacion$FDR_FN < 0.05 &
  abs(comparacion$logFC_FN) >= 1
)

comparacion$criterio_FOXO1FP <- (
  comparacion$FDR_FOXO1FP < 0.05 &
  abs(comparacion$logFC_FOXO1FP) >= 1
)

ambos <- (
  comparacion$criterio_FN &
  comparacion$criterio_FOXO1FP
)

solo_FN <- (
  comparacion$criterio_FN &
  !comparacion$criterio_FOXO1FP
)

solo_FOXO1FP <- (
  !comparacion$criterio_FN &
  comparacion$criterio_FOXO1FP
)

ninguno <- (
  !comparacion$criterio_FN &
  !comparacion$criterio_FOXO1FP
)

cat("\nCumplen criterio en FN vs Normal:\n")
print(sum(comparacion$criterio_FN))

cat("\nCumplen criterio en FOXO1-FP vs Normal:\n")
print(sum(comparacion$criterio_FOXO1FP))

cat("\nCumplen criterio en ambos:\n")
print(sum(ambos))

cat("\nCumplen criterio sólo en FN:\n")
print(sum(solo_FN))

cat("\nCumplen criterio sólo en FOXO1-FP:\n")
print(sum(solo_FOXO1FP))

cat("\nNo cumplen criterio en ninguno:\n")
print(sum(ninguno))


# ------------------------------------------------------------
# 8. Dirección entre genes que cumplen en ambos
# ------------------------------------------------------------

comp_ambos <- comparacion[
  ambos,
  ,
  drop = FALSE
]

if (nrow(comp_ambos) > 0) {

  misma_dir_ambos <- (
    sign(comp_ambos$logFC_FN) ==
    sign(comp_ambos$logFC_FOXO1FP)
  )

  cat(
    "\nGenes que cumplen criterio en ambos y conservan dirección:\n"
  )
  print(sum(misma_dir_ambos))

  cat(
    "\nGenes que cumplen criterio en ambos pero tienen dirección opuesta:\n"
  )
  print(sum(!misma_dir_ambos))

  cat(
    "\nProporción direccional entre los que cumplen en ambos:\n"
  )
  print(mean(misma_dir_ambos))
}


# ------------------------------------------------------------
# 9. Separar señal compartida positiva y negativa
#
# Positiva:
# mayor expresión que Normal en ambos grupos RMS
#
# Negativa:
# menor expresión que Normal en ambos grupos RMS
# ------------------------------------------------------------

compartidos_arriba <- (
  ambos &
  comparacion$logFC_FN >= 1 &
  comparacion$logFC_FOXO1FP >= 1
)

compartidos_abajo <- (
  ambos &
  comparacion$logFC_FN <= -1 &
  comparacion$logFC_FOXO1FP <= -1
)

cat(
  "\nCompartidos con logFC >= 1 en FN y FOXO1-FP:\n"
)
print(sum(compartidos_arriba))

cat(
  "\nCompartidos con logFC <= -1 en FN y FOXO1-FP:\n"
)
print(sum(compartidos_abajo))


# ------------------------------------------------------------
# 10. Índice de Jaccard
# ------------------------------------------------------------

union_criterio <- (
  comparacion$criterio_FN |
  comparacion$criterio_FOXO1FP
)

jaccard <- sum(ambos) / sum(union_criterio)

cat("\nÍndice de Jaccard:\n")
print(jaccard)


# ------------------------------------------------------------
# 11. Guardar clasificación
# ------------------------------------------------------------

comparacion$clasificacion <- "ninguno"

comparacion$clasificacion[solo_FN] <- "solo_FN"

comparacion$clasificacion[solo_FOXO1FP] <- "solo_FOXO1FP"

comparacion$clasificacion[
  compartidos_arriba
] <- "compartido_mayor_en_RMS"

comparacion$clasificacion[
  compartidos_abajo
] <- "compartido_mayor_en_Normal"

# Control: ningún gen que cumpla en ambos debería quedar
# sin clasificar si conserva dirección.
sin_clasificar_ambos <- (
  ambos &
  comparacion$clasificacion == "ninguno"
)

cat(
  "\nGenes que cumplen criterio en ambos pero no entran en una categoría compartida:\n"
)
print(sum(sin_clasificar_ambos))


# ------------------------------------------------------------
# 12. Guardar tabla completa
# ------------------------------------------------------------

comparacion <- comparacion[
  order(
    comparacion$diferencia_abs_logFC,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

rownames(comparacion) <- NULL

write.csv(
  comparacion,
  "concordancia_FN_vs_FOXO1FP_respecto_Normal.csv",
  row.names = FALSE
)

cat("\nArchivo generado:\n")
cat("- concordancia_FN_vs_FOXO1FP_respecto_Normal.csv\n")

cat(
  "\nRESULTADO: concordancia FN y FOXO1-FP respecto de Normal completada.\n"
)