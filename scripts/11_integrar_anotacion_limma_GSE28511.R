# ============================================================
# 11_integrar_anotacion_limma_GSE28511.R
#
# Objetivo:
# - recuperar la identidad exacta de las 48.802 sondas
# - integrar el ajuste limma con la anotación oficial GPL6947
# - verificar correspondencia sonda -> anotación
# - conservar todas las sondas individualmente
#
# IMPORTANTE:
# - NO colapsa sondas por Symbol
# - NO selecciona sondas significativas
# - NO aplica umbrales de FDR o logFC
# ============================================================


# ------------------------------------------------------------
# 1. Archivos necesarios
# ------------------------------------------------------------

archivo_fit <- "GSE28511_primario_limma_ebayes.rds"
archivo_expr <- "GSE28511_primario_expression_log2.rds"
archivo_anotacion <- "GPL6947_anotacion_GEO.csv"

archivos <- c(
  archivo_fit,
  archivo_expr,
  archivo_anotacion
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
# 2. Cargar objetos
# ------------------------------------------------------------

fit <- readRDS(
  archivo_fit
)

expr <- readRDS(
  archivo_expr
)

anotacion <- read.csv(
  archivo_anotacion,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("Dimensiones del ajuste limma:\n")
print(dim(fit$coefficients))

cat("\nDimensiones de expresión primaria:\n")
print(dim(expr))

cat("\nDimensiones de anotación GPL6947:\n")
print(dim(anotacion))


# ------------------------------------------------------------
# 3. Recuperar IDs de las sondas
# ------------------------------------------------------------

probe_expr <- rownames(expr)
probe_fit <- rownames(fit$coefficients)

if (is.null(probe_expr)) {
  stop(
    "ERROR: la matriz de expresión no tiene IDs de sonda como rownames."
  )
}

if (is.null(probe_fit)) {
  stop(
    "ERROR: el ajuste limma no tiene IDs de sonda como rownames."
  )
}

if (anyNA(probe_expr) || any(probe_expr == "")) {
  stop(
    "ERROR: existen IDs de sonda vacíos o NA en expresión."
  )
}

if (anyNA(probe_fit) || any(probe_fit == "")) {
  stop(
    "ERROR: existen IDs de sonda vacíos o NA en limma."
  )
}


cat("\nNúmero de IDs de sonda en expresión:\n")
cat(length(probe_expr), "\n")

cat("\nNúmero de IDs de sonda en limma:\n")
cat(length(probe_fit), "\n")


# ------------------------------------------------------------
# 4. Verificar identidad y orden expresión <-> limma
# ------------------------------------------------------------

if (!identical(
  probe_expr,
  probe_fit
)) {
  stop(
    "ERROR: los IDs o el orden de las sondas difieren entre expresión y limma."
  )
}

cat(
  "\nVERIFICACIÓN: identidad y orden de sondas",
  "preservados entre expresión y limma.\n"
)


# ------------------------------------------------------------
# 5. Verificar unicidad de sondas
# ------------------------------------------------------------

duplicadas_expr <- sum(
  duplicated(probe_expr)
)

cat("\nIDs de sonda duplicados en expresión:\n")
cat(duplicadas_expr, "\n")

if (duplicadas_expr != 0) {
  stop(
    "ERROR: existen IDs de sonda duplicados en expresión."
  )
}


# ------------------------------------------------------------
# 6. Identificar columna ID en GPL6947
# ------------------------------------------------------------

if (!("ID" %in% colnames(anotacion))) {
  stop(
    "ERROR: la anotación GPL6947 no contiene la columna ID."
  )
}

probe_annot <- as.character(
  anotacion$ID
)

if (anyNA(probe_annot) || any(probe_annot == "")) {
  stop(
    "ERROR: existen IDs vacíos o NA en la anotación GPL6947."
  )
}


duplicadas_annot <- sum(
  duplicated(probe_annot)
)

cat("\nIDs duplicados en anotación GPL6947:\n")
cat(duplicadas_annot, "\n")

if (duplicadas_annot != 0) {
  stop(
    "ERROR: existen IDs duplicados en la anotación GPL6947."
  )
}


# ------------------------------------------------------------
# 7. Emparejar las 48.802 sondas con GPL6947
# ------------------------------------------------------------

indice <- match(
  probe_expr,
  probe_annot
)

n_sin_anotacion <- sum(
  is.na(indice)
)

cat("\nSondas de expresión sin correspondencia en GPL6947:\n")
cat(n_sin_anotacion, "\n")

if (n_sin_anotacion != 0) {
  stop(
    "ERROR: existen sondas de expresión sin correspondencia en GPL6947."
  )
}

anotacion_expr <- anotacion[
  indice,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 8. Verificar correspondencia exacta
# ------------------------------------------------------------

if (!identical(
  as.character(anotacion_expr$ID),
  probe_expr
)) {
  stop(
    "ERROR: la anotación no quedó alineada exactamente con las sondas."
  )
}

cat(
  "\nVERIFICACIÓN: las 48.802 sondas tienen",
  "correspondencia exacta con GPL6947.\n"
)


# ------------------------------------------------------------
# 9. Verificar columna Symbol
# ------------------------------------------------------------

if (!("Symbol" %in% colnames(anotacion_expr))) {
  stop(
    "ERROR: la anotación no contiene la columna Symbol."
  )
}

symbol <- trimws(
  as.character(anotacion_expr$Symbol)
)

symbol[is.na(symbol)] <- ""


cat("\n====================================================\n")
cat("AUDITORÍA DE SYMBOL\n")
cat("====================================================\n")

cat("\nSondas totales:\n")
cat(length(symbol), "\n")

cat("\nSondas con Symbol no vacío:\n")
cat(sum(symbol != ""), "\n")

cat("\nSondas sin Symbol:\n")
cat(sum(symbol == ""), "\n")

cat("\nSymbols distintos no vacíos:\n")
cat(
  length(unique(symbol[symbol != ""])),
  "\n"
)


# ------------------------------------------------------------
# 10. Extraer estadísticos limma SIN ordenar
#
# Se mantiene exactamente el orden original de las sondas.
# ------------------------------------------------------------

coef_nombre <- "RMS_vs_Skeletal_muscle"

if (!(coef_nombre %in% colnames(fit$coefficients))) {
  stop(
    paste(
      "ERROR: no existe el coeficiente",
      coef_nombre
    )
  )
}


logFC <- fit$coefficients[
  ,
  coef_nombre
]

t_stat <- fit$t[
  ,
  coef_nombre
]

p_value <- fit$p.value[
  ,
  coef_nombre
]


if (
  length(logFC) != length(probe_expr) ||
  length(t_stat) != length(probe_expr) ||
  length(p_value) != length(probe_expr)
) {
  stop(
    "ERROR: longitud incompatible de los estadísticos limma."
  )
}


# ------------------------------------------------------------
# 11. Calcular FDR sobre las 48.802 sondas
# ------------------------------------------------------------

FDR <- p.adjust(
  p_value,
  method = "BH"
)


# ------------------------------------------------------------
# 12. Construir tabla integrada
# ------------------------------------------------------------

tabla <- data.frame(
  Probe_ID = probe_expr,
  Symbol = symbol,
  logFC = as.numeric(logFC),
  t = as.numeric(t_stat),
  PValue = as.numeric(p_value),
  FDR = as.numeric(FDR),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Añadir identificadores oficiales útiles si existen,
# sin modificar la identidad de las sondas.

columnas_extra <- c(
  "ILMN_Gene",
  "Entrez_Gene_ID",
  "RefSeq_ID"
)

for (col in columnas_extra) {

  if (col %in% colnames(anotacion_expr)) {

    tabla[[col]] <- as.character(
      anotacion_expr[[col]]
    )
  }
}


# ------------------------------------------------------------
# 13. Verificaciones finales de identidad
# ------------------------------------------------------------

if (!identical(
  tabla$Probe_ID,
  probe_expr
)) {
  stop(
    "ERROR: se perdió la identidad u orden de las sondas."
  )
}

if (anyNA(tabla$logFC)) {
  stop(
    "ERROR: existen NA en logFC."
  )
}

if (anyNA(tabla$t)) {
  stop(
    "ERROR: existen NA en t."
  )
}

if (anyNA(tabla$PValue)) {
  stop(
    "ERROR: existen NA en PValue."
  )
}

if (anyNA(tabla$FDR)) {
  stop(
    "ERROR: existen NA en FDR."
  )
}


cat("\n====================================================\n")
cat("VERIFICACIÓN DE ESTADÍSTICOS\n")
cat("====================================================\n")

cat("\nFilas de la tabla integrada:\n")
cat(nrow(tabla), "\n")

cat("\nResumen de logFC:\n")
print(
  summary(tabla$logFC)
)

cat("\nResumen de PValue:\n")
print(
  summary(tabla$PValue)
)

cat("\nResumen de FDR:\n")
print(
  summary(tabla$FDR)
)


# ------------------------------------------------------------
# 14. Verificación independiente con topTable
#
# topTable(sort.by = "none") debe devolver los mismos
# estadísticos en el mismo orden.
# ------------------------------------------------------------

suppressPackageStartupMessages(
  library(limma)
)

control <- topTable(
  fit,
  coef = coef_nombre,
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)


if (nrow(control) != nrow(tabla)) {
  stop(
    "ERROR: topTable tiene un número inesperado de filas."
  )
}

if (!identical(
  rownames(control),
  tabla$Probe_ID
)) {
  stop(
    "ERROR: topTable no conserva la misma identidad/orden de sondas."
  )
}


error_logFC <- max(
  abs(
    control$logFC -
      tabla$logFC
  )
)

error_PValue <- max(
  abs(
    control$P.Value -
      tabla$PValue
  )
)

error_FDR <- max(
  abs(
    control$adj.P.Val -
      tabla$FDR
  )
)


cat("\nMáximo error logFC vs topTable:\n")
cat(
  format(error_logFC, scientific = TRUE),
  "\n"
)

cat("\nMáximo error PValue vs topTable:\n")
cat(
  format(error_PValue, scientific = TRUE),
  "\n"
)

cat("\nMáximo error FDR vs topTable:\n")
cat(
  format(error_FDR, scientific = TRUE),
  "\n"
)


tolerancia <- 1e-12

if (
  error_logFC > tolerancia ||
  error_PValue > tolerancia ||
  error_FDR > tolerancia
) {
  stop(
    "ERROR: los estadísticos no coinciden con topTable."
  )
}

cat(
  "\nVERIFICACIÓN: identidad de sonda y estadísticos",
  "coinciden exactamente con topTable.\n"
)


# ------------------------------------------------------------
# 15. Guardar tabla completa
#
# TODAS las sondas.
# Sin filtro de significación.
# ------------------------------------------------------------

write.csv(
  tabla,
  "GSE28511_primario_limma_todas_las_sondas.csv",
  row.names = FALSE
)

saveRDS(
  tabla,
  "GSE28511_primario_limma_todas_las_sondas.rds"
)


cat("\nArchivos generados:\n")
cat(
  "- GSE28511_primario_limma_todas_las_sondas.csv\n"
)
cat(
  "- GSE28511_primario_limma_todas_las_sondas.rds\n"
)


# ------------------------------------------------------------
# 16. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("INTEGRACIÓN LIMMA + GPL6947 COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: todavía NO se colapsaron sondas",
  "por gen y NO se seleccionaron genes diferenciales.\n"
)