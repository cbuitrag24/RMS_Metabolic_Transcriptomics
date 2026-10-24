# ============================================================
# 03_preparar_log2_GSE28511.R
# Preparación de la matriz procesada de GSE28511
#
# Objetivo:
# - leer Series Matrix
# - comprobar integridad
# - aplicar log2 a los valores ya normalizados por cuantiles
# - NO volver a normalizar
# - NO hacer expresión diferencial
# ============================================================


# ------------------------------------------------------------
# 1. Localizar Series Matrix
# ------------------------------------------------------------

archivos <- list.files(
  pattern = "^GSE28511_series_matrix.*\\.txt\\.gz$",
  full.names = TRUE
)

if (length(archivos) != 1) {
  stop(
    paste(
      "ERROR: se esperaba exactamente una Series Matrix;",
      "se encontraron",
      length(archivos)
    )
  )
}

archivo <- archivos[1]

cat("Archivo encontrado:\n")
cat(archivo, "\n")


# ------------------------------------------------------------
# 2. Leer únicamente la tabla de expresión
# ------------------------------------------------------------

lineas <- readLines(
  gzfile(archivo),
  warn = FALSE
)

inicio <- grep(
  "^!series_matrix_table_begin",
  lineas
)

fin <- grep(
  "^!series_matrix_table_end",
  lineas
)

if (length(inicio) != 1 || length(fin) != 1) {
  stop("ERROR: no se localizaron correctamente los límites de la tabla.")
}

if (fin <= inicio + 1) {
  stop("ERROR: la tabla de expresión está vacía.")
}

tabla_lineas <- lineas[
  (inicio + 1):(fin - 1)
]

conexion <- textConnection(tabla_lineas)

tabla <- read.delim(
  conexion,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)

close(conexion)


# ------------------------------------------------------------
# 3. Verificar dimensiones
# ------------------------------------------------------------

cat("\nDimensiones de la tabla:\n")
print(dim(tabla))

if (nrow(tabla) != 48802) {
  stop(
    paste(
      "ERROR: se esperaban 48802 filas y se encontraron",
      nrow(tabla)
    )
  )
}

if (ncol(tabla) != 25) {
  stop(
    paste(
      "ERROR: se esperaban 25 columnas y se encontraron",
      ncol(tabla)
    )
  )
}

if (colnames(tabla)[1] != "ID_REF") {
  stop("ERROR: la primera columna no es ID_REF.")
}

if (anyDuplicated(tabla$ID_REF) > 0) {
  stop("ERROR: existen ID_REF duplicados.")
}

cat(
  "\nVERIFICACIÓN: dimensiones e identificadores correctos.\n"
)


# ------------------------------------------------------------
# 4. Construir matriz numérica
# ------------------------------------------------------------

id_ref <- tabla$ID_REF

expr <- as.matrix(
  tabla[, -1, drop = FALSE]
)

storage.mode(expr) <- "numeric"

rownames(expr) <- id_ref

if (anyNA(expr)) {
  stop("ERROR: existen NA en la matriz de expresión.")
}

if (any(!is.finite(expr))) {
  stop("ERROR: existen valores no finitos.")
}

if (any(expr <= 0)) {
  stop(
    "ERROR: existen valores <= 0; no se aplicará log2 automáticamente."
  )
}

cat("\nDimensiones de la matriz numérica:\n")
print(dim(expr))

cat("\nRango antes de log2:\n")
print(range(expr))


# ------------------------------------------------------------
# 5. Transformación log2
#
# Todos los valores son estrictamente positivos.
# Por eso usamos log2(x), sin sumar una constante.
# ------------------------------------------------------------

expr_log2 <- log2(expr)

if (anyNA(expr_log2) ||
    any(!is.finite(expr_log2))) {
  stop("ERROR: la transformación log2 produjo valores inválidos.")
}

cat("\nRango después de log2:\n")
print(range(expr_log2))

cat("\nResumen global después de log2:\n")
print(summary(as.vector(expr_log2)))


# ------------------------------------------------------------
# 6. Resumen por muestra después de log2
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RESUMEN POR MUESTRA DESPUÉS DE LOG2\n")
cat("====================================================\n")

resumen <- data.frame(
  muestra = colnames(expr_log2),
  minimo = apply(expr_log2, 2, min),
  Q1 = apply(expr_log2, 2, quantile, probs = 0.25),
  mediana = apply(expr_log2, 2, median),
  media = colMeans(expr_log2),
  Q3 = apply(expr_log2, 2, quantile, probs = 0.75),
  maximo = apply(expr_log2, 2, max),
  stringsAsFactors = FALSE
)

print(
  resumen,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Guardar matriz transformada
# ------------------------------------------------------------

saveRDS(
  expr_log2,
  "GSE28511_expression_quantileNormalized_log2.rds"
)

write.csv(
  resumen,
  "GSE28511_resumen_muestras_log2.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Verificación final
# ------------------------------------------------------------

cat("\nArchivos generados:\n")
cat(
  "- GSE28511_expression_quantileNormalized_log2.rds\n"
)
cat(
  "- GSE28511_resumen_muestras_log2.csv\n"
)

cat(
  "\nIMPORTANTE: no se realizó una nueva normalización.\n"
)

cat(
  "RESULTADO FINAL: matriz procesada transformada a log2 y verificada.\n"
)