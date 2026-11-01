# ============================================================
# 06_auditar_sondas_repetidas_GSE28511.R
#
# Objetivo:
# - evaluar la concordancia entre sondas asignadas al mismo Symbol
# - trabajar sobre la matriz log2 ya preparada
# - usar correlación de Pearson entre perfiles de las 24 muestras
#
# NO colapsa sondas.
# NO realiza expresión diferencial.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos necesarios
# ------------------------------------------------------------

archivo_expr <- "GSE28511_expression_quantileNormalized_log2.rds"
archivo_anotacion <- "GPL6947_anotacion_GEO.csv"

if (!file.exists(archivo_expr)) {
  stop(paste("ERROR: no se encontró", archivo_expr))
}

if (!file.exists(archivo_anotacion)) {
  stop(paste("ERROR: no se encontró", archivo_anotacion))
}


# ------------------------------------------------------------
# 2. Leer datos
# ------------------------------------------------------------

expr <- readRDS(archivo_expr)

anotacion <- read.csv(
  archivo_anotacion,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensiones de expresión:\n")
print(dim(expr))

cat("\nDimensiones de anotación:\n")
print(dim(anotacion))


# ------------------------------------------------------------
# 3. Alinear anotación con matriz
# ------------------------------------------------------------

if (!all(c("ID", "Symbol") %in% colnames(anotacion))) {
  stop("ERROR: faltan las columnas ID o Symbol.")
}

ids_expr <- rownames(expr)

indice <- match(
  ids_expr,
  anotacion$ID
)

if (anyNA(indice)) {
  stop("ERROR: existen sondas sin correspondencia en GPL6947.")
}

ann <- anotacion[indice, , drop = FALSE]

if (!identical(
  as.character(ann$ID),
  as.character(ids_expr)
)) {
  stop("ERROR: expresión y anotación no quedaron alineadas.")
}

symbol <- trimws(as.character(ann$Symbol))
symbol[is.na(symbol)] <- ""

cat(
  "\nVERIFICACIÓN: expresión y anotación correctamente alineadas.\n"
)


# ------------------------------------------------------------
# 4. Seleccionar símbolos representados por >1 sonda
# ------------------------------------------------------------

indices_con_symbol <- which(symbol != "")

grupos <- split(
  indices_con_symbol,
  symbol[indices_con_symbol]
)

grupos_repetidos <- grupos[
  lengths(grupos) > 1
]

cat("\n====================================================\n")
cat("SÍMBOLOS CON MÚLTIPLES SONDAS\n")
cat("====================================================\n")

cat("\nNúmero de símbolos con >1 sonda:\n")
cat(length(grupos_repetidos), "\n")

cat("\nNúmero total de sondas involucradas:\n")
cat(sum(lengths(grupos_repetidos)), "\n")


# ------------------------------------------------------------
# 5. Calcular correlaciones entre sondas del mismo símbolo
#
# Para cada símbolo:
# - se calculan todas las correlaciones por pares
# - se resume mediana, mínimo y máximo
# ------------------------------------------------------------

resultados <- vector(
  "list",
  length(grupos_repetidos)
)

nombres_genes <- names(grupos_repetidos)


for (i in seq_along(grupos_repetidos)) {

  idx <- grupos_repetidos[[i]]

  x <- expr[
    idx,
    ,
    drop = FALSE
  ]

  # Correlación entre FILAS (sondas) usando las 24 muestras
  correlaciones <- cor(
    t(x),
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  valores <- correlaciones[
    upper.tri(correlaciones)
  ]

  valores_finitos <- valores[
    is.finite(valores)
  ]

  n_pares_total <- length(valores)

  n_pares_validos <- length(valores_finitos)

  if (n_pares_validos > 0) {

    mediana_cor <- median(valores_finitos)

    media_cor <- mean(valores_finitos)

    minimo_cor <- min(valores_finitos)

    maximo_cor <- max(valores_finitos)

  } else {

    mediana_cor <- NA_real_
    media_cor <- NA_real_
    minimo_cor <- NA_real_
    maximo_cor <- NA_real_
  }

  resultados[[i]] <- data.frame(
    Symbol = nombres_genes[i],
    n_sondas = length(idx),
    n_pares_total = n_pares_total,
    n_pares_validos = n_pares_validos,
    correlacion_mediana = mediana_cor,
    correlacion_media = media_cor,
    correlacion_minima = minimo_cor,
    correlacion_maxima = maximo_cor,
    stringsAsFactors = FALSE
  )
}


resultado <- do.call(
  rbind,
  resultados
)


# ------------------------------------------------------------
# 6. Verificar correlaciones no calculables
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("CONTROL DE CORRELACIONES\n")
cat("====================================================\n")

cat("\nSímbolos sin ninguna correlación válida:\n")
cat(
  sum(resultado$n_pares_validos == 0),
  "\n"
)

cat("\nSímbolos con al menos una correlación válida:\n")
cat(
  sum(resultado$n_pares_validos > 0),
  "\n"
)


# ------------------------------------------------------------
# 7. Distribución de la correlación mediana por símbolo
# ------------------------------------------------------------

validos <- resultado$correlacion_mediana[
  is.finite(resultado$correlacion_mediana)
]

cat("\n====================================================\n")
cat("DISTRIBUCIÓN DE CORRELACIÓN MEDIANA\n")
cat("====================================================\n\n")

print(
  summary(validos)
)


# ------------------------------------------------------------
# 8. Umbrales descriptivos
#
# Son únicamente descriptivos.
# NO se usan todavía para excluir genes.
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RESUMEN DESCRIPTIVO DE CONCORDANCIA\n")
cat("====================================================\n")

umbrales <- c(
  0,
  0.25,
  0.50,
  0.75,
  0.90
)

for (u in umbrales) {

  cat(
    "\nCorrelación mediana >=",
    u,
    ":",
    sum(validos >= u),
    "de",
    length(validos),
    "(",
    round(
      100 * mean(validos >= u),
      2
    ),
    "%)\n"
  )
}

cat(
  "\nCorrelación mediana < 0:",
  sum(validos < 0),
  "de",
  length(validos),
  "(",
  round(
    100 * mean(validos < 0),
    2
  ),
  "%)\n"
)


# ------------------------------------------------------------
# 9. Genes con menor concordancia
# ------------------------------------------------------------

resultado_ordenado <- resultado[
  order(
    resultado$correlacion_mediana,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

cat("\n====================================================\n")
cat("30 SÍMBOLOS CON MENOR CORRELACIÓN MEDIANA\n")
cat("====================================================\n\n")

print(
  head(
    resultado_ordenado,
    30
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Genes con mayor concordancia
# ------------------------------------------------------------

resultado_mayor <- resultado[
  order(
    -resultado$correlacion_mediana,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

cat("\n====================================================\n")
cat("30 SÍMBOLOS CON MAYOR CORRELACIÓN MEDIANA\n")
cat("====================================================\n\n")

print(
  head(
    resultado_mayor,
    30
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Guardar resultados
# ------------------------------------------------------------

write.csv(
  resultado_ordenado,
  "GSE28511_concordancia_sondas_mismo_symbol.csv",
  row.names = FALSE
)

cat("\nArchivo generado:\n")
cat(
  "- GSE28511_concordancia_sondas_mismo_symbol.csv\n"
)


# ------------------------------------------------------------
# 12. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA DE SONDAS REPETIDAS COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: este análisis es diagnóstico;",
  "todavía NO se colapsaron sondas a genes.\n"
)