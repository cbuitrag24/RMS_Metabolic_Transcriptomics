# ============================================================
# 12_auditar_seleccion_sonda_representativa_GSE28511.R
#
# Objetivo:
# - auditar una regla para obtener una sonda representativa
#   por Symbol
#
# Regla:
# - entre sondas que comparten Symbol, elegir la sonda con
#   mayor expresión media log2 en las 21 muestras primarias
#
# IMPORTANTE:
# - la selección NO usa logFC
# - NO usa t
# - NO usa PValue
# - NO usa FDR
# - NO se realiza todavía validación contra GSE108022
# ============================================================


# ------------------------------------------------------------
# 1. Archivos necesarios
# ------------------------------------------------------------

archivo_expr <-
  "GSE28511_primario_expression_log2.rds"

archivo_resultados <-
  "GSE28511_primario_limma_todas_las_sondas.rds"


archivos <- c(
  archivo_expr,
  archivo_resultados
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

expr <- readRDS(
  archivo_expr
)

resultados <- readRDS(
  archivo_resultados
)


cat("Dimensiones de expresión:\n")
print(dim(expr))

cat("\nDimensiones de resultados limma:\n")
print(dim(resultados))


# ------------------------------------------------------------
# 3. Verificar correspondencia exacta
# ------------------------------------------------------------

if (!("Probe_ID" %in% colnames(resultados))) {
  stop(
    "ERROR: resultados no contiene Probe_ID."
  )
}

if (!("Symbol" %in% colnames(resultados))) {
  stop(
    "ERROR: resultados no contiene Symbol."
  )
}


if (!identical(
  rownames(expr),
  as.character(resultados$Probe_ID)
)) {
  stop(
    paste(
      "ERROR: expresión y resultados limma",
      "no tienen la misma identidad/orden de sondas."
    )
  )
}

cat(
  "\nVERIFICACIÓN: expresión y resultados limma",
  "perfectamente alineados.\n"
)


# ------------------------------------------------------------
# 4. Limpiar únicamente espacios externos de Symbol
# ------------------------------------------------------------

symbol <- trimws(
  as.character(resultados$Symbol)
)

symbol[is.na(symbol)] <- ""


# ------------------------------------------------------------
# 5. Calcular expresión media de cada sonda
#
# IMPORTANTE:
# esta media utiliza las 21 muestras conjuntamente,
# sin utilizar las etiquetas de condición.
# ------------------------------------------------------------

media_expresion <- rowMeans(expr)


if (anyNA(media_expresion)) {
  stop(
    "ERROR: existen NA en la expresión media."
  )
}

if (any(!is.finite(media_expresion))) {
  stop(
    "ERROR: existen valores no finitos en la expresión media."
  )
}


cat("\n====================================================\n")
cat("EXPRESIÓN MEDIA POR SONDA\n")
cat("====================================================\n")

print(
  summary(media_expresion)
)


# ------------------------------------------------------------
# 6. Considerar solamente sondas con Symbol explícito
# ------------------------------------------------------------

usar <- symbol != ""

probe_symbol <- as.character(
  resultados$Probe_ID[usar]
)

symbol_validos <- symbol[usar]

media_validos <- media_expresion[usar]


cat("\n====================================================\n")
cat("SONDAS CON SYMBOL\n")
cat("====================================================\n")

cat("\nSondas totales:\n")
cat(nrow(expr), "\n")

cat("\nSondas con Symbol:\n")
cat(length(symbol_validos), "\n")

cat("\nSondas sin Symbol excluidas de esta auditoría:\n")
cat(sum(!usar), "\n")

cat("\nSymbols distintos:\n")
cat(length(unique(symbol_validos)), "\n")


# ------------------------------------------------------------
# 7. Multiplicidad de sondas por Symbol
# ------------------------------------------------------------

multiplicidad <- table(symbol_validos)

cat("\n====================================================\n")
cat("MULTIPLICIDAD SONDA -> SYMBOL\n")
cat("====================================================\n")

cat("\nSymbols con exactamente 1 sonda:\n")
cat(sum(multiplicidad == 1), "\n")

cat("\nSymbols con más de 1 sonda:\n")
cat(sum(multiplicidad > 1), "\n")

cat("\nMáximo número de sondas para un Symbol:\n")
cat(max(multiplicidad), "\n")


# ------------------------------------------------------------
# 8. Auditar empates exactos en máxima expresión media
#
# Para cada Symbol se identifica el máximo.
# Si más de una sonda tiene exactamente ese mismo valor,
# se registra como empate.
# ------------------------------------------------------------

indices_por_symbol <- split(
  seq_along(symbol_validos),
  symbol_validos
)

n_maximos <- integer(
  length(indices_por_symbol)
)

names(n_maximos) <- names(indices_por_symbol)


for (s in names(indices_por_symbol)) {

  idx <- indices_por_symbol[[s]]

  valores <- media_validos[idx]

  maximo <- max(valores)

  n_maximos[s] <- sum(
    valores == maximo
  )
}


symbols_con_empate <- names(
  n_maximos[n_maximos > 1]
)


cat("\n====================================================\n")
cat("AUDITORÍA DE EMPATES EXACTOS\n")
cat("====================================================\n")

cat("\nSymbols evaluados:\n")
cat(length(n_maximos), "\n")

cat("\nSymbols con máximo único:\n")
cat(sum(n_maximos == 1), "\n")

cat("\nSymbols con empate exacto en el máximo:\n")
cat(length(symbols_con_empate), "\n")


if (length(symbols_con_empate) > 0) {

  cat("\nPrimeros Symbols con empate exacto:\n")

  print(
    head(
      symbols_con_empate,
      20
    )
  )
}


# ------------------------------------------------------------
# 9. Seleccionar provisionalmente la sonda de mayor media
#
# Si hubiera empate exacto, se conserva provisionalmente
# la primera en el orden original SOLO para poder auditarla.
#
# La presencia de empates se informará explícitamente.
# ------------------------------------------------------------

seleccion <- integer(
  length(indices_por_symbol)
)

names(seleccion) <- names(indices_por_symbol)


for (s in names(indices_por_symbol)) {

  idx <- indices_por_symbol[[s]]

  valores <- media_validos[idx]

  seleccion[s] <- idx[
    which.max(valores)
  ]
}


probe_seleccionada <- probe_symbol[
  seleccion
]

media_seleccionada <- media_validos[
  seleccion
]


# ------------------------------------------------------------
# 10. Construir tabla de auditoría por Symbol
# ------------------------------------------------------------

tabla_auditoria <- data.frame(
  Symbol = names(seleccion),
  n_sondas = as.integer(
    multiplicidad[names(seleccion)]
  ),
  Probe_ID_seleccionada = probe_seleccionada,
  media_log2_seleccionada =
    as.numeric(media_seleccionada),
  empate_exacto_maximo =
    n_maximos[names(seleccion)] > 1,
  n_sondas_empatadas_maximo =
    as.integer(
      n_maximos[names(seleccion)]
    ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 11. Verificaciones
# ------------------------------------------------------------

if (
  nrow(tabla_auditoria) !=
    length(unique(symbol_validos))
) {
  stop(
    "ERROR: número inesperado de Symbols en la auditoría."
  )
}


if (anyDuplicated(tabla_auditoria$Symbol)) {
  stop(
    "ERROR: hay Symbols duplicados en la tabla de auditoría."
  )
}


if (anyDuplicated(
  tabla_auditoria$Probe_ID_seleccionada
)) {
  stop(
    paste(
      "ERROR: una misma sonda fue seleccionada",
      "para más de un Symbol."
    )
  )
}


# Verificación explícita:
# para cada Symbol, la sonda elegida debe tener
# exactamente la máxima media observada.

errores_seleccion <- 0L

for (i in seq_len(nrow(tabla_auditoria))) {

  s <- tabla_auditoria$Symbol[i]

  idx <- indices_por_symbol[[s]]

  maximo_real <- max(
    media_validos[idx]
  )

  media_elegida <-
    tabla_auditoria$media_log2_seleccionada[i]

  if (!isTRUE(
    all.equal(
      media_elegida,
      maximo_real,
      tolerance = 0
    )
  )) {
    errores_seleccion <-
      errores_seleccion + 1L
  }
}


cat("\n====================================================\n")
cat("VERIFICACIÓN DE LA REGLA\n")
cat("====================================================\n")

cat("\nNúmero de selecciones incorrectas:\n")
cat(errores_seleccion, "\n")

if (errores_seleccion != 0) {
  stop(
    "ERROR: alguna sonda seleccionada no tiene la máxima media."
  )
}

cat(
  "\nVERIFICACIÓN: para todos los Symbols,",
  "la sonda provisional seleccionada posee",
  "la máxima expresión media.\n"
)


# ------------------------------------------------------------
# 12. Resumen de las sondas seleccionadas
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RESUMEN DE SELECCIÓN PROVISIONAL\n")
cat("====================================================\n")

cat("\nSymbols representados:\n")
cat(nrow(tabla_auditoria), "\n")

cat("\nSondas seleccionadas:\n")
cat(length(probe_seleccionada), "\n")

cat("\nSondas descartadas por redundancia de Symbol:\n")
cat(
  length(symbol_validos) -
    length(probe_seleccionada),
  "\n"
)

cat("\nResumen de expresión media de sondas seleccionadas:\n")

print(
  summary(
    tabla_auditoria$media_log2_seleccionada
  )
)


# ------------------------------------------------------------
# 13. Guardar auditoría
#
# IMPORTANTE:
# esto NO es todavía la tabla génica definitiva.
# ------------------------------------------------------------

write.csv(
  tabla_auditoria,
  "GSE28511_auditoria_sonda_representativa.csv",
  row.names = FALSE
)

saveRDS(
  tabla_auditoria,
  "GSE28511_auditoria_sonda_representativa.rds"
)


cat("\nArchivos generados:\n")
cat(
  "- GSE28511_auditoria_sonda_representativa.csv\n"
)
cat(
  "- GSE28511_auditoria_sonda_representativa.rds\n"
)


# ------------------------------------------------------------
# 14. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: todavía NO se generó la tabla",
  "génica definitiva y NO se realizó la comparación",
  "con GSE108022.\n"
)