# ============================================================
# 05_auditar_simbolos_GPL6947.R
#
# Objetivo:
# - auditar la correspondencia ID_REF -> Symbol de GPL6947
# - cuantificar símbolos ausentes
# - cuantificar símbolos repetidos en múltiples sondas
# - detectar anotaciones potencialmente ambiguas
#
# NO colapsa sondas a genes.
# NO modifica la matriz de expresión.
# NO realiza expresión diferencial.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos necesarios
# ------------------------------------------------------------

archivo_anotacion <- "GPL6947_anotacion_GEO.csv"
archivo_expr <- "GSE28511_expression_quantileNormalized_log2.rds"

if (!file.exists(archivo_anotacion)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_anotacion
    )
  )
}

if (!file.exists(archivo_expr)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_expr
    )
  )
}


# ------------------------------------------------------------
# 2. Leer anotación y matriz
# ------------------------------------------------------------

anotacion <- read.csv(
  archivo_anotacion,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expr <- readRDS(archivo_expr)

cat("Dimensiones de la anotación:\n")
print(dim(anotacion))

cat("\nDimensiones de la matriz de expresión:\n")
print(dim(expr))


# ------------------------------------------------------------
# 3. Verificaciones básicas
# ------------------------------------------------------------

columnas_necesarias <- c(
  "ID",
  "Symbol",
  "ILMN_Gene",
  "Entrez_Gene_ID",
  "RefSeq_ID"
)

faltantes <- setdiff(
  columnas_necesarias,
  colnames(anotacion)
)

if (length(faltantes) > 0) {
  stop(
    paste(
      "ERROR: faltan columnas necesarias:",
      paste(faltantes, collapse = ", ")
    )
  )
}

ids_expr <- rownames(expr)

if (is.null(ids_expr)) {
  stop("ERROR: la matriz no tiene ID_REF como rownames.")
}

if (anyDuplicated(ids_expr) > 0) {
  stop("ERROR: existen ID_REF duplicados en la matriz.")
}

if (anyDuplicated(anotacion$ID) > 0) {
  stop("ERROR: existen ID duplicados en la anotación.")
}

indice <- match(
  ids_expr,
  anotacion$ID
)

if (anyNA(indice)) {
  stop(
    "ERROR: existen ID_REF de expresión sin correspondencia en GPL6947."
  )
}

# Mantener exactamente el orden de la matriz
ann <- anotacion[indice, , drop = FALSE]

if (!identical(
  as.character(ann$ID),
  as.character(ids_expr)
)) {
  stop("ERROR: falló la alineación entre expresión y anotación.")
}

cat(
  "\nVERIFICACIÓN: las 48802 sondas están alineadas",
  "con la anotación oficial.\n"
)


# ------------------------------------------------------------
# 4. Limpiar campos SOLO para auditoría
#
# No se reemplazan ni actualizan símbolos.
# Únicamente se eliminan espacios externos.
# ------------------------------------------------------------

symbol <- trimws(as.character(ann$Symbol))
ilmn_gene <- trimws(as.character(ann$ILMN_Gene))
entrez <- trimws(as.character(ann$Entrez_Gene_ID))
refseq <- trimws(as.character(ann$RefSeq_ID))

# read.csv puede representar campos vacíos como NA
symbol[is.na(symbol)] <- ""
ilmn_gene[is.na(ilmn_gene)] <- ""
entrez[is.na(entrez)] <- ""
refseq[is.na(refseq)] <- ""


# ------------------------------------------------------------
# 5. Auditoría del campo Symbol
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA DEL CAMPO Symbol\n")
cat("====================================================\n")

n_total <- length(symbol)
n_con_symbol <- sum(symbol != "")
n_sin_symbol <- sum(symbol == "")

cat("\nSondas totales:\n")
cat(n_total, "\n")

cat("\nSondas con Symbol no vacío:\n")
cat(n_con_symbol, "\n")

cat("\nSondas sin Symbol:\n")
cat(n_sin_symbol, "\n")

cat("\nProporción con Symbol:\n")
cat(n_con_symbol / n_total, "\n")


# ------------------------------------------------------------
# 6. Símbolos únicos y repetidos
# ------------------------------------------------------------

symbol_no_vacio <- symbol[symbol != ""]

tabla_symbol <- sort(
  table(symbol_no_vacio),
  decreasing = TRUE
)

n_symbols_unicos <- length(tabla_symbol)

n_symbols_una_sonda <- sum(tabla_symbol == 1)

n_symbols_multiples_sondas <- sum(tabla_symbol > 1)

n_sondas_en_symbols_repetidos <- sum(
  tabla_symbol[tabla_symbol > 1]
)

cat("\nNúmero de símbolos distintos no vacíos:\n")
cat(n_symbols_unicos, "\n")

cat("\nSímbolos representados por exactamente 1 sonda:\n")
cat(n_symbols_una_sonda, "\n")

cat("\nSímbolos representados por >1 sonda:\n")
cat(n_symbols_multiples_sondas, "\n")

cat("\nSondas pertenecientes a símbolos repetidos:\n")
cat(n_sondas_en_symbols_repetidos, "\n")

cat("\nMáximo número de sondas para un mismo Symbol:\n")
cat(max(tabla_symbol), "\n")


# ------------------------------------------------------------
# 7. Mostrar los símbolos con más sondas
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("SÍMBOLOS CON MAYOR NÚMERO DE SONDAS\n")
cat("====================================================\n\n")

print(
  head(
    tabla_symbol,
    30
  )
)


# ------------------------------------------------------------
# 8. Auditar posibles separadores/ambigüedad textual
#
# Esto NO supone que sean realmente genes múltiples.
# Sólo detecta patrones de texto que debemos inspeccionar.
# ------------------------------------------------------------

patrones <- c(
  "///",
  "//",
  ";",
  ",",
  "|"
)

cat("\n====================================================\n")
cat("PATRONES POTENCIALMENTE AMBIGUOS EN Symbol\n")
cat("====================================================\n")

for (p in patrones) {

  n <- sum(
    grepl(
      p,
      symbol,
      fixed = TRUE
    )
  )

  cat(
    "\nPatrón",
    shQuote(p),
    ":",
    n,
    "\n"
  )
}


# ------------------------------------------------------------
# 9. Comparar disponibilidad de otros identificadores
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("DISPONIBILIDAD DE OTROS IDENTIFICADORES\n")
cat("====================================================\n")

cat("\nILMN_Gene no vacío:\n")
cat(sum(ilmn_gene != ""), "\n")

cat("\nEntrez_Gene_ID no vacío:\n")
cat(sum(entrez != ""), "\n")

cat("\nRefSeq_ID no vacío:\n")
cat(sum(refseq != ""), "\n")


# ------------------------------------------------------------
# 10. Sondas sin Symbol
# ------------------------------------------------------------

sin_symbol <- data.frame(
  ID = ann$ID[symbol == ""],
  ILMN_Gene = ilmn_gene[symbol == ""],
  Entrez_Gene_ID = entrez[symbol == ""],
  RefSeq_ID = refseq[symbol == ""],
  stringsAsFactors = FALSE
)

cat("\n====================================================\n")
cat("PRIMERAS SONDAS SIN Symbol\n")
cat("====================================================\n\n")

print(
  head(
    sin_symbol,
    20
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Tabla de multiplicidad por símbolo
# ------------------------------------------------------------

multiplicidad <- data.frame(
  Symbol = names(tabla_symbol),
  n_sondas = as.integer(tabla_symbol),
  stringsAsFactors = FALSE
)

multiplicidad <- multiplicidad[
  order(
    -multiplicidad$n_sondas,
    multiplicidad$Symbol
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Guardar auditorías
# ------------------------------------------------------------

write.csv(
  multiplicidad,
  "GPL6947_multiplicidad_sondas_por_symbol.csv",
  row.names = FALSE
)

write.csv(
  sin_symbol,
  "GPL6947_sondas_sin_symbol.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("- GPL6947_multiplicidad_sondas_por_symbol.csv\n")
cat("- GPL6947_sondas_sin_symbol.csv\n")


# ------------------------------------------------------------
# 13. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AUDITORÍA DE SÍMBOLOS COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: todavía NO se colapsaron sondas a genes.\n"
)