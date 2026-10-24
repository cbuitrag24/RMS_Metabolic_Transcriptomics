# ============================================================
# 07_preparar_diseno_GSE28511.R
#
# Objetivo:
# - recuperar la clasificación de las 24 muestras desde GEO
# - alinear metadata y matriz de expresión
# - definir cuatro grupos:
#     Skeletal_muscle
#     Tumor_adjacent_muscle
#     ERMS
#     ARMS
# - preparar el diseño para análisis posteriores
#
# IMPORTANTE:
# - NO combina músculo normal con músculo adyacente a tumor
# - NO interpreta ERMS/ARMS como estado de fusión
# - NO realiza expresión diferencial
# ============================================================


# ------------------------------------------------------------
# 1. Archivos necesarios
# ------------------------------------------------------------

archivo_series <- "GSE28511_series_matrix.txt.gz"
archivo_expr <- "GSE28511_expression_quantileNormalized_log2.rds"

if (!file.exists(archivo_series)) {
  stop(paste("ERROR: no se encontró", archivo_series))
}

if (!file.exists(archivo_expr)) {
  stop(paste("ERROR: no se encontró", archivo_expr))
}


# ------------------------------------------------------------
# 2. Leer matriz de expresión
# ------------------------------------------------------------

expr <- readRDS(archivo_expr)

cat("Dimensiones de la matriz de expresión:\n")
print(dim(expr))

if (ncol(expr) != 24) {
  stop("ERROR: se esperaban 24 muestras.")
}

muestras_expr <- colnames(expr)

if (is.null(muestras_expr)) {
  stop("ERROR: la matriz no tiene nombres de muestras.")
}

if (anyDuplicated(muestras_expr) > 0) {
  stop("ERROR: existen nombres de muestra duplicados.")
}


# ------------------------------------------------------------
# 3. Leer líneas de metadata de la Series Matrix
# ------------------------------------------------------------

lineas <- readLines(
  gzfile(archivo_series),
  warn = FALSE
)

extraer_campo <- function(prefijo) {

  x <- grep(
    paste0("^", prefijo),
    lineas,
    value = TRUE
  )

  if (length(x) != 1) {
    stop(
      paste(
        "ERROR: se esperaba una sola línea para",
        prefijo,
        "y se encontraron",
        length(x)
      )
    )
  }

  campos <- strsplit(
    x,
    "\t",
    fixed = TRUE
  )[[1]]

  campos <- campos[-1]

  campos <- gsub(
    '^"|"$',
    "",
    campos
  )

  campos
}


gsm <- extraer_campo("!Sample_geo_accession")
titulo <- extraer_campo("!Sample_title")
source <- extraer_campo("!Sample_source_name_ch1")


# ------------------------------------------------------------
# 4. Verificar longitudes
# ------------------------------------------------------------

cat("\nNúmero de GSM recuperados:\n")
cat(length(gsm), "\n")

cat("\nNúmero de títulos recuperados:\n")
cat(length(titulo), "\n")

cat("\nNúmero de source_name recuperados:\n")
cat(length(source), "\n")

if (
  length(gsm) != 24 ||
  length(titulo) != 24 ||
  length(source) != 24
) {
  stop("ERROR: metadata incompleta o inesperada.")
}

if (anyDuplicated(gsm) > 0) {
  stop("ERROR: existen GSM duplicados en metadata.")
}


# ------------------------------------------------------------
# 5. Construir metadata inicial
# ------------------------------------------------------------

metadata <- data.frame(
  samplename = gsm,
  title = titulo,
  source_name = source,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 6. Alinear metadata con expresión
# ------------------------------------------------------------

indice <- match(
  muestras_expr,
  metadata$samplename
)

if (anyNA(indice)) {
  stop(
    "ERROR: existen muestras de expresión sin metadata."
  )
}

metadata <- metadata[
  indice,
  ,
  drop = FALSE
]

if (!identical(
  as.character(metadata$samplename),
  as.character(muestras_expr)
)) {
  stop("ERROR: metadata y expresión no quedaron alineadas.")
}

cat(
  "\nVERIFICACIÓN: metadata y expresión correctamente alineadas.\n"
)


# ------------------------------------------------------------
# 7. Mostrar valores exactos de source_name
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("VALORES EXACTOS DE source_name\n")
cat("====================================================\n\n")

print(
  table(metadata$source_name)
)


# ------------------------------------------------------------
# 8. Definir grupos conservando la clasificación de GEO
# ------------------------------------------------------------

metadata$grupo <- NA_character_

metadata$grupo[
  metadata$source_name == "skeletal muscle"
] <- "Skeletal_muscle"

metadata$grupo[
  metadata$source_name == "tumor adjacent skeletal muscle"
] <- "Tumor_adjacent_muscle"

metadata$grupo[
  metadata$source_name == "ERMS"
] <- "ERMS"

metadata$grupo[
  metadata$source_name == "ARMS"
] <- "ARMS"


# ------------------------------------------------------------
# 9. Verificar clasificación completa
# ------------------------------------------------------------

if (anyNA(metadata$grupo)) {

  cat("\nSource_name no clasificados:\n")

  print(
    unique(
      metadata$source_name[
        is.na(metadata$grupo)
      ]
    )
  )

  stop(
    "ERROR: existen muestras que no pudieron clasificarse."
  )
}


metadata$grupo <- factor(
  metadata$grupo,
  levels = c(
    "Skeletal_muscle",
    "Tumor_adjacent_muscle",
    "ERMS",
    "ARMS"
  )
)


cat("\n====================================================\n")
cat("NÚMERO DE MUESTRAS POR GRUPO\n")
cat("====================================================\n\n")

print(
  table(metadata$grupo)
)


# ------------------------------------------------------------
# 10. Mostrar clasificación muestra por muestra
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("CLASIFICACIÓN DE LAS 24 MUESTRAS\n")
cat("====================================================\n\n")

print(
  metadata,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Crear diseño de cuatro grupos
#
# Sin intercepto:
# cada columna representa directamente un grupo.
# ------------------------------------------------------------

design <- model.matrix(
  ~ 0 + grupo,
  data = metadata
)

colnames(design) <- levels(
  metadata$grupo
)

rownames(design) <- metadata$samplename


# ------------------------------------------------------------
# 12. Verificaciones del diseño
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("MATRIZ DE DISEÑO\n")
cat("====================================================\n")

cat("\nDimensiones:\n")
print(dim(design))

cat("\nColumnas:\n")
print(colnames(design))

cat("\nRango de la matriz:\n")
cat(
  qr(design)$rank,
  "de",
  ncol(design),
  "\n"
)

if (qr(design)$rank != ncol(design)) {
  stop("ERROR: la matriz de diseño no tiene rango completo.")
}

if (!identical(
  rownames(design),
  colnames(expr)
)) {
  stop(
    "ERROR: el orden del diseño no coincide con la matriz de expresión."
  )
}

cat(
  "\nVERIFICACIÓN: diseño de rango completo y",
  "perfectamente alineado con expresión.\n"
)


# ------------------------------------------------------------
# 13. Guardar
# ------------------------------------------------------------

write.csv(
  metadata,
  "GSE28511_diseno_muestras.csv",
  row.names = FALSE
)

saveRDS(
  design,
  "GSE28511_design_matrix.rds"
)

cat("\nArchivos generados:\n")
cat("- GSE28511_diseno_muestras.csv\n")
cat("- GSE28511_design_matrix.rds\n")


# ------------------------------------------------------------
# 14. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("PREPARACIÓN DEL DISEÑO COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: todavía NO se realizó",
  "expresión diferencial.\n"
)