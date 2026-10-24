# ============================================================
# 09_preparar_limma_primario_GSE28511.R
#
# Objetivo:
# - preparar el subconjunto primario de validación:
#       Skeletal_muscle = 3
#       RMS = 18 (8 ERMS + 10 ARMS)
# - excluir del contraste primario:
#       Tumor_adjacent_muscle = 3
# - mantener las sondas individualmente
# - construir y verificar el diseño para limma
#
# IMPORTANTE:
# - NO colapsa sondas a genes
# - NO vuelve a normalizar
# - NO realiza expresión diferencial
# ============================================================


# ------------------------------------------------------------
# 1. Archivos necesarios
# ------------------------------------------------------------

archivo_expr <- "GSE28511_expression_quantileNormalized_log2.rds"
archivo_metadata <- "GSE28511_diseno_muestras.csv"

if (!file.exists(archivo_expr)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_expr
    )
  )
}

if (!file.exists(archivo_metadata)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_metadata
    )
  )
}


# ------------------------------------------------------------
# 2. Leer matriz y metadata
# ------------------------------------------------------------

expr <- readRDS(
  archivo_expr
)

metadata <- read.csv(
  archivo_metadata,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensiones de expresión completa:\n")
print(dim(expr))

cat("\nDimensiones de metadata completa:\n")
print(dim(metadata))


# ------------------------------------------------------------
# 3. Verificar columnas necesarias
# ------------------------------------------------------------

columnas_necesarias <- c(
  "samplename",
  "title",
  "source_name",
  "grupo"
)

faltantes <- setdiff(
  columnas_necesarias,
  colnames(metadata)
)

if (length(faltantes) > 0) {
  stop(
    paste(
      "ERROR: faltan columnas en metadata:",
      paste(faltantes, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 4. Alinear metadata con expresión
# ------------------------------------------------------------

indice <- match(
  colnames(expr),
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
  as.character(colnames(expr))
)) {
  stop(
    "ERROR: metadata y expresión no están alineadas."
  )
}

cat(
  "\nVERIFICACIÓN: metadata completa y expresión",
  "correctamente alineadas.\n"
)


# ------------------------------------------------------------
# 5. Definir muestras del análisis primario
#
# Incluimos:
# - Skeletal_muscle
# - ERMS
# - ARMS
#
# Excluimos:
# - Tumor_adjacent_muscle
# ------------------------------------------------------------

usar <- metadata$grupo %in% c(
  "Skeletal_muscle",
  "ERMS",
  "ARMS"
)

metadata_primario <- metadata[
  usar,
  ,
  drop = FALSE
]

expr_primario <- expr[
  ,
  usar,
  drop = FALSE
]


# ------------------------------------------------------------
# 6. Verificar alineación del subconjunto
# ------------------------------------------------------------

if (!identical(
  as.character(metadata_primario$samplename),
  as.character(colnames(expr_primario))
)) {
  stop(
    "ERROR: el subconjunto de metadata y expresión no está alineado."
  )
}

cat("\n====================================================\n")
cat("SUBCONJUNTO PRIMARIO\n")
cat("====================================================\n")

cat("\nDimensiones de expresión primaria:\n")
print(dim(expr_primario))

cat("\nNúmero de muestras primarias:\n")
cat(ncol(expr_primario), "\n")

if (ncol(expr_primario) != 21) {
  stop(
    paste(
      "ERROR: se esperaban 21 muestras primarias y se encontraron",
      ncol(expr_primario)
    )
  )
}


# ------------------------------------------------------------
# 7. Crear condición binaria para validación primaria
#
# Skeletal_muscle = control
# RMS = ERMS + ARMS
#
# Esto NO implica que ERMS y ARMS sean biológicamente
# equivalentes. Se agrupan aquí exclusivamente para validar
# el contraste global RMS vs músculo esquelético.
# ------------------------------------------------------------

metadata_primario$condicion <- NA_character_

metadata_primario$condicion[
  metadata_primario$grupo == "Skeletal_muscle"
] <- "Skeletal_muscle"

metadata_primario$condicion[
  metadata_primario$grupo %in% c(
    "ERMS",
    "ARMS"
  )
] <- "RMS"


if (anyNA(metadata_primario$condicion)) {
  stop(
    "ERROR: existen muestras primarias sin condición asignada."
  )
}

metadata_primario$condicion <- factor(
  metadata_primario$condicion,
  levels = c(
    "Skeletal_muscle",
    "RMS"
  )
)


# ------------------------------------------------------------
# 8. Verificar composición
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("COMPOSICIÓN DEL ANÁLISIS PRIMARIO\n")
cat("====================================================\n")

cat("\nCondición:\n")
print(
  table(metadata_primario$condicion)
)

cat("\nGrupo original:\n")
print(
  table(metadata_primario$grupo)
)


n_normal <- sum(
  metadata_primario$condicion == "Skeletal_muscle"
)

n_rms <- sum(
  metadata_primario$condicion == "RMS"
)

n_erms <- sum(
  metadata_primario$grupo == "ERMS"
)

n_arms <- sum(
  metadata_primario$grupo == "ARMS"
)


if (
  n_normal != 3 ||
  n_rms != 18 ||
  n_erms != 8 ||
  n_arms != 10
) {
  stop(
    paste(
      "ERROR: composición inesperada.",
      "Skeletal_muscle =", n_normal,
      "RMS =", n_rms,
      "ERMS =", n_erms,
      "ARMS =", n_arms
    )
  )
}

cat(
  "\nVERIFICACIÓN: composición primaria exacta:",
  "3 Skeletal_muscle y 18 RMS",
  "(8 ERMS + 10 ARMS).\n"
)


# ------------------------------------------------------------
# 9. Verificar que Tumor_adjacent_muscle quedó fuera
# ------------------------------------------------------------

if (
  any(
    metadata_primario$grupo ==
      "Tumor_adjacent_muscle"
  )
) {
  stop(
    "ERROR: Tumor_adjacent_muscle entró en el análisis primario."
  )
}

cat(
  "\nVERIFICACIÓN:",
  "Tumor_adjacent_muscle está excluido",
  "del análisis primario.\n"
)


# ------------------------------------------------------------
# 10. Construir diseño para limma
#
# Sin intercepto:
# cada columna representa una condición.
# ------------------------------------------------------------

design_primario <- model.matrix(
  ~ 0 + condicion,
  data = metadata_primario
)

colnames(design_primario) <- levels(
  metadata_primario$condicion
)

rownames(design_primario) <- metadata_primario$samplename


# ------------------------------------------------------------
# 11. Verificar diseño
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("DISEÑO PRIMARIO PARA LIMMA\n")
cat("====================================================\n")

cat("\nDimensiones:\n")
print(dim(design_primario))

cat("\nColumnas:\n")
print(
  colnames(design_primario)
)

cat("\nRango:\n")
cat(
  qr(design_primario)$rank,
  "de",
  ncol(design_primario),
  "\n"
)


if (
  qr(design_primario)$rank !=
    ncol(design_primario)
) {
  stop(
    "ERROR: el diseño primario no tiene rango completo."
  )
}


if (!identical(
  rownames(design_primario),
  colnames(expr_primario)
)) {
  stop(
    "ERROR: diseño y expresión primaria no están alineados."
  )
}

cat(
  "\nVERIFICACIÓN: diseño primario de rango completo",
  "y perfectamente alineado con expresión.\n"
)


# ------------------------------------------------------------
# 12. Verificar matriz de expresión primaria
# ------------------------------------------------------------

if (anyNA(expr_primario)) {
  stop(
    "ERROR: existen NA en la matriz primaria."
  )
}

if (any(!is.finite(expr_primario))) {
  stop(
    "ERROR: existen valores no finitos en la matriz primaria."
  )
}


cat("\nRango de expresión log2 en subconjunto primario:\n")
print(
  range(expr_primario)
)


# ------------------------------------------------------------
# 13. Mostrar muestras incluidas
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("MUESTRAS INCLUIDAS EN EL ANÁLISIS PRIMARIO\n")
cat("====================================================\n\n")

print(
  metadata_primario[
    ,
    c(
      "samplename",
      "title",
      "source_name",
      "grupo",
      "condicion"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Guardar objetos
# ------------------------------------------------------------

saveRDS(
  expr_primario,
  "GSE28511_primario_expression_log2.rds"
)

saveRDS(
  design_primario,
  "GSE28511_primario_design_matrix.rds"
)

write.csv(
  metadata_primario,
  "GSE28511_primario_diseno_muestras.csv",
  row.names = FALSE
)


cat("\nArchivos generados:\n")

cat(
  "- GSE28511_primario_expression_log2.rds\n"
)

cat(
  "- GSE28511_primario_design_matrix.rds\n"
)

cat(
  "- GSE28511_primario_diseno_muestras.csv\n"
)


# ------------------------------------------------------------
# 15. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("PREPARACIÓN PRIMARIA PARA LIMMA COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: todavía NO se realizó",
  "expresión diferencial.\n"
)