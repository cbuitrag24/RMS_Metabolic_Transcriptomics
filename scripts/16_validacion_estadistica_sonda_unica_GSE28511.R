# ============================================================
# 16_validacion_estadistica_sonda_unica_GSE28511.R
#
# OBJETIVO:
# Evaluar validación estadística independiente de los genes
# fuertes definidos EXCLUSIVAMENTE en GSE108022.
#
# Universo analizado:
# genes con una única sonda en GSE28511.
#
# Criterio de descubrimiento (GSE108022):
#   FDR < 0.05
#   |logFC| >= 1
#
# Criterio de validación (GSE28511):
#   misma dirección del efecto
#   FDR < 0.05
#
# IMPORTANTE:
# - NO se exige |logFC| >= 1 en GSE28511.
# - NO se modifica el criterio según los resultados.
# - NO se analizan aquí genes con múltiples sondas.
# ============================================================


# ------------------------------------------------------------
# 1. Archivo de concordancia previamente verificado
# ------------------------------------------------------------

archivo <-
  "concordancia_GSE108022_GSE28511_sonda_unica.rds"

if (!file.exists(archivo)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo
    )
  )
}


# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

tabla <- readRDS(
  archivo
)

cat("Dimensiones de la tabla integrada:\n")
print(dim(tabla))


# ------------------------------------------------------------
# 3. Verificar columnas necesarias
# ------------------------------------------------------------

necesarias <- c(
  "gene",
  "logFC_GSE108022",
  "FDR_GSE108022",
  "logFC_GSE28511",
  "FDR_GSE28511"
)

faltantes <- setdiff(
  necesarias,
  colnames(tabla)
)

if (length(faltantes) > 0) {
  stop(
    paste(
      "ERROR: faltan columnas:",
      paste(faltantes, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 4. Controles básicos
# ------------------------------------------------------------

if (nrow(tabla) != 9428) {
  stop(
    paste(
      "ERROR: se esperaban 9428 genes y se encontraron",
      nrow(tabla)
    )
  )
}

if (anyDuplicated(tabla$gene)) {
  stop(
    "ERROR: existen genes duplicados."
  )
}

columnas_numericas <- c(
  "logFC_GSE108022",
  "FDR_GSE108022",
  "logFC_GSE28511",
  "FDR_GSE28511"
)

for (col in columnas_numericas) {

  x <- tabla[[col]]

  if (anyNA(x)) {
    stop(
      paste(
        "ERROR: existen NA en",
        col
      )
    )
  }

  if (any(!is.finite(x))) {
    stop(
      paste(
        "ERROR: existen valores no finitos en",
        col
      )
    )
  }
}

cat(
  "\nVERIFICACIÓN: tabla de 9428 genes",
  "sin duplicados ni valores inválidos.\n"
)


# ------------------------------------------------------------
# 5. Recalcular conjunto fuerte de descubrimiento
#
# Definido exclusivamente en GSE108022.
# ------------------------------------------------------------

fuerte_descubrimiento <- (
  tabla$FDR_GSE108022 < 0.05 &
  abs(tabla$logFC_GSE108022) >= 1
)

n_fuertes <- sum(
  fuerte_descubrimiento
)


cat("\n====================================================\n")
cat("CONJUNTO FUERTE DE DESCUBRIMIENTO\n")
cat("====================================================\n")

cat("\nCriterio:\n")
cat("GSE108022 FDR < 0.05 y |logFC| >= 1\n")

cat("\nGenes fuertes:\n")
cat(n_fuertes, "\n")


if (n_fuertes != 3715) {
  stop(
    paste(
      "ERROR: se esperaban 3715 genes fuertes",
      "según el script 15 y se encontraron",
      n_fuertes
    )
  )
}


# ------------------------------------------------------------
# 6. Extraer genes fuertes
# ------------------------------------------------------------

fuertes <- tabla[
  fuerte_descubrimiento,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Determinar concordancia direccional
# ------------------------------------------------------------

fuertes$misma_direccion <- (
  sign(fuertes$logFC_GSE108022) ==
  sign(fuertes$logFC_GSE28511)
)


cat("\n====================================================\n")
cat("REPLICACIÓN DIRECCIONAL\n")
cat("====================================================\n")

cat("\nGenes fuertes evaluados:\n")
cat(nrow(fuertes), "\n")

cat("\nMisma dirección:\n")
cat(
  sum(fuertes$misma_direccion),
  "\n"
)

cat("\nDirección opuesta:\n")
cat(
  sum(!fuertes$misma_direccion),
  "\n"
)

cat("\nProporción misma dirección:\n")
cat(
  mean(fuertes$misma_direccion),
  "\n"
)


# ------------------------------------------------------------
# 8. Evidencia estadística en GSE28511
#
# FDR ya fue calculado sobre TODAS las 48.802 sondas
# en el análisis limma.
#
# No se recalcula FDR solamente sobre los genes seleccionados.
# ------------------------------------------------------------

fuertes$significativo_GSE28511 <- (
  fuertes$FDR_GSE28511 < 0.05
)


cat("\n====================================================\n")
cat("EVIDENCIA ESTADÍSTICA EN GSE28511\n")
cat("====================================================\n")

cat("\nGenes fuertes con FDR < 0.05 en GSE28511:\n")
cat(
  sum(fuertes$significativo_GSE28511),
  "\n"
)

cat("\nProporción:\n")
cat(
  mean(fuertes$significativo_GSE28511),
  "\n"
)


# ------------------------------------------------------------
# 9. Criterio predefinido de validación
#
# Validado =
# - fuerte en descubrimiento
# - misma dirección
# - FDR < 0.05 en GSE28511
# ------------------------------------------------------------

fuertes$validado <- (
  fuertes$misma_direccion &
  fuertes$significativo_GSE28511
)


n_validados <- sum(
  fuertes$validado
)


cat("\n====================================================\n")
cat("VALIDACIÓN ESTADÍSTICA PREDEFINIDA\n")
cat("====================================================\n")

cat("\nCriterio:\n")
cat(
  "fuerte en GSE108022 + misma dirección",
  "+ FDR < 0.05 en GSE28511\n"
)

cat("\nGenes que cumplen el criterio:\n")
cat(n_validados, "\n")

cat("\nProporción sobre los 3715 genes fuertes:\n")
cat(
  n_validados / nrow(fuertes),
  "\n"
)


# ------------------------------------------------------------
# 10. Separar según dirección original del descubrimiento
# ------------------------------------------------------------

fuertes$direccion_GSE108022 <- ifelse(
  fuertes$logFC_GSE108022 > 0,
  "Mayor_en_RMS",
  "Mayor_en_musculo"
)


cat("\n====================================================\n")
cat("VALIDACIÓN SEGÚN DIRECCIÓN DEL DESCUBRIMIENTO\n")
cat("====================================================\n")

cat("\nGenes fuertes por dirección:\n")
print(
  table(
    fuertes$direccion_GSE108022
  )
)

cat("\nGenes validados por dirección:\n")
print(
  table(
    fuertes$direccion_GSE108022[
      fuertes$validado
    ]
  )
)


# ------------------------------------------------------------
# 11. Tabla cruzada completa
#
# Permite distinguir:
# - misma dirección + significativo
# - misma dirección + no significativo
# - dirección opuesta + significativo
# - dirección opuesta + no significativo
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("TABLA CRUZADA DIRECCIÓN x SIGNIFICACIÓN\n")
cat("====================================================\n\n")

print(
  table(
    Misma_direccion =
      fuertes$misma_direccion,
    FDR005_GSE28511 =
      fuertes$significativo_GSE28511
  )
)


# ------------------------------------------------------------
# 12. Verificación interna
# ------------------------------------------------------------

n_misma <- sum(
  fuertes$misma_direccion
)

n_signif <- sum(
  fuertes$significativo_GSE28511
)

n_validado_control <- sum(
  fuertes$misma_direccion &
  fuertes$FDR_GSE28511 < 0.05
)


if (n_misma != 3197) {
  stop(
    paste(
      "ERROR: la concordancia direccional",
      "no reproduce los 3197 genes del script 15."
    )
  )
}

if (n_validados != n_validado_control) {
  stop(
    "ERROR: inconsistencia en el número de genes validados."
  )
}


cat(
  "\nVERIFICACIÓN: los 3197 genes con misma",
  "dirección reproducen exactamente el script 15.\n"
)

cat(
  "VERIFICACIÓN: criterio de validación",
  "aplicado de forma consistente.\n"
)


# ------------------------------------------------------------
# 13. Guardar tabla de los 3715 genes fuertes
# ------------------------------------------------------------

write.csv(
  fuertes,
  "GSE28511_validacion_genes_fuertes_sonda_unica.csv",
  row.names = FALSE
)

saveRDS(
  fuertes,
  "GSE28511_validacion_genes_fuertes_sonda_unica.rds"
)


# ------------------------------------------------------------
# 14. Guardar solamente los genes que cumplen
# el criterio predefinido
# ------------------------------------------------------------

validados <- fuertes[
  fuertes$validado,
  ,
  drop = FALSE
]

write.csv(
  validados,
  "GSE28511_genes_validacion_estadistica_sonda_unica.csv",
  row.names = FALSE
)


cat("\nArchivos generados:\n")

cat(
  "- GSE28511_validacion_genes_fuertes_sonda_unica.csv\n"
)

cat(
  "- GSE28511_validacion_genes_fuertes_sonda_unica.rds\n"
)

cat(
  "- GSE28511_genes_validacion_estadistica_sonda_unica.csv\n"
)


# ------------------------------------------------------------
# 15. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("VALIDACIÓN ESTADÍSTICA COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nLos genes informados como validados cumplen",
  "el criterio definido antes de observar este resultado:\n"
)

cat(
  "GSE108022: FDR < 0.05 y |logFC| >= 1\n"
)

cat(
  "GSE28511: misma dirección y FDR < 0.05\n"
)

cat(
  "\nEste análisis está restringido a genes",
  "representados por una única sonda en GSE28511.\n"
)