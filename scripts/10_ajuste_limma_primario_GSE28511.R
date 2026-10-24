# ============================================================
# 10_ajuste_limma_primario_GSE28511.R
#
# Objetivo:
# - ajustar el modelo lineal con limma a nivel de sonda
# - contraste primario:
#       RMS - Skeletal_muscle
# - aplicar moderación empírica robusta
# - realizar controles del ajuste
#
# IMPORTANTE:
# - NO vuelve a normalizar
# - NO colapsa sondas a genes
# - NO selecciona genes diferencialmente expresados
# ============================================================


# ------------------------------------------------------------
# 1. Cargar limma
# ------------------------------------------------------------

suppressPackageStartupMessages(
  library(limma)
)

cat("Versión de limma:\n")
cat(as.character(packageVersion("limma")), "\n")


# ------------------------------------------------------------
# 2. Archivos necesarios
# ------------------------------------------------------------

archivo_expr <- "GSE28511_primario_expression_log2.rds"
archivo_design <- "GSE28511_primario_design_matrix.rds"
archivo_metadata <- "GSE28511_primario_diseno_muestras.csv"

archivos <- c(
  archivo_expr,
  archivo_design,
  archivo_metadata
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
# 3. Cargar objetos
# ------------------------------------------------------------

expr <- readRDS(
  archivo_expr
)

design <- readRDS(
  archivo_design
)

metadata <- read.csv(
  archivo_metadata,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nDimensiones de expresión:\n")
print(dim(expr))

cat("\nDimensiones del diseño:\n")
print(dim(design))

cat("\nDimensiones de metadata:\n")
print(dim(metadata))


# ------------------------------------------------------------
# 4. Controles de alineación
# ------------------------------------------------------------

if (!identical(
  colnames(expr),
  rownames(design)
)) {
  stop(
    "ERROR: expresión y diseño no están alineados."
  )
}

if (!identical(
  colnames(expr),
  as.character(metadata$samplename)
)) {
  stop(
    "ERROR: expresión y metadata no están alineadas."
  )
}

if (ncol(expr) != nrow(design)) {
  stop(
    "ERROR: número de muestras incompatible entre expresión y diseño."
  )
}

if (anyNA(expr)) {
  stop(
    "ERROR: existen NA en expresión."
  )
}

if (any(!is.finite(expr))) {
  stop(
    "ERROR: existen valores no finitos en expresión."
  )
}

cat(
  "\nVERIFICACIÓN: expresión, diseño y metadata",
  "perfectamente alineados.\n"
)


# ------------------------------------------------------------
# 5. Verificar diseño esperado
# ------------------------------------------------------------

columnas_esperadas <- c(
  "Skeletal_muscle",
  "RMS"
)

if (!identical(
  colnames(design),
  columnas_esperadas
)) {
  stop(
    paste(
      "ERROR: columnas inesperadas en el diseño:",
      paste(colnames(design), collapse = ", ")
    )
  )
}

if (qr(design)$rank != ncol(design)) {
  stop(
    "ERROR: la matriz de diseño no tiene rango completo."
  )
}

cat("\nDiseño:\n")
print(colnames(design))

cat(
  "Rango:",
  qr(design)$rank,
  "de",
  ncol(design),
  "\n"
)


# ------------------------------------------------------------
# 6. Ajuste lineal inicial
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AJUSTE LINEAL INICIAL\n")
cat("====================================================\n")

fit <- lmFit(
  expr,
  design
)

cat("\nDimensiones de coeficientes:\n")
print(dim(fit$coefficients))

cat("\nGrados de libertad residuales:\n")
print(
  summary(fit$df.residual)
)

if (anyNA(fit$coefficients)) {
  stop(
    "ERROR: existen NA en los coeficientes del modelo."
  )
}

cat(
  "\nVERIFICACIÓN: ajuste lineal inicial completado.\n"
)


# ------------------------------------------------------------
# 7. Definir contraste primario
#
# logFC positivo:
# mayor expresión en RMS
#
# logFC negativo:
# mayor expresión en Skeletal_muscle
# ------------------------------------------------------------

contraste <- makeContrasts(
  RMS_vs_Skeletal_muscle =
    RMS - Skeletal_muscle,
  levels = design
)

cat("\n====================================================\n")
cat("CONTRASTE PRIMARIO\n")
cat("====================================================\n\n")

print(contraste)


# ------------------------------------------------------------
# 8. Aplicar contraste
# ------------------------------------------------------------

fit_contraste <- contrasts.fit(
  fit,
  contraste
)


# ------------------------------------------------------------
# 9. Moderación empírica robusta
#
# robust = TRUE:
# limita la influencia de sondas con varianzas residuales
# extremas sobre la estimación empírica Bayesiana.
#
# No usamos trend=TRUE en este paso.
# ------------------------------------------------------------

fit_ebayes <- eBayes(
  fit_contraste,
  robust = TRUE
)


# ------------------------------------------------------------
# 10. Controles del ajuste moderado
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("CONTROL DEL AJUSTE MODERADO\n")
cat("====================================================\n")

cat("\nNúmero de sondas:\n")
cat(
  nrow(fit_ebayes$coefficients),
  "\n"
)

cat("\nResumen de logFC RMS - Skeletal_muscle:\n")
print(
  summary(
    fit_ebayes$coefficients[
      ,
      "RMS_vs_Skeletal_muscle"
    ]
  )
)

cat("\nResumen de sigma residual:\n")
print(
  summary(fit_ebayes$sigma)
)

cat("\nResumen de s2.post:\n")
print(
  summary(fit_ebayes$s2.post)
)

cat("\nResumen de df.prior:\n")
print(
  summary(fit_ebayes$df.prior)
)

cat("\nResumen de df.total:\n")
print(
  summary(fit_ebayes$df.total)
)


# ------------------------------------------------------------
# 11. Verificaciones de valores
# ------------------------------------------------------------

objetos_control <- list(
  logFC =
    fit_ebayes$coefficients[
      ,
      "RMS_vs_Skeletal_muscle"
    ],
  sigma = fit_ebayes$sigma,
  s2_post = fit_ebayes$s2.post,
  t = fit_ebayes$t[
    ,
    "RMS_vs_Skeletal_muscle"
  ],
  p_value = fit_ebayes$p.value[
    ,
    "RMS_vs_Skeletal_muscle"
  ]
)

for (nombre in names(objetos_control)) {

  x <- objetos_control[[nombre]]

  if (anyNA(x)) {
    stop(
      paste(
        "ERROR: existen NA en",
        nombre
      )
    )
  }

  if (any(!is.finite(x))) {
    stop(
      paste(
        "ERROR: existen valores no finitos en",
        nombre
      )
    )
  }
}

cat(
  "\nVERIFICACIÓN: no hay NA ni valores",
  "no finitos en los estadísticos principales.\n"
)


# ------------------------------------------------------------
# 12. Verificar dirección del contraste
# ------------------------------------------------------------

media_normal <- rowMeans(
  expr[
    ,
    metadata$condicion == "Skeletal_muscle",
    drop = FALSE
  ]
)

media_rms <- rowMeans(
  expr[
    ,
    metadata$condicion == "RMS",
    drop = FALSE
  ]
)

diferencia_directa <- media_rms - media_normal

error_contraste <- max(
  abs(
    diferencia_directa -
      fit_ebayes$coefficients[
        ,
        "RMS_vs_Skeletal_muscle"
      ]
  )
)

cat("\nMáximo error entre:")
cat("\nmedia(RMS) - media(Skeletal_muscle)")
cat("\ny el coeficiente limma:\n")
cat(
  format(
    error_contraste,
    scientific = TRUE
  ),
  "\n"
)

if (
  !is.finite(error_contraste) ||
  error_contraste > 1e-10
) {
  stop(
    "ERROR: la dirección del contraste no coincide con las medias directas."
  )
}

cat(
  "\nVERIFICACIÓN: dirección del contraste correcta.\n"
)

cat(
  "logFC > 0 = mayor expresión en RMS.\n"
)

cat(
  "logFC < 0 = mayor expresión en Skeletal_muscle.\n"
)


# ------------------------------------------------------------
# 13. Guardar ajuste
# ------------------------------------------------------------

saveRDS(
  fit_ebayes,
  "GSE28511_primario_limma_ebayes.rds"
)

saveRDS(
  contraste,
  "GSE28511_primario_contraste.rds"
)

cat("\nArchivos generados:\n")
cat(
  "- GSE28511_primario_limma_ebayes.rds\n"
)
cat(
  "- GSE28511_primario_contraste.rds\n"
)


# ------------------------------------------------------------
# 14. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("AJUSTE LIMMA PRIMARIO COMPLETADO\n")
cat("====================================================\n")

cat(
  "\nIMPORTANTE: todavía NO se aplicó ningún",
  "umbral de significación ni se seleccionaron sondas.\n"
)