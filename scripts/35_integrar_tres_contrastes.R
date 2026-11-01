# ============================================================
# 35_integrar_tres_contrastes.R
# Integración de tres contrastes:
#   1) FN vs Normal
#   2) FOXO1-FP vs Normal
#   3) FOXO1-FP vs FN
#
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

directo <- read.csv(
  "GPL16791_DE_FOXO1FP_vs_FN_todos.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Tablas cargadas.\n")

cat("\nNúmero de genes por contraste:\n")
cat("FN vs Normal:", nrow(fn), "\n")
cat("FOXO1-FP vs Normal:", nrow(fp), "\n")
cat("FOXO1-FP vs FN:", nrow(directo), "\n")


# ------------------------------------------------------------
# 2. Verificaciones básicas
# ------------------------------------------------------------

columnas_necesarias <- c(
  "gene",
  "logFC",
  "PValue",
  "FDR"
)

verificar_tabla <- function(x, nombre) {

  if (!all(columnas_necesarias %in% colnames(x))) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- faltan columnas necesarias."
      )
    )
  }

  if (anyNA(x$gene) || any(x$gene == "")) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- existen genes inválidos."
      )
    )
  }

  if (anyDuplicated(x$gene) > 0) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- existen genes duplicados."
      )
    )
  }

  if (anyNA(x$logFC) ||
      anyNA(x$PValue) ||
      anyNA(x$FDR)) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- existen NA en estadísticas."
      )
    )
  }
}

verificar_tabla(fn, "FN vs Normal")
verificar_tabla(fp, "FOXO1-FP vs Normal")
verificar_tabla(directo, "FOXO1-FP vs FN")

cat(
  "\nVERIFICACIÓN: las tres tablas tienen identidad génica válida.\n"
)


# ------------------------------------------------------------
# 3. Integrar estrictamente por símbolo génico
# ------------------------------------------------------------

a <- fn[
  ,
  c("gene", "logFC", "PValue", "FDR")
]

colnames(a) <- c(
  "gene",
  "logFC_FN_vs_Normal",
  "PValue_FN_vs_Normal",
  "FDR_FN_vs_Normal"
)

b <- fp[
  ,
  c("gene", "logFC", "PValue", "FDR")
]

colnames(b) <- c(
  "gene",
  "logFC_FOXO1FP_vs_Normal",
  "PValue_FOXO1FP_vs_Normal",
  "FDR_FOXO1FP_vs_Normal"
)

c <- directo[
  ,
  c("gene", "logFC", "PValue", "FDR")
]

colnames(c) <- c(
  "gene",
  "logFC_FOXO1FP_vs_FN",
  "PValue_FOXO1FP_vs_FN",
  "FDR_FOXO1FP_vs_FN"
)

integrado <- merge(
  a,
  b,
  by = "gene",
  all = FALSE,
  sort = FALSE
)

integrado <- merge(
  integrado,
  c,
  by = "gene",
  all = FALSE,
  sort = FALSE
)

cat("\nGenes presentes en los tres contrastes:\n")
print(nrow(integrado))

if (nrow(integrado) != nrow(fn) ||
    nrow(integrado) != nrow(fp) ||
    nrow(integrado) != nrow(directo)) {
  stop(
    "ERROR: los tres contrastes no contienen exactamente el mismo conjunto de genes."
  )
}

if (anyDuplicated(integrado$gene) > 0) {
  stop("ERROR: genes duplicados después de integrar.")
}


# ------------------------------------------------------------
# 4. Comprobar coherencia algebraica de los logFC
#
# Por definición:
#
# FOXO1FP - FN =
# (FOXO1FP - Normal) - (FN - Normal)
# ------------------------------------------------------------

integrado$logFC_directo_esperado <- (
  integrado$logFC_FOXO1FP_vs_Normal -
  integrado$logFC_FN_vs_Normal
)

integrado$error_algebraico <- abs(
  integrado$logFC_FOXO1FP_vs_FN -
  integrado$logFC_directo_esperado
)

cat("\nResumen del error algebraico de logFC:\n")
print(summary(integrado$error_algebraico))

cat("\nMáximo error algebraico:\n")
print(max(integrado$error_algebraico))

if (max(integrado$error_algebraico) > 1e-8) {
  stop(
    "ERROR: los logFC de los tres contrastes no son algebraicamente coherentes."
  )
}

cat(
  "\nVERIFICACIÓN: los logFC de los tres contrastes son algebraicamente coherentes.\n"
)


# ------------------------------------------------------------
# 5. Criterio de efecto fuerte usado previamente
#
# FDR < 0.05 y |logFC| >= 1
# ------------------------------------------------------------

integrado$criterio_FN_vs_Normal <- (
  integrado$FDR_FN_vs_Normal < 0.05 &
  abs(integrado$logFC_FN_vs_Normal) >= 1
)

integrado$criterio_FOXO1FP_vs_Normal <- (
  integrado$FDR_FOXO1FP_vs_Normal < 0.05 &
  abs(integrado$logFC_FOXO1FP_vs_Normal) >= 1
)

integrado$criterio_FOXO1FP_vs_FN <- (
  integrado$FDR_FOXO1FP_vs_FN < 0.05 &
  abs(integrado$logFC_FOXO1FP_vs_FN) >= 1
)


# ------------------------------------------------------------
# 6. Señal compartida RMS vs Normal
#
# Requiere:
# - criterio en FN vs Normal
# - criterio en FOXO1-FP vs Normal
# - misma dirección frente a Normal
#
# NO exige ausencia de diferencia directa entre FN y FP.
# ------------------------------------------------------------

integrado$compartido_RMS_vs_Normal <- (
  integrado$criterio_FN_vs_Normal &
  integrado$criterio_FOXO1FP_vs_Normal &
  sign(integrado$logFC_FN_vs_Normal) ==
    sign(integrado$logFC_FOXO1FP_vs_Normal)
)

integrado$compartido_mayor_RMS <- (
  integrado$compartido_RMS_vs_Normal &
  integrado$logFC_FN_vs_Normal > 0 &
  integrado$logFC_FOXO1FP_vs_Normal > 0
)

integrado$compartido_mayor_Normal <- (
  integrado$compartido_RMS_vs_Normal &
  integrado$logFC_FN_vs_Normal < 0 &
  integrado$logFC_FOXO1FP_vs_Normal < 0
)


# ------------------------------------------------------------
# 7. Diferencia directa fuerte FOXO1-FP vs FN
# ------------------------------------------------------------

integrado$mayor_FOXO1FP_vs_FN <- (
  integrado$criterio_FOXO1FP_vs_FN &
  integrado$logFC_FOXO1FP_vs_FN > 0
)

integrado$mayor_FN_vs_FOXO1FP <- (
  integrado$criterio_FOXO1FP_vs_FN &
  integrado$logFC_FOXO1FP_vs_FN < 0
)


# ------------------------------------------------------------
# 8. Cruce de ambos componentes
#
# Esto permite saber cuántos genes compartidos frente a Normal
# también presentan una diferencia fuerte entre los dos grupos RMS.
# ------------------------------------------------------------

compartido_y_directo <- (
  integrado$compartido_RMS_vs_Normal &
  integrado$criterio_FOXO1FP_vs_FN
)

compartido_sin_directo_fuerte <- (
  integrado$compartido_RMS_vs_Normal &
  !integrado$criterio_FOXO1FP_vs_FN
)


# ------------------------------------------------------------
# 9. Resumen cuantitativo
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("SEÑAL COMPARTIDA RMS vs NORMAL\n")
cat("====================================================\n")

cat("\nCompartidos RMS vs Normal, misma dirección:\n")
print(sum(integrado$compartido_RMS_vs_Normal))

cat("\nCompartidos con mayor expresión en ambos RMS:\n")
print(sum(integrado$compartido_mayor_RMS))

cat("\nCompartidos con mayor expresión en Normal:\n")
print(sum(integrado$compartido_mayor_Normal))


cat("\n====================================================\n")
cat("DIFERENCIA DIRECTA FOXO1-FP vs FN\n")
cat("====================================================\n")

cat("\nGenes con diferencia directa fuerte:\n")
print(sum(integrado$criterio_FOXO1FP_vs_FN))

cat("\nMayor expresión en FOXO1-FP:\n")
print(sum(integrado$mayor_FOXO1FP_vs_FN))

cat("\nMayor expresión en FN:\n")
print(sum(integrado$mayor_FN_vs_FOXO1FP))


cat("\n====================================================\n")
cat("CRUCE ENTRE AMBOS COMPONENTES\n")
cat("====================================================\n")

cat(
  "\nGenes compartidos RMS vs Normal que además difieren fuertemente entre FOXO1-FP y FN:\n"
)
print(sum(compartido_y_directo))

cat(
  "\nGenes compartidos RMS vs Normal sin diferencia directa fuerte entre FOXO1-FP y FN:\n"
)
print(sum(compartido_sin_directo_fuerte))


# ------------------------------------------------------------
# 10. Categoría descriptiva
#
# Las categorías son estadísticas/descriptivas.
# No implican todavía mecanismo biológico.
# ------------------------------------------------------------

integrado$categoria <- "otros"

integrado$categoria[
  integrado$compartido_mayor_RMS &
  !integrado$criterio_FOXO1FP_vs_FN
] <- "compartido_mayor_RMS_sin_diferencia_directa_fuerte"

integrado$categoria[
  integrado$compartido_mayor_Normal &
  !integrado$criterio_FOXO1FP_vs_FN
] <- "compartido_mayor_Normal_sin_diferencia_directa_fuerte"

integrado$categoria[
  integrado$compartido_RMS_vs_Normal &
  integrado$criterio_FOXO1FP_vs_FN
] <- "compartido_RMS_Normal_y_diferente_FN_FOXO1FP"

integrado$categoria[
  !integrado$compartido_RMS_vs_Normal &
  integrado$mayor_FOXO1FP_vs_FN
] <- "diferencia_directa_mayor_FOXO1FP"

integrado$categoria[
  !integrado$compartido_RMS_vs_Normal &
  integrado$mayor_FN_vs_FOXO1FP
] <- "diferencia_directa_mayor_FN"


cat("\nNúmero de genes por categoría descriptiva:\n")
print(
  sort(
    table(integrado$categoria),
    decreasing = TRUE
  )
)


# ------------------------------------------------------------
# 11. Guardar tabla integrada
# ------------------------------------------------------------

integrado <- integrado[
  order(integrado$gene),
  ,
  drop = FALSE
]

rownames(integrado) <- NULL

write.csv(
  integrado,
  "integracion_tres_contrastes_FN_FOXO1FP_Normal.csv",
  row.names = FALSE
)

cat("\nArchivo generado:\n")
cat(
  "- integracion_tres_contrastes_FN_FOXO1FP_Normal.csv\n"
)

cat(
  "\nRESULTADO FINAL: integración de los tres contrastes completada y verificada.\n"
)