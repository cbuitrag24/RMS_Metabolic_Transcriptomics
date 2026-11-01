# ============================================================
# 27_DE_sensibilidad_106_muestras.R
# Análisis de sensibilidad:
# RMS vs músculo normal en las 106 muestras
# ajustado por plataforma
# GSE108022
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar ajuste QL del conjunto completo
# ------------------------------------------------------------

fit <- readRDS("GSE108022_edgeR_QLfit.rds")

cat("Ajuste completo cargado.\n")

cat("\nNúmero de genes:\n")
print(nrow(fit$coefficients))

cat("\nColumnas del diseño:\n")
print(colnames(fit$design))

if (!"conditionRMS" %in% colnames(fit$design)) {
  stop("ERROR: no existe el coeficiente conditionRMS.")
}

if (is.null(fit$genes)) {
  stop("ERROR: fit$genes no está disponible.")
}

if (!"gene" %in% colnames(fit$genes)) {
  stop("ERROR: fit$genes no contiene la columna gene.")
}

if (nrow(fit$genes) != nrow(fit$coefficients)) {
  stop("ERROR: el número de genes no coincide con los coeficientes.")
}


# ------------------------------------------------------------
# 2. Test QL
#
# conditionRMS = RMS - Normal
#
# logFC > 0 -> mayor expresión en RMS
# logFC < 0 -> mayor expresión en músculo normal
# ------------------------------------------------------------

qlf <- glmQLFTest(
  fit,
  coef = "conditionRMS"
)

cat("\nTest QL completado.\n")


# ------------------------------------------------------------
# 3. Construir tabla preservando identidad génica
# ------------------------------------------------------------

if (nrow(qlf$table) != nrow(qlf$genes)) {
  stop("ERROR: qlf$table y qlf$genes tienen distinto número de filas.")
}

resultado_base <- data.frame(
  gene = qlf$genes$gene,
  qlf$table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!identical(
  as.character(resultado_base$gene),
  as.character(qlf$genes$gene)
)) {
  stop("ERROR: se perdió la correspondencia gen-fila.")
}

cat(
  "\nVERIFICACIÓN: identidad gen-fila preservada antes del ordenamiento.\n"
)


# ------------------------------------------------------------
# 4. Calcular FDR
# ------------------------------------------------------------

resultado_base$FDR <- p.adjust(
  resultado_base$PValue,
  method = "BH"
)

if (anyNA(resultado_base$logFC) ||
    anyNA(resultado_base$PValue) ||
    anyNA(resultado_base$FDR)) {
  stop("ERROR: existen NA en logFC, PValue o FDR.")
}


# ------------------------------------------------------------
# 5. Ordenar por PValue
# ------------------------------------------------------------

orden <- order(
  resultado_base$PValue,
  decreasing = FALSE
)

resultado <- resultado_base[
  orden,
  ,
  drop = FALSE
]

rownames(resultado) <- NULL


# ------------------------------------------------------------
# 6. Verificación independiente con topTags
# ------------------------------------------------------------

tt <- topTags(
  qlf,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

if (!isTRUE(all.equal(
  resultado$logFC,
  tt$logFC,
  tolerance = 1e-12
))) {
  stop("ERROR: logFC no coincide con topTags.")
}

if (!isTRUE(all.equal(
  resultado$PValue,
  tt$PValue,
  tolerance = 1e-12
))) {
  stop("ERROR: PValue no coincide con topTags.")
}

if (!isTRUE(all.equal(
  resultado$FDR,
  tt$FDR,
  tolerance = 1e-12
))) {
  stop("ERROR: FDR no coincide con topTags.")
}

cat(
  "VERIFICACIÓN: logFC, PValue y FDR coinciden con topTags.\n"
)


# ------------------------------------------------------------
# 7. Resumen
# ------------------------------------------------------------

cat("\nNúmero total de genes:\n")
print(nrow(resultado))

cat("\nResumen de logFC:\n")
print(summary(resultado$logFC))

cat("\nGenes con FDR < 0.05:\n")
print(sum(resultado$FDR < 0.05))

cat("\nGenes con FDR < 0.01:\n")
print(sum(resultado$FDR < 0.01))

cat("\nFDR < 0.05 y |logFC| >= 1:\n")
print(
  sum(
    resultado$FDR < 0.05 &
    abs(resultado$logFC) >= 1
  )
)

cat("\nFDR < 0.05 y logFC >= 1 (mayor en RMS):\n")
print(
  sum(
    resultado$FDR < 0.05 &
    resultado$logFC >= 1
  )
)

cat("\nFDR < 0.05 y logFC <= -1 (mayor en normal):\n")
print(
  sum(
    resultado$FDR < 0.05 &
    resultado$logFC <= -1
  )
)


# ------------------------------------------------------------
# 8. Mostrar 20 genes con menor PValue
# ------------------------------------------------------------

cat("\n20 genes con menor PValue:\n")

print(
  resultado[
    seq_len(min(20, nrow(resultado))),
    c("gene", "logFC", "logCPM", "F", "PValue", "FDR")
  ]
)


# ------------------------------------------------------------
# 9. Guardar resultados
# ------------------------------------------------------------

write.csv(
  resultado,
  "ALL106_DE_RMS_vs_Normal_ajustado_plataforma_todos.csv",
  row.names = FALSE
)

significativos <- resultado[
  resultado$FDR < 0.05,
  ,
  drop = FALSE
]

write.csv(
  significativos,
  "ALL106_DE_RMS_vs_Normal_ajustado_plataforma_FDR005.csv",
  row.names = FALSE
)

significativos_fc1 <- resultado[
  resultado$FDR < 0.05 &
  abs(resultado$logFC) >= 1,
  ,
  drop = FALSE
]

write.csv(
  significativos_fc1,
  "ALL106_DE_RMS_vs_Normal_ajustado_plataforma_FDR005_logFC1.csv",
  row.names = FALSE
)

saveRDS(
  qlf,
  "GSE108022_ALL106_QLFTest_RMS_vs_Normal.rds"
)

cat("\nArchivos generados:\n")
cat("- ALL106_DE_RMS_vs_Normal_ajustado_plataforma_todos.csv\n")
cat("- ALL106_DE_RMS_vs_Normal_ajustado_plataforma_FDR005.csv\n")
cat("- ALL106_DE_RMS_vs_Normal_ajustado_plataforma_FDR005_logFC1.csv\n")
cat("- GSE108022_ALL106_QLFTest_RMS_vs_Normal.rds\n")

cat("\nRESULTADO: análisis de sensibilidad completado con identidad génica verificada.\n")