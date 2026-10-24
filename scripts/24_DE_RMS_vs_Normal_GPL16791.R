# ============================================================
# 24_DE_RMS_vs_Normal_GPL16791.R
# Expresión diferencial RMS vs músculo normal
# Análisis primario: GPL16791
# GSE108022
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar ajuste QL
# ------------------------------------------------------------

fit <- readRDS("GSE108022_GPL16791_QLfit.rds")

cat("Ajuste cargado:\n")
cat("Genes:", nrow(fit), "\n")
cat("Muestras:", ncol(fit), "\n")

cat("\nColumnas del diseño:\n")
print(colnames(fit$design))

if (!"conditionRMS" %in% colnames(fit$design)) {
  stop("ERROR: no existe el coeficiente conditionRMS.")
}


# ------------------------------------------------------------
# 2. Test quasi-likelihood
#
# conditionRMS = RMS - Normal
#
# logFC > 0  -> mayor expresión en RMS
# logFC < 0  -> mayor expresión en músculo normal
# ------------------------------------------------------------

qlf <- glmQLFTest(
  fit,
  coef = "conditionRMS"
)

cat("\nTest QL completado para conditionRMS.\n")


# ------------------------------------------------------------
# 3. Obtener TODOS los genes
# ------------------------------------------------------------

resultado <- topTags(
  qlf,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

resultado$gene <- rownames(resultado)

# Colocar gene como primera columna
resultado <- resultado[
  ,
  c(
    "gene",
    setdiff(colnames(resultado), "gene")
  )
]

cat("\nNúmero total de genes en el resultado:\n")
print(nrow(resultado))

cat("\nColumnas del resultado:\n")
print(colnames(resultado))


# ------------------------------------------------------------
# 4. Controles básicos
# ------------------------------------------------------------

if (anyNA(resultado$PValue)) {
  stop("ERROR: existen PValue NA.")
}

if (anyNA(resultado$FDR)) {
  stop("ERROR: existen FDR NA.")
}

if (anyNA(resultado$logFC)) {
  stop("ERROR: existen logFC NA.")
}

cat("\nRESULTADO: no se detectaron NA en logFC, PValue ni FDR.\n")


# ------------------------------------------------------------
# 5. Resumen estadístico
# ------------------------------------------------------------

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
# 6. Mostrar los genes con menor PValue
#
# Solo inspección estadística inicial.
# No constituye todavía interpretación biológica.
# ------------------------------------------------------------

cat("\n20 genes con menor PValue:\n")

print(
  resultado[
    seq_len(min(20, nrow(resultado))),
    c("gene", "logFC", "logCPM", "F", "PValue", "FDR")
  ]
)


# ------------------------------------------------------------
# 7. Guardar resultados completos
# ------------------------------------------------------------

write.csv(
  resultado,
  "GPL16791_DE_RMS_vs_Normal_todos_los_genes.csv",
  row.names = FALSE
)

significativos <- resultado[
  resultado$FDR < 0.05,
  ,
  drop = FALSE
]

write.csv(
  significativos,
  "GPL16791_DE_RMS_vs_Normal_FDR005.csv",
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
  "GPL16791_DE_RMS_vs_Normal_FDR005_logFC1.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Guardar objeto del test
# ------------------------------------------------------------

saveRDS(
  qlf,
  "GSE108022_GPL16791_QLFTest_RMS_vs_Normal.rds"
)

cat("\nArchivos generados:\n")
cat("- GPL16791_DE_RMS_vs_Normal_todos_los_genes.csv\n")
cat("- GPL16791_DE_RMS_vs_Normal_FDR005.csv\n")
cat("- GPL16791_DE_RMS_vs_Normal_FDR005_logFC1.csv\n")
cat("- GSE108022_GPL16791_QLFTest_RMS_vs_Normal.rds\n")