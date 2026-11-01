# ============================================================
# 34_DE_FOXO1FP_vs_FN.R
# Expresión diferencial directa:
#   FOXO1-FP - FN
# GPL16791 - GSE108022
# ============================================================

library(edgeR)
library(limma)

# ------------------------------------------------------------
# 1. Cargar ajuste QL validado
# ------------------------------------------------------------

fit <- readRDS(
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_QLfit.rds"
)

cat("Ajuste QL cargado.\n")

cat("\nColumnas del diseño:\n")
print(colnames(fit$design))

if (!identical(
  colnames(fit$design),
  c("Normal", "FN", "FOXO1_FP")
)) {
  stop("ERROR: columnas inesperadas en la matriz de diseño.")
}

if (is.null(fit$genes) ||
    !"gene" %in% colnames(fit$genes)) {
  stop("ERROR: no se dispone de identidad génica válida.")
}


# ------------------------------------------------------------
# 2. Definir contraste directo
#
# logFC > 0 = mayor expresión en FOXO1-FP
# logFC < 0 = mayor expresión en FN
# ------------------------------------------------------------

contraste <- makeContrasts(
  FOXO1FP_vs_FN = FOXO1_FP - FN,
  levels = fit$design
)

cat("\nContraste:\n")
print(contraste)


# ------------------------------------------------------------
# 3. Test quasi-likelihood
# ------------------------------------------------------------

qlf <- glmQLFTest(
  fit,
  contrast = contraste[, "FOXO1FP_vs_FN"]
)

if (nrow(qlf$table) != nrow(qlf$genes)) {
  stop("ERROR: qlf$table y qlf$genes no coinciden.")
}


# ------------------------------------------------------------
# 4. Reconstruir tabla preservando identidad génica
# ------------------------------------------------------------

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

resultado_base$FDR <- p.adjust(
  resultado_base$PValue,
  method = "BH"
)

if (anyNA(resultado_base$logFC) ||
    anyNA(resultado_base$PValue) ||
    anyNA(resultado_base$FDR)) {
  stop("ERROR: existen NA en los resultados.")
}

resultado <- resultado_base[
  order(resultado_base$PValue),
  ,
  drop = FALSE
]

rownames(resultado) <- NULL


# ------------------------------------------------------------
# 5. Verificación independiente con topTags
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
  "\nVERIFICACIÓN: identidad génica y estadísticas verificadas.\n"
)


# ------------------------------------------------------------
# 6. Resumen
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

cat(
  "\nFDR < 0.05 y logFC >= 1 (mayor en FOXO1-FP):\n"
)
print(
  sum(
    resultado$FDR < 0.05 &
    resultado$logFC >= 1
  )
)

cat(
  "\nFDR < 0.05 y logFC <= -1 (mayor en FN):\n"
)
print(
  sum(
    resultado$FDR < 0.05 &
    resultado$logFC <= -1
  )
)


# ------------------------------------------------------------
# 7. Mostrar 20 genes con menor PValue
# ------------------------------------------------------------

cat("\n20 genes con menor PValue:\n")

print(
  resultado[
    seq_len(min(20, nrow(resultado))),
    c(
      "gene",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR"
    )
  ]
)


# ------------------------------------------------------------
# 8. Guardar resultados
# ------------------------------------------------------------

write.csv(
  resultado,
  "GPL16791_DE_FOXO1FP_vs_FN_todos.csv",
  row.names = FALSE
)

write.csv(
  resultado[
    resultado$FDR < 0.05,
    ,
    drop = FALSE
  ],
  "GPL16791_DE_FOXO1FP_vs_FN_FDR005.csv",
  row.names = FALSE
)

write.csv(
  resultado[
    resultado$FDR < 0.05 &
    abs(resultado$logFC) >= 1,
    ,
    drop = FALSE
  ],
  "GPL16791_DE_FOXO1FP_vs_FN_FDR005_logFC1.csv",
  row.names = FALSE
)

saveRDS(
  qlf,
  "GSE108022_GPL16791_FOXO1FP_vs_FN_QLFTest.rds"
)


# ------------------------------------------------------------
# 9. Confirmación final
# ------------------------------------------------------------

cat("\nArchivos generados:\n")
cat("- GPL16791_DE_FOXO1FP_vs_FN_todos.csv\n")
cat("- GPL16791_DE_FOXO1FP_vs_FN_FDR005.csv\n")
cat("- GPL16791_DE_FOXO1FP_vs_FN_FDR005_logFC1.csv\n")
cat("- GSE108022_GPL16791_FOXO1FP_vs_FN_QLFTest.rds\n")

cat(
  "\nRESULTADO FINAL: contraste directo FOXO1-FP vs FN completado y verificado.\n"
)