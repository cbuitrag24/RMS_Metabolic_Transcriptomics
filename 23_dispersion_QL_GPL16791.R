# ============================================================
# 23_dispersion_QL_GPL16791.R
# Estimación de dispersiones y ajuste QL robusto
# Análisis primario GPL16791
# GSE108022
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar objeto y diseño
# ------------------------------------------------------------

dge <- readRDS("GSE108022_GPL16791_filtrado_TMM.rds")
design <- readRDS("GSE108022_GPL16791_design_matrix.rds")

cat("Objeto GPL16791 cargado:\n")
cat("Genes:", nrow(dge), "\n")
cat("Muestras:", ncol(dge), "\n")

cat("\nColumnas del diseño:\n")
print(colnames(design))

cat("\nDimensiones del diseño:\n")
print(dim(design))

cat("\nRango del diseño:", qr(design)$rank, "\n")

if (nrow(design) != ncol(dge)) {
  stop("ERROR: diseño y objeto DGE tienen distinto número de muestras.")
}

if (qr(design)$rank != ncol(design)) {
  stop("ERROR: la matriz de diseño no es de rango completo.")
}


# ------------------------------------------------------------
# 2. Estimar dispersiones
# ------------------------------------------------------------

dge <- estimateDisp(
  dge,
  design = design,
  robust = TRUE
)

cat("\nDispersión común:\n")
print(dge$common.dispersion)

cat("\nBCV común:\n")
print(sqrt(dge$common.dispersion))

cat("\nResumen de dispersión tagwise:\n")
print(summary(dge$tagwise.dispersion))

cat("\nResumen de BCV tagwise:\n")
print(summary(sqrt(dge$tagwise.dispersion)))

if (anyNA(dge$tagwise.dispersion)) {
  stop("ERROR: existen dispersiones tagwise NA.")
}

if (any(dge$tagwise.dispersion < 0)) {
  stop("ERROR: existen dispersiones negativas.")
}


# ------------------------------------------------------------
# 3. Ajuste quasi-likelihood robusto
# ------------------------------------------------------------

fit <- glmQLFit(
  dge,
  design = design,
  robust = TRUE
)

cat("\nRESULTADO: ajuste quasi-likelihood completado.\n")

cat("\nGrados de libertad residuales:\n")
print(summary(fit$df.residual))

cat("\nGrados de libertad residuales ajustados:\n")
print(summary(fit$df.residual.adj))

cat("\nResumen de varianza QL posterior (s2.post):\n")
print(summary(fit$s2.post))

cat("\nResumen de df.prior:\n")
print(summary(fit$df.prior))

cat("\nDispersión QL promedio:\n")
print(fit$average.ql.dispersion)


# ------------------------------------------------------------
# 4. Comprobaciones de integridad
# ------------------------------------------------------------

if (anyNA(fit$coefficients)) {
  stop("ERROR: existen coeficientes NA.")
}

if (anyNA(fit$s2.post)) {
  stop("ERROR: existen valores NA en s2.post.")
}

if (any(fit$s2.post <= 0)) {
  stop("ERROR: existen valores no positivos en s2.post.")
}

cat("\nRESULTADO: ajuste QL válido, sin NA ni s2.post no positivos.\n")


# ------------------------------------------------------------
# 5. Diagnóstico de genes con mayor varianza QL
#    NO se elimina ningún gen
# ------------------------------------------------------------

ord <- order(
  fit$s2.post,
  decreasing = TRUE
)

n_show <- min(20, length(ord))

diagnostico <- data.frame(
  gene = dge$genes$gene[ord[seq_len(n_show)]],
  dispersion = dge$tagwise.dispersion[ord[seq_len(n_show)]],
  QL_s2_post = fit$s2.post[ord[seq_len(n_show)]],
  df_prior = fit$df.prior[ord[seq_len(n_show)]],
  stringsAsFactors = FALSE
)

cat("\n20 genes con mayor varianza QL posterior:\n")
print(diagnostico)


# ------------------------------------------------------------
# 6. Guardar objetos
# ------------------------------------------------------------

saveRDS(
  dge,
  "GSE108022_GPL16791_dispersion.rds"
)

saveRDS(
  fit,
  "GSE108022_GPL16791_QLfit.rds"
)

write.csv(
  diagnostico,
  "GPL16791_QL_diagnostico_genes_mayor_varianza.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("- GSE108022_GPL16791_dispersion.rds\n")
cat("- GSE108022_GPL16791_QLfit.rds\n")
cat("- GPL16791_QL_diagnostico_genes_mayor_varianza.csv\n")