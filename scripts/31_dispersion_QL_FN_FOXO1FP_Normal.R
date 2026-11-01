# ============================================================
# 31_dispersion_QL_FN_FOXO1FP_Normal.R
# Estimación de dispersión y ajuste quasi-likelihood robusto
# Normal / FN / FOXO1-FP
# GPL16791 - GSE108022
#
# Este script NO realiza todavía los contrastes DE.
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar objetos preparados
# ------------------------------------------------------------

dge <- readRDS(
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_filtrado_TMM.rds"
)

design <- readRDS(
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_design.rds"
)

cat("Objetos cargados.\n")

cat("\nDimensiones del DGEList:\n")
print(dim(dge))

cat("\nDimensiones del diseño:\n")
print(dim(design))

cat("\nColumnas del diseño:\n")
print(colnames(design))


# ------------------------------------------------------------
# 2. Controles previos
# ------------------------------------------------------------

if (ncol(dge) != nrow(design)) {
  stop("ERROR: número de muestras incompatible entre DGEList y diseño.")
}

if (qr(design)$rank != ncol(design)) {
  stop("ERROR: la matriz de diseño no tiene rango completo.")
}

if (!identical(
  colnames(dge$counts),
  rownames(design)
)) {
  stop("ERROR: el orden de las muestras no coincide entre counts y diseño.")
}

cat("\nVERIFICACIÓN: DGEList y diseño correctamente alineados.\n")


# ------------------------------------------------------------
# 3. Estimar dispersiones
# ------------------------------------------------------------

dge <- estimateDisp(
  dge,
  design = design,
  robust = TRUE
)

cat("\nEstimación de dispersión completada.\n")

cat("\nDispersión común:\n")
print(dge$common.dispersion)

cat("\nBCV común:\n")
print(sqrt(dge$common.dispersion))

cat("\nResumen de dispersiones tagwise:\n")
print(summary(dge$tagwise.dispersion))

cat("\nResumen de BCV tagwise:\n")
print(summary(sqrt(dge$tagwise.dispersion)))

if (anyNA(dge$tagwise.dispersion) ||
    any(dge$tagwise.dispersion <= 0)) {
  stop("ERROR: dispersiones tagwise inválidas.")
}


# ------------------------------------------------------------
# 4. Ajuste quasi-likelihood robusto
# ------------------------------------------------------------

fit <- glmQLFit(
  dge,
  design = design,
  robust = TRUE
)

cat("\nAjuste quasi-likelihood completado.\n")


# ------------------------------------------------------------
# 5. Diagnóstico de grados de libertad
# ------------------------------------------------------------

cat("\nResumen de df.residual:\n")
print(summary(fit$df.residual))

if (!is.null(fit$df.residual.adj)) {

  cat("\nResumen de df.residual.adj:\n")
  print(summary(fit$df.residual.adj))

}


# ------------------------------------------------------------
# 6. Diagnóstico de dispersión QL posterior
# ------------------------------------------------------------

if (is.null(fit$s2.post)) {
  stop("ERROR: fit$s2.post no está disponible.")
}

cat("\nResumen de s2.post:\n")
print(summary(fit$s2.post))

if (anyNA(fit$s2.post) ||
    any(fit$s2.post <= 0)) {
  stop("ERROR: valores inválidos en s2.post.")
}


# ------------------------------------------------------------
# 7. Prior degrees of freedom
# ------------------------------------------------------------

if (!is.null(fit$df.prior)) {

  cat("\nResumen de df.prior:\n")
  print(summary(fit$df.prior))

}


# ------------------------------------------------------------
# 8. Dispersión QL promedio
# ------------------------------------------------------------

if (!is.null(fit$average.ql.dispersion)) {

  cat("\nDispersión QL promedio:\n")
  print(fit$average.ql.dispersion)

}


# ------------------------------------------------------------
# 9. Genes con mayor s2.post
#    Sólo diagnóstico; NO se excluyen.
# ------------------------------------------------------------

orden_var <- order(
  fit$s2.post,
  decreasing = TRUE
)

n_top <- min(20, length(orden_var))

diagnostico <- data.frame(
  gene = fit$genes$gene[orden_var[seq_len(n_top)]],
  s2.post = fit$s2.post[orden_var[seq_len(n_top)]],
  stringsAsFactors = FALSE
)

cat("\n20 genes con mayor s2.post:\n")
print(diagnostico, row.names = FALSE)


# ------------------------------------------------------------
# 10. Guardar objetos
# ------------------------------------------------------------

saveRDS(
  dge,
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_dispersion.rds"
)

saveRDS(
  fit,
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_QLfit.rds"
)

write.csv(
  diagnostico,
  "GPL16791_FN_FOXO1FP_Normal_QL_diagnostico_mayor_varianza.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("- GSE108022_GPL16791_FN_FOXO1FP_Normal_dispersion.rds\n")
cat("- GSE108022_GPL16791_FN_FOXO1FP_Normal_QLfit.rds\n")
cat("- GPL16791_FN_FOXO1FP_Normal_QL_diagnostico_mayor_varianza.csv\n")

cat(
  "\nRESULTADO: dispersiones y ajuste QL completados; todavía no se realizaron contrastes DE.\n"
)