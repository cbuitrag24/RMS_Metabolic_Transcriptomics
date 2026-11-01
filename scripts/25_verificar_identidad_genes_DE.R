# ============================================================
# 25_verificar_identidad_genes_DE.R
# Verificación de identidad de genes en los objetos edgeR
# GSE108022 - GPL16791
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar objetos
# ------------------------------------------------------------

dge <- readRDS("GSE108022_GPL16791_dispersion.rds")
fit <- readRDS("GSE108022_GPL16791_QLfit.rds")
qlf <- readRDS("GSE108022_GPL16791_QLFTest_RMS_vs_Normal.rds")


# ------------------------------------------------------------
# 2. Inspeccionar dge$genes
# ------------------------------------------------------------

cat("\n========== dge$genes ==========\n")

cat("Clase:\n")
print(class(dge$genes))

cat("\nDimensiones:\n")
print(dim(dge$genes))

cat("\nNombres de columnas:\n")
print(colnames(dge$genes))

cat("\nPrimeras 10 filas:\n")
print(head(dge$genes, 10))


# ------------------------------------------------------------
# 3. Inspeccionar fit$genes
# ------------------------------------------------------------

cat("\n========== fit$genes ==========\n")

cat("Clase:\n")
print(class(fit$genes))

cat("\nDimensiones:\n")
print(dim(fit$genes))

cat("\nNombres de columnas:\n")
print(colnames(fit$genes))

cat("\nPrimeras 10 filas:\n")
print(head(fit$genes, 10))


# ------------------------------------------------------------
# 4. Inspeccionar qlf$genes
# ------------------------------------------------------------

cat("\n========== qlf$genes ==========\n")

cat("Clase:\n")
print(class(qlf$genes))

cat("\nDimensiones:\n")
print(dim(qlf$genes))

cat("\nNombres de columnas:\n")
print(colnames(qlf$genes))

cat("\nPrimeras 10 filas:\n")
print(head(qlf$genes, 10))


# ------------------------------------------------------------
# 5. Inspeccionar nombres de fila
# ------------------------------------------------------------

cat("\n========== ROWNAMES ==========\n")

cat("\nPrimeros 10 rownames de counts:\n")
print(head(rownames(dge$counts), 10))

cat("\nPrimeros 10 rownames de coefficients:\n")
print(head(rownames(fit$coefficients), 10))

cat("\nPrimeros 10 rownames de qlf$table:\n")
print(head(rownames(qlf$table), 10))


# ------------------------------------------------------------
# 6. Comparar identidad y orden
# ------------------------------------------------------------

cat("\n========== COMPARACIONES ==========\n")

cat("\n¿rownames(counts) idénticos a rownames(coefficients)?\n")
print(
  identical(
    rownames(dge$counts),
    rownames(fit$coefficients)
  )
)

cat("\n¿rownames(coefficients) idénticos a rownames(qlf$table)?\n")
print(
  identical(
    rownames(fit$coefficients),
    rownames(qlf$table)
  )
)

cat("\nInspección terminada. No se modificó ningún archivo.\n")