# ============================================================
# 15_preparar_counts_edgeR.R
# Verificación de conteos y preparación del objeto edgeR
# GSE108022
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Leer tabla de diseño ya verificada
# ------------------------------------------------------------

design_data <- read.csv(
  "DE_diseno_muestras.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Muestras en tabla de diseño:", nrow(design_data), "\n")


# ------------------------------------------------------------
# 2. Leer matriz de conteos
# ------------------------------------------------------------

counts_raw <- read.delim(
  "GSE108022_counts_matrix.txt.gz",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensiones de la matriz original:\n")
print(dim(counts_raw))

cat("\nPrimera columna:\n")
print(colnames(counts_raw)[1])


# ------------------------------------------------------------
# 3. Separar genes y matriz numérica
# ------------------------------------------------------------

genes <- counts_raw[[1]]

counts <- counts_raw[, -1, drop = FALSE]

cat("\nNúmero de genes:", length(genes), "\n")
cat("Número de columnas de muestras:", ncol(counts), "\n")


# ------------------------------------------------------------
# 4. Comprobar identificadores de muestras
# ------------------------------------------------------------

samples_counts <- colnames(counts)
samples_design <- design_data$samplename

cat("\nMuestras de counts ausentes en diseño:\n")
print(setdiff(samples_counts, samples_design))

cat("\nMuestras del diseño ausentes en counts:\n")
print(setdiff(samples_design, samples_counts))

if (!setequal(samples_counts, samples_design)) {
  stop("ERROR: los nombres de muestra de counts y diseño no coinciden.")
}

cat("\nRESULTADO: los dos archivos contienen exactamente las mismas muestras.\n")


# ------------------------------------------------------------
# 5. Reordenar counts según la tabla de diseño
# ------------------------------------------------------------

counts <- counts[, samples_design, drop = FALSE]

if (!identical(colnames(counts), design_data$samplename)) {
  stop("ERROR: el orden final de las muestras no coincide.")
}

cat("\nRESULTADO: orden de muestras verificado.\n")


# ------------------------------------------------------------
# 6. Comprobar integridad de los conteos
# ------------------------------------------------------------

counts_matrix <- as.matrix(counts)

if (anyNA(counts_matrix)) {
  stop("ERROR: existen valores NA en la matriz de conteos.")
}

if (any(counts_matrix < 0)) {
  stop("ERROR: existen conteos negativos.")
}

if (any(counts_matrix != floor(counts_matrix))) {
  stop("ERROR: existen valores no enteros en la matriz de conteos.")
}

cat("\nRESULTADO: conteos enteros, no negativos y sin NA.\n")


# ------------------------------------------------------------
# 7. Crear objeto DGEList
#    Todavía NO filtramos genes ni normalizamos
# ------------------------------------------------------------

dge <- DGEList(
  counts = counts_matrix,
  genes = data.frame(
    gene = genes,
    stringsAsFactors = FALSE
  )
)

cat("\nObjeto DGEList creado:\n")
cat("Genes:", nrow(dge), "\n")
cat("Muestras:", ncol(dge), "\n")

cat("\nResumen de tamaños de biblioteca calculados por edgeR:\n")
print(summary(dge$samples$lib.size))


# ------------------------------------------------------------
# 8. Guardar objeto para los próximos pasos
# ------------------------------------------------------------

saveRDS(
  dge,
  "GSE108022_edgeR_sin_filtrar.rds"
)

cat("\nArchivo generado: GSE108022_edgeR_sin_filtrar.rds\n")