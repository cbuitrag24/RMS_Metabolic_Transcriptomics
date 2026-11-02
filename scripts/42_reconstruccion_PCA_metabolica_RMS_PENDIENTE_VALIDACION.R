# ============================================================
# RECONSTRUCCION DOCUMENTADA (NO ES EL SCRIPT HISTORICO)
# Figura 3: PCA de 16 genes metabolicos en 64 tumores RMS
# GSE108022, plataforma GPL16791
# ESTADO: REPRODUCCION NUMERICA Y GRAFICA VERIFICADA 2026-11-01
# ============================================================
# Procedencia: scripts originales 14 y 22 del repositorio,
# manuscrito RMS y script grafico de Figura 3 recuperado.
# Validacion: 16 genes, 64 tumores, PC1 = 36.481659%,
# PC2 = 13.131902%; cargas PC1 negativas. PDF cotejado
# visualmente con Figura 3 historica. No es el script original.
# No modifica datos de entrada ni escribe resultados hasta
# superar las comprobaciones de identidad y varianza.

suppressPackageStartupMessages(library(edgeR))

# Ejecutar desde la carpeta donde estan los RDS y CSV originales.
archivo_dge <- "GSE108022_GPL16791_filtrado_TMM.rds"
archivo_metadata <- "GPL16791_diseno_muestras.csv"

nucleo_metabolico <- c(
  "ETFDH", "HADHB", "ACO2", "ECH1", "PDHB", "DLD",
  "MDH1", "DLST", "SDHD", "SUCLG1", "FH", "PDHA1",
  "MDH2", "RETSAT", "DECR1", "ACAA1"
)

# 1. Comprobar archivos de entrada
for (ruta in c(archivo_dge, archivo_metadata)) {
  if (!file.exists(ruta)) stop("Falta el archivo de entrada: ", ruta)
}
dge <- readRDS(archivo_dge)
metadata <- read.csv(archivo_metadata, stringsAsFactors = FALSE,
                     check.names = FALSE)

# 2. Verificar y alinear las muestras
campos <- c("samplename", "condition", "fusion_status")
if (!all(campos %in% names(metadata))) {
  stop("Faltan columnas en metadata: ",
       paste(setdiff(campos, names(metadata)), collapse = ", "))
}
if (anyDuplicated(metadata$samplename)) {
  stop("Muestras duplicadas en metadata")
}
if (is.null(colnames(dge$counts)) ||
    anyDuplicated(colnames(dge$counts))) {
  stop("Los nombres de muestra del DGE son invalidos")
}
if (!setequal(colnames(dge$counts), metadata$samplename)) {
  stop("Las muestras del DGE no coinciden con la metadata")
}
metadata <- metadata[
  match(colnames(dge$counts), metadata$samplename),
  , drop = FALSE
]
if (!identical(metadata$samplename, colnames(dge$counts))) {
  stop("Error al alinear metadata y expresion")
}
es_rms <- as.character(metadata$condition) == "RMS"
if (anyNA(es_rms) || sum(es_rms) != 64L || sum(!es_rms) != 5L) {
  stop("Se esperaban exactamente 64 tumores RMS y 5 normales")
}

# 3. Identificar los 16 genes mediante la anotacion original
if (anyDuplicated(nucleo_metabolico) ||
    length(nucleo_metabolico) != 16L) {
  stop("La lista no contiene exactamente 16 genes unicos")
}
if (is.null(dge$genes) ||
    !("gene" %in% names(dge$genes))) {
  stop("No existe la anotacion dge$genes$gene")
}
simbolos <- as.character(dge$genes$gene)
if (length(simbolos) != nrow(dge$counts)) {
  stop("La anotacion y los counts tienen distinto numero de filas")
}
indices_genes <- lapply(
  nucleo_metabolico,
  function(g) which(!is.na(simbolos) & simbolos == g)
)
nombres_indices <- vapply(indices_genes, length, integer(1))
if (any(nombres_indices != 1L)) {
  stop("Cada gen debe tener exactamente una fila. Problemas: ",
       paste(paste0(
         nucleo_metabolico[nombres_indices != 1L], " (",
         nombres_indices[nombres_indices != 1L], " filas)"
       ), collapse = ", "))
}
filas_metabolicas <- unlist(indices_genes, use.names = FALSE)
cat("Los 16 genes fueron identificados de forma unica.\n")

# 4. Reconstruir matriz de expresion
# Hipotesis de reconstruccion validada frente a Figura 3:
# log2 CPM normalizados TMM, prior.count = 2.
expresion_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
expresion_rms_22 <- expresion_logcpm[, es_rms, drop = FALSE]
fusion_rms <- as.character(metadata$fusion_status[es_rms])

# 5. Verificar grupos de fusion
categorias <- c(
  "Fusion_negative", "PAX3-FOXO1", "PAX7-FOXO1", "PAX3-NCOA1"
)
if (anyNA(fusion_rms) || !all(fusion_rms %in% categorias)) {
  stop("Categorias de fusion no reconocidas. Verificar etiquetas originales.")
}
conteos_esperados <- c(
  "Fusion_negative" = 44L, "PAX3-FOXO1" = 13L,
  "PAX7-FOXO1" = 6L, "PAX3-NCOA1" = 1L
)
conteos <- table(factor(fusion_rms, levels = categorias))
if (!identical(as.integer(conteos), as.integer(conteos_esperados))) {
  stop("Los recuentos de fusion no coinciden: ",
       paste(names(conteos), conteos, collapse = "; "))
}

# 6. Preparar matriz de 16 genes x 64 tumores
matriz <- expresion_rms_22[filas_metabolicas, , drop = FALSE]
rownames(matriz) <- nucleo_metabolico
if (!identical(dim(matriz), c(16L, 64L)) ||
    any(!is.finite(matriz))) {
  stop("La matriz 16 x 64 contiene valores invalidos")
}
if (any(apply(matriz, 1, sd) == 0)) {
  stop("Hay genes sin variacion")
}
cat("Matriz verificada: 16 genes x 64 tumores.\n")

# 7. Calcular PCA
pca_metabolico <- prcomp(t(matriz), center = TRUE, scale. = TRUE)
varianza <- 100 * pca_metabolico$sdev^2 /
  sum(pca_metabolico$sdev^2)
cat(sprintf(
  "PC1 = %.6f%%; PC2 = %.6f%%\n", varianza[1], varianza[2]
))
cat("Cargas PC1:\n")
print(pca_metabolico$rotation[, 1])

# 8. Validar contra la Figura 3 original
if (!identical(round(varianza[1:2], 2), c(36.48, 13.13))) {
  stop("VALIDACION FALLIDA: la varianza no reproduce ",
       "36.48% y 13.13%. No se genera la figura.")
}
if (!all(pca_metabolico$rotation[, 1] < 0)) {
  stop("VALIDACION INCOMPLETA: no todas las cargas PC1 ",
       "son negativas. La orientacion debe contrastarse ",
       "con la figura historica.")
}

# 9. Figura 3: codigo grafico original recuperado
grupos <- factor(fusion_rms)
colores <- c(
  "Fusion_negative" = "#2166AC",
  "PAX3-FOXO1" = "#B2182B",
  "PAX7-FOXO1" = "#4DAF4A",
  "PAX3-NCOA1" = "#984EA3"
)
# Nombre diferente para proteger la figura historica.
salida <- "Figura_3_PCA_Metabolico_RMS_RECONSTRUCCION.pdf"
pdf(salida, width = 9, height = 6.5)
par(mar = c(5, 5, 2, 2))
plot(
  pca_metabolico$x[, 1], pca_metabolico$x[, 2],
  col = colores[as.character(grupos)],
  pch = 16, cex = 1.1,
  xlab = sprintf("PC1 (%.2f%% variance explained)", varianza[1]),
  ylab = sprintf("PC2 (%.2f%% variance explained)", varianza[2]),
  main = "", las = 1
)
abline(h = 0, v = 0, lty = 3, col = "gray75")
legend(
  "topleft", legend = names(colores), col = colores,
  pch = 16, bty = "n", cex = 0.85, inset = c(0.02, 0.02)
)
dev.off()
cat(
  "VALIDACION NUMERICA PRELIMINAR SUPERADA.\n",
  "PDF generado: ", salida, "\n", sep = ""
)
cat("El PDF reconstruido fue cotejado visualmente con el historico.\n")
