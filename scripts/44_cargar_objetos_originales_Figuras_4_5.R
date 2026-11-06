# Load retained original R objects supporting manuscript Figures 4 and 5.
# This script does not rerun statistical analyses or regenerate figures.
# Run from the repository root (or adjust the relative path).

archivo <- file.path("results", "RMS_Figuras_4_5_objetos_originales.RData")
if (!file.exists(archivo)) {
  stop("Archivo no encontrado: ", archivo)
}

objetos_originales <- new.env(parent = baseenv())
nombres <- load(archivo, envir = objetos_originales)

cat("Objetos originales cargados desde:", archivo, "\n")
print(nombres)

# To inspect an original object without recalculating anything, e.g.:
# print(objetos_originales$tabla_cor_CKM)
# print(summary(objetos_originales$modelo_final))
# print(objetos_originales$resultado_PC1_miogenesis)

# Original numerical output is also preserved in:
# results/resultados_originales_Figuras_4_5.txt
