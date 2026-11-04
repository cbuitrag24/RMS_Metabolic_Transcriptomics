# ============================================================
# 43_verificar_seleccion_16_genes_metabolicos.R
#
# Objetivo: reproducir la interseccion de leading edges
# compartidos entre GSE108022 y GSE28511 para dos Hallmarks.
# No se recalcula enriquecimiento ni se aplican nuevos umbrales.
#
# Entrada: Hallmark_leadingEdge_genes_compartidos.csv
# Generado por: 25_leading_edge_Hallmarks_replicados.R
# Ejecutar desde la raiz del repositorio (Rscript scripts/43_...).
# ============================================================

entrada <- file.path("data", "Hallmark_leadingEdge_genes_compartidos.csv")
if (!file.exists(entrada)) stop("No se encontro el CSV de entrada: ", entrada)

datos <- read.csv(entrada, stringsAsFactors = FALSE, check.names = FALSE)
columnas <- c("pathway", "gene", "leadingEdge_GSE108022",
              "leadingEdge_GSE28511", "compartido")
if (!all(columnas %in% names(datos))) stop("Columnas esperadas ausentes.")
if (anyNA(datos[, columnas]) ||
    any(!nzchar(trimws(datos$gene))) ||
    anyDuplicated(datos[c("pathway", "gene")])) {
  stop("Valores faltantes, genes vacios o asociaciones duplicadas.")
}
if (!all(datos$compartido & datos$leadingEdge_GSE108022 &
         datos$leadingEdge_GSE28511)) {
  stop("El archivo no contiene exclusivamente leading edges compartidos.")
}

oxphos <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
fatty <- "HALLMARK_FATTY_ACID_METABOLISM"
if (!all(c(oxphos, fatty) %in% datos$pathway)) {
  stop("Falta uno de los dos Hallmarks metabolicos.")
}

genes_oxphos <- unique(datos$gene[datos$pathway == oxphos])
genes_fatty <- unique(datos$gene[datos$pathway == fatty])
seleccion <- sort(intersect(genes_oxphos, genes_fatty))

# Lista independiente transcrita de la Tabla 2 del manuscrito.
tabla2 <- sort(c("ETFDH", "ACO2", "HADHB", "ECH1", "MDH1", "SDHD",
                 "DLD", "MDH2", "RETSAT", "PDHB", "FH", "PDHA1",
                 "DECR1", "SUCLG1", "DLST", "ACAA1"))

cat("Asociaciones totales:", nrow(datos), "\n")
cat("Genes OXPHOS:", length(genes_oxphos), "\n")
cat("Genes fatty acid metabolism:", length(genes_fatty), "\n")
cat("Interseccion:", length(seleccion), "\n")
cat("Genes seleccionados:", paste(seleccion, collapse = ", "), "\n")

if (!identical(seleccion, tabla2)) {
  cat("Faltantes respecto de Tabla 2:",
      paste(setdiff(tabla2, seleccion), collapse = ", "), "\n")
  cat("Adicionales respecto de Tabla 2:",
      paste(setdiff(seleccion, tabla2), collapse = ", "), "\n")
  stop("La interseccion no coincide con la Tabla 2.")
}

write.csv(data.frame(gene = seleccion),
          "seleccion_16_genes_metabolicos_verificada.csv",
          row.names = FALSE)
cat("VALIDADO: coincidencia exacta con los 16 genes de Tabla 2.\n")
