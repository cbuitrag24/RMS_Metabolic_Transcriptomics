# ============================================================
# 25_leading_edge_Hallmarks_replicados.R
#
# OBJETIVO:
# Comparar los leading edges de los Hallmarks que cumplieron
# el criterio predefinido de replicación funcional entre:
#
#   GSE108022
#   GSE28511
#
# Ambos enriquecimientos fueron realizados sobre exactamente
# el mismo universo de 9.428 genes.
#
# Para cada Hallmark replicado se calculará:
# - genes del leading edge en GSE108022
# - genes del leading edge en GSE28511
# - intersección
# - unión
# - proporción compartida respecto de cada leading edge
# - índice de Jaccard
#
# IMPORTANTE:
# - No se redefine qué Hallmarks están replicados.
# - No se seleccionan genes por PValue/FDR.
# - No se modifica ningún resultado de fgsea.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_108022 <-
  "Hallmark_fgsea_GSE108022_universo_comun_9428.rds"

archivo_28511 <-
  "Hallmark_fgsea_GSE28511_universo_comun_9428.rds"

archivo_replicados <-
  "Hallmark_replicados_GSE108022_GSE28511.csv"


archivos <- c(
  archivo_108022,
  archivo_28511,
  archivo_replicados
)


faltantes <- archivos[
  !file.exists(archivos)
]


if (length(faltantes) > 0) {

  cat("\nArchivos faltantes:\n")
  print(faltantes)

  stop(
    "ERROR: faltan archivos necesarios."
  )
}


# ------------------------------------------------------------
# 2. Cargar
# ------------------------------------------------------------

gse108022 <- readRDS(
  archivo_108022
)

gse28511 <- readRDS(
  archivo_28511
)

replicados <- read.csv(
  archivo_replicados,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\nDimensiones GSE108022:\n")
print(
  dim(gse108022)
)

cat("\nDimensiones GSE28511:\n")
print(
  dim(gse28511)
)

cat("\nDimensiones tabla de Hallmarks replicados:\n")
print(
  dim(replicados)
)


# ------------------------------------------------------------
# 3. Verificar estructuras
# ------------------------------------------------------------

necesarias_fg <- c(
  "pathway",
  "NES",
  "padj",
  "leadingEdge"
)


if (!all(
  necesarias_fg %in% colnames(gse108022)
)) {
  stop(
    "ERROR: faltan columnas necesarias en GSE108022."
  )
}


if (!all(
  necesarias_fg %in% colnames(gse28511)
)) {
  stop(
    "ERROR: faltan columnas necesarias en GSE28511."
  )
}


if (!("pathway" %in% colnames(replicados))) {
  stop(
    "ERROR: la tabla de replicados no contiene pathway."
  )
}


if (anyDuplicated(gse108022$pathway)) {
  stop(
    "ERROR: pathways duplicados en GSE108022."
  )
}


if (anyDuplicated(gse28511$pathway)) {
  stop(
    "ERROR: pathways duplicados en GSE28511."
  )
}


if (anyDuplicated(replicados$pathway)) {
  stop(
    "ERROR: pathways duplicados en tabla de replicados."
  )
}


# ------------------------------------------------------------
# 4. Verificar los 13 Hallmarks replicados
# ------------------------------------------------------------

cat("\nNúmero de Hallmarks replicados:\n")
cat(
  nrow(replicados),
  "\n"
)


if (nrow(replicados) != 13) {
  stop(
    paste(
      "ERROR: se esperaban 13 Hallmarks replicados y se encontraron",
      nrow(replicados)
    )
  )
}


faltan_108022 <- setdiff(
  replicados$pathway,
  gse108022$pathway
)

faltan_28511 <- setdiff(
  replicados$pathway,
  gse28511$pathway
)


cat("\nReplicados ausentes de GSE108022:\n")
cat(
  length(faltan_108022),
  "\n"
)

cat("\nReplicados ausentes de GSE28511:\n")
cat(
  length(faltan_28511),
  "\n"
)


if (
  length(faltan_108022) > 0 ||
  length(faltan_28511) > 0
) {
  stop(
    "ERROR: algún Hallmark replicado falta en los resultados fgsea."
  )
}


# ------------------------------------------------------------
# 5. Verificar que leadingEdge sea realmente una lista
# ------------------------------------------------------------

if (!is.list(gse108022$leadingEdge)) {
  stop(
    "ERROR: leadingEdge de GSE108022 no es una lista."
  )
}

if (!is.list(gse28511$leadingEdge)) {
  stop(
    "ERROR: leadingEdge de GSE28511 no es una lista."
  )
}


# ------------------------------------------------------------
# 6. Preparar almacenamiento
# ------------------------------------------------------------

resumen_lista <- vector(
  "list",
  nrow(replicados)
)

genes_lista <- vector(
  "list",
  nrow(replicados)
)


# ------------------------------------------------------------
# 7. Comparar cada Hallmark replicado
# ------------------------------------------------------------

for (i in seq_len(nrow(replicados))) {

  pathway_actual <- replicados$pathway[i]


  idx_108022 <- match(
    pathway_actual,
    gse108022$pathway
  )

  idx_28511 <- match(
    pathway_actual,
    gse28511$pathway
  )


  if (
    is.na(idx_108022) ||
    is.na(idx_28511)
  ) {
    stop(
      paste(
        "ERROR: no pudo alinearse",
        pathway_actual
      )
    )
  }


  genes_108022 <- unique(
    as.character(
      gse108022$leadingEdge[[idx_108022]]
    )
  )

  genes_28511 <- unique(
    as.character(
      gse28511$leadingEdge[[idx_28511]]
    )
  )


  genes_108022 <- genes_108022[
    !is.na(genes_108022) &
    genes_108022 != ""
  ]

  genes_28511 <- genes_28511[
    !is.na(genes_28511) &
    genes_28511 != ""
  ]


  compartidos <- intersect(
    genes_108022,
    genes_28511
  )

  union_genes <- union(
    genes_108022,
    genes_28511
  )


  n_108022 <- length(
    genes_108022
  )

  n_28511 <- length(
    genes_28511
  )

  n_compartidos <- length(
    compartidos
  )

  n_union <- length(
    union_genes
  )


  prop_108022_compartida <- (
    n_compartidos /
    n_108022
  )

  prop_28511_compartida <- (
    n_compartidos /
    n_28511
  )

  jaccard <- (
    n_compartidos /
    n_union
  )


  # ----------------------------------------------------------
  # Resumen por Hallmark
  # ----------------------------------------------------------

  resumen_lista[[i]] <- data.frame(

    pathway =
      pathway_actual,

    NES_GSE108022 =
      gse108022$NES[idx_108022],

    FDR_GSE108022 =
      gse108022$padj[idx_108022],

    NES_GSE28511 =
      gse28511$NES[idx_28511],

    FDR_GSE28511 =
      gse28511$padj[idx_28511],

    n_leadingEdge_GSE108022 =
      n_108022,

    n_leadingEdge_GSE28511 =
      n_28511,

    n_compartidos =
      n_compartidos,

    n_union =
      n_union,

    proporcion_GSE108022_compartida =
      prop_108022_compartida,

    proporcion_GSE28511_compartida =
      prop_28511_compartida,

    jaccard =
      jaccard,

    stringsAsFactors = FALSE
  )


  # ----------------------------------------------------------
  # Tabla larga de genes
  # ----------------------------------------------------------

  genes_lista[[i]] <- data.frame(

    pathway =
      pathway_actual,

    gene =
      union_genes,

    leadingEdge_GSE108022 =
      union_genes %in% genes_108022,

    leadingEdge_GSE28511 =
      union_genes %in% genes_28511,

    compartido =
      union_genes %in% compartidos,

    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 8. Combinar resultados
# ------------------------------------------------------------

resumen <- do.call(
  rbind,
  resumen_lista
)

genes_detalle <- do.call(
  rbind,
  genes_lista
)


rownames(resumen) <- NULL
rownames(genes_detalle) <- NULL


# ------------------------------------------------------------
# 9. Verificaciones
# ------------------------------------------------------------

if (nrow(resumen) != 13) {
  stop(
    "ERROR: el resumen no contiene 13 Hallmarks."
  )
}


if (anyDuplicated(resumen$pathway)) {
  stop(
    "ERROR: pathways duplicados en el resumen."
  )
}


if (
  any(resumen$n_leadingEdge_GSE108022 <= 0) ||
  any(resumen$n_leadingEdge_GSE28511 <= 0)
) {
  stop(
    "ERROR: algún leading edge está vacío."
  )
}


if (
  any(resumen$n_compartidos >
      resumen$n_leadingEdge_GSE108022) ||
  any(resumen$n_compartidos >
      resumen$n_leadingEdge_GSE28511)
) {
  stop(
    "ERROR: intersección imposible detectada."
  )
}


if (
  any(resumen$jaccard < 0) ||
  any(resumen$jaccard > 1)
) {
  stop(
    "ERROR: Jaccard fuera del rango [0,1]."
  )
}


# ------------------------------------------------------------
# 10. Resumen global
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("CONCORDANCIA DE LEADING EDGES\n")
cat("====================================================\n")


cat("\nHallmarks evaluados:\n")
cat(
  nrow(resumen),
  "\n"
)


cat("\nResumen del número de genes compartidos:\n")
print(
  summary(
    resumen$n_compartidos
  )
)


cat("\nResumen de proporción del leading edge de GSE108022")
cat(" compartida con GSE28511:\n")
print(
  summary(
    resumen$proporcion_GSE108022_compartida
  )
)


cat("\nResumen de proporción del leading edge de GSE28511")
cat(" compartida con GSE108022:\n")
print(
  summary(
    resumen$proporcion_GSE28511_compartida
  )
)


cat("\nResumen del índice de Jaccard:\n")
print(
  summary(
    resumen$jaccard
  )
)


# ------------------------------------------------------------
# 11. Mostrar los 13 Hallmarks
# ------------------------------------------------------------

resumen_ordenado <- resumen[
  order(
    -resumen$jaccard,
    -resumen$n_compartidos
  ),
  ,
  drop = FALSE
]

rownames(resumen_ordenado) <- NULL


cat("\n====================================================\n")
cat("LEADING EDGE POR HALLMARK REPLICADO\n")
cat("====================================================\n\n")


print(
  resumen_ordenado[
    ,
    c(
      "pathway",
      "n_leadingEdge_GSE108022",
      "n_leadingEdge_GSE28511",
      "n_compartidos",
      "proporcion_GSE108022_compartida",
      "proporcion_GSE28511_compartida",
      "jaccard"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Genes compartidos
# ------------------------------------------------------------

genes_compartidos <- genes_detalle[
  genes_detalle$compartido,
  ,
  drop = FALSE
]


cat("\nNúmero total de asociaciones Hallmark-gen")
cat(" compartidas en leading edge:\n")
cat(
  nrow(genes_compartidos),
  "\n"
)


cat("\nGenes únicos presentes en al menos un leading edge")
cat(" compartido:\n")
cat(
  length(
    unique(
      genes_compartidos$gene
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 13. Frecuencia de genes en leading edges compartidos
#
# Esto NO constituye una prueba estadística adicional.
# Sólo cuantifica en cuántos Hallmarks replicados aparece
# cada gen dentro del leading edge de ambos estudios.
# ------------------------------------------------------------

frecuencia_genes <- as.data.frame(
  table(
    genes_compartidos$gene
  ),
  stringsAsFactors = FALSE
)

colnames(frecuencia_genes) <- c(
  "gene",
  "n_Hallmarks_leadingEdge_compartido"
)


frecuencia_genes <- frecuencia_genes[
  order(
    -frecuencia_genes$n_Hallmarks_leadingEdge_compartido,
    frecuencia_genes$gene
  ),
  ,
  drop = FALSE
]

rownames(frecuencia_genes) <- NULL


cat("\nGenes con mayor frecuencia en leading edges compartidos:\n\n")

print(
  head(
    frecuencia_genes,
    30
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Guardar
# ------------------------------------------------------------

write.csv(
  resumen_ordenado,
  "Hallmark_leadingEdge_concordancia_resumen.csv",
  row.names = FALSE
)


write.csv(
  genes_detalle,
  "Hallmark_leadingEdge_genes_detalle.csv",
  row.names = FALSE
)


write.csv(
  genes_compartidos,
  "Hallmark_leadingEdge_genes_compartidos.csv",
  row.names = FALSE
)


write.csv(
  frecuencia_genes,
  "Hallmark_leadingEdge_genes_frecuencia.csv",
  row.names = FALSE
)


saveRDS(
  list(
    resumen = resumen_ordenado,
    genes_detalle = genes_detalle,
    genes_compartidos = genes_compartidos,
    frecuencia_genes = frecuencia_genes
  ),
  "Hallmark_leadingEdge_concordancia.rds"
)


cat("\nArchivos generados:\n")

cat(
  "- Hallmark_leadingEdge_concordancia_resumen.csv\n"
)

cat(
  "- Hallmark_leadingEdge_genes_detalle.csv\n"
)

cat(
  "- Hallmark_leadingEdge_genes_compartidos.csv\n"
)

cat(
  "- Hallmark_leadingEdge_genes_frecuencia.csv\n"
)

cat(
  "- Hallmark_leadingEdge_concordancia.rds\n"
)


cat("\n====================================================\n")
cat("ANÁLISIS DE LEADING EDGES COMPLETADO\n")
cat("====================================================\n")