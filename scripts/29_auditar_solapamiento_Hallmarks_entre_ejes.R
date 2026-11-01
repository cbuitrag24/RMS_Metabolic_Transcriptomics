# ============================================================
# 29_auditar_solapamiento_Hallmarks_entre_ejes.R
#
# OBJETIVO:
# Determinar hasta qué punto la separación observada entre
# genes de leading edge RMS y músculo puede estar condicionada
# por la estructura de los conjuntos Hallmark originales.
#
# Se comparan:
#   A) genes Hallmark presentes en el universo común de 9428
#   B) genes observados en los leading edges compartidos
#
# No se realiza selección de genes ni nueva inferencia
# funcional.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_hallmark <-
  "MSigDB_Hallmark_Homo_sapiens.rds"

archivo_ranking <-
  "ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.rds"

archivo_replicados <-
  "Hallmark_replicados_GSE108022_GSE28511.csv"

archivo_mapa <-
  "Hallmark_leadingEdge_mapa_genes.rds"


archivos <- c(
  archivo_hallmark,
  archivo_ranking,
  archivo_replicados,
  archivo_mapa
)


faltantes <- archivos[
  !file.exists(archivos)
]


if (length(faltantes) > 0) {
  cat("\nArchivos faltantes:\n")
  print(faltantes)
  stop("ERROR: faltan archivos necesarios.")
}


# ------------------------------------------------------------
# 2. Cargar
# ------------------------------------------------------------

hallmark <- readRDS(
  archivo_hallmark
)

ranking <- readRDS(
  archivo_ranking
)

replicados <- read.csv(
  archivo_replicados,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapa <- readRDS(
  archivo_mapa
)


# ------------------------------------------------------------
# 3. Verificaciones
# ------------------------------------------------------------

if (!all(
  c("gene_set", "gene") %in% colnames(hallmark)
)) {
  stop("ERROR: estructura Hallmark inesperada.")
}


if (!("gene" %in% colnames(ranking))) {
  stop("ERROR: el ranking no contiene gene.")
}


if (!all(
  c(
    "pathway",
    "NES_GSE108022",
    "NES_GSE28511"
  ) %in% colnames(replicados)
)) {
  stop("ERROR: estructura de replicados inesperada.")
}


if (nrow(ranking) != 9428) {
  stop("ERROR: el universo común no contiene 9428 genes.")
}


if (nrow(replicados) != 13) {
  stop("ERROR: se esperaban 13 Hallmarks replicados.")
}


if (!("resumen_genes" %in% names(mapa))) {
  stop("ERROR: el objeto mapa no contiene resumen_genes.")
}


# ------------------------------------------------------------
# 4. Universo común
# ------------------------------------------------------------

universo <- unique(
  trimws(
    as.character(
      ranking$gene
    )
  )
)


if (length(universo) != 9428) {
  stop(
    "ERROR: el universo contiene genes duplicados o inválidos."
  )
}


# ------------------------------------------------------------
# 5. Dirección de los Hallmarks replicados
# ------------------------------------------------------------

signo1 <- sign(
  replicados$NES_GSE108022
)

signo2 <- sign(
  replicados$NES_GSE28511
)


if (any(signo1 != signo2)) {
  stop(
    "ERROR: existe un Hallmark replicado con signo discordante."
  )
}


replicados$eje <- ifelse(
  signo1 > 0,
  "RMS",
  "Musculo"
)


hallmarks_RMS <- replicados$pathway[
  replicados$eje == "RMS"
]

hallmarks_musculo <- replicados$pathway[
  replicados$eje == "Musculo"
]


cat("\nHallmarks RMS:\n")
cat(
  length(hallmarks_RMS),
  "\n"
)

cat("\nHallmarks músculo:\n")
cat(
  length(hallmarks_musculo),
  "\n"
)


if (
  length(hallmarks_RMS) != 5 ||
  length(hallmarks_musculo) != 8
) {
  stop(
    "ERROR: se esperaba distribución 5 RMS / 8 músculo."
  )
}


# ------------------------------------------------------------
# 6. Preparar colección Hallmark
# ------------------------------------------------------------

hallmark$gene_set <- as.character(
  hallmark$gene_set
)

hallmark$gene <- trimws(
  as.character(
    hallmark$gene
  )
)


pathways <- split(
  hallmark$gene,
  hallmark$gene_set
)

pathways <- lapply(
  pathways,
  unique
)


if (!all(
  c(
    hallmarks_RMS,
    hallmarks_musculo
  ) %in% names(pathways)
)) {
  stop(
    "ERROR: algún Hallmark replicado falta en MSigDB."
  )
}


# ------------------------------------------------------------
# 7. Genes disponibles de cada eje en el universo común
#
# Se toma la unión de todos los genes pertenecientes a los
# Hallmarks de cada eje, restringida a los 9428 genes.
# ------------------------------------------------------------

genes_Hallmark_RMS <- unique(
  unlist(
    pathways[
      hallmarks_RMS
    ],
    use.names = FALSE
  )
)

genes_Hallmark_musculo <- unique(
  unlist(
    pathways[
      hallmarks_musculo
    ],
    use.names = FALSE
  )
)


genes_Hallmark_RMS <- intersect(
  genes_Hallmark_RMS,
  universo
)

genes_Hallmark_musculo <- intersect(
  genes_Hallmark_musculo,
  universo
)


compartidos_Hallmark <- intersect(
  genes_Hallmark_RMS,
  genes_Hallmark_musculo
)

union_Hallmark <- union(
  genes_Hallmark_RMS,
  genes_Hallmark_musculo
)


jaccard_Hallmark <- (
  length(compartidos_Hallmark) /
  length(union_Hallmark)
)


# ------------------------------------------------------------
# 8. Genes observados en leading edges compartidos
# ------------------------------------------------------------

resumen_genes <- mapa$resumen_genes


if (!all(
  c(
    "gene",
    "n_Hallmarks_RMS",
    "n_Hallmarks_musculo"
  ) %in% colnames(resumen_genes)
)) {
  stop(
    "ERROR: resumen_genes tiene estructura inesperada."
  )
}


genes_LE_RMS <- resumen_genes$gene[
  resumen_genes$n_Hallmarks_RMS > 0
]

genes_LE_musculo <- resumen_genes$gene[
  resumen_genes$n_Hallmarks_musculo > 0
]


compartidos_LE <- intersect(
  genes_LE_RMS,
  genes_LE_musculo
)

union_LE <- union(
  genes_LE_RMS,
  genes_LE_musculo
)


jaccard_LE <- (
  length(compartidos_LE) /
  length(union_LE)
)


# ------------------------------------------------------------
# 9. De los genes estructuralmente disponibles en ambos ejes,
#    ¿cuántos aparecen en leading edge de RMS, músculo o ambos?
# ------------------------------------------------------------

detalle_compartidos <- data.frame(
  gene = sort(compartidos_Hallmark),
  stringsAsFactors = FALSE
)


detalle_compartidos$leadingEdge_RMS <- (
  detalle_compartidos$gene %in%
  genes_LE_RMS
)


detalle_compartidos$leadingEdge_musculo <- (
  detalle_compartidos$gene %in%
  genes_LE_musculo
)


detalle_compartidos$leadingEdge_ambos <- (
  detalle_compartidos$leadingEdge_RMS &
  detalle_compartidos$leadingEdge_musculo
)


detalle_compartidos$clasificacion <- ifelse(

  detalle_compartidos$leadingEdge_ambos,

  "Ambos",

  ifelse(
    detalle_compartidos$leadingEdge_RMS,
    "Solo_RMS",

    ifelse(
      detalle_compartidos$leadingEdge_musculo,
      "Solo_musculo",
      "Ninguno"
    )
  )
)


# ------------------------------------------------------------
# 10. Resumen
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("ESTRUCTURA DE LOS HALLMARKS REPLICADOS\n")
cat("====================================================\n")


cat("\nGenes disponibles en los 5 Hallmarks RMS:\n")
cat(
  length(genes_Hallmark_RMS),
  "\n"
)


cat("\nGenes disponibles en los 8 Hallmarks músculo:\n")
cat(
  length(genes_Hallmark_musculo),
  "\n"
)


cat("\nGenes que pertenecen estructuralmente a Hallmarks")
cat(" de ambos ejes:\n")
cat(
  length(compartidos_Hallmark),
  "\n"
)


cat("\nJaccard entre las uniones Hallmark RMS y músculo:\n")
cat(
  jaccard_Hallmark,
  "\n"
)


cat("\n====================================================\n")
cat("LEADING EDGES OBSERVADOS\n")
cat("====================================================\n")


cat("\nGenes únicos de leading edge hacia RMS:\n")
cat(
  length(genes_LE_RMS),
  "\n"
)


cat("\nGenes únicos de leading edge hacia músculo:\n")
cat(
  length(genes_LE_musculo),
  "\n"
)


cat("\nGenes presentes en leading edges de ambos ejes:\n")
cat(
  length(compartidos_LE),
  "\n"
)


cat("\nJaccard entre genes de leading edge RMS y músculo:\n")
cat(
  jaccard_LE,
  "\n"
)


cat("\n====================================================\n")
cat("GENES ESTRUCTURALMENTE DISPONIBLES EN AMBOS EJES\n")
cat("====================================================\n")


if (nrow(detalle_compartidos) == 0) {

  cat(
    "\nNo existen genes pertenecientes a Hallmarks de ambos ejes.\n"
  )

} else {

  cat("\nClasificación según leading edge:\n")

  print(
    table(
      detalle_compartidos$clasificacion
    )
  )


  cat("\nDetalle:\n\n")

  print(
    detalle_compartidos,
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 11. Tabla resumen
# ------------------------------------------------------------

resumen <- data.frame(

  medida = c(
    "Genes_Hallmark_RMS",
    "Genes_Hallmark_musculo",
    "Genes_Hallmark_ambos_ejes",
    "Jaccard_Hallmark_entre_ejes",
    "Genes_leadingEdge_RMS",
    "Genes_leadingEdge_musculo",
    "Genes_leadingEdge_ambos_ejes",
    "Jaccard_leadingEdge_entre_ejes"
  ),

  valor = c(
    length(genes_Hallmark_RMS),
    length(genes_Hallmark_musculo),
    length(compartidos_Hallmark),
    jaccard_Hallmark,
    length(genes_LE_RMS),
    length(genes_LE_musculo),
    length(compartidos_LE),
    jaccard_LE
  ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12. Guardar
# ------------------------------------------------------------

write.csv(
  resumen,
  "Hallmark_auditoria_solapamiento_entre_ejes.csv",
  row.names = FALSE
)


write.csv(
  detalle_compartidos,
  "Hallmark_genes_estructuralmente_compartidos_entre_ejes.csv",
  row.names = FALSE
)


saveRDS(
  list(
    resumen = resumen,
    genes_Hallmark_RMS = genes_Hallmark_RMS,
    genes_Hallmark_musculo = genes_Hallmark_musculo,
    genes_Hallmark_compartidos = compartidos_Hallmark,
    genes_LE_RMS = genes_LE_RMS,
    genes_LE_musculo = genes_LE_musculo,
    genes_LE_compartidos = compartidos_LE,
    detalle_compartidos = detalle_compartidos
  ),
  "Hallmark_auditoria_solapamiento_entre_ejes.rds"
)


cat("\nArchivos generados:\n")
cat("- Hallmark_auditoria_solapamiento_entre_ejes.csv\n")
cat("- Hallmark_genes_estructuralmente_compartidos_entre_ejes.csv\n")
cat("- Hallmark_auditoria_solapamiento_entre_ejes.rds\n")


cat("\n====================================================\n")
cat("AUDITORÍA ENTRE EJES COMPLETADA\n")
cat("====================================================\n")