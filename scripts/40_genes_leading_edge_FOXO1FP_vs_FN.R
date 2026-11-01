
# ============================================================
# SCRIPT 40
# Genes leading-edge de dos Hallmarks significativos
# Contraste directo FOXO1-FP vs FN
#
# No recalcula edgeR ni fgsea.
# ============================================================

archivo_fgsea <-
  "Hallmark_fgsea_GSE108022_FOXO1FP_vs_FN_directo_9427.rds"

archivo_ranking <-
  "ranking_GSE108022_FOXO1FP_vs_FN_directo_9427.rds"

if (!file.exists(archivo_fgsea) ||
    !file.exists(archivo_ranking)) {
  stop("ERROR: faltan archivos de entrada.")
}

fgsea <- readRDS(archivo_fgsea)
ranking <- readRDS(archivo_ranking)

programas <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

# 1. Verificar estructura

columnas_ranking <- c(
  "gene", "logFC", "F",
  "PValue", "FDR", "ranking_stat"
)

if (!all(columnas_ranking %in% names(ranking))) {
  stop("ERROR: faltan columnas en ranking.")
}

if (!all(c("pathway", "NES", "padj", "leadingEdge") %in%
         names(fgsea))) {
  stop("ERROR: estructura inesperada de fgsea.")
}

if (anyDuplicated(ranking$gene)) {
  stop("ERROR: genes duplicados.")
}

if (nrow(ranking) != 9427) {
  stop("ERROR: universo inesperado.")
}

if (anyDuplicated(fgsea$pathway)) {
  stop("ERROR: Hallmarks duplicados.")
}

# 2. Extraer cada leading-edge

resumen <- data.frame(
  pathway = character(),
  NES = numeric(),
  FDR_Hallmark = numeric(),
  n_leading_edge = integer(),
  n_genes_FDR005 = integer(),
  n_logFC_positivo = integer(),
  n_logFC_negativo = integer()
)

for (programa in programas) {

  indice <- match(programa, fgsea$pathway)

  if (is.na(indice)) {
    stop(paste("ERROR: Hallmark ausente:", programa))
  }

  genes <- as.character(
    fgsea$leadingEdge[[indice]]
  )

  if (anyNA(genes) ||
      anyDuplicated(genes) ||
      length(genes) == 0) {
    stop("ERROR: leading-edge invalido.")
  }

  posiciones <- match(genes, ranking$gene)

  if (anyNA(posiciones)) {
    stop("ERROR: genes leading-edge ausentes del ranking.")
  }

  tabla <- ranking[
    posiciones,
    columnas_ranking,
    drop = FALSE
  ]

  rownames(tabla) <- NULL

  if (!identical(tabla$gene, genes)) {
    stop("ERROR: identidad de genes inconsistente.")
  }

  tabla$Hallmark <- programa

  tabla$FDR005 <- tabla$FDR < 0.05

  tabla$direccion <- ifelse(
    tabla$logFC > 0,
    "FOXO1-FP",
    ifelse(tabla$logFC < 0, "FN", "Sin cambio")
  )

  nombre <- if (
    programa == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  ) {
    "LeadingEdge_TNFA_NFKB_FOXO1FP_vs_FN.csv"
  } else {
    "LeadingEdge_OXPHOS_FOXO1FP_vs_FN.csv"
  }

  write.csv(
    tabla,
    nombre,
    row.names = FALSE
  )

  resumen <- rbind(
    resumen,
    data.frame(
      pathway = programa,
      NES = fgsea$NES[indice],
      FDR_Hallmark = fgsea$padj[indice],
      n_leading_edge = nrow(tabla),
      n_genes_FDR005 = sum(tabla$FDR005),
      n_logFC_positivo = sum(tabla$logFC > 0),
      n_logFC_negativo = sum(tabla$logFC < 0)
    )
  )

  cat("\nArchivo generado:", nombre, "\n")
}

# 3. Guardar resumen

write.csv(
  resumen,
  "Resumen_leading_edge_FOXO1FP_vs_FN.csv",
  row.names = FALSE
)

cat("\nRESUMEN:\n")
print(resumen, row.names = FALSE)

cat("\nAnalisis leading-edge completado.\n")