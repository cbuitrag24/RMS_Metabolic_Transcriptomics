
# ============================================================
# 38_Hallmark_FOXO1FP_vs_FN_directo.R
#
# Enriquecimiento funcional del contraste directo
# FOXO1-FP vs FN, mediante fgseaMultilevel.
#
# Universo comparable: 9427 genes.
# MSigDB Hallmark: misma referencia del script 32.
#
# NES > 0: enriquecimiento hacia FOXO1-FP.
# NES < 0: enriquecimiento hacia FN.
#
# No se filtra el ranking por FDR.
# ============================================================

suppressPackageStartupMessages({
  library(fgsea)
})

# ------------------------------------------------------------
# 1. Archivos originales
# ------------------------------------------------------------

archivo_ranking <-
  "ranking_GSE108022_FOXO1FP_vs_FN_directo_9427.rds"

archivo_universo <-
  "ranking_GSE108022_FN_vs_Normal_universo_comparable.rds"

archivo_hallmark <- file.path(
  "validacion_GSE28511",
  "MSigDB_Hallmark_Homo_sapiens.rds"
)

archivos <- c(
  archivo_ranking,
  archivo_universo,
  archivo_hallmark
)

if (any(!file.exists(archivos))) {
  print(archivos[!file.exists(archivos)])
  stop("ERROR: faltan archivos de entrada.")
}

# ------------------------------------------------------------
# 2. Cargar y verificar datos
# ------------------------------------------------------------

ranking <- readRDS(archivo_ranking)
referencia <- readRDS(archivo_universo)
hallmark <- readRDS(archivo_hallmark)

if (!all(c("gene", "ranking_stat") %in% names(ranking))) {
  stop("ERROR: estructura inesperada del ranking.")
}

if (!("gene" %in% names(referencia))) {
  stop("ERROR: universo de referencia invalido.")
}

if (!all(c("gene_set", "gene") %in% names(hallmark))) {
  stop("ERROR: estructura Hallmark invalida.")
}

ranking$gene <- trimws(as.character(ranking$gene))
referencia$gene <- trimws(as.character(referencia$gene))

if (nrow(ranking) != 9427 ||
    nrow(referencia) != 9427) {
  stop("ERROR: universo distinto de 9427 genes.")
}

if (anyNA(ranking$gene) ||
    any(ranking$gene == "") ||
    anyDuplicated(ranking$gene)) {
  stop("ERROR: genes invalidos o duplicados.")
}

if (!identical(ranking$gene, referencia$gene)) {
  stop("ERROR: el universo no coincide con la referencia.")
}

if (!is.numeric(ranking$ranking_stat) ||
    any(!is.finite(ranking$ranking_stat))) {
  stop("ERROR: estadisticos invalidos.")
}

cat("\nGenes del ranking:", nrow(ranking), "\n")

# ------------------------------------------------------------
# 3. Preparar los 50 conjuntos Hallmark
# ------------------------------------------------------------

hallmark$gene_set <- trimws(
  as.character(hallmark$gene_set)
)

hallmark$gene <- trimws(
  as.character(hallmark$gene)
)

if (anyNA(hallmark$gene_set) ||
    anyNA(hallmark$gene) ||
    any(hallmark$gene_set == "") ||
    any(hallmark$gene == "")) {
  stop("ERROR: identificadores Hallmark invalidos.")
}

hallmark <- unique(hallmark[, c("gene_set", "gene")])

pathways <- split(
  hallmark$gene,
  hallmark$gene_set
)

pathways <- lapply(pathways, unique)

if (length(pathways) != 50) {
  stop("ERROR: se esperaban 50 conjuntos Hallmark.")
}

cat("Conjuntos Hallmark:", length(pathways), "\n")

# ------------------------------------------------------------
# 4. Auditar cobertura
# ------------------------------------------------------------

universo <- ranking$gene

cobertura <- data.frame(
  pathway = names(pathways),
  n_original = vapply(
    pathways,
    length,
    integer(1)
  ),
  n_universo = vapply(
    pathways,
    function(g) length(intersect(g, universo)),
    integer(1)
  ),
  row.names = NULL
)

if (any(cobertura$n_universo < 15) ||
    any(cobertura$n_universo > 500)) {
  stop(
    "ERROR: algun Hallmark queda fuera de minSize/maxSize."
  )
}

# ------------------------------------------------------------
# 5. Construir vector de estadisticos
# ------------------------------------------------------------

stats <- ranking$ranking_stat
names(stats) <- ranking$gene

stats <- sort(stats, decreasing = TRUE)

if (length(stats) != 9427 ||
    anyDuplicated(names(stats)) ||
    any(!is.finite(stats))) {
  stop("ERROR: vector de ranking invalido.")
}

cat("Genes en stats:", length(stats), "\n")

# ------------------------------------------------------------
# 6. Ejecutar fgseaMultilevel
# ------------------------------------------------------------

set.seed(20261007)

resultado <- fgseaMultilevel(
  pathways = pathways,
  stats = stats,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

resultado <- as.data.frame(resultado)

resultado <- resultado[
  order(resultado$padj, -abs(resultado$NES)),
  ,
  drop = FALSE
]

rownames(resultado) <- NULL

if (nrow(resultado) != 50) {
  stop("ERROR: no se evaluaron los 50 Hallmarks.")
}

# ------------------------------------------------------------
# 7. Auditoria de resultados
# ------------------------------------------------------------

cat("\nHallmarks evaluados:", nrow(resultado), "\n")

cat(
  "Hallmarks con FDR < 0.05:",
  sum(resultado$padj < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Significativos hacia FOXO1-FP:",
  sum(
    resultado$padj < 0.05 &
      resultado$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significativos hacia FN:",
  sum(
    resultado$padj < 0.05 &
      resultado$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)

cat("\nResultados significativos:\n")

print(
  resultado[
    !is.na(resultado$padj) &
      resultado$padj < 0.05,
    c("pathway", "NES", "pval", "padj", "size"),
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 8. Guardar resultados completos
# ------------------------------------------------------------

archivo_rds <-
  "Hallmark_fgsea_GSE108022_FOXO1FP_vs_FN_directo_9427.rds"

archivo_csv <-
  "Hallmark_fgsea_GSE108022_FOXO1FP_vs_FN_directo_9427.csv"

saveRDS(resultado, archivo_rds)

resultado_csv <- resultado

# Conservar la informacion de leadingEdge en CSV.
# El RDS mantiene la columna original tipo lista.

resultado_csv$leadingEdge <- vapply(
  resultado$leadingEdge,
  function(g) paste(g, collapse = ";"),
  character(1)
)

write.csv(
  resultado_csv,
  archivo_csv,
  row.names = FALSE
)

write.csv(
  cobertura,
  "Hallmark_cobertura_FOXO1FP_vs_FN_directo_9427.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("-", archivo_rds, "\n")
cat("-", archivo_csv, "\n")
cat("- Hallmark_cobertura_FOXO1FP_vs_FN_directo_9427.csv\n")

cat("\n=============================================\n")
cat("ANALISIS HALLMARK DIRECTO FINALIZADO\n")
cat("=============================================\n")