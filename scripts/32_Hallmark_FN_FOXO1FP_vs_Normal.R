# ============================================================
# 32_Hallmark_FN_FOXO1FP_vs_Normal.R
#
# Enriquecimiento Hallmark:
#   FN vs Normal
#   FOXO1-FP vs Normal
#
# Universo comparable: intersección auditada de 9427 genes.
#
# MSigDB Hallmark 2026.1.Hs
# fgseaMultilevel
# minSize = 15
# maxSize = 500
# eps = 0
# Semilla = 20261007
#
# NES > 0: enriquecimiento hacia el subgrupo RMS
# NES < 0: enriquecimiento hacia músculo normal
#
# CORRECCIÓN:
# Exportación segura de columnas tipo lista a CSV.
# Los RDS conservan los objetos originales completos.
# ============================================================

suppressPackageStartupMessages({
  library(fgsea)
})

# ------------------------------------------------------------
# 1. Definir archivos
# ------------------------------------------------------------

archivo_FN <-
  "ranking_GSE108022_FN_vs_Normal_universo_comparable.rds"

archivo_FP <-
  "ranking_GSE108022_FOXO1FP_vs_Normal_universo_comparable.rds"

archivo_Hallmark <-
  file.path(
    "validacion_GSE28511",
    "MSigDB_Hallmark_Homo_sapiens.rds"
  )

archivos_entrada <- c(
  archivo_FN,
  archivo_FP,
  archivo_Hallmark
)

faltantes <- archivos_entrada[
  !file.exists(archivos_entrada)
]

if (length(faltantes) > 0) {
  print(faltantes)
  stop("ERROR: faltan archivos de entrada.")
}

# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

FN <- readRDS(archivo_FN)
FP <- readRDS(archivo_FP)
hallmark <- readRDS(archivo_Hallmark)

cat("\nDimensiones ranking FN:\n")
print(dim(FN))

cat("\nDimensiones ranking FOXO1-FP:\n")
print(dim(FP))

cat("\nDimensiones colección Hallmark:\n")
print(dim(hallmark))

# ------------------------------------------------------------
# 3. Auditar rankings
# ------------------------------------------------------------

columnas_necesarias <- c(
  "gene",
  "ranking_stat"
)

if (!all(columnas_necesarias %in% names(FN))) {
  stop("ERROR: faltan columnas en FN.")
}

if (!all(columnas_necesarias %in% names(FP))) {
  stop("ERROR: faltan columnas en FOXO1-FP.")
}

if (!all(c("gene_set", "gene") %in% names(hallmark))) {
  stop("ERROR: estructura inesperada de Hallmark.")
}

FN$gene <- trimws(as.character(FN$gene))
FP$gene <- trimws(as.character(FP$gene))

if (anyNA(FN$gene) || any(FN$gene == "")) {
  stop("ERROR: genes inválidos en FN.")
}

if (anyNA(FP$gene) || any(FP$gene == "")) {
  stop("ERROR: genes inválidos en FOXO1-FP.")
}

if (anyDuplicated(FN$gene)) {
  stop("ERROR: genes duplicados en FN.")
}

if (anyDuplicated(FP$gene)) {
  stop("ERROR: genes duplicados en FOXO1-FP.")
}

if (nrow(FN) != 9427 || nrow(FP) != 9427) {
  stop("ERROR: tamaño inesperado del universo comparable.")
}

if (!identical(FN$gene, FP$gene)) {
  stop(
    "ERROR: FN y FOXO1-FP no tienen idénticos genes y orden."
  )
}

if (
  !is.numeric(FN$ranking_stat) ||
  !is.numeric(FP$ranking_stat)
) {
  stop("ERROR: ranking_stat no es numérico.")
}

if (
  any(!is.finite(FN$ranking_stat)) ||
  any(!is.finite(FP$ranking_stat))
) {
  stop("ERROR: hay estadísticos no finitos.")
}

cat("\nGenes del universo comparable:\n")
cat(nrow(FN), "\n")

# ------------------------------------------------------------
# 4. Preparar Hallmarks
# ------------------------------------------------------------

hallmark$gene_set <- as.character(hallmark$gene_set)
hallmark$gene <- trimws(as.character(hallmark$gene))

if (
  anyNA(hallmark$gene_set) ||
  anyNA(hallmark$gene) ||
  any(hallmark$gene_set == "") ||
  any(hallmark$gene == "")
) {
  stop("ERROR: identificadores inválidos en Hallmark.")
}

hallmark <- unique(
  hallmark[, c("gene_set", "gene")]
)

pathways <- split(
  hallmark$gene,
  hallmark$gene_set
)

pathways <- lapply(pathways, unique)

cat("\nHallmarks originales:\n")
cat(length(pathways), "\n")

if (length(pathways) != 50) {
  stop("ERROR: se esperaban 50 Hallmarks.")
}

# ------------------------------------------------------------
# 5. Auditar cobertura
# ------------------------------------------------------------

universo <- FN$gene

cobertura <- data.frame(
  pathway = names(pathways),

  n_original = vapply(
    pathways,
    length,
    integer(1)
  ),

  n_universo = vapply(
    pathways,
    function(x) {
      length(intersect(x, universo))
    },
    integer(1)
  ),

  stringsAsFactors = FALSE
)

cobertura$proporcion_cubierta <-
  cobertura$n_universo / cobertura$n_original

cat("\nCobertura de Hallmarks:\n")
print(summary(cobertura$n_universo))

cat("\nProporción cubierta:\n")
print(summary(cobertura$proporcion_cubierta))

# ------------------------------------------------------------
# 6. Preparar rankings
# ------------------------------------------------------------

stats_FN <- FN$ranking_stat
names(stats_FN) <- FN$gene

stats_FP <- FP$ranking_stat
names(stats_FP) <- FP$gene

stats_FN <- sort(stats_FN, decreasing = TRUE)
stats_FP <- sort(stats_FP, decreasing = TRUE)

if (
  anyDuplicated(names(stats_FN)) ||
  anyDuplicated(names(stats_FP))
) {
  stop("ERROR: genes duplicados en rankings.")
}

if (
  !setequal(names(stats_FN), names(stats_FP))
) {
  stop("ERROR: universos distintos.")
}

# ------------------------------------------------------------
# 7. fgsea FN vs Normal
# ------------------------------------------------------------

set.seed(20261007)

fgsea_FN <- fgseaMultilevel(
  pathways = pathways,
  stats = stats_FN,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

fgsea_FN <- as.data.frame(fgsea_FN)

# ------------------------------------------------------------
# 8. fgsea FOXO1-FP vs Normal
# ------------------------------------------------------------

set.seed(20261007)

fgsea_FP <- fgseaMultilevel(
  pathways = pathways,
  stats = stats_FP,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

fgsea_FP <- as.data.frame(fgsea_FP)

# ------------------------------------------------------------
# 9. Auditar Hallmarks analizados
# ------------------------------------------------------------

cat("\nHallmarks analizados FN:\n")
cat(nrow(fgsea_FN), "\n")

cat("\nHallmarks analizados FOXO1-FP:\n")
cat(nrow(fgsea_FP), "\n")

if (
  nrow(fgsea_FN) != 50 ||
  nrow(fgsea_FP) != 50
) {
  stop("ERROR: no se analizaron los 50 Hallmarks.")
}

if (!setequal(fgsea_FN$pathway, fgsea_FP$pathway)) {
  stop("ERROR: los Hallmarks analizados no coinciden.")
}

if (
  anyNA(fgsea_FN$pval) ||
  anyNA(fgsea_FP$pval) ||
  anyNA(fgsea_FN$padj) ||
  anyNA(fgsea_FP$padj)
) {
  stop("ERROR: hay valores p o FDR ausentes.")
}

# ------------------------------------------------------------
# 10. Verificar corrección BH
# ------------------------------------------------------------

BH_FN <- p.adjust(
  fgsea_FN$pval,
  method = "BH"
)

BH_FP <- p.adjust(
  fgsea_FP$pval,
  method = "BH"
)

error_BH_FN <- max(
  abs(BH_FN - fgsea_FN$padj)
)

error_BH_FP <- max(
  abs(BH_FP - fgsea_FP$padj)
)

cat("\nError máximo padj vs BH FN:\n")
cat(error_BH_FN, "\n")

cat("\nError máximo padj vs BH FOXO1-FP:\n")
cat(error_BH_FP, "\n")

if (
  !isTRUE(all.equal(BH_FN, fgsea_FN$padj)) ||
  !isTRUE(all.equal(BH_FP, fgsea_FP$padj))
) {
  stop("ERROR: padj no coincide con BH.")
}

# ------------------------------------------------------------
# 11. Ordenar resultados
# ------------------------------------------------------------

resultado_FN <- fgsea_FN[
  order(
    fgsea_FN$padj,
    -abs(fgsea_FN$NES)
  ),
  ,
  drop = FALSE
]

resultado_FP <- fgsea_FP[
  order(
    fgsea_FP$padj,
    -abs(fgsea_FP$NES)
  ),
  ,
  drop = FALSE
]

rownames(resultado_FN) <- NULL
rownames(resultado_FP) <- NULL

# ------------------------------------------------------------
# 12. Resumir resultados
# ------------------------------------------------------------

sig_FN <- resultado_FN$padj < 0.05
sig_FP <- resultado_FP$padj < 0.05

cat("\n============================================\n")
cat("HALLMARK FN vs NORMAL\n")
cat("============================================\n")

cat("\nSignificativos FDR < 0.05:\n")
cat(sum(sig_FN), "\n")

cat("\nSignificativos hacia FN:\n")
cat(sum(sig_FN & resultado_FN$NES > 0), "\n")

cat("\nSignificativos hacia músculo normal:\n")
cat(sum(sig_FN & resultado_FN$NES < 0), "\n")

cat("\nResultados FN:\n")

print(
  resultado_FN[
    ,
    c("pathway", "size", "NES", "pval", "padj")
  ],
  row.names = FALSE
)

cat("\n============================================\n")
cat("HALLMARK FOXO1-FP vs NORMAL\n")
cat("============================================\n")

cat("\nSignificativos FDR < 0.05:\n")
cat(sum(sig_FP), "\n")

cat("\nSignificativos hacia FOXO1-FP:\n")
cat(sum(sig_FP & resultado_FP$NES > 0), "\n")

cat("\nSignificativos hacia músculo normal:\n")
cat(sum(sig_FP & resultado_FP$NES < 0), "\n")

cat("\nResultados FOXO1-FP:\n")

print(
  resultado_FP[
    ,
    c("pathway", "size", "NES", "pval", "padj")
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Comparación descriptiva
# ------------------------------------------------------------

comparacion <- merge(
  resultado_FN[
    ,
    c("pathway", "size", "NES", "pval", "padj")
  ],

  resultado_FP[
    ,
    c("pathway", "size", "NES", "pval", "padj")
  ],

  by = "pathway",

  suffixes = c("_FN", "_FOXO1FP")
)

if (
  nrow(comparacion) != 50 ||
  any(comparacion$size_FN != comparacion$size_FOXO1FP)
) {
  stop("ERROR: comparación Hallmark inconsistente.")
}

pearson_NES <- cor(
  comparacion$NES_FN,
  comparacion$NES_FOXO1FP,
  method = "pearson"
)

spearman_NES <- cor(
  comparacion$NES_FN,
  comparacion$NES_FOXO1FP,
  method = "spearman"
)

misma_direccion_NES <- mean(
  sign(comparacion$NES_FN) ==
    sign(comparacion$NES_FOXO1FP)
)

cat("\n============================================\n")
cat("CONCORDANCIA GLOBAL DE HALLMARKS\n")
cat("============================================\n")

cat("\nPearson NES:", pearson_NES, "\n")
cat("Spearman NES:", spearman_NES, "\n")
cat("Misma dirección NES:", misma_direccion_NES, "\n")

# ------------------------------------------------------------
# 14. Convertir columnas lista SOLO para CSV
#
# Los RDS conservan las columnas originales,
# especialmente leadingEdge.
# ------------------------------------------------------------

preparar_csv <- function(tabla) {

  salida <- tabla

  columnas_lista <- names(salida)[
    vapply(salida, is.list, logical(1))
  ]

  for (columna in columnas_lista) {

    salida[[columna]] <- vapply(
      salida[[columna]],
      function(x) {
        paste(as.character(x), collapse = ";")
      },
      character(1)
    )
  }

  return(salida)
}

resultado_FN_csv <- preparar_csv(resultado_FN)
resultado_FP_csv <- preparar_csv(resultado_FP)

if (
  any(vapply(resultado_FN_csv, is.list, logical(1))) ||
  any(vapply(resultado_FP_csv, is.list, logical(1)))
) {
  stop("ERROR: persisten columnas lista en CSV.")
}

# ------------------------------------------------------------
# 15. Guardar resultados
# ------------------------------------------------------------

archivo_FN_csv <-
  "Hallmark_fgsea_GSE108022_FN_vs_Normal_universo_comparable.csv"

archivo_FN_rds <-
  "Hallmark_fgsea_GSE108022_FN_vs_Normal_universo_comparable.rds"

archivo_FP_csv <-
  "Hallmark_fgsea_GSE108022_FOXO1FP_vs_Normal_universo_comparable.csv"

archivo_FP_rds <-
  "Hallmark_fgsea_GSE108022_FOXO1FP_vs_Normal_universo_comparable.rds"

archivo_cobertura <-
  "Hallmark_cobertura_universo_comparable_subgrupos.csv"

archivo_comparacion <-
  "Hallmark_comparacion_FN_FOXO1FP_vs_Normal.csv"

write.csv(
  resultado_FN_csv,
  archivo_FN_csv,
  row.names = FALSE
)

saveRDS(
  resultado_FN,
  archivo_FN_rds
)

write.csv(
  resultado_FP_csv,
  archivo_FP_csv,
  row.names = FALSE
)

saveRDS(
  resultado_FP,
  archivo_FP_rds
)

write.csv(
  cobertura,
  archivo_cobertura,
  row.names = FALSE
)

write.csv(
  comparacion,
  archivo_comparacion,
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Verificar archivos generados
# ------------------------------------------------------------

archivos_salida <- c(
  archivo_FN_csv,
  archivo_FN_rds,
  archivo_FP_csv,
  archivo_FP_rds,
  archivo_cobertura,
  archivo_comparacion
)

estado <- file.exists(archivos_salida)

cat("\n============================================\n")
cat("VERIFICACIÓN DE ARCHIVOS DE SALIDA\n")
cat("============================================\n\n")

print(
  data.frame(
    archivo = archivos_salida,
    existe = estado
  ),
  row.names = FALSE
)

if (!all(estado)) {
  stop("ERROR: faltan archivos de salida.")
}

# ------------------------------------------------------------
# 17. Verificación de integridad de los RDS
# ------------------------------------------------------------

FN_guardado <- readRDS(archivo_FN_rds)
FP_guardado <- readRDS(archivo_FP_rds)

if (
  !identical(resultado_FN, FN_guardado) ||
  !identical(resultado_FP, FP_guardado)
) {
  stop("ERROR: los RDS no conservaron los resultados.")
}

cat("\nIntegridad de RDS FN: OK\n")
cat("Integridad de RDS FOXO1-FP: OK\n")

cat("\n============================================\n")
cat("ENRIQUECIMIENTO HALLMARK COMPLETADO\n")
cat("============================================\n")