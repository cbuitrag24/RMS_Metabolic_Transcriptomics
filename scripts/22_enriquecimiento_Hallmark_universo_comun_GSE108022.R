# ============================================================
# 22_enriquecimiento_Hallmark_universo_comun_GSE108022.R
#
# OBJETIVO:
# Realizar enriquecimiento prerankeado Hallmark en GSE108022
# restringido exactamente al universo común de 9.428 genes
# representados por una única sonda en GSE28511.
#
# Este análisis permitirá posteriormente una comparación
# funcional GSE108022 vs GSE28511 sobre el mismo universo.
#
# Dirección:
#   NES > 0 = enriquecimiento hacia mayor expresión en RMS
#   NES < 0 = enriquecimiento hacia mayor expresión en músculo
#
# IMPORTANTE:
# - No se utiliza ningún resultado de enriquecimiento de GSE28511.
# - No se seleccionan genes por significación.
# - Se utiliza el ranking continuo completo de los 9.428 genes.
# ============================================================


# ------------------------------------------------------------
# 1. Cargar fgsea
# ------------------------------------------------------------

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop(
    "ERROR: fgsea no está instalado."
  )
}

suppressPackageStartupMessages(
  library(fgsea)
)

cat("\nVersión de fgsea:\n")
cat(
  as.character(packageVersion("fgsea")),
  "\n"
)


# ------------------------------------------------------------
# 2. Archivos
# ------------------------------------------------------------

archivo_ranking <-
  "ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.rds"

archivo_hallmark <-
  "MSigDB_Hallmark_Homo_sapiens.rds"


if (!file.exists(archivo_ranking)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_ranking
    )
  )
}

if (!file.exists(archivo_hallmark)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_hallmark
    )
  )
}


# ------------------------------------------------------------
# 3. Cargar objetos
# ------------------------------------------------------------

ranking <- readRDS(
  archivo_ranking
)

hallmark <- readRDS(
  archivo_hallmark
)


cat("\nDimensiones ranking GSE108022 universo común:\n")
print(
  dim(ranking)
)

cat("\nDimensiones Hallmark:\n")
print(
  dim(hallmark)
)


# ------------------------------------------------------------
# 4. Verificar columnas
# ------------------------------------------------------------

if (!all(
  c("gene", "ranking_stat") %in% colnames(ranking)
)) {
  stop(
    "ERROR: faltan columnas gene/ranking_stat en el ranking."
  )
}


if (!all(
  c("gene_set", "gene") %in% colnames(hallmark)
)) {
  stop(
    "ERROR: faltan columnas gene_set/gene en Hallmark."
  )
}


# ------------------------------------------------------------
# 5. Auditar ranking
# ------------------------------------------------------------

ranking$gene <- trimws(
  as.character(ranking$gene)
)

ranking$ranking_stat <- as.numeric(
  ranking$ranking_stat
)


if (nrow(ranking) != 9428) {
  stop(
    paste(
      "ERROR: se esperaban 9428 genes y se encontraron",
      nrow(ranking)
    )
  )
}


if (anyDuplicated(ranking$gene)) {
  stop(
    "ERROR: existen genes duplicados en el ranking."
  )
}


if (
  anyNA(ranking$gene) ||
  any(ranking$gene == "")
) {
  stop(
    "ERROR: existen genes vacíos."
  )
}


if (
  anyNA(ranking$ranking_stat) ||
  any(!is.finite(ranking$ranking_stat))
) {
  stop(
    "ERROR: existen estadísticas inválidas."
  )
}


if (is.unsorted(
  -ranking$ranking_stat
)) {
  stop(
    "ERROR: el ranking no está ordenado de mayor a menor."
  )
}


# ------------------------------------------------------------
# 6. Vector nombrado para fgsea
# ------------------------------------------------------------

stats <- ranking$ranking_stat

names(stats) <- ranking$gene


if (anyDuplicated(names(stats))) {
  stop(
    "ERROR: nombres duplicados en stats."
  )
}


cat("\nGenes del ranking:\n")
cat(
  length(stats),
  "\n"
)

cat("\nEstadísticos positivos:\n")
cat(
  sum(stats > 0),
  "\n"
)

cat("\nEstadísticos negativos:\n")
cat(
  sum(stats < 0),
  "\n"
)

cat("\nEstadísticos iguales a cero:\n")
cat(
  sum(stats == 0),
  "\n"
)


# ------------------------------------------------------------
# 7. Preparar Hallmarks
# ------------------------------------------------------------

hallmark$gene_set <- as.character(
  hallmark$gene_set
)

hallmark$gene <- trimws(
  as.character(hallmark$gene)
)


pathways <- split(
  hallmark$gene,
  hallmark$gene_set
)

pathways <- lapply(
  pathways,
  unique
)


cat("\nNúmero de Hallmarks originales:\n")
cat(
  length(pathways),
  "\n"
)


if (length(pathways) != 50) {
  stop(
    paste(
      "ERROR: se esperaban 50 Hallmarks y se encontraron",
      length(pathways)
    )
  )
}


# ------------------------------------------------------------
# 8. Auditar cobertura en el universo común
# ------------------------------------------------------------

cobertura <- data.frame(
  pathway = names(pathways),

  n_Hallmark = vapply(
    pathways,
    length,
    integer(1)
  ),

  n_en_universo_9428 = vapply(
    pathways,
    function(x) {
      sum(x %in% names(stats))
    },
    integer(1)
  ),

  stringsAsFactors = FALSE
)


cobertura$proporcion_en_universo_9428 <- (
  cobertura$n_en_universo_9428 /
  cobertura$n_Hallmark
)


cobertura <- cobertura[
  order(cobertura$pathway),
  ,
  drop = FALSE
]

rownames(cobertura) <- NULL


cat("\n====================================================\n")
cat("COBERTURA HALLMARK - UNIVERSO COMÚN 9428\n")
cat("====================================================\n")


cat("\nGenes Hallmark presentes por conjunto:\n")
print(
  summary(
    cobertura$n_en_universo_9428
  )
)


cat("\nProporción de cada Hallmark representada:\n")
print(
  summary(
    cobertura$proporcion_en_universo_9428
  )
)


# ------------------------------------------------------------
# 9. Criterios de tamaño
#
# Se mantienen exactamente los mismos parámetros utilizados
# en el análisis Hallmark completo de GSE108022.
# ------------------------------------------------------------

min_size <- 15
max_size <- 500


n_elegibles <- sum(
  cobertura$n_en_universo_9428 >= min_size &
  cobertura$n_en_universo_9428 <= max_size
)


cat("\nHallmarks elegibles con minSize = 15 y maxSize = 500:\n")
cat(
  n_elegibles,
  "\n"
)


# ------------------------------------------------------------
# 10. Ejecutar fgseaMultilevel
# ------------------------------------------------------------

set.seed(20261007)


fg <- fgseaMultilevel(
  pathways = pathways,
  stats = stats,
  minSize = min_size,
  maxSize = max_size,
  eps = 0
)


fg <- as.data.frame(
  fg
)


# ------------------------------------------------------------
# 11. Verificar resultado
# ------------------------------------------------------------

columnas_necesarias <- c(
  "pathway",
  "pval",
  "padj",
  "ES",
  "NES",
  "size",
  "leadingEdge"
)


if (!all(
  columnas_necesarias %in% colnames(fg)
)) {
  stop(
    paste(
      "ERROR: faltan columnas en fgsea:",
      paste(
        setdiff(
          columnas_necesarias,
          colnames(fg)
        ),
        collapse = ", "
      )
    )
  )
}


if (anyDuplicated(fg$pathway)) {
  stop(
    "ERROR: pathways duplicados."
  )
}


if (
  anyNA(fg$NES) ||
  any(!is.finite(fg$NES))
) {
  stop(
    "ERROR: existen NES inválidos."
  )
}


if (
  anyNA(fg$pval) ||
  anyNA(fg$padj)
) {
  stop(
    "ERROR: existen pval/padj inválidos."
  )
}


# ------------------------------------------------------------
# 12. Dirección
# ------------------------------------------------------------

fg$direccion <- ifelse(
  fg$NES > 0,
  "Mayor_en_RMS",
  "Mayor_en_musculo"
)


# ------------------------------------------------------------
# 13. Leading edge a texto
# ------------------------------------------------------------

fg$leadingEdge_genes <- vapply(
  fg$leadingEdge,
  function(x) {
    paste(
      x,
      collapse = ";"
    )
  },
  character(1)
)


# ------------------------------------------------------------
# 14. Ordenar
# ------------------------------------------------------------

fg_ordenado <- fg[
  order(
    fg$padj,
    fg$pval,
    -abs(fg$NES)
  ),
  ,
  drop = FALSE
]

rownames(fg_ordenado) <- NULL


# ------------------------------------------------------------
# 15. Resumen
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RESULTADOS HALLMARK - GSE108022 - UNIVERSO 9428\n")
cat("====================================================\n")


cat("\nHallmarks analizados:\n")
cat(
  nrow(fg_ordenado),
  "\n"
)


cat("\nHallmarks con FDR < 0.05:\n")
cat(
  sum(fg_ordenado$padj < 0.05),
  "\n"
)


cat("\nFDR < 0.05 y NES > 0 (hacia RMS):\n")
cat(
  sum(
    fg_ordenado$padj < 0.05 &
    fg_ordenado$NES > 0
  ),
  "\n"
)


cat("\nFDR < 0.05 y NES < 0 (hacia músculo):\n")
cat(
  sum(
    fg_ordenado$padj < 0.05 &
    fg_ordenado$NES < 0
  ),
  "\n"
)


# ------------------------------------------------------------
# 16. Mostrar los 50 Hallmarks
# ------------------------------------------------------------

tabla_mostrar <- fg_ordenado[
  ,
  c(
    "pathway",
    "NES",
    "pval",
    "padj",
    "size",
    "direccion"
  )
]


cat("\nResultados ordenados:\n\n")

print(
  tabla_mostrar,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Guardar
# ------------------------------------------------------------

fg_csv <- fg_ordenado

fg_csv$leadingEdge <- NULL


write.csv(
  fg_csv,
  "Hallmark_fgsea_GSE108022_universo_comun_9428.csv",
  row.names = FALSE
)


saveRDS(
  fg_ordenado,
  "Hallmark_fgsea_GSE108022_universo_comun_9428.rds"
)


write.csv(
  cobertura,
  "Hallmark_cobertura_GSE108022_universo_comun_9428.csv",
  row.names = FALSE
)


cat("\nArchivos generados:\n")

cat(
  "- Hallmark_fgsea_GSE108022_universo_comun_9428.csv\n"
)

cat(
  "- Hallmark_fgsea_GSE108022_universo_comun_9428.rds\n"
)

cat(
  "- Hallmark_cobertura_GSE108022_universo_comun_9428.csv\n"
)


# ------------------------------------------------------------
# 18. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("ENRIQUECIMIENTO GSE108022 UNIVERSO 9428 COMPLETADO\n")
cat("====================================================\n")


cat(
  "\nNES positivo = enriquecimiento hacia RMS.\n"
)

cat(
  "NES negativo = enriquecimiento hacia músculo.\n"
)

cat(
  "\nEste análisis utiliza exactamente el universo génico",
  "definido para la validación primaria de GSE28511.\n"
)

cat(
  "\nTodavía NO se realizó el enriquecimiento Hallmark",
  "de GSE28511 ni la comparación entre estudios.\n"
)