# ============================================================
# 23_enriquecimiento_Hallmark_GSE28511_universo_comun_9428.R
#
# OBJETIVO:
# Realizar enriquecimiento prerankeado Hallmark en la cohorte
# independiente GSE28511 utilizando exactamente el universo
# común de 9.428 genes representados por una única sonda.
#
# Ranking:
#   ranking_stat = estadístico t moderado de limma
#
# Dirección:
#   NES > 0 = enriquecimiento hacia mayor expresión en RMS
#   NES < 0 = enriquecimiento hacia mayor expresión en músculo
#
# IMPORTANTE:
# - No se seleccionan genes por FDR ni por logFC.
# - Se utiliza el ranking continuo completo de 9.428 genes.
# - Se mantienen los mismos Hallmarks y parámetros utilizados
#   para GSE108022 en el universo común.
# - Todavía NO se realiza la comparación entre estudios.
# ============================================================


# ------------------------------------------------------------
# 1. Cargar fgsea
# ------------------------------------------------------------

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("ERROR: fgsea no está instalado.")
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
  "ranking_GSE28511_RMS_vs_musculo_sonda_unica.rds"

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


cat("\nDimensiones ranking GSE28511:\n")
print(
  dim(ranking)
)

cat("\nColumnas ranking GSE28511:\n")
print(
  colnames(ranking)
)

cat("\nDimensiones Hallmark:\n")
print(
  dim(hallmark)
)


# ------------------------------------------------------------
# 4. Verificar estructura
# ------------------------------------------------------------

if (!("gene" %in% colnames(ranking))) {
  stop(
    "ERROR: el ranking no contiene la columna gene."
  )
}

if (!("ranking_stat" %in% colnames(ranking))) {
  stop(
    "ERROR: el ranking no contiene la columna ranking_stat."
  )
}


cat("\nColumna utilizada como estadístico de ranking:\n")
cat(
  "ranking_stat\n"
)


# ------------------------------------------------------------
# 5. Auditar ranking
# ------------------------------------------------------------

ranking$gene <- trimws(
  as.character(ranking$gene)
)

ranking_stat <- as.numeric(
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
    "ERROR: existen identificadores génicos vacíos."
  )
}


if (
  anyNA(ranking_stat) ||
  any(!is.finite(ranking_stat))
) {
  stop(
    "ERROR: existen estadísticas de ranking inválidas."
  )
}


# ------------------------------------------------------------
# 6. Verificar que ranking_stat coincide con t
#
# En el script 18 se definió ranking_stat a partir del
# estadístico t moderado de limma.
# Si la columna t está disponible, verificamos que ambas
# columnas sean exactamente equivalentes numéricamente.
# ------------------------------------------------------------

if ("t" %in% colnames(ranking)) {

  t_original <- as.numeric(
    ranking$t
  )

  if (
    anyNA(t_original) ||
    any(!is.finite(t_original))
  ) {
    stop(
      "ERROR: existen valores inválidos en la columna t."
    )
  }

  error_t <- max(
    abs(
      ranking_stat - t_original
    )
  )

  cat("\nError máximo entre ranking_stat y t:\n")
  cat(
    format(
      error_t,
      scientific = TRUE
    ),
    "\n"
  )

  if (error_t != 0) {
    stop(
      paste(
        "ERROR: ranking_stat no coincide exactamente",
        "con el estadístico t."
      )
    )
  }
}


# ------------------------------------------------------------
# 7. Ordenar explícitamente por estadístico decreciente
#
# No se modifica ningún valor.
# ------------------------------------------------------------

orden <- order(
  ranking_stat,
  decreasing = TRUE
)

ranking <- ranking[
  orden,
  ,
  drop = FALSE
]

ranking_stat <- ranking_stat[
  orden
]

rownames(ranking) <- NULL


if (is.unsorted(
  -ranking_stat
)) {
  stop(
    "ERROR: el ranking no quedó ordenado correctamente."
  )
}


# ------------------------------------------------------------
# 8. Crear vector nombrado para fgsea
# ------------------------------------------------------------

stats <- ranking_stat

names(stats) <- ranking$gene


if (anyDuplicated(names(stats))) {
  stop(
    "ERROR: existen nombres duplicados en stats."
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
# 9. Preparar Hallmarks
# ------------------------------------------------------------

if (!all(
  c("gene_set", "gene") %in% colnames(hallmark)
)) {
  stop(
    "ERROR: faltan columnas gene_set/gene en Hallmark."
  )
}


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
# 10. Auditar cobertura
#
# Al utilizar exactamente el mismo universo de 9.428 genes,
# la cobertura debe coincidir con la obtenida para GSE108022
# restringido al universo común.
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
cat("COBERTURA HALLMARK - GSE28511 - UNIVERSO 9428\n")
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
# 11. Parámetros predefinidos
#
# Exactamente los mismos utilizados para GSE108022.
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


if (n_elegibles != 50) {
  stop(
    paste(
      "ERROR: se esperaban 50 Hallmarks elegibles y se encontraron",
      n_elegibles
    )
  )
}


# ------------------------------------------------------------
# 12. Ejecutar fgseaMultilevel
#
# Misma semilla y mismos parámetros que GSE108022.
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
# 13. Verificar resultado
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


if (nrow(fg) != 50) {
  stop(
    paste(
      "ERROR: se esperaban 50 resultados y se obtuvieron",
      nrow(fg)
    )
  )
}


if (anyDuplicated(fg$pathway)) {
  stop(
    "ERROR: existen pathways duplicados."
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
# 14. Dirección
# ------------------------------------------------------------

fg$direccion <- ifelse(
  fg$NES > 0,
  "Mayor_en_RMS",
  "Mayor_en_musculo"
)


# ------------------------------------------------------------
# 15. Leading edge a texto
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
# 16. Ordenar resultados
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
# 17. Resumen
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("RESULTADOS HALLMARK - GSE28511 - UNIVERSO 9428\n")
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
# 18. Mostrar los 50 Hallmarks
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
# 19. Guardar
# ------------------------------------------------------------

fg_csv <- fg_ordenado

fg_csv$leadingEdge <- NULL


write.csv(
  fg_csv,
  "Hallmark_fgsea_GSE28511_universo_comun_9428.csv",
  row.names = FALSE
)


saveRDS(
  fg_ordenado,
  "Hallmark_fgsea_GSE28511_universo_comun_9428.rds"
)


write.csv(
  cobertura,
  "Hallmark_cobertura_GSE28511_universo_comun_9428.csv",
  row.names = FALSE
)


cat("\nArchivos generados:\n")

cat(
  "- Hallmark_fgsea_GSE28511_universo_comun_9428.csv\n"
)

cat(
  "- Hallmark_fgsea_GSE28511_universo_comun_9428.rds\n"
)

cat(
  "- Hallmark_cobertura_GSE28511_universo_comun_9428.csv\n"
)


# ------------------------------------------------------------
# 20. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("ENRIQUECIMIENTO GSE28511 UNIVERSO 9428 COMPLETADO\n")
cat("====================================================\n")


cat(
  "\nNES positivo = enriquecimiento hacia RMS.\n"
)

cat(
  "NES negativo = enriquecimiento hacia músculo.\n"
)

cat(
  "\nTodavía NO se realizó la comparación formal",
  "entre GSE108022 y GSE28511.\n"
)