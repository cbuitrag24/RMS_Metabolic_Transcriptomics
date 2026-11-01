# ============================================================
# 28_mapear_genes_leading_edges_replicados.R
#
# OBJETIVO:
# Construir un mapa gen x Hallmark para los genes presentes
# en los leading edges compartidos de los 13 Hallmarks
# funcionalmente replicados entre GSE108022 y GSE28511.
#
# Se cuantificará para cada gen:
# - número total de Hallmarks compartidos en los que aparece
# - número de Hallmarks dirigidos hacia RMS
# - número de Hallmarks dirigidos hacia músculo
# - identidad exacta de esos Hallmarks
#
# IMPORTANTE:
# - No se seleccionan genes.
# - No se aplica ningún umbral nuevo.
# - No se asigna importancia biológica a la recurrencia.
# - El análisis es descriptivo y estructural.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_genes <-
  "Hallmark_leadingEdge_genes_compartidos.csv"

archivo_replicados <-
  "Hallmark_replicados_GSE108022_GSE28511.csv"


archivos <- c(
  archivo_genes,
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

genes <- read.csv(
  archivo_genes,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

replicados <- read.csv(
  archivo_replicados,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\nDimensiones asociaciones Hallmark-gen compartidas:\n")
print(
  dim(genes)
)

cat("\nDimensiones Hallmarks replicados:\n")
print(
  dim(replicados)
)


# ------------------------------------------------------------
# 3. Verificar estructura
# ------------------------------------------------------------

if (!all(
  c(
    "pathway",
    "gene",
    "compartido"
  ) %in% colnames(genes)
)) {
  stop(
    "ERROR: faltan columnas necesarias en la tabla de genes."
  )
}


if (!all(
  c(
    "pathway",
    "NES_GSE108022",
    "NES_GSE28511"
  ) %in% colnames(replicados)
)) {
  stop(
    "ERROR: faltan columnas necesarias en la tabla de Hallmarks replicados."
  )
}


if (nrow(replicados) != 13) {
  stop(
    paste(
      "ERROR: se esperaban 13 Hallmarks replicados y se encontraron",
      nrow(replicados)
    )
  )
}


if (anyDuplicated(replicados$pathway)) {
  stop(
    "ERROR: existen Hallmarks duplicados."
  )
}


# ------------------------------------------------------------
# 4. Confirmar que esta tabla contiene sólo asociaciones
#    compartidas
# ------------------------------------------------------------

if (!all(genes$compartido)) {
  stop(
    paste(
      "ERROR: la tabla contiene asociaciones que no están",
      "compartidas entre ambos estudios."
    )
  )
}


if (anyNA(genes$gene) || any(genes$gene == "")) {
  stop(
    "ERROR: existen genes vacíos."
  )
}


if (
  anyNA(genes$pathway) ||
  any(genes$pathway == "")
) {
  stop(
    "ERROR: existen pathways vacíos."
  )
}


# ------------------------------------------------------------
# 5. Verificar asociaciones únicas
# ------------------------------------------------------------

clave <- paste(
  genes$pathway,
  genes$gene,
  sep = "|||"
)


if (anyDuplicated(clave)) {
  stop(
    "ERROR: existen asociaciones Hallmark-gen duplicadas."
  )
}


cat("\nAsociaciones Hallmark-gen compartidas:\n")
cat(
  nrow(genes),
  "\n"
)


cat("\nGenes únicos:\n")
cat(
  length(unique(genes$gene)),
  "\n"
)


# Los valores esperados provienen directamente del script 25.

if (nrow(genes) != 495) {
  stop(
    paste(
      "ERROR: se esperaban 495 asociaciones y se encontraron",
      nrow(genes)
    )
  )
}


if (length(unique(genes$gene)) != 374) {
  stop(
    paste(
      "ERROR: se esperaban 374 genes únicos y se encontraron",
      length(unique(genes$gene))
    )
  )
}


# ------------------------------------------------------------
# 6. Determinar dirección de cada Hallmark
#
# Como los 13 Hallmarks fueron replicados, el signo de NES
# debe coincidir entre los dos estudios.
# ------------------------------------------------------------

replicados$signo_108022 <- sign(
  replicados$NES_GSE108022
)

replicados$signo_28511 <- sign(
  replicados$NES_GSE28511
)


if (any(
  replicados$signo_108022 !=
  replicados$signo_28511
)) {
  stop(
    "ERROR: existe un Hallmark replicado con dirección discordante."
  )
}


if (any(replicados$signo_108022 == 0)) {
  stop(
    "ERROR: existe un Hallmark con NES igual a cero."
  )
}


replicados$eje <- ifelse(
  replicados$signo_108022 > 0,
  "RMS",
  "Musculo"
)


cat("\nHallmarks hacia RMS:\n")
cat(
  sum(replicados$eje == "RMS"),
  "\n"
)

cat("\nHallmarks hacia músculo:\n")
cat(
  sum(replicados$eje == "Musculo"),
  "\n"
)


if (
  sum(replicados$eje == "RMS") != 5 ||
  sum(replicados$eje == "Musculo") != 8
) {
  stop(
    "ERROR: la distribución esperada 5 RMS / 8 músculo no coincide."
  )
}


# ------------------------------------------------------------
# 7. Incorporar dirección a cada asociación Hallmark-gen
# ------------------------------------------------------------

idx <- match(
  genes$pathway,
  replicados$pathway
)


if (anyNA(idx)) {

  faltantes_pathway <- unique(
    genes$pathway[
      is.na(idx)
    ]
  )

  cat("\nHallmarks sin correspondencia:\n")
  print(faltantes_pathway)

  stop(
    "ERROR: existen asociaciones pertenecientes a Hallmarks no replicados."
  )
}


genes$eje <- replicados$eje[
  idx
]


# ------------------------------------------------------------
# 8. Crear matriz gen x Hallmark
# ------------------------------------------------------------

genes_unicos <- sort(
  unique(
    genes$gene
  )
)

hallmarks_ordenados <- c(
  sort(
    replicados$pathway[
      replicados$eje == "RMS"
    ]
  ),
  sort(
    replicados$pathway[
      replicados$eje == "Musculo"
    ]
  )
)


matriz <- matrix(
  0L,
  nrow = length(genes_unicos),
  ncol = length(hallmarks_ordenados),
  dimnames = list(
    genes_unicos,
    hallmarks_ordenados
  )
)


for (i in seq_len(nrow(genes))) {

  matriz[
    genes$gene[i],
    genes$pathway[i]
  ] <- 1L
}


# ------------------------------------------------------------
# 9. Verificar matriz
# ------------------------------------------------------------

if (sum(matriz) != nrow(genes)) {
  stop(
    paste(
      "ERROR: la suma de la matriz no coincide con",
      "el número de asociaciones."
    )
  )
}


if (nrow(matriz) != 374) {
  stop(
    "ERROR: la matriz no contiene 374 genes."
  )
}


if (ncol(matriz) != 13) {
  stop(
    "ERROR: la matriz no contiene 13 Hallmarks."
  )
}


# ------------------------------------------------------------
# 10. Hallmarks RMS y músculo
# ------------------------------------------------------------

hallmarks_RMS <- replicados$pathway[
  replicados$eje == "RMS"
]

hallmarks_musculo <- replicados$pathway[
  replicados$eje == "Musculo"
]


# ------------------------------------------------------------
# 11. Resumen por gen
# ------------------------------------------------------------

n_total <- rowSums(
  matriz
)

n_RMS <- rowSums(
  matriz[
    ,
    hallmarks_RMS,
    drop = FALSE
  ]
)

n_musculo <- rowSums(
  matriz[
    ,
    hallmarks_musculo,
    drop = FALSE
  ]
)


lista_hallmarks <- vapply(

  genes_unicos,

  function(g) {

    presentes <- colnames(matriz)[
      matriz[g, ] == 1
    ]

    paste(
      presentes,
      collapse = ";"
    )
  },

  character(1)
)


lista_RMS <- vapply(

  genes_unicos,

  function(g) {

    presentes <- hallmarks_RMS[
      matriz[g, hallmarks_RMS] == 1
    ]

    paste(
      sort(presentes),
      collapse = ";"
    )
  },

  character(1)
)


lista_musculo <- vapply(

  genes_unicos,

  function(g) {

    presentes <- hallmarks_musculo[
      matriz[g, hallmarks_musculo] == 1
    ]

    paste(
      sort(presentes),
      collapse = ";"
    )
  },

  character(1)
)


resumen_genes <- data.frame(

  gene =
    genes_unicos,

  n_Hallmarks_total =
    n_total,

  n_Hallmarks_RMS =
    n_RMS,

  n_Hallmarks_musculo =
    n_musculo,

  Hallmarks_total =
    lista_hallmarks,

  Hallmarks_RMS =
    lista_RMS,

  Hallmarks_musculo =
    lista_musculo,

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12. Clasificación estructural
#
# Esta clasificación NO implica función biológica.
# Sólo describe en qué eje(s) aparece el gen dentro de
# los leading edges compartidos.
# ------------------------------------------------------------

resumen_genes$patron <- ifelse(

  resumen_genes$n_Hallmarks_RMS > 0 &
  resumen_genes$n_Hallmarks_musculo == 0,

  "Solo_RMS",

  ifelse(

    resumen_genes$n_Hallmarks_RMS == 0 &
    resumen_genes$n_Hallmarks_musculo > 0,

    "Solo_musculo",

    "Ambos_ejes"
  )
)


# ------------------------------------------------------------
# 13. Verificaciones del resumen
# ------------------------------------------------------------

if (nrow(resumen_genes) != 374) {
  stop(
    "ERROR: el resumen no contiene 374 genes."
  )
}


if (any(
  resumen_genes$n_Hallmarks_total !=
  resumen_genes$n_Hallmarks_RMS +
  resumen_genes$n_Hallmarks_musculo
)) {
  stop(
    "ERROR: inconsistencia en los conteos por gen."
  )
}


if (any(resumen_genes$n_Hallmarks_total < 1)) {
  stop(
    "ERROR: existe un gen sin Hallmark asociado."
  )
}


# ------------------------------------------------------------
# 14. Ordenar resumen
# ------------------------------------------------------------

resumen_genes <- resumen_genes[
  order(
    -resumen_genes$n_Hallmarks_total,
    -resumen_genes$n_Hallmarks_RMS,
    -resumen_genes$n_Hallmarks_musculo,
    resumen_genes$gene
  ),
  ,
  drop = FALSE
]

rownames(resumen_genes) <- NULL


# ------------------------------------------------------------
# 15. Resumen global
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("MAPA GEN x HALLMARK REPLICADO\n")
cat("====================================================\n")


cat("\nGenes evaluados:\n")
cat(
  nrow(resumen_genes),
  "\n"
)


cat("\nDistribución del número total de Hallmarks por gen:\n")
print(
  table(
    resumen_genes$n_Hallmarks_total
  )
)


cat("\nDistribución estructural por eje:\n")
print(
  table(
    resumen_genes$patron
  )
)


cat("\nMáximo número de Hallmarks compartidos por un gen:\n")
cat(
  max(resumen_genes$n_Hallmarks_total),
  "\n"
)


cat("\nGenes que alcanzan ese máximo:\n")

print(
  resumen_genes[
    resumen_genes$n_Hallmarks_total ==
      max(resumen_genes$n_Hallmarks_total),
    c(
      "gene",
      "n_Hallmarks_total",
      "n_Hallmarks_RMS",
      "n_Hallmarks_musculo",
      "patron",
      "Hallmarks_total"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Genes presentes en >= 2 Hallmarks
#
# Sólo resumen descriptivo.
# No constituye selección de candidatos.
# ------------------------------------------------------------

recurrentes <- resumen_genes[
  resumen_genes$n_Hallmarks_total >= 2,
  ,
  drop = FALSE
]


cat("\nGenes presentes en >= 2 Hallmarks compartidos:\n")
cat(
  nrow(recurrentes),
  "\n"
)


cat("\nPrimeros genes del resumen de recurrencia:\n\n")

print(
  head(
    recurrentes[
      ,
      c(
        "gene",
        "n_Hallmarks_total",
        "n_Hallmarks_RMS",
        "n_Hallmarks_musculo",
        "patron"
      )
    ],
    40
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Resumen por Hallmark
# ------------------------------------------------------------

resumen_hallmarks <- data.frame(

  pathway =
    hallmarks_ordenados,

  eje =
    replicados$eje[
      match(
        hallmarks_ordenados,
        replicados$pathway
      )
    ],

  n_genes_leadingEdge_compartido =
    colSums(
      matriz[
        ,
        hallmarks_ordenados,
        drop = FALSE
      ]
    ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 18. Preparar matriz para exportación
# ------------------------------------------------------------

matriz_exportar <- data.frame(
  gene = rownames(matriz),
  matriz,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(matriz_exportar) <- NULL


# ------------------------------------------------------------
# 19. Guardar
# ------------------------------------------------------------

write.csv(
  resumen_genes,
  "Hallmark_leadingEdge_mapa_por_gen.csv",
  row.names = FALSE
)


write.csv(
  recurrentes,
  "Hallmark_leadingEdge_genes_recurrentes_descriptivo.csv",
  row.names = FALSE
)


write.csv(
  resumen_hallmarks,
  "Hallmark_leadingEdge_resumen_por_Hallmark.csv",
  row.names = FALSE
)


write.csv(
  matriz_exportar,
  "Hallmark_leadingEdge_matriz_gen_x_Hallmark.csv",
  row.names = FALSE
)


saveRDS(
  list(
    resumen_genes = resumen_genes,
    recurrentes = recurrentes,
    resumen_hallmarks = resumen_hallmarks,
    matriz = matriz
  ),
  "Hallmark_leadingEdge_mapa_genes.rds"
)


cat("\nArchivos generados:\n")

cat(
  "- Hallmark_leadingEdge_mapa_por_gen.csv\n"
)

cat(
  "- Hallmark_leadingEdge_genes_recurrentes_descriptivo.csv\n"
)

cat(
  "- Hallmark_leadingEdge_resumen_por_Hallmark.csv\n"
)

cat(
  "- Hallmark_leadingEdge_matriz_gen_x_Hallmark.csv\n"
)

cat(
  "- Hallmark_leadingEdge_mapa_genes.rds\n"
)


cat("\n====================================================\n")
cat("MAPEO DE GENES COMPLETADO\n")
cat("====================================================\n")