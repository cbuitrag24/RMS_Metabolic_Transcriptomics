# ============================================================
# 27_significancia_leading_edges_universo_Hallmark.R
#
# OBJETIVO:
# Evaluar el solapamiento de los leading edges de los
# 13 Hallmarks replicados usando como universo estadístico
# específico los genes de cada Hallmark presentes en el
# universo común de 9.428 genes.
#
# Prueba:
#   hipergeométrica unilateral
#
# Corrección:
#   Benjamini-Hochberg sobre las 13 pruebas
#
# NOTA:
# Este análisis reemplaza, para inferencia estadística,
# al script 26 que utilizó 9.428 como universo para cada
# Hallmark.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_resumen <-
  "Hallmark_leadingEdge_concordancia_resumen.csv"

archivo_hallmark <-
  "MSigDB_Hallmark_Homo_sapiens.rds"

archivo_ranking <-
  "ranking_GSE108022_RMS_vs_musculo_universo_comun_9428.rds"


archivos <- c(
  archivo_resumen,
  archivo_hallmark,
  archivo_ranking
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

datos <- read.csv(
  archivo_resumen,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hallmark <- readRDS(
  archivo_hallmark
)

ranking <- readRDS(
  archivo_ranking
)


cat("\nHallmarks a evaluar:\n")
cat(
  nrow(datos),
  "\n"
)

cat("\nGenes del universo común:\n")
cat(
  nrow(ranking),
  "\n"
)


# ------------------------------------------------------------
# 3. Verificaciones
# ------------------------------------------------------------

if (nrow(datos) != 13) {
  stop(
    "ERROR: se esperaban 13 Hallmarks replicados."
  )
}


if (!all(
  c(
    "pathway",
    "n_leadingEdge_GSE108022",
    "n_leadingEdge_GSE28511",
    "n_compartidos",
    "jaccard"
  ) %in% colnames(datos)
)) {
  stop(
    "ERROR: faltan columnas necesarias en el resumen."
  )
}


if (!all(
  c("gene_set", "gene") %in% colnames(hallmark)
)) {
  stop(
    "ERROR: estructura Hallmark inesperada."
  )
}


if (!("gene" %in% colnames(ranking))) {
  stop(
    "ERROR: el ranking no contiene la columna gene."
  )
}


if (nrow(ranking) != 9428) {
  stop(
    paste(
      "ERROR: se esperaban 9428 genes y se encontraron",
      nrow(ranking)
    )
  )
}


ranking$gene <- trimws(
  as.character(ranking$gene)
)


if (anyDuplicated(ranking$gene)) {
  stop(
    "ERROR: genes duplicados en el universo común."
  )
}


universo_comun <- ranking$gene


# ------------------------------------------------------------
# 4. Preparar Hallmarks
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


# ------------------------------------------------------------
# 5. Calcular tamaño del universo específico de cada Hallmark
# ------------------------------------------------------------

datos$n_Hallmark_universo_comun <- vapply(

  datos$pathway,

  function(p) {

    if (!(p %in% names(pathways))) {
      stop(
        paste(
          "ERROR: no se encontró el Hallmark",
          p
        )
      )
    }

    length(
      intersect(
        pathways[[p]],
        universo_comun
      )
    )
  },

  integer(1)
)


# ------------------------------------------------------------
# 6. Verificaciones de tamaños
# ------------------------------------------------------------

if (
  any(
    datos$n_leadingEdge_GSE108022 >
    datos$n_Hallmark_universo_comun
  )
) {
  stop(
    "ERROR: leading edge GSE108022 mayor que su Hallmark."
  )
}


if (
  any(
    datos$n_leadingEdge_GSE28511 >
    datos$n_Hallmark_universo_comun
  )
) {
  stop(
    "ERROR: leading edge GSE28511 mayor que su Hallmark."
  )
}


if (
  any(
    datos$n_compartidos >
    pmin(
      datos$n_leadingEdge_GSE108022,
      datos$n_leadingEdge_GSE28511
    )
  )
) {
  stop(
    "ERROR: intersección imposible."
  )
}


# ------------------------------------------------------------
# 7. Solapamiento esperado bajo el universo de cada Hallmark
# ------------------------------------------------------------

datos$esperados_azar_Hallmark <- (
  datos$n_leadingEdge_GSE108022 *
  datos$n_leadingEdge_GSE28511 /
  datos$n_Hallmark_universo_comun
)


datos$fold_sobre_azar_Hallmark <- (
  datos$n_compartidos /
  datos$esperados_azar_Hallmark
)


# ------------------------------------------------------------
# 8. Prueba hipergeométrica unilateral
#
# Para cada pathway:
#
# N = genes de ese Hallmark presentes en los 9.428
# K = leading edge GSE108022
# n = leading edge GSE28511
# x = intersección observada
#
# P(X >= x)
# ------------------------------------------------------------

datos$p_hipergeometrica_Hallmark <- mapply(

  FUN = function(N, K, n, x) {

    phyper(
      q = x - 1,
      m = K,
      n = N - K,
      k = n,
      lower.tail = FALSE
    )
  },

  N = datos$n_Hallmark_universo_comun,
  K = datos$n_leadingEdge_GSE108022,
  n = datos$n_leadingEdge_GSE28511,
  x = datos$n_compartidos
)


# ------------------------------------------------------------
# 9. Corrección BH sobre las 13 pruebas
# ------------------------------------------------------------

datos$FDR_BH_Hallmark <- p.adjust(
  datos$p_hipergeometrica_Hallmark,
  method = "BH"
)


datos$significativo_FDR005 <- (
  datos$FDR_BH_Hallmark < 0.05
)


# ------------------------------------------------------------
# 10. Verificaciones
# ------------------------------------------------------------

if (
  anyNA(datos$p_hipergeometrica_Hallmark) ||
  anyNA(datos$FDR_BH_Hallmark)
) {
  stop(
    "ERROR: se produjeron valores NA."
  )
}


if (
  any(datos$p_hipergeometrica_Hallmark < 0) ||
  any(datos$p_hipergeometrica_Hallmark > 1)
) {
  stop(
    "ERROR: p-value fuera del rango [0,1]."
  )
}


# ------------------------------------------------------------
# 11. Ordenar
# ------------------------------------------------------------

datos <- datos[
  order(
    datos$FDR_BH_Hallmark,
    -datos$jaccard
  ),
  ,
  drop = FALSE
]

rownames(datos) <- NULL


# ------------------------------------------------------------
# 12. Resumen
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("LEADING EDGES: UNIVERSO ESPECÍFICO DE CADA HALLMARK\n")
cat("====================================================\n")


cat("\nHallmarks evaluados:\n")
cat(
  nrow(datos),
  "\n"
)


cat("\nHallmarks con FDR BH < 0.05:\n")
cat(
  sum(datos$significativo_FDR005),
  "\n"
)


cat("\nTamaño de los universos Hallmark:\n")
print(
  summary(
    datos$n_Hallmark_universo_comun
  )
)


cat("\nGenes compartidos observados:\n")
print(
  summary(
    datos$n_compartidos
  )
)


cat("\nGenes compartidos esperados bajo cada Hallmark:\n")
print(
  summary(
    datos$esperados_azar_Hallmark
  )
)


cat("\nFold observado/esperado:\n")
print(
  summary(
    datos$fold_sobre_azar_Hallmark
  )
)


# ------------------------------------------------------------
# 13. Resultados completos
# ------------------------------------------------------------

cat("\nResultados por Hallmark:\n\n")


print(
  datos[
    ,
    c(
      "pathway",
      "n_Hallmark_universo_comun",
      "n_leadingEdge_GSE108022",
      "n_leadingEdge_GSE28511",
      "n_compartidos",
      "esperados_azar_Hallmark",
      "fold_sobre_azar_Hallmark",
      "jaccard",
      "p_hipergeometrica_Hallmark",
      "FDR_BH_Hallmark",
      "significativo_FDR005"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Guardar
# ------------------------------------------------------------

write.csv(
  datos,
  "Hallmark_leadingEdge_significancia_universo_Hallmark.csv",
  row.names = FALSE
)


saveRDS(
  datos,
  "Hallmark_leadingEdge_significancia_universo_Hallmark.rds"
)


cat("\nArchivos generados:\n")

cat(
  "- Hallmark_leadingEdge_significancia_universo_Hallmark.csv\n"
)

cat(
  "- Hallmark_leadingEdge_significancia_universo_Hallmark.rds\n"
)


cat("\n====================================================\n")
cat("ANÁLISIS CORREGIDO COMPLETADO\n")
cat("====================================================\n")