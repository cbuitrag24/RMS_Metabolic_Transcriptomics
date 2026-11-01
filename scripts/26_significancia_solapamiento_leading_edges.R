# ============================================================
# 26_significancia_solapamiento_leading_edges.R
#
# OBJETIVO:
# Evaluar si el solapamiento observado entre los leading edges
# de GSE108022 y GSE28511 para los 13 Hallmarks replicados es
# mayor que el esperado por azar.
#
# Prueba:
#   hipergeométrica unilateral (enriquecimiento)
#
# Universo:
#   9.428 genes comunes utilizados en ambos análisis fgsea
#
# Corrección múltiple:
#   Benjamini-Hochberg sobre las 13 pruebas
#
# IMPORTANTE:
# Esta prueba evalúa exclusivamente el solapamiento entre
# leading edges ya definidos. No redefine pathways, genes,
# rankings ni criterios de replicación.
# ============================================================


# ------------------------------------------------------------
# 1. Archivo de entrada
# ------------------------------------------------------------

archivo <-
  "Hallmark_leadingEdge_concordancia_resumen.csv"


if (!file.exists(archivo)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo
    )
  )
}


# ------------------------------------------------------------
# 2. Cargar
# ------------------------------------------------------------

datos <- read.csv(
  archivo,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\nDimensiones del archivo:\n")
print(
  dim(datos)
)

cat("\nColumnas:\n")
print(
  colnames(datos)
)


# ------------------------------------------------------------
# 3. Verificar estructura
# ------------------------------------------------------------

columnas_necesarias <- c(
  "pathway",
  "n_leadingEdge_GSE108022",
  "n_leadingEdge_GSE28511",
  "n_compartidos",
  "jaccard"
)


if (!all(
  columnas_necesarias %in% colnames(datos)
)) {
  stop(
    paste(
      "ERROR: faltan columnas:",
      paste(
        setdiff(
          columnas_necesarias,
          colnames(datos)
        ),
        collapse = ", "
      )
    )
  )
}


if (nrow(datos) != 13) {
  stop(
    paste(
      "ERROR: se esperaban 13 Hallmarks y se encontraron",
      nrow(datos)
    )
  )
}


if (anyDuplicated(datos$pathway)) {
  stop(
    "ERROR: existen pathways duplicados."
  )
}


# ------------------------------------------------------------
# 4. Universo común
# ------------------------------------------------------------

N <- 9428


cat("\nUniverso utilizado:\n")
cat(
  N,
  "genes\n"
)


# ------------------------------------------------------------
# 5. Verificar tamaños
# ------------------------------------------------------------

if (
  any(datos$n_leadingEdge_GSE108022 <= 0) ||
  any(datos$n_leadingEdge_GSE28511 <= 0)
) {
  stop(
    "ERROR: existe algún leading edge vacío."
  )
}


if (
  any(datos$n_leadingEdge_GSE108022 > N) ||
  any(datos$n_leadingEdge_GSE28511 > N)
) {
  stop(
    "ERROR: algún leading edge supera el universo."
  )
}


if (
  any(datos$n_compartidos >
      datos$n_leadingEdge_GSE108022) ||
  any(datos$n_compartidos >
      datos$n_leadingEdge_GSE28511)
) {
  stop(
    "ERROR: intersección imposible."
  )
}


# ------------------------------------------------------------
# 6. Prueba hipergeométrica unilateral
#
# Para cada Hallmark:
#
# N = universo total
# K = genes leading edge GSE108022
# n = genes leading edge GSE28511
# x = genes compartidos observados
#
# P(X >= x)
# ------------------------------------------------------------

datos$esperados_azar <- (
  datos$n_leadingEdge_GSE108022 *
  datos$n_leadingEdge_GSE28511 /
  N
)


datos$fold_sobre_azar <- (
  datos$n_compartidos /
  datos$esperados_azar
)


datos$p_hipergeometrica <- mapply(

  FUN = function(K, n, x) {

    phyper(
      q = x - 1,
      m = K,
      n = N - K,
      k = n,
      lower.tail = FALSE
    )
  },

  K = datos$n_leadingEdge_GSE108022,
  n = datos$n_leadingEdge_GSE28511,
  x = datos$n_compartidos
)


# ------------------------------------------------------------
# 7. Corrección Benjamini-Hochberg
# ------------------------------------------------------------

datos$FDR_BH <- p.adjust(
  datos$p_hipergeometrica,
  method = "BH"
)


datos$significativo_FDR005 <- (
  datos$FDR_BH < 0.05
)


# ------------------------------------------------------------
# 8. Verificaciones numéricas
# ------------------------------------------------------------

if (
  anyNA(datos$p_hipergeometrica) ||
  anyNA(datos$FDR_BH)
) {
  stop(
    "ERROR: se generaron valores NA."
  )
}


if (
  any(datos$p_hipergeometrica < 0) ||
  any(datos$p_hipergeometrica > 1)
) {
  stop(
    "ERROR: p-value fuera del rango [0,1]."
  )
}


if (
  any(datos$FDR_BH < 0) ||
  any(datos$FDR_BH > 1)
) {
  stop(
    "ERROR: FDR fuera del rango [0,1]."
  )
}


# ------------------------------------------------------------
# 9. Ordenar
# ------------------------------------------------------------

datos <- datos[
  order(
    datos$FDR_BH,
    -datos$jaccard
  ),
  ,
  drop = FALSE
]

rownames(datos) <- NULL


# ------------------------------------------------------------
# 10. Resumen
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("SIGNIFICANCIA DEL SOLAPAMIENTO DE LEADING EDGES\n")
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


cat("\nResumen de genes compartidos observados:\n")
print(
  summary(
    datos$n_compartidos
  )
)


cat("\nResumen de genes compartidos esperados por azar:\n")
print(
  summary(
    datos$esperados_azar
  )
)


cat("\nResumen del fold observado/esperado:\n")
print(
  summary(
    datos$fold_sobre_azar
  )
)


# ------------------------------------------------------------
# 11. Mostrar resultados
# ------------------------------------------------------------

cat("\nResultados por Hallmark:\n\n")


print(
  datos[
    ,
    c(
      "pathway",
      "n_leadingEdge_GSE108022",
      "n_leadingEdge_GSE28511",
      "n_compartidos",
      "esperados_azar",
      "fold_sobre_azar",
      "jaccard",
      "p_hipergeometrica",
      "FDR_BH",
      "significativo_FDR005"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Guardar
# ------------------------------------------------------------

write.csv(
  datos,
  "Hallmark_leadingEdge_significancia_solapamiento.csv",
  row.names = FALSE
)


saveRDS(
  datos,
  "Hallmark_leadingEdge_significancia_solapamiento.rds"
)


cat("\nArchivos generados:\n")

cat(
  "- Hallmark_leadingEdge_significancia_solapamiento.csv\n"
)

cat(
  "- Hallmark_leadingEdge_significancia_solapamiento.rds\n"
)


cat("\n====================================================\n")
cat("PRUEBA DE SOLAPAMIENTO COMPLETADA\n")
cat("====================================================\n")