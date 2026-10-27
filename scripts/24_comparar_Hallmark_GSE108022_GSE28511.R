# ============================================================
# 24_comparar_Hallmark_GSE108022_GSE28511.R
#
# OBJETIVO:
# Comparar formalmente los resultados Hallmark obtenidos en
# GSE108022 y en la cohorte independiente GSE28511 utilizando
# exactamente el mismo universo de 9.428 genes.
#
# Criterio predefinido de replicación funcional:
#
#   1. FDR < 0.05 en GSE108022
#   2. FDR < 0.05 en GSE28511
#   3. mismo signo de NES en ambos estudios
#
# IMPORTANTE:
# - No se modifica ningún resultado de fgsea.
# - No se seleccionan pathways después de observar GSE28511.
# - Se comparan los 50 Hallmarks.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_descubrimiento <-
  "Hallmark_fgsea_GSE108022_universo_comun_9428.rds"

archivo_validacion <-
  "Hallmark_fgsea_GSE28511_universo_comun_9428.rds"


if (!file.exists(archivo_descubrimiento)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_descubrimiento
    )
  )
}

if (!file.exists(archivo_validacion)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_validacion
    )
  )
}


# ------------------------------------------------------------
# 2. Cargar resultados
# ------------------------------------------------------------

gse108022 <- readRDS(
  archivo_descubrimiento
)

gse28511 <- readRDS(
  archivo_validacion
)


cat("\nDimensiones GSE108022:\n")
print(
  dim(gse108022)
)

cat("\nDimensiones GSE28511:\n")
print(
  dim(gse28511)
)


# ------------------------------------------------------------
# 3. Verificaciones básicas
# ------------------------------------------------------------

columnas_necesarias <- c(
  "pathway",
  "NES",
  "pval",
  "padj",
  "size"
)


if (!all(
  columnas_necesarias %in% colnames(gse108022)
)) {
  stop(
    "ERROR: faltan columnas necesarias en GSE108022."
  )
}


if (!all(
  columnas_necesarias %in% colnames(gse28511)
)) {
  stop(
    "ERROR: faltan columnas necesarias en GSE28511."
  )
}


if (nrow(gse108022) != 50) {
  stop(
    "ERROR: GSE108022 no contiene exactamente 50 Hallmarks."
  )
}


if (nrow(gse28511) != 50) {
  stop(
    "ERROR: GSE28511 no contiene exactamente 50 Hallmarks."
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


# ------------------------------------------------------------
# 4. Verificar que sean exactamente los mismos Hallmarks
# ------------------------------------------------------------

solo_108022 <- setdiff(
  gse108022$pathway,
  gse28511$pathway
)

solo_28511 <- setdiff(
  gse28511$pathway,
  gse108022$pathway
)


cat("\nHallmarks sólo en GSE108022:\n")
cat(
  length(solo_108022),
  "\n"
)

cat("\nHallmarks sólo en GSE28511:\n")
cat(
  length(solo_28511),
  "\n"
)


if (
  length(solo_108022) != 0 ||
  length(solo_28511) != 0
) {
  stop(
    "ERROR: los estudios no contienen exactamente los mismos Hallmarks."
  )
}


# ------------------------------------------------------------
# 5. Verificar que el tamaño efectivo de cada Hallmark
#    sea idéntico entre estudios
# ------------------------------------------------------------

idx_size <- match(
  gse108022$pathway,
  gse28511$pathway
)


if (anyNA(idx_size)) {
  stop(
    "ERROR: no pudieron alinearse los Hallmarks."
  )
}


size_identico <- (
  gse108022$size ==
  gse28511$size[idx_size]
)


cat("\nHallmarks con tamaño efectivo idéntico:\n")
cat(
  sum(size_identico),
  "/",
  length(size_identico),
  "\n"
)


if (!all(size_identico)) {
  stop(
    "ERROR: los tamaños efectivos de los Hallmarks difieren entre estudios."
  )
}


# ------------------------------------------------------------
# 6. Construir tabla integrada
# ------------------------------------------------------------

tabla <- data.frame(
  pathway = gse108022$pathway,

  NES_GSE108022 = gse108022$NES,

  pval_GSE108022 = gse108022$pval,

  FDR_GSE108022 = gse108022$padj,

  NES_GSE28511 =
    gse28511$NES[
      match(
        gse108022$pathway,
        gse28511$pathway
      )
    ],

  pval_GSE28511 =
    gse28511$pval[
      match(
        gse108022$pathway,
        gse28511$pathway
      )
    ],

  FDR_GSE28511 =
    gse28511$padj[
      match(
        gse108022$pathway,
        gse28511$pathway
      )
    ],

  size =
    gse108022$size,

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 7. Dirección y significación
# ------------------------------------------------------------

tabla$misma_direccion <- (
  sign(tabla$NES_GSE108022) ==
  sign(tabla$NES_GSE28511)
)


tabla$significativo_GSE108022 <- (
  tabla$FDR_GSE108022 < 0.05
)


tabla$significativo_GSE28511 <- (
  tabla$FDR_GSE28511 < 0.05
)


tabla$significativo_ambos <- (
  tabla$significativo_GSE108022 &
  tabla$significativo_GSE28511
)


# ------------------------------------------------------------
# 8. Criterio predefinido de replicación funcional
# ------------------------------------------------------------

tabla$replicado <- (
  tabla$significativo_ambos &
  tabla$misma_direccion
)


tabla$direccion_replicada <- NA_character_


tabla$direccion_replicada[
  tabla$replicado &
  tabla$NES_GSE108022 > 0
] <- "Mayor_en_RMS"


tabla$direccion_replicada[
  tabla$replicado &
  tabla$NES_GSE108022 < 0
] <- "Mayor_en_musculo"


# ------------------------------------------------------------
# 9. Concordancia global de NES
# ------------------------------------------------------------

pearson <- cor(
  tabla$NES_GSE108022,
  tabla$NES_GSE28511,
  method = "pearson"
)


spearman <- cor(
  tabla$NES_GSE108022,
  tabla$NES_GSE28511,
  method = "spearman"
)


n_misma_direccion <- sum(
  tabla$misma_direccion
)


prop_misma_direccion <- mean(
  tabla$misma_direccion
)


# ------------------------------------------------------------
# 10. Resumen de significación
# ------------------------------------------------------------

n_sig_108022 <- sum(
  tabla$significativo_GSE108022
)

n_sig_28511 <- sum(
  tabla$significativo_GSE28511
)

n_sig_ambos <- sum(
  tabla$significativo_ambos
)

n_replicados <- sum(
  tabla$replicado
)


n_replicados_RMS <- sum(
  tabla$replicado &
  tabla$NES_GSE108022 > 0
)


n_replicados_musculo <- sum(
  tabla$replicado &
  tabla$NES_GSE108022 < 0
)


# ------------------------------------------------------------
# 11. Verificar si algún Hallmark significativo en ambos
#     presenta dirección opuesta
# ------------------------------------------------------------

n_sig_ambos_opuestos <- sum(
  tabla$significativo_ambos &
  !tabla$misma_direccion
)


# ------------------------------------------------------------
# 12. Mostrar resumen
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("CONCORDANCIA HALLMARK ENTRE ESTUDIOS\n")
cat("====================================================\n")


cat("\nHallmarks comparados:\n")
cat(
  nrow(tabla),
  "\n"
)


cat("\nCorrelación Pearson de NES:\n")
cat(
  pearson,
  "\n"
)


cat("\nCorrelación Spearman de NES:\n")
cat(
  spearman,
  "\n"
)


cat("\nHallmarks con la misma dirección:\n")
cat(
  n_misma_direccion,
  "/",
  nrow(tabla),
  "=",
  prop_misma_direccion,
  "\n"
)


cat("\nHallmarks FDR < 0.05 en GSE108022:\n")
cat(
  n_sig_108022,
  "\n"
)


cat("\nHallmarks FDR < 0.05 en GSE28511:\n")
cat(
  n_sig_28511,
  "\n"
)


cat("\nHallmarks FDR < 0.05 en ambos estudios:\n")
cat(
  n_sig_ambos,
  "\n"
)


cat("\nHallmarks significativos en ambos pero con dirección opuesta:\n")
cat(
  n_sig_ambos_opuestos,
  "\n"
)


cat("\nHallmarks que cumplen el criterio de replicación funcional:\n")
cat(
  n_replicados,
  "\n"
)


cat("\nReplicados hacia RMS:\n")
cat(
  n_replicados_RMS,
  "\n"
)


cat("\nReplicados hacia músculo:\n")
cat(
  n_replicados_musculo,
  "\n"
)


# ------------------------------------------------------------
# 13. Tabla de Hallmarks replicados
# ------------------------------------------------------------

replicados <- tabla[
  tabla$replicado,
  ,
  drop = FALSE
]


replicados <- replicados[
  order(
    replicados$direccion_replicada,
    -abs(replicados$NES_GSE108022)
  ),
  ,
  drop = FALSE
]

rownames(replicados) <- NULL


cat("\n====================================================\n")
cat("HALLMARKS REPLICADOS\n")
cat("====================================================\n\n")


if (nrow(replicados) == 0) {

  cat(
    "Ningún Hallmark cumplió el criterio predefinido.\n"
  )

} else {

  print(
    replicados[
      ,
      c(
        "pathway",
        "NES_GSE108022",
        "FDR_GSE108022",
        "NES_GSE28511",
        "FDR_GSE28511",
        "direccion_replicada"
      )
    ],
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 14. Ordenar tabla completa
# ------------------------------------------------------------

tabla <- tabla[
  order(
    !tabla$replicado,
    tabla$FDR_GSE108022,
    tabla$FDR_GSE28511
  ),
  ,
  drop = FALSE
]

rownames(tabla) <- NULL


# ------------------------------------------------------------
# 15. Guardar
# ------------------------------------------------------------

write.csv(
  tabla,
  "Hallmark_comparacion_GSE108022_GSE28511.csv",
  row.names = FALSE
)


saveRDS(
  tabla,
  "Hallmark_comparacion_GSE108022_GSE28511.rds"
)


write.csv(
  replicados,
  "Hallmark_replicados_GSE108022_GSE28511.csv",
  row.names = FALSE
)


cat("\nArchivos generados:\n")

cat(
  "- Hallmark_comparacion_GSE108022_GSE28511.csv\n"
)

cat(
  "- Hallmark_comparacion_GSE108022_GSE28511.rds\n"
)

cat(
  "- Hallmark_replicados_GSE108022_GSE28511.csv\n"
)


cat("\n====================================================\n")
cat("COMPARACIÓN HALLMARK COMPLETADA\n")
cat("====================================================\n")