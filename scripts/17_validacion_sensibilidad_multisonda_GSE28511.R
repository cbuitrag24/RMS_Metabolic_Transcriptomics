# ============================================================
# 17_validacion_sensibilidad_multisonda_GSE28511.R
#
# OBJETIVO:
# Sensibilidad de la concordancia RMS vs músculo para los
# genes comunes representados por >1 sonda en GSE28511.
#
# REGLA PREDEFINIDA:
# El efecto del gen en GSE28511 se resume como la MEDIANA
# de los logFC de TODAS las sondas asignadas a ese Symbol.
#
# Esta regla:
# - no usa GSE108022 para escoger sondas;
# - no usa PValue/FDR para escoger sondas;
# - no selecciona la sonda con mayor efecto;
# - conserva explícitamente la heterogeneidad entre sondas.
#
# IMPORTANTE:
# Este es un análisis de SENSIBILIDAD.
# NO reemplaza la validación primaria de genes de sonda única.
# NO se calcula una FDR génica a partir de las sondas.
# ============================================================


# ------------------------------------------------------------
# 1. Archivos
# ------------------------------------------------------------

archivo_rnaseq <- file.path(
  "..",
  "GPL16791_DE_RMS_vs_Normal_todos_los_genes_CORREGIDO.csv"
)

archivo_micro <-
  "GSE28511_primario_limma_todas_las_sondas.rds"

archivo_universo <-
  "universo_comun_GSE108022_GSE28511.csv"


archivos <- c(
  archivo_rnaseq,
  archivo_micro,
  archivo_universo
)

faltantes <- archivos[!file.exists(archivos)]

if (length(faltantes) > 0) {
  stop(
    paste(
      "ERROR: faltan archivos:",
      paste(faltantes, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

rnaseq <- read.csv(
  archivo_rnaseq,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

micro <- readRDS(
  archivo_micro
)

universo <- read.csv(
  archivo_universo,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("Dimensiones GSE108022:\n")
print(dim(rnaseq))

cat("\nDimensiones GSE28511:\n")
print(dim(micro))

cat("\nDimensiones universo común:\n")
print(dim(universo))


# ------------------------------------------------------------
# 3. Verificar columnas
# ------------------------------------------------------------

necesarias_rna <- c(
  "gene",
  "logFC",
  "PValue",
  "FDR"
)

necesarias_micro <- c(
  "Probe_ID",
  "Symbol",
  "logFC",
  "PValue",
  "FDR"
)

necesarias_universo <- c(
  "gene",
  "n_sondas_GSE28511"
)


if (!all(necesarias_rna %in% colnames(rnaseq))) {
  stop("ERROR: faltan columnas necesarias en GSE108022.")
}

if (!all(necesarias_micro %in% colnames(micro))) {
  stop("ERROR: faltan columnas necesarias en GSE28511.")
}

if (!all(necesarias_universo %in% colnames(universo))) {
  stop("ERROR: faltan columnas necesarias en el universo común.")
}


# ------------------------------------------------------------
# 4. Limpiar identificadores
# ------------------------------------------------------------

rnaseq$gene <- trimws(
  as.character(rnaseq$gene)
)

micro$Symbol <- trimws(
  as.character(micro$Symbol)
)

micro$Symbol[is.na(micro$Symbol)] <- ""

universo$gene <- trimws(
  as.character(universo$gene)
)


# ------------------------------------------------------------
# 5. Seleccionar genes multisondeados
# ------------------------------------------------------------

universo_multi <- universo[
  universo$n_sondas_GSE28511 > 1,
  ,
  drop = FALSE
]


cat("\n====================================================\n")
cat("UNIVERSO MULTISONDA\n")
cat("====================================================\n")

cat("\nGenes comunes totales:\n")
cat(nrow(universo), "\n")

cat("\nGenes con múltiples sondas:\n")
cat(nrow(universo_multi), "\n")

cat("\nNúmero total esperado de sondas para esos genes:\n")
cat(sum(universo_multi$n_sondas_GSE28511), "\n")


if (nrow(universo_multi) != 5125) {
  stop(
    paste(
      "ERROR: se esperaban 5125 genes multisondeados y se encontraron",
      nrow(universo_multi)
    )
  )
}


# ------------------------------------------------------------
# 6. Extraer las sondas correspondientes
# ------------------------------------------------------------

micro_multi <- micro[
  micro$Symbol %in% universo_multi$gene,
  ,
  drop = FALSE
]


cat("\nSondas encontradas:\n")
cat(nrow(micro_multi), "\n")


# Verificar multiplicidad exacta por gen

conteos <- table(
  micro_multi$Symbol
)

conteos_ordenados <- as.integer(
  conteos[
    universo_multi$gene
  ]
)

if (anyNA(conteos_ordenados)) {
  stop(
    "ERROR: faltan genes multisondeados en GSE28511."
  )
}

if (!identical(
  conteos_ordenados,
  as.integer(universo_multi$n_sondas_GSE28511)
)) {
  stop(
    "ERROR: la multiplicidad de sondas no coincide con el universo común."
  )
}


cat(
  "\nVERIFICACIÓN: multiplicidad de sondas",
  "coincide exactamente con el universo común.\n"
)


# ------------------------------------------------------------
# 7. Resumir GSE28511 por gen
#
# Efecto central = mediana de logFC.
#
# También se conservan:
# - número de sondas
# - positivas
# - negativas
# - cero
# - mínimo
# - máximo
# - clasificación direccional
# ------------------------------------------------------------

resumen_lista <- lapply(
  universo_multi$gene,
  function(g) {

    x <- micro_multi[
      micro_multi$Symbol == g,
      ,
      drop = FALSE
    ]

    efectos <- as.numeric(
      x$logFC
    )

    n_pos <- sum(
      efectos > 0
    )

    n_neg <- sum(
      efectos < 0
    )

    n_cero <- sum(
      efectos == 0
    )

    clasificacion <- if (
      n_pos == nrow(x)
    ) {
      "todas_positivas"

    } else if (
      n_neg == nrow(x)
    ) {
      "todas_negativas"

    } else if (
      n_cero == nrow(x)
    ) {
      "todas_cero"

    } else if (
      n_pos > 0 &&
      n_neg > 0
    ) {
      "direcciones_mixtas"

    } else if (
      n_pos > 0 &&
      n_cero > 0
    ) {
      "positivas_y_cero"

    } else if (
      n_neg > 0 &&
      n_cero > 0
    ) {
      "negativas_y_cero"

    } else {
      "otra"
    }

    data.frame(
      gene = g,

      n_sondas_GSE28511 =
        nrow(x),

      n_positivas =
        n_pos,

      n_negativas =
        n_neg,

      n_cero =
        n_cero,

      logFC_micro_min =
        min(efectos),

      logFC_micro_mediana =
        median(efectos),

      logFC_micro_max =
        max(efectos),

      clasificacion_sondas =
        clasificacion,

      stringsAsFactors = FALSE
    )
  }
)


resumen_micro <- do.call(
  rbind,
  resumen_lista
)

rownames(resumen_micro) <- NULL


# ------------------------------------------------------------
# 8. Verificaciones del resumen
# ------------------------------------------------------------

if (nrow(resumen_micro) != 5125) {
  stop(
    "ERROR: número incorrecto de genes en resumen multisondeado."
  )
}

if (anyDuplicated(resumen_micro$gene)) {
  stop(
    "ERROR: genes duplicados en resumen multisondeado."
  )
}

if (!identical(
  resumen_micro$gene,
  universo_multi$gene
)) {
  stop(
    "ERROR: se perdió el orden del universo multisondeado."
  )
}

if (!identical(
  as.integer(resumen_micro$n_sondas_GSE28511),
  as.integer(universo_multi$n_sondas_GSE28511)
)) {
  stop(
    "ERROR: multiplicidad incorrecta después del resumen."
  )
}


cat("\n====================================================\n")
cat("AUDITORÍA DE DIRECCIÓN ENTRE SONDAS\n")
cat("====================================================\n\n")

print(
  table(
    resumen_micro$clasificacion_sondas
  )
)


# ------------------------------------------------------------
# 9. Integrar GSE108022
# ------------------------------------------------------------

idx_rna <- match(
  resumen_micro$gene,
  rnaseq$gene
)

if (anyNA(idx_rna)) {
  stop(
    "ERROR: faltan genes multisondeados en GSE108022."
  )
}

rna_multi <- rnaseq[
  idx_rna,
  ,
  drop = FALSE
]


if (!identical(
  rna_multi$gene,
  resumen_micro$gene
)) {
  stop(
    "ERROR: identidad gen-fila incorrecta al integrar GSE108022."
  )
}


tabla <- data.frame(
  gene =
    resumen_micro$gene,

  n_sondas_GSE28511 =
    resumen_micro$n_sondas_GSE28511,

  n_positivas =
    resumen_micro$n_positivas,

  n_negativas =
    resumen_micro$n_negativas,

  n_cero =
    resumen_micro$n_cero,

  clasificacion_sondas =
    resumen_micro$clasificacion_sondas,

  logFC_GSE108022 =
    as.numeric(rna_multi$logFC),

  PValue_GSE108022 =
    as.numeric(rna_multi$PValue),

  FDR_GSE108022 =
    as.numeric(rna_multi$FDR),

  logFC_GSE28511_mediana =
    resumen_micro$logFC_micro_mediana,

  logFC_GSE28511_min =
    resumen_micro$logFC_micro_min,

  logFC_GSE28511_max =
    resumen_micro$logFC_micro_max,

  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 10. Control de valores
# ------------------------------------------------------------

cols_num <- c(
  "logFC_GSE108022",
  "PValue_GSE108022",
  "FDR_GSE108022",
  "logFC_GSE28511_mediana",
  "logFC_GSE28511_min",
  "logFC_GSE28511_max"
)

for (col in cols_num) {

  if (anyNA(tabla[[col]])) {
    stop(
      paste(
        "ERROR: existen NA en",
        col
      )
    )
  }

  if (any(!is.finite(tabla[[col]]))) {
    stop(
      paste(
        "ERROR: existen valores no finitos en",
        col
      )
    )
  }
}


# ------------------------------------------------------------
# 11. Concordancia global
# ------------------------------------------------------------

pearson <- cor(
  tabla$logFC_GSE108022,
  tabla$logFC_GSE28511_mediana,
  method = "pearson"
)

spearman <- cor(
  tabla$logFC_GSE108022,
  tabla$logFC_GSE28511_mediana,
  method = "spearman"
)


evaluables <- (
  tabla$logFC_GSE108022 != 0 &
  tabla$logFC_GSE28511_mediana != 0
)

misma_direccion <- (
  sign(
    tabla$logFC_GSE108022[evaluables]
  ) ==
  sign(
    tabla$logFC_GSE28511_mediana[evaluables]
  )
)


cat("\n====================================================\n")
cat("CONCORDANCIA GLOBAL MULTISONDA\n")
cat("MEDIANA DE logFC DE TODAS LAS SONDAS\n")
cat("====================================================\n")

cat("\nPearson:\n")
cat(pearson, "\n")

cat("\nSpearman:\n")
cat(spearman, "\n")

cat("\nGenes evaluables direccionalmente:\n")
cat(sum(evaluables), "\n")

cat("\nGenes con misma dirección:\n")
cat(sum(misma_direccion), "\n")

cat("\nProporción con misma dirección:\n")
cat(mean(misma_direccion), "\n")


# ------------------------------------------------------------
# 12. Genes fuertes definidos SOLO en GSE108022
# ------------------------------------------------------------

tabla$fuerte_GSE108022 <- (
  tabla$FDR_GSE108022 < 0.05 &
  abs(tabla$logFC_GSE108022) >= 1
)

fuertes <- tabla[
  tabla$fuerte_GSE108022,
  ,
  drop = FALSE
]


cat("\n====================================================\n")
cat("GENES MULTISONDA FUERTES EN GSE108022\n")
cat("====================================================\n")

cat("\nCriterio:\n")
cat("FDR < 0.05 y |logFC| >= 1\n")

cat("\nNúmero de genes:\n")
cat(nrow(fuertes), "\n")


# ------------------------------------------------------------
# 13. Concordancia de genes fuertes
# ------------------------------------------------------------

evaluables_fuertes <- (
  fuertes$logFC_GSE108022 != 0 &
  fuertes$logFC_GSE28511_mediana != 0
)

misma_fuertes <- (
  sign(
    fuertes$logFC_GSE108022[evaluables_fuertes]
  ) ==
  sign(
    fuertes$logFC_GSE28511_mediana[evaluables_fuertes]
  )
)


cat("\nGenes fuertes evaluables:\n")
cat(sum(evaluables_fuertes), "\n")

cat("\nGenes fuertes con misma dirección:\n")
cat(sum(misma_fuertes), "\n")

cat("\nProporción con misma dirección:\n")
cat(mean(misma_fuertes), "\n")


if (nrow(fuertes) >= 3) {

  pearson_fuertes <- cor(
    fuertes$logFC_GSE108022,
    fuertes$logFC_GSE28511_mediana,
    method = "pearson"
  )

  spearman_fuertes <- cor(
    fuertes$logFC_GSE108022,
    fuertes$logFC_GSE28511_mediana,
    method = "spearman"
  )

  cat("\nPearson en genes fuertes:\n")
  cat(pearson_fuertes, "\n")

  cat("\nSpearman en genes fuertes:\n")
  cat(spearman_fuertes, "\n")
}


# ------------------------------------------------------------
# 14. Separar concordancia según consistencia de las sondas
#
# Esto es descriptivo.
# No redefine la regla de mediana.
# ------------------------------------------------------------

fuertes$sondas_consistentes <- (
  fuertes$clasificacion_sondas %in%
    c(
      "todas_positivas",
      "todas_negativas"
    )
)

fuertes$misma_direccion_mediana <- (
  sign(fuertes$logFC_GSE108022) ==
  sign(fuertes$logFC_GSE28511_mediana)
)


cat("\n====================================================\n")
cat("SENSIBILIDAD SEGÚN CONSISTENCIA ENTRE SONDAS\n")
cat("====================================================\n")

cat("\nGenes fuertes con todas las sondas consistentes:\n")
cat(
  sum(fuertes$sondas_consistentes),
  "\n"
)

cat("\nGenes fuertes con direcciones mixtas:\n")
cat(
  sum(
    fuertes$clasificacion_sondas ==
      "direcciones_mixtas"
  ),
  "\n"
)


if (sum(fuertes$sondas_consistentes) > 0) {

  cat(
    "\nConcordancia con GSE108022 entre genes",
    "con todas las sondas consistentes:\n"
  )

  cat(
    mean(
      fuertes$misma_direccion_mediana[
        fuertes$sondas_consistentes
      ]
    ),
    "\n"
  )
}


if (
  sum(
    fuertes$clasificacion_sondas ==
      "direcciones_mixtas"
  ) > 0
) {

  idx_mixtos <- (
    fuertes$clasificacion_sondas ==
      "direcciones_mixtas"
  )

  cat(
    "\nConcordancia con GSE108022 entre genes",
    "con sondas de direcciones mixtas:\n"
  )

  cat(
    mean(
      fuertes$misma_direccion_mediana[
        idx_mixtos
      ]
    ),
    "\n"
  )
}


# ------------------------------------------------------------
# 15. Guardar variables finales
# ------------------------------------------------------------

tabla$misma_direccion_mediana <- (
  sign(tabla$logFC_GSE108022) ==
  sign(tabla$logFC_GSE28511_mediana)
)


# ------------------------------------------------------------
# 16. Controles finales
# ------------------------------------------------------------

if (nrow(tabla) != 5125) {
  stop(
    "ERROR: número final incorrecto de genes."
  )
}

if (anyDuplicated(tabla$gene)) {
  stop(
    "ERROR: genes duplicados en tabla final."
  )
}


cat("\n====================================================\n")
cat("VERIFICACIÓN FINAL\n")
cat("====================================================\n")

cat("\nGenes multisondeados integrados:\n")
cat(nrow(tabla), "\n")

cat("\nGenes duplicados:\n")
cat(sum(duplicated(tabla$gene)), "\n")

cat(
  "\nVERIFICACIÓN: integración completada",
  "sin pérdida ni duplicación de genes.\n"
)


# ------------------------------------------------------------
# 17. Guardar resultados
# ------------------------------------------------------------

write.csv(
  tabla,
  "sensibilidad_multisonda_GSE108022_GSE28511.csv",
  row.names = FALSE
)

saveRDS(
  tabla,
  "sensibilidad_multisonda_GSE108022_GSE28511.rds"
)


cat("\nArchivos generados:\n")
cat(
  "- sensibilidad_multisonda_GSE108022_GSE28511.csv\n"
)
cat(
  "- sensibilidad_multisonda_GSE108022_GSE28511.rds\n"
)


# ------------------------------------------------------------
# 18. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("SENSIBILIDAD MULTISONDA COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nEste análisis resume cada gen mediante",
  "la mediana de los logFC de todas sus sondas.\n"
)

cat(
  "\nNo se seleccionó ninguna sonda utilizando",
  "los resultados de GSE108022.\n"
)

cat(
  "\nNo se asignó significación estadística génica",
  "a la mediana de las sondas.\n"
)