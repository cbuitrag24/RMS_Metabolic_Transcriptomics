
# ============================================================
# 36_integrar_Hallmarks_replicados_FN_FOXO1FP.R
#
# Objetivo:
# Evaluar la consistencia en FN y FOXO1-FP de los
# 13 Hallmarks previamente replicados entre
# GSE108022 y GSE28511.
#
# No recalcula fgsea.
# No modifica FDR ni NES.
# No selecciona Hallmarks nuevos por los subgrupos.
# ============================================================

# 1. Archivos de entrada

archivo_replicados <- file.path(
  "validacion_GSE28511",
  "Hallmark_replicados_GSE108022_GSE28511.csv"
)

archivo_FN <-
  "Hallmark_fgsea_GSE108022_FN_vs_Normal_universo_comparable.rds"

archivo_FP <-
  "Hallmark_fgsea_GSE108022_FOXO1FP_vs_Normal_universo_comparable.rds"

archivos <- c(
  archivo_replicados,
  archivo_FN,
  archivo_FP
)

if (!all(file.exists(archivos))) {
  print(archivos[!file.exists(archivos)])
  stop("ERROR: faltan archivos de entrada.")
}

# 2. Cargar resultados

replicados <- read.csv(
  archivo_replicados,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

FN <- as.data.frame(readRDS(archivo_FN))
FP <- as.data.frame(readRDS(archivo_FP))

cat("\nDimensiones de los archivos:\n")
print(dim(replicados))
print(dim(FN))
print(dim(FP))

cat("\nColumnas del archivo de replicados:\n")
print(names(replicados))

# 3. Auditar identificadores

if (!("pathway" %in% names(replicados))) {
  stop("ERROR: falta pathway en replicados.")
}

necesarias <- c("pathway", "NES", "padj")

if (
  !all(necesarias %in% names(FN)) ||
  !all(necesarias %in% names(FP))
) {
  stop("ERROR: faltan columnas en resultados fgsea.")
}

if (
  anyDuplicated(replicados$pathway) ||
  anyDuplicated(FN$pathway) ||
  anyDuplicated(FP$pathway)
) {
  stop("ERROR: hay Hallmarks duplicados.")
}

if (nrow(replicados) != 13) {
  stop(
    "ERROR: se esperaban 13 Hallmarks replicados."
  )
}

if (nrow(FN) != 50 || nrow(FP) != 50) {
  stop("ERROR: se esperaban 50 Hallmarks por subgrupo.")
}

if (
  !all(replicados$pathway %in% FN$pathway) ||
  !all(replicados$pathway %in% FP$pathway)
) {
  stop(
    "ERROR: algún Hallmark replicado falta en los subgrupos."
  )
}

# 4. Identificar columnas NES de los estudios originales
#
# Se inspeccionan nombres para evitar suponer
# cómo fueron denominados en el archivo previo.

columnas_NES <- grep(
  "^NES",
  names(replicados),
  value = TRUE
)

cat("\nColumnas NES encontradas:\n")
print(columnas_NES)

if (length(columnas_NES) != 2) {
  stop(
    "ERROR: se esperaban exactamente dos columnas NES
     en la tabla de replicación independiente."
  )
}

# Ambas direcciones deben coincidir, porque los
# Hallmarks ya habían cumplido el criterio de replicación.

NES_original_1 <- replicados[[columnas_NES[1]]]
NES_original_2 <- replicados[[columnas_NES[2]]]

if (
  anyNA(NES_original_1) ||
  anyNA(NES_original_2) ||
  !is.numeric(NES_original_1) ||
  !is.numeric(NES_original_2)
) {
  stop("ERROR: NES originales inválidos.")
}

if (
  any(sign(NES_original_1) != sign(NES_original_2)) ||
  any(NES_original_1 == 0)
) {
  stop(
    "ERROR: los Hallmarks originales no tienen
     dirección concordante."
  )
}

direccion_original <- sign(NES_original_1)

# 5. Recuperar resultados de subgrupos
# preservando el orden de los 13 Hallmarks

idx_FN <- match(
  replicados$pathway,
  FN$pathway
)

idx_FP <- match(
  replicados$pathway,
  FP$pathway
)

FN_13 <- FN[idx_FN, , drop = FALSE]
FP_13 <- FP[idx_FP, , drop = FALSE]

if (
  !identical(
    as.character(FN_13$pathway),
    as.character(replicados$pathway)
  ) ||
  !identical(
    as.character(FP_13$pathway),
    as.character(replicados$pathway)
  )
) {
  stop("ERROR: falló la alineación por pathway.")
}

# 6. Construir tabla integradora

integracion <- data.frame(
  pathway = replicados$pathway,

  NES_replicacion_1 = NES_original_1,
  NES_replicacion_2 = NES_original_2,

  NES_FN = FN_13$NES,
  FDR_FN = FN_13$padj,

  NES_FOXO1FP = FP_13$NES,
  FDR_FOXO1FP = FP_13$padj,

  stringsAsFactors = FALSE
)

# 7. Aplicar criterio previamente establecido

integracion$significativo_FN <-
  !is.na(integracion$FDR_FN) &
  integracion$FDR_FN < 0.05

integracion$significativo_FOXO1FP <-
  !is.na(integracion$FDR_FOXO1FP) &
  integracion$FDR_FOXO1FP < 0.05

integracion$direccion_FN_correcta <-
  sign(integracion$NES_FN) == direccion_original

integracion$direccion_FOXO1FP_correcta <-
  sign(integracion$NES_FOXO1FP) == direccion_original

integracion$consistente_ambos <-
  integracion$significativo_FN &
  integracion$significativo_FOXO1FP &
  integracion$direccion_FN_correcta &
  integracion$direccion_FOXO1FP_correcta

integracion$eje_original <- ifelse(
  direccion_original > 0,
  "RMS",
  "Musculo_normal"
)

# 8. Auditoría

cat("\n============================================\n")
cat("CONSISTENCIA DE LOS 13 HALLMARKS REPLICADOS\n")
cat("============================================\n")

cat("\nHallmarks previamente replicados:\n")
cat(nrow(integracion), "\n")

cat("\nSignificativos en FN:\n")
cat(sum(integracion$significativo_FN), "\n")

cat("\nSignificativos en FOXO1-FP:\n")
cat(sum(integracion$significativo_FOXO1FP), "\n")

cat("\nConsistentes en ambos subgrupos:\n")
cat(sum(integracion$consistente_ambos), "\n")

cat("\nConsistentes hacia RMS:\n")
cat(sum(
  integracion$consistente_ambos &
  integracion$eje_original == "RMS"
), "\n")

cat("\nConsistentes hacia músculo normal:\n")
cat(sum(
  integracion$consistente_ambos &
  integracion$eje_original == "Musculo_normal"
), "\n")

cat("\nTabla completa:\n")
print(integracion, row.names = FALSE)

# 9. Guardar

write.csv(
  integracion,
  "Hallmark_13_replicados_consistencia_FN_FOXO1FP.csv",
  row.names = FALSE
)

saveRDS(
  integracion,
  "Hallmark_13_replicados_consistencia_FN_FOXO1FP.rds"
)

cat("\nArchivos generados:\n")
cat(
  "- Hallmark_13_replicados_consistencia_FN_FOXO1FP.csv\n"
)
cat(
  "- Hallmark_13_replicados_consistencia_FN_FOXO1FP.rds\n"
)

cat("\nANÁLISIS DE CONSISTENCIA COMPLETADO\n")