# ============================================================
# 30_preparar_DE_FN_FOXO1FP_vs_Normal.R
# Preparación del análisis:
# Normal vs FN vs FOXO1-FP
# restringido a GPL16791
# GSE108022
#
# Este script prepara, filtra y normaliza los datos.
# TODAVÍA NO realiza expresión diferencial.
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar DGEList original sin filtrar
#    y diseño de muestras
# ------------------------------------------------------------

dge0 <- readRDS("GSE108022_edgeR_sin_filtrar.rds")

diseno <- read.csv(
  "DE_diseno_muestras.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Datos cargados.\n")


# ------------------------------------------------------------
# 2. Seleccionar GPL16791
# ------------------------------------------------------------

d <- diseno[
  diseno$platform == "GPL16791",
  ,
  drop = FALSE
]

cat("\nMuestras GPL16791 antes de selección específica:\n")
print(nrow(d))


# ------------------------------------------------------------
# 3. Definir grupos biológicos
#
# Normal:
#   Normal
#
# FN:
#   Fusion_negative
#
# FOXO1_FP:
#   PAX3-FOXO1 o PAX7-FOXO1
#
# PAX3-NCOA1 queda fuera de este análisis específico.
# ------------------------------------------------------------

d$grupo <- NA_character_

d$grupo[
  d$fusion_status == "Normal"
] <- "Normal"

d$grupo[
  d$fusion_status == "Fusion_negative"
] <- "FN"

d$grupo[
  d$fusion_status %in% c(
    "PAX3-FOXO1",
    "PAX7-FOXO1"
  )
] <- "FOXO1_FP"


# ------------------------------------------------------------
# 4. Mostrar muestras excluidas
# ------------------------------------------------------------

excluidas <- d[
  is.na(d$grupo),
  c("samplename", "fusion_status"),
  drop = FALSE
]

cat("\nMuestras excluidas de este análisis específico:\n")
print(excluidas, row.names = FALSE)


# ------------------------------------------------------------
# 5. Conservar sólo los tres grupos definidos
# ------------------------------------------------------------

d <- d[
  !is.na(d$grupo),
  ,
  drop = FALSE
]

d$grupo <- factor(
  d$grupo,
  levels = c(
    "Normal",
    "FN",
    "FOXO1_FP"
  )
)

cat("\nNúmero de muestras por grupo:\n")
print(table(d$grupo))

cat("\nTabla grupo x fusion_status:\n")
print(
  table(
    grupo = d$grupo,
    fusion_status = d$fusion_status
  )
)

if (nrow(d) != 68) {
  stop("ERROR: se esperaban 68 muestras.")
}

if (sum(d$grupo == "Normal") != 5) {
  stop("ERROR: se esperaban 5 muestras normales.")
}

if (sum(d$grupo == "FN") != 44) {
  stop("ERROR: se esperaban 44 muestras FN.")
}

if (sum(d$grupo == "FOXO1_FP") != 19) {
  stop("ERROR: se esperaban 19 muestras FOXO1-FP.")
}


# ------------------------------------------------------------
# 6. Verificar correspondencia con counts
# ------------------------------------------------------------

if (!all(d$samplename %in% colnames(dge0$counts))) {
  stop("ERROR: alguna muestra no está presente en los counts.")
}

dge <- dge0[
  ,
  d$samplename,
  keep.lib.sizes = FALSE
]

if (!identical(
  colnames(dge$counts),
  d$samplename
)) {
  stop("ERROR: el orden de las muestras no coincide.")
}

cat("\nVERIFICACIÓN: muestras y counts correctamente alineados.\n")


# ------------------------------------------------------------
# 7. Construir matriz de diseño
#
# Sin intercepto para facilitar posteriormente
# contrastes explícitos:
#
# FN - Normal
# FOXO1_FP - Normal
# ------------------------------------------------------------

design <- model.matrix(
  ~ 0 + grupo,
  data = d
)

colnames(design) <- levels(d$grupo)
rownames(design) <- d$samplename

cat("\nColumnas de la matriz de diseño:\n")
print(colnames(design))

cat("\nDimensiones de la matriz de diseño:\n")
print(dim(design))

cat("\nRango de la matriz de diseño:\n")
print(qr(design)$rank)

if (qr(design)$rank != ncol(design)) {
  stop("ERROR: la matriz de diseño no tiene rango completo.")
}


# ------------------------------------------------------------
# 8. Filtrado específico con filterByExpr
# ------------------------------------------------------------

n_antes <- nrow(dge)

keep <- filterByExpr(
  dge,
  design = design
)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

n_despues <- nrow(dge)

cat("\nGenes antes del filtrado:\n")
print(n_antes)

cat("\nGenes conservados:\n")
print(n_despues)

cat("\nGenes eliminados:\n")
print(n_antes - n_despues)

cat("\nPorcentaje conservado:\n")
print(100 * n_despues / n_antes)


# ------------------------------------------------------------
# 9. Normalización TMM
# ------------------------------------------------------------

dge <- normLibSizes(
  dge,
  method = "TMM"
)

cat("\nResumen de factores TMM:\n")
print(summary(dge$samples$norm.factors))

if (anyNA(dge$samples$norm.factors) ||
    any(dge$samples$norm.factors <= 0)) {
  stop("ERROR: factores TMM inválidos.")
}


# ------------------------------------------------------------
# 10. Control TMM por grupo
# ------------------------------------------------------------

control_tmm <- data.frame(
  samplename = colnames(dge$counts),
  grupo = d$grupo,
  fusion_status = d$fusion_status,
  lib.size = dge$samples$lib.size,
  norm.factor = dge$samples$norm.factors,
  effective.lib.size =
    dge$samples$lib.size *
    dge$samples$norm.factors,
  stringsAsFactors = FALSE
)

cat("\nFactores TMM por grupo:\n")

print(
  aggregate(
    norm.factor ~ grupo,
    data = control_tmm,
    FUN = function(x) {
      c(
        n = length(x),
        min = min(x),
        median = median(x),
        mean = mean(x),
        max = max(x)
      )
    }
  )
)


# ------------------------------------------------------------
# 11. Guardar objetos preparados
# ------------------------------------------------------------

saveRDS(
  dge,
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_filtrado_TMM.rds"
)

saveRDS(
  design,
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_design.rds"
)

write.csv(
  d,
  "GPL16791_FN_FOXO1FP_Normal_diseno_muestras.csv",
  row.names = FALSE
)

write.csv(
  control_tmm,
  "GPL16791_FN_FOXO1FP_Normal_control_TMM.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("- GSE108022_GPL16791_FN_FOXO1FP_Normal_filtrado_TMM.rds\n")
cat("- GSE108022_GPL16791_FN_FOXO1FP_Normal_design.rds\n")
cat("- GPL16791_FN_FOXO1FP_Normal_diseno_muestras.csv\n")
cat("- GPL16791_FN_FOXO1FP_Normal_control_TMM.csv\n")

cat(
  "\nRESULTADO: datos preparados para FN vs Normal y FOXO1-FP vs Normal; todavía no se realizó DE.\n"
)