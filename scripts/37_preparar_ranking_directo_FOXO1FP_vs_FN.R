
# ============================================================
# 37_preparar_ranking_directo_FOXO1FP_vs_FN.R
#
# OBJETIVO:
# Preparar el ranking del contraste directo:
# FOXO1-FP vs Fusion-negative (FN)
#
# Universo: interseccion previamente auditada de 9427 genes.
#
# Estadistico:
# sign(logFC) * sqrt(F)
#
# Positivo: mayor expresion en FOXO1-FP
# Negativo: mayor expresion en FN
#
# No se filtra por FDR ni por magnitud de logFC.
# No se recalcula el analisis diferencial.
# ============================================================

archivo_DE <- "GPL16791_DE_FOXO1FP_vs_FN_todos.csv"

archivo_universo <-
  "ranking_GSE108022_FN_vs_Normal_universo_comparable.rds"

archivo_FP <-
  "ranking_GSE108022_FOXO1FP_vs_Normal_universo_comparable.rds"

archivos <- c(
  archivo_DE,
  archivo_universo,
  archivo_FP
)

faltantes <- archivos[!file.exists(archivos)]

if (length(faltantes) > 0) {
  print(faltantes)
  stop("ERROR: faltan archivos de entrada.")
}

# ------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------

DE <- read.csv(
  archivo_DE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

FN_normal <- readRDS(archivo_universo)
FP_normal <- readRDS(archivo_FP)

cat("\nDimensiones contraste directo:\n")
print(dim(DE))

cat("\nDimensiones universo FN:\n")
print(dim(FN_normal))

cat("\nDimensiones universo FOXO1-FP:\n")
print(dim(FP_normal))

# ------------------------------------------------------------
# 2. Verificar columnas e identificadores
# ------------------------------------------------------------

columnas <- c(
  "gene", "logFC", "F", "PValue", "FDR"
)

if (!all(columnas %in% names(DE))) {
  stop("ERROR: faltan columnas en la tabla DE.")
}

if (!("gene" %in% names(FN_normal)) ||
    !("gene" %in% names(FP_normal))) {
  stop("ERROR: faltan identificadores en los universos.")
}

DE$gene <- as.character(DE$gene)
FN_normal$gene <- as.character(FN_normal$gene)
FP_normal$gene <- as.character(FP_normal$gene)

if (anyNA(DE$gene) ||
    any(!nzchar(trimws(DE$gene))) ||
    anyDuplicated(DE$gene)) {
  stop("ERROR: identificadores invalidos o duplicados en DE.")
}

if (anyNA(FN_normal$gene) ||
    anyNA(FP_normal$gene) ||
    anyDuplicated(FN_normal$gene) ||
    anyDuplicated(FP_normal$gene)) {
  stop("ERROR: identificadores invalidos en los universos.")
}

if (!identical(FN_normal$gene, FP_normal$gene)) {
  stop("ERROR: FN y FOXO1-FP no tienen identico universo.")
}

universo <- FN_normal$gene

cat("\nGenes del universo comparable:\n")
cat(length(universo), "\n")

# ------------------------------------------------------------
# 3. Verificar estadisticos originales
# ------------------------------------------------------------

for (col in c("logFC", "F", "PValue", "FDR")) {
  if (!is.numeric(DE[[col]]) ||
      any(!is.finite(DE[[col]]))) {
    stop(paste("ERROR: valores invalidos en", col))
  }
}

if (any(DE$F < 0)) {
  stop("ERROR: existen estadisticos F negativos.")
}

if (any(DE$PValue < 0 | DE$PValue > 1) ||
    any(DE$FDR < 0 | DE$FDR > 1)) {
  stop("ERROR: probabilidades fuera del intervalo [0,1].")
}

# ------------------------------------------------------------
# 4. Auditar cobertura
# ------------------------------------------------------------

ausentes <- setdiff(universo, DE$gene)

cat("\nGenes del universo ausentes en DE:\n")
cat(length(ausentes), "\n")

if (length(ausentes) > 0) {
  print(ausentes)
  stop("ERROR: el universo no esta completamente representado.")
}

# ------------------------------------------------------------
# 5. Restringir al universo comparable
# ------------------------------------------------------------

indices <- match(universo, DE$gene)

directo <- DE[
  indices,
  c("gene", "logFC", "F", "PValue", "FDR"),
  drop = FALSE
]

rownames(directo) <- NULL

if (!identical(directo$gene, universo)) {
  stop("ERROR: no se conserva el orden del universo.")
}

# ------------------------------------------------------------
# 6. Construir ranking
# ------------------------------------------------------------

directo$ranking_stat <-
  sign(directo$logFC) * sqrt(directo$F)

if (any(!is.finite(directo$ranking_stat))) {
  stop("ERROR: ranking con valores no finitos.")
}

error_formula <- max(abs(
  directo$ranking_stat -
    sign(directo$logFC) * sqrt(directo$F)
))

cat("\nError maximo de la formula:\n")
cat(error_formula, "\n")

# ------------------------------------------------------------
# 7. Verificar sentido del contraste
#
# IMPORTANTE:
# El signo positivo se interpreta como FOXO1-FP
# solamente si la tabla original fue exportada con
# el contraste FOXO1FP - FN.
# ------------------------------------------------------------

cat("\nPrimeros genes del contraste directo:\n")
print(head(directo, 5))

cat("\nResumen ranking_stat:\n")
print(summary(directo$ranking_stat))

cat("\nGenes con ranking positivo:\n")
cat(sum(directo$ranking_stat > 0), "\n")

cat("\nGenes con ranking negativo:\n")
cat(sum(directo$ranking_stat < 0), "\n")

cat("\nGenes con ranking cero:\n")
cat(sum(directo$ranking_stat == 0), "\n")

# ------------------------------------------------------------
# 8. Auditoria de identidad y recuperacion
# ------------------------------------------------------------

error_logFC <- max(abs(
  directo$logFC - DE$logFC[indices]
))

error_F <- max(abs(
  directo$F - DE$F[indices]
))

cat("\nError recuperacion logFC:\n")
cat(error_logFC, "\n")

cat("\nError recuperacion F:\n")
cat(error_F, "\n")

if (error_logFC != 0 ||
    error_F != 0 ||
    error_formula != 0) {
  stop("ERROR: inconsistencia en los estadisticos.")
}

# ------------------------------------------------------------
# 9. Guardar archivos
# ------------------------------------------------------------

write.csv(
  directo,
  "ranking_GSE108022_FOXO1FP_vs_FN_directo_9427.csv",
  row.names = FALSE
)

saveRDS(
  directo,
  "ranking_GSE108022_FOXO1FP_vs_FN_directo_9427.rds"
)

auditoria <- data.frame(
  medida = c(
    "genes_DE_original",
    "genes_universo",
    "genes_ausentes",
    "genes_ranking_final",
    "ranking_positivos",
    "ranking_negativos",
    "ranking_ceros",
    "error_formula",
    "error_logFC",
    "error_F"
  ),
  valor = c(
    nrow(DE),
    length(universo),
    length(ausentes),
    nrow(directo),
    sum(directo$ranking_stat > 0),
    sum(directo$ranking_stat < 0),
    sum(directo$ranking_stat == 0),
    error_formula,
    error_logFC,
    error_F
  )
)

write.csv(
  auditoria,
  "auditoria_ranking_directo_FOXO1FP_vs_FN.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("- ranking_GSE108022_FOXO1FP_vs_FN_directo_9427.csv\n")
cat("- ranking_GSE108022_FOXO1FP_vs_FN_directo_9427.rds\n")
cat("- auditoria_ranking_directo_FOXO1FP_vs_FN.csv\n")

cat("\n====================================================\n")
cat("RANKING DIRECTO CONSTRUIDO Y AUDITADO\n")
cat("====================================================\n")