
# ============================================================
# SCRIPT 41
# Integracion de genes leading-edge en tres contrastes:
#
# 1. FN vs Normal
# 2. FOXO1-FP vs Normal
# 3. FOXO1-FP vs FN
#
# No recalcula edgeR ni fgsea.
# No filtra genes por significacion.
# ============================================================

archivos <- c(
  FN = "GPL16791_DE_FN_vs_Normal_todos.csv",
  FP = "GPL16791_DE_FOXO1FP_vs_Normal_todos.csv",
  Directo = "GPL16791_DE_FOXO1FP_vs_FN_todos.csv",
  LE_TNFA = "LeadingEdge_TNFA_NFKB_FOXO1FP_vs_FN.csv",
  LE_OXPHOS = "LeadingEdge_OXPHOS_FOXO1FP_vs_FN.csv"
)

faltantes <- archivos[!file.exists(archivos)]

if (length(faltantes) > 0) {
  print(faltantes)
  stop("ERROR: faltan archivos de entrada.")
}

# ------------------------------------------------------------
# 1. Leer tablas
# ------------------------------------------------------------

FN <- read.csv(archivos["FN"], stringsAsFactors = FALSE)
FP <- read.csv(archivos["FP"], stringsAsFactors = FALSE)
DIR <- read.csv(archivos["Directo"], stringsAsFactors = FALSE)

LE_TNFA <- read.csv(
  archivos["LE_TNFA"],
  stringsAsFactors = FALSE
)

LE_OXPHOS <- read.csv(
  archivos["LE_OXPHOS"],
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Control de calidad
# ------------------------------------------------------------

tablas <- list(
  FN = FN,
  FP = FP,
  DIR = DIR
)

for (nombre in names(tablas)) {

  x <- tablas[[nombre]]

  if (!all(c("gene", "logFC", "FDR") %in% names(x))) {
    stop(paste("ERROR: columnas faltantes en", nombre))
  }

  if (anyNA(x$gene) ||
      anyDuplicated(x$gene) ||
      any(!nzchar(trimws(x$gene)))) {
    stop(paste("ERROR: genes invalidos en", nombre))
  }

  if (!is.numeric(x$logFC) ||
      !is.numeric(x$FDR) ||
      any(!is.finite(x$logFC)) ||
      any(!is.finite(x$FDR))) {
    stop(paste("ERROR: estadisticos invalidos en", nombre))
  }
}

# ------------------------------------------------------------
# 3. Funcion de integracion
# ------------------------------------------------------------

integrar <- function(leading_edge, nombre_programa) {

  genes <- as.character(leading_edge$gene)

  if (anyNA(genes) ||
      anyDuplicated(genes) ||
      any(!nzchar(trimws(genes)))) {
    stop("ERROR: leading-edge invalido.")
  }

  i_FN <- match(genes, FN$gene)
  i_FP <- match(genes, FP$gene)
  i_DIR <- match(genes, DIR$gene)

  if (anyNA(i_FN) || anyNA(i_FP) || anyNA(i_DIR)) {
    stop(paste(
      "ERROR: genes ausentes en algun contraste:",
      nombre_programa
    ))
  }

  salida <- data.frame(
    gene = genes,

    logFC_FN_vs_Normal = FN$logFC[i_FN],
    FDR_FN_vs_Normal = FN$FDR[i_FN],

    logFC_FOXO1FP_vs_Normal = FP$logFC[i_FP],
    FDR_FOXO1FP_vs_Normal = FP$FDR[i_FP],

    logFC_FOXO1FP_vs_FN = DIR$logFC[i_DIR],
    FDR_FOXO1FP_vs_FN = DIR$FDR[i_DIR],

    stringsAsFactors = FALSE
  )

  # Verificar correspondencia con leading-edge original

  if (any(abs(
    salida$logFC_FOXO1FP_vs_FN -
      leading_edge$logFC
  ) > 1e-10)) {
    stop("ERROR: discrepancia con leading-edge original.")
  }

  salida$significativo_FN_vs_Normal <-
    salida$FDR_FN_vs_Normal < 0.05

  salida$significativo_FOXO1FP_vs_Normal <-
    salida$FDR_FOXO1FP_vs_Normal < 0.05

  salida$significativo_directo <-
    salida$FDR_FOXO1FP_vs_FN < 0.05

  salida$direccion_FN_vs_Normal <- ifelse(
    salida$logFC_FN_vs_Normal > 0,
    "Mayor_FN",
    ifelse(
      salida$logFC_FN_vs_Normal < 0,
      "Mayor_Normal",
      "Igual"
    )
  )

  salida$direccion_FOXO1FP_vs_Normal <- ifelse(
    salida$logFC_FOXO1FP_vs_Normal > 0,
    "Mayor_FOXO1FP",
    ifelse(
      salida$logFC_FOXO1FP_vs_Normal < 0,
      "Mayor_Normal",
      "Igual"
    )
  )

  salida$programa <- nombre_programa

  salida
}

# ------------------------------------------------------------
# 4. Integrar los dos programas
# ------------------------------------------------------------

TNFA <- integrar(
  LE_TNFA,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

OXPHOS <- integrar(
  LE_OXPHOS,
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

# ------------------------------------------------------------
# 5. Guardar resultados
# ------------------------------------------------------------

write.csv(
  TNFA,
  "Integracion_LE_TNFA_NFKB_tres_contrastes.csv",
  row.names = FALSE
)

write.csv(
  OXPHOS,
  "Integracion_LE_OXPHOS_tres_contrastes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Resumen descriptivo
# ------------------------------------------------------------

resumir <- function(x) {

  data.frame(
    programa = unique(x$programa),
    genes = nrow(x),

    FDR005_FN_vs_Normal =
      sum(x$significativo_FN_vs_Normal),

    FDR005_FOXO1FP_vs_Normal =
      sum(x$significativo_FOXO1FP_vs_Normal),

    FDR005_directo =
      sum(x$significativo_directo),

    logFC_positivo_FN_vs_Normal =
      sum(x$logFC_FN_vs_Normal > 0),

    logFC_positivo_FOXO1FP_vs_Normal =
      sum(x$logFC_FOXO1FP_vs_Normal > 0)
  )
}

resumen <- rbind(
  resumir(TNFA),
  resumir(OXPHOS)
)

write.csv(
  resumen,
  "Resumen_LE_tres_contrastes.csv",
  row.names = FALSE
)

cat("\nRESUMEN DE INTEGRACION:\n")
print(resumen, row.names = FALSE)

cat("\nArchivos generados correctamente.\n")