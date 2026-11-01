# ============================================================
# 29_preparar_contrastes_FN_FP_normal.R
# Verificación de grupos para análisis:
# Normal vs Fusion_negative vs Fusion_positive
# restringido a GPL16791
# GSE108022
#
# IMPORTANTE:
# Este script NO realiza todavía expresión diferencial.
# Sólo define y audita la composición de los grupos.
# ============================================================


# ------------------------------------------------------------
# 1. Cargar diseño previamente verificado
# ------------------------------------------------------------

diseno <- read.csv(
  "DE_diseno_muestras.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Diseño cargado.\n")

cat("\nDimensiones:\n")
print(dim(diseno))

cat("\nColumnas disponibles:\n")
print(colnames(diseno))


# ------------------------------------------------------------
# 2. Verificar columnas necesarias
# ------------------------------------------------------------

columnas_necesarias <- c(
  "samplename",
  "fusion_status",
  "platform"
)

if (!all(columnas_necesarias %in% colnames(diseno))) {
  stop(
    paste(
      "ERROR: faltan columnas necesarias.",
      "Columnas encontradas:",
      paste(colnames(diseno), collapse = ", ")
    )
  )
}

if (anyNA(diseno$samplename) ||
    any(diseno$samplename == "")) {
  stop("ERROR: existen nombres de muestra vacíos o NA.")
}

if (anyDuplicated(diseno$samplename) > 0) {
  stop("ERROR: existen nombres de muestra duplicados.")
}


# ------------------------------------------------------------
# 3. Restringir a GPL16791
# ------------------------------------------------------------

d <- diseno[
  diseno$platform == "GPL16791",
  ,
  drop = FALSE
]

cat("\nMuestras GPL16791:\n")
print(nrow(d))

cat("\nFusion_status dentro de GPL16791:\n")
print(table(d$fusion_status, useNA = "ifany"))


# ------------------------------------------------------------
# 4. Definir categorías de forma explícita
#
# Normal:
#   fusion_status == "Normal"
#
# FN:
#   fusion_status == "Fusion_negative"
#
# FP:
#   cualquier categoría distinta de Normal y Fusion_negative
#
# Todavía NO se utilizará esta definición para DE.
# Primero imprimiremos exactamente qué categorías entran.
# ------------------------------------------------------------

d$grupo <- NA_character_

d$grupo[d$fusion_status == "Normal"] <- "Normal"

d$grupo[
  d$fusion_status == "Fusion_negative"
] <- "FN"

d$grupo[
  !d$fusion_status %in% c(
    "Normal",
    "Fusion_negative"
  )
] <- "FP"


# ------------------------------------------------------------
# 5. Verificar que ninguna muestra quede sin clasificar
# ------------------------------------------------------------

if (anyNA(d$grupo)) {

  cat("\nEstados no clasificados:\n")

  print(
    unique(
      d$fusion_status[is.na(d$grupo)]
    )
  )

  stop("ERROR: existen muestras sin grupo asignado.")
}


# ------------------------------------------------------------
# 6. Mostrar composición exacta
# ------------------------------------------------------------

cat("\nNúmero de muestras por grupo propuesto:\n")
print(table(d$grupo))

cat("\nTabla grupo x fusion_status:\n")
print(
  table(
    grupo = d$grupo,
    fusion_status = d$fusion_status
  )
)


# ------------------------------------------------------------
# 7. Listar todas las muestras FP
# para comprobar qué fusiones contiene el grupo
# ------------------------------------------------------------

fp <- d[
  d$grupo == "FP",
  c("samplename", "fusion_status"),
  drop = FALSE
]

fp <- fp[
  order(fp$fusion_status, fp$samplename),
  ,
  drop = FALSE
]

cat("\nMuestras propuestas como Fusion_positive:\n")
print(fp, row.names = FALSE)


# ------------------------------------------------------------
# 8. Verificación de totales
# ------------------------------------------------------------

n_normal <- sum(d$grupo == "Normal")
n_fn     <- sum(d$grupo == "FN")
n_fp     <- sum(d$grupo == "FP")

cat("\nResumen final propuesto:\n")
cat("Normal:", n_normal, "\n")
cat("FN:", n_fn, "\n")
cat("FP:", n_fp, "\n")
cat("Total:", n_normal + n_fn + n_fp, "\n")

if ((n_normal + n_fn + n_fp) != nrow(d)) {
  stop("ERROR: el total de grupos no coincide con GPL16791.")
}

cat(
  "\nRESULTADO: composición de grupos preparada para revisión; todavía no se realizó DE.\n"
)