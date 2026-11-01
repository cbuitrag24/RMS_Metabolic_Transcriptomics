# ============================================================
# 04_inspeccionar_anotacion_GPL6947.R
#
# Objetivo:
# - consultar la anotación oficial de la plataforma GPL6947
# - inspeccionar sus columnas
# - comprobar correspondencia con los ID_REF de GSE28511
#
# NO modifica la matriz de expresión.
# NO realiza expresión diferencial.
# ============================================================


# ------------------------------------------------------------
# 1. Verificar GEOquery
# ------------------------------------------------------------

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop(
    paste(
      "ERROR: el paquete GEOquery no está instalado.",
      "No se realizará ninguna instalación automática."
    )
  )
}

library(GEOquery)


# ------------------------------------------------------------
# 2. Descargar la plataforma GPL6947 desde GEO
# ------------------------------------------------------------

cat("Descargando anotación oficial de GPL6947 desde GEO...\n")

gpl <- getGEO(
  "GPL6947",
  AnnotGPL = FALSE
)

cat("\nPlataforma obtenida:\n")
cat(Meta(gpl)$title, "\n")


# ------------------------------------------------------------
# 3. Extraer tabla de anotación
# ------------------------------------------------------------

anotacion <- Table(gpl)

cat("\n====================================================\n")
cat("DIMENSIONES DE LA ANOTACIÓN\n")
cat("====================================================\n")

print(dim(anotacion))


# ------------------------------------------------------------
# 4. Inspeccionar nombres de columnas
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("COLUMNAS DISPONIBLES\n")
cat("====================================================\n")

for (i in seq_along(colnames(anotacion))) {
  cat(
    i,
    ":",
    colnames(anotacion)[i],
    "\n"
  )
}


# ------------------------------------------------------------
# 5. Mostrar primeras filas sin interpretar
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("PRIMERAS 5 FILAS DE LA ANOTACIÓN\n")
cat("====================================================\n")

print(
  head(anotacion, 5)
)


# ------------------------------------------------------------
# 6. Leer IDs de nuestra matriz log2
# ------------------------------------------------------------

archivo_expr <- "GSE28511_expression_quantileNormalized_log2.rds"

if (!file.exists(archivo_expr)) {
  stop(
    paste(
      "ERROR: no se encontró",
      archivo_expr
    )
  )
}

expr <- readRDS(archivo_expr)

ids_expr <- rownames(expr)

if (is.null(ids_expr)) {
  stop(
    "ERROR: la matriz de expresión no tiene rownames."
  )
}

cat("\n====================================================\n")
cat("CORRESPONDENCIA CON GSE28511\n")
cat("====================================================\n")

cat("\nIDs en matriz de expresión:\n")
cat(length(ids_expr), "\n")


# ------------------------------------------------------------
# 7. Identificar columna ID de la plataforma
# ------------------------------------------------------------

if (!("ID" %in% colnames(anotacion))) {
  stop(
    "ERROR: la tabla GPL6947 no contiene una columna llamada ID."
  )
}

ids_gpl <- as.character(anotacion$ID)

cat("\nIDs en GPL6947:\n")
cat(length(ids_gpl), "\n")

cat("\nIDs GPL6947 duplicados:\n")
cat(
  length(ids_gpl) - length(unique(ids_gpl)),
  "\n"
)


# ------------------------------------------------------------
# 8. Cuantificar coincidencia
# ------------------------------------------------------------

coinciden <- ids_expr %in% ids_gpl

cat("\nIDs de expresión presentes en GPL6947:\n")
cat(sum(coinciden), "\n")

cat("\nIDs de expresión NO presentes en GPL6947:\n")
cat(sum(!coinciden), "\n")

cat("\nProporción de IDs de expresión anotados:\n")
cat(mean(coinciden), "\n")


# ------------------------------------------------------------
# 9. Mostrar algunos IDs no encontrados, si existen
# ------------------------------------------------------------

if (any(!coinciden)) {

  cat("\nPrimeros ID_REF no encontrados:\n")

  print(
    head(
      ids_expr[!coinciden],
      20
    )
  )
}


# ------------------------------------------------------------
# 10. Guardar copia local de la tabla oficial
# ------------------------------------------------------------

write.csv(
  anotacion,
  "GPL6947_anotacion_GEO.csv",
  row.names = FALSE
)

cat("\nArchivo generado:\n")
cat("- GPL6947_anotacion_GEO.csv\n")

cat("\n====================================================\n")
cat("INSPECCIÓN DE ANOTACIÓN COMPLETADA\n")
cat("====================================================\n")