# Export archived results supporting manuscript Figures 4 and 5.
# Uses the existing saved R objects; no correlations or models are refitted.
# Run from the repository root after downloading the repository.

archivo <- file.path("results", "RMS_Figuras_4_5_objetos_originales.RData")
if (!file.exists(archivo)) stop("Archivo no encontrado: ", archivo)
original <- new.env(parent = globalenv())
load(archivo, envir = original)

salida <- file.path("results", "exportacion_objetos_originales")
dir.create(salida, recursive = TRUE, showWarnings = FALSE)

# Export only already-computed tables when present.
tablas <- c("resultado_PC1_miogenesis", "tabla_cor_CKM", "correlaciones_CKM")
for (nombre in tablas) {
  if (!exists(nombre, envir = original, inherits = FALSE)) next
  objeto <- get(nombre, envir = original, inherits = FALSE)
  if (is.data.frame(objeto) || is.matrix(objeto)) {
    write.csv(as.data.frame(objeto), file.path(salida, paste0(nombre, ".csv")), row.names = TRUE)
  } else {
    saveRDS(objeto, file.path(salida, paste0(nombre, ".rds")))
  }
}

# Preserve original fitted model summaries, without refitting models.
modelos <- c("modelo_final", "modelo_sensibilidad")
for (nombre in modelos) {
  if (!exists(nombre, envir = original, inherits = FALSE)) next
  objeto <- get(nombre, envir = original, inherits = FALSE)
  capture.output(summary(objeto), file = file.path(salida, paste0(nombre, "_summary.txt")))
}

# Save the original differentiation score and leave-one-out results as objects.
for (nombre in c("score_diferenciacion", "rho_LOO_CKM")) {
  if (exists(nombre, envir = original, inherits = FALSE)) {
    saveRDS(get(nombre, envir = original, inherits = FALSE),
            file.path(salida, paste0(nombre, ".rds")))
  }
}

message("Exportación terminada: ", salida)
# These exports are supporting data, not scripts for regenerating figure artwork.
