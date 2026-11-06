# Recreate the investigator-supplied Figure 5 coefficient plot.
# Uses archived fitted model; does not refit any regression.
# Run from repository root.
archivo <- file.path("results", "RMS_Figuras_4_5_objetos_originales.RData")
if (!file.exists(archivo)) stop("Missing original RData: ", archivo)
o <- new.env(parent = globalenv())
load(archivo, envir = o)
if (!exists("modelo_final", envir = o, inherits = FALSE)) stop("modelo_final not found")
modelo <- o$modelo_final
if (!inherits(modelo, "lm")) stop("Archived modelo_final is not an lm object")
b <- coef(modelo)
# Standard 95% confidence intervals from the existing fit; no model refitting.
ci <- confint.default(modelo, level = 0.95)
terminos <- c("CKM", "diferenciacion", "fusionPAX3-FOXO1",
              "fusionPAX3-NCOA1", "fusionPAX7-FOXO1")
if (!all(terminos %in% names(b))) stop("Archived model terms differ from original Figure 5")
etiquetas <- c("CKM", "Differentiation signature",
               "PAX3-FOXO1 vs fusion-negative",
               "PAX3-NCOA1 vs fusion-negative",
               "PAX7-FOXO1 vs fusion-negative")
y <- rev(seq_along(terminos))
lim <- range(c(ci[terminos, ], 0), finite = TRUE)
dir.create("figures", showWarnings = FALSE)
salida <- file.path("figures", "Figura_5_Modelo_Ajustado_RMS_recreada.pdf")
pdf(salida, width = 10.8, height = 5.2, useDingbats = FALSE)
par(mar = c(4.6, 14.5, 0.9, 1.1), xpd = NA)
plot(NA, xlim = lim + c(-0.2, 0.2), ylim = c(0.3, 5.7), yaxt = "n",
     ylab = "", xlab = "Regression coefficient (95% CI)", bty = "o")
axis(2, at = y, labels = etiquetas, las = 1, tick = FALSE)
abline(v = 0, lty = 3, col = "grey55")
segments(ci[terminos, 1], y, ci[terminos, 2], y, col = "#2877B6", lwd = 1.6)
points(b[terminos], y, pch = 16, col = "#2877B6", cex = 1)
dev.off()
cat("Created:", salida, "\n")
# Uses normal-approximation intervals (confint.default) matching a coefficient/SE plot.
# Confirm exact CI method against the historical figure-generating script if found.
