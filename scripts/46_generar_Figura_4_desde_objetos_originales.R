# Recreate the three-panel layout of the investigator-supplied Figure 4.
# This script draws the existing observations; it does not recalculate correlations.
# Run from repository root. Requires the archived RData in results/.
archivo <- file.path("results", "RMS_Figuras_4_5_objetos_originales.RData")
if (!file.exists(archivo)) stop("Missing original RData: ", archivo)
o <- new.env(parent = globalenv())
load(archivo, envir = o)
necesarios <- c("pca_metabolico", "expresion_rms_22", "score_diferenciacion")
faltan <- necesarios[!vapply(necesarios, exists, logical(1), envir = o, inherits = FALSE)]
if (length(faltan)) stop("Missing archived objects: ", paste(faltan, collapse = ", "))

# Model formula in original archived summary defines the three variables.
pc1 <- as.numeric(o$pca_metabolico$x[, 1])
ckm <- as.numeric(o$expresion_rms_22["CKM", ])
dif <- as.numeric(o$score_diferenciacion)
if (!(length(pc1) == length(ckm) && length(ckm) == length(dif)))
  stop("Original vectors have different lengths; cannot plot without alignment.")
if (!all(is.finite(pc1) & is.finite(ckm) & is.finite(dif)))
  stop("Missing/non-finite observations; inspect original objects.")
# Preserve source order; if names are supplied, require matching sample identities.
ids <- list(rownames(o$pca_metabolico$x), colnames(o$expresion_rms_22), names(o$score_diferenciacion))
ids <- Filter(function(x) !is.null(x), ids)
if (length(ids) > 1L && !all(vapply(ids[-1L], identical, logical(1), ids[[1L]])))
  stop("Sample identifiers differ across archived objects; do not plot misaligned data.")

dir.create("figures", showWarnings = FALSE)
salida <- file.path("figures", "Figura_4_Metabolismo_Diferenciacion_RMS_recreada.pdf")
pdf(salida, width = 12.5, height = 5.2, useDingbats = FALSE)
op <- par(mfrow = c(1, 3), mar = c(4.6, 4.5, 3.2, 0.8), mgp = c(2.6, 0.8, 0))
on.exit <- function() { par(op); dev.off() }
azul <- "#2877B6"
plot(ckm, pc1, pch = 16, cex = 0.55, col = azul,
     xlab = "CKM expression (log2 CPM)", ylab = "Metabolic PC1",
     main = "Spearman rho = -0.413; p = 0.00069")
mtext("A", side = 3, adj = -0.09, line = 1, font = 2)
plot(ckm, dif, pch = 16, cex = 0.55, col = azul,
     xlab = "CKM expression (log2 CPM)", ylab = "Differentiation signature score",
     main = "Spearman rho = 0.539; p = 4.31e-06")
mtext("B", side = 3, adj = -0.09, line = 1, font = 2)
plot(dif, pc1, pch = 16, cex = 0.55, col = azul,
     xlab = "Differentiation signature score", ylab = "Metabolic PC1",
     main = "Spearman rho = -0.139; p = 0.2733")
mtext("C", side = 3, adj = -0.09, line = 1, font = 2)
dev.off()
cat("Created:", salida, "\n")
# Correlation labels are transcribed from the original supplied Figure 4 PDF.
# Plot is a layout reconstruction, not a verified pixel-identical reproduction.
