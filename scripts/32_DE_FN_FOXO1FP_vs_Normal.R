# ============================================================
# 32_DE_FN_FOXO1FP_vs_Normal.R
# Expresión diferencial:
#   1) FN - Normal
#   2) FOXO1-FP - Normal
# GPL16791 - GSE108022
# ============================================================

library(edgeR)
library(limma)

# ------------------------------------------------------------
# 1. Cargar ajuste QL validado
# ------------------------------------------------------------

fit <- readRDS(
  "GSE108022_GPL16791_FN_FOXO1FP_Normal_QLfit.rds"
)

cat("Ajuste QL cargado.\n")

cat("\nColumnas del diseño:\n")
print(colnames(fit$design))

if (!identical(
  colnames(fit$design),
  c("Normal", "FN", "FOXO1_FP")
)) {
  stop("ERROR: columnas inesperadas en la matriz de diseño.")
}

if (is.null(fit$genes) ||
    !"gene" %in% colnames(fit$genes)) {
  stop("ERROR: no se dispone de identidad génica válida.")
}


# ------------------------------------------------------------
# 2. Definir contrastes explícitos
# ------------------------------------------------------------

contrastes <- makeContrasts(
  FN_vs_Normal = FN - Normal,
  FOXO1FP_vs_Normal = FOXO1_FP - Normal,
  levels = fit$design
)

cat("\nMatriz de contrastes:\n")
print(contrastes)


# ------------------------------------------------------------
# 3. Función para ejecutar y verificar cada contraste
# ------------------------------------------------------------

analizar_contraste <- function(
  fit,
  contraste,
  nombre
) {

  cat("\n====================================================\n")
  cat("CONTRASTE:", nombre, "\n")
  cat("====================================================\n")

  qlf <- glmQLFTest(
    fit,
    contrast = contraste
  )

  if (nrow(qlf$table) != nrow(qlf$genes)) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- qlf$table y qlf$genes no coinciden."
      )
    )
  }

  resultado_base <- data.frame(
    gene = qlf$genes$gene,
    qlf$table,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!identical(
    as.character(resultado_base$gene),
    as.character(qlf$genes$gene)
  )) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- se perdió la correspondencia gen-fila."
      )
    )
  }

  resultado_base$FDR <- p.adjust(
    resultado_base$PValue,
    method = "BH"
  )

  if (anyNA(resultado_base$logFC) ||
      anyNA(resultado_base$PValue) ||
      anyNA(resultado_base$FDR)) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- existen NA en los resultados."
      )
    )
  }

  resultado <- resultado_base[
    order(resultado_base$PValue),
    ,
    drop = FALSE
  ]

  rownames(resultado) <- NULL

  # ----------------------------------------------------------
  # Verificación independiente con topTags
  # ----------------------------------------------------------

  tt <- topTags(
    qlf,
    n = Inf,
    adjust.method = "BH",
    sort.by = "PValue"
  )$table

  if (!isTRUE(all.equal(
    resultado$logFC,
    tt$logFC,
    tolerance = 1e-12
  ))) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- logFC no coincide con topTags."
      )
    )
  }

  if (!isTRUE(all.equal(
    resultado$PValue,
    tt$PValue,
    tolerance = 1e-12
  ))) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- PValue no coincide con topTags."
      )
    )
  }

  if (!isTRUE(all.equal(
    resultado$FDR,
    tt$FDR,
    tolerance = 1e-12
  ))) {
    stop(
      paste(
        "ERROR:",
        nombre,
        "- FDR no coincide con topTags."
      )
    )
  }

  cat(
    "\nVERIFICACIÓN:",
    nombre,
    "- identidad génica y estadísticas verificadas.\n"
  )

  # ----------------------------------------------------------
  # Resumen
  # ----------------------------------------------------------

  cat("\nNúmero total de genes:\n")
  print(nrow(resultado))

  cat("\nResumen de logFC:\n")
  print(summary(resultado$logFC))

  cat("\nGenes con FDR < 0.05:\n")
  print(sum(resultado$FDR < 0.05))

  cat("\nGenes con FDR < 0.01:\n")
  print(sum(resultado$FDR < 0.01))

  cat("\nFDR < 0.05 y |logFC| >= 1:\n")
  print(
    sum(
      resultado$FDR < 0.05 &
      abs(resultado$logFC) >= 1
    )
  )

  cat("\nFDR < 0.05 y logFC >= 1:\n")
  print(
    sum(
      resultado$FDR < 0.05 &
      resultado$logFC >= 1
    )
  )

  cat("\nFDR < 0.05 y logFC <= -1:\n")
  print(
    sum(
      resultado$FDR < 0.05 &
      resultado$logFC <= -1
    )
  )

  cat("\n20 genes con menor PValue:\n")

  print(
    resultado[
      seq_len(min(20, nrow(resultado))),
      c(
        "gene",
        "logFC",
        "logCPM",
        "F",
        "PValue",
        "FDR"
      )
    ]
  )

  # ----------------------------------------------------------
  # Guardar tablas
  # ----------------------------------------------------------

  archivo_base <- paste0(
    "GPL16791_DE_",
    nombre
  )

  write.csv(
    resultado,
    paste0(archivo_base, "_todos.csv"),
    row.names = FALSE
  )

  write.csv(
    resultado[
      resultado$FDR < 0.05,
      ,
      drop = FALSE
    ],
    paste0(archivo_base, "_FDR005.csv"),
    row.names = FALSE
  )

  write.csv(
    resultado[
      resultado$FDR < 0.05 &
      abs(resultado$logFC) >= 1,
      ,
      drop = FALSE
    ],
    paste0(
      archivo_base,
      "_FDR005_logFC1.csv"
    ),
    row.names = FALSE
  )

  saveRDS(
    qlf,
    paste0(
      "GSE108022_GPL16791_",
      nombre,
      "_QLFTest.rds"
    )
  )

  return(resultado)
}


# ------------------------------------------------------------
# 4. Ejecutar FN - Normal
#
# logFC > 0 = mayor en FN
# logFC < 0 = mayor en Normal
# ------------------------------------------------------------

resultado_FN <- analizar_contraste(
  fit = fit,
  contraste = contrastes[, "FN_vs_Normal"],
  nombre = "FN_vs_Normal"
)


# ------------------------------------------------------------
# 5. Ejecutar FOXO1-FP - Normal
#
# logFC > 0 = mayor en FOXO1-FP
# logFC < 0 = mayor en Normal
# ------------------------------------------------------------

resultado_FOXO1FP <- analizar_contraste(
  fit = fit,
  contraste = contrastes[, "FOXO1FP_vs_Normal"],
  nombre = "FOXO1FP_vs_Normal"
)


# ------------------------------------------------------------
# 6. Confirmar archivos
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("ARCHIVOS GENERADOS\n")
cat("====================================================\n")

cat("- GPL16791_DE_FN_vs_Normal_todos.csv\n")
cat("- GPL16791_DE_FN_vs_Normal_FDR005.csv\n")
cat("- GPL16791_DE_FN_vs_Normal_FDR005_logFC1.csv\n")

cat("- GPL16791_DE_FOXO1FP_vs_Normal_todos.csv\n")
cat("- GPL16791_DE_FOXO1FP_vs_Normal_FDR005.csv\n")
cat("- GPL16791_DE_FOXO1FP_vs_Normal_FDR005_logFC1.csv\n")

cat("- GSE108022_GPL16791_FN_vs_Normal_QLFTest.rds\n")
cat("- GSE108022_GPL16791_FOXO1FP_vs_Normal_QLFTest.rds\n")

cat(
  "\nRESULTADO FINAL: ambos contrastes completados con identidad génica y estadísticas verificadas.\n"
)