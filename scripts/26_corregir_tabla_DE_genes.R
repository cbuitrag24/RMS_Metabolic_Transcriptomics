# ============================================================
# 26_corregir_tabla_DE_genes.R
# Reconstrucción correcta de la tabla de expresión diferencial
# preservando la identidad de los genes
# GSE108022 - GPL16791
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Cargar el objeto del test QL ya calculado
# ------------------------------------------------------------

qlf <- readRDS(
  "GSE108022_GPL16791_QLFTest_RMS_vs_Normal.rds"
)

cat("Objeto QL cargado.\n")

cat("\nNúmero de filas en qlf$table:\n")
print(nrow(qlf$table))

cat("\nNúmero de genes en qlf$genes:\n")
print(nrow(qlf$genes))


# ------------------------------------------------------------
# 2. Verificaciones antes de construir la tabla
# ------------------------------------------------------------

if (nrow(qlf$table) != nrow(qlf$genes)) {
  stop("ERROR: qlf$table y qlf$genes tienen distinto número de filas.")
}

if (!"gene" %in% colnames(qlf$genes)) {
  stop("ERROR: qlf$genes no contiene la columna gene.")
}

if (anyNA(qlf$genes$gene)) {
  stop("ERROR: existen símbolos génicos NA.")
}

if (any(qlf$genes$gene == "")) {
  stop("ERROR: existen símbolos génicos vacíos.")
}

cat("\nPrimeros 10 símbolos almacenados en qlf$genes:\n")
print(head(qlf$genes$gene, 10))


# ------------------------------------------------------------
# 3. Construir tabla SIN ordenar
#
# La columna gene se incorpora antes de cualquier ordenamiento.
# Así cada símbolo queda unido a su fila estadística original.
# ------------------------------------------------------------

resultado_base <- data.frame(
  gene = qlf$genes$gene,
  qlf$table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nDimensiones de la tabla base:\n")
print(dim(resultado_base))

cat("\nPrimeras 10 filas de la tabla base:\n")
print(
  head(
    resultado_base[
      ,
      c("gene", "logFC", "logCPM", "F", "PValue")
    ],
    10
  )
)


# ------------------------------------------------------------
# 4. Verificación explícita de correspondencia
# ------------------------------------------------------------

if (!identical(
  as.character(resultado_base$gene),
  as.character(qlf$genes$gene)
)) {
  stop("ERROR: se perdió la correspondencia entre genes y estadísticas.")
}

cat(
  "\nVERIFICACIÓN: identidad gen-fila preservada antes del ordenamiento.\n"
)


# ------------------------------------------------------------
# 5. Calcular FDR por Benjamini-Hochberg
# ------------------------------------------------------------

resultado_base$FDR <- p.adjust(
  resultado_base$PValue,
  method = "BH"
)

if (anyNA(resultado_base$FDR)) {
  stop("ERROR: existen FDR NA.")
}


# ------------------------------------------------------------
# 6. Ordenar por PValue
#
# Ahora se ordena la fila COMPLETA, incluyendo el símbolo génico.
# ------------------------------------------------------------

orden <- order(
  resultado_base$PValue,
  decreasing = FALSE
)

resultado <- resultado_base[
  orden,
  ,
  drop = FALSE
]

rownames(resultado) <- NULL


# ------------------------------------------------------------
# 7. Comprobación independiente con topTags
# ------------------------------------------------------------

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
  stop("ERROR: los logFC no coinciden con topTags.")
}

if (!isTRUE(all.equal(
  resultado$PValue,
  tt$PValue,
  tolerance = 1e-12
))) {
  stop("ERROR: los PValue no coinciden con topTags.")
}

if (!isTRUE(all.equal(
  resultado$FDR,
  tt$FDR,
  tolerance = 1e-12
))) {
  stop("ERROR: los FDR no coinciden con topTags.")
}

cat(
  "VERIFICACIÓN: logFC, PValue y FDR coinciden con topTags.\n"
)


# ------------------------------------------------------------
# 8. Resumen del resultado corregido
# ------------------------------------------------------------

cat("\nNúmero total de genes:\n")
print(nrow(resultado))

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

cat("\nFDR < 0.05 y logFC >= 1 (mayor en RMS):\n")
print(
  sum(
    resultado$FDR < 0.05 &
    resultado$logFC >= 1
  )
)

cat("\nFDR < 0.05 y logFC <= -1 (mayor en normal):\n")
print(
  sum(
    resultado$FDR < 0.05 &
    resultado$logFC <= -1
  )
)


# ------------------------------------------------------------
# 9. Mostrar los 20 genes con menor PValue
# ------------------------------------------------------------

cat("\n20 genes con menor PValue, ahora con identidad verificada:\n")

print(
  resultado[
    seq_len(min(20, nrow(resultado))),
    c("gene", "logFC", "logCPM", "F", "PValue", "FDR")
  ]
)


# ------------------------------------------------------------
# 10. Guardar tablas CORREGIDAS
# ------------------------------------------------------------

write.csv(
  resultado,
  "GPL16791_DE_RMS_vs_Normal_todos_los_genes_CORREGIDO.csv",
  row.names = FALSE
)

significativos <- resultado[
  resultado$FDR < 0.05,
  ,
  drop = FALSE
]

write.csv(
  significativos,
  "GPL16791_DE_RMS_vs_Normal_FDR005_CORREGIDO.csv",
  row.names = FALSE
)

significativos_fc1 <- resultado[
  resultado$FDR < 0.05 &
  abs(resultado$logFC) >= 1,
  ,
  drop = FALSE
]

write.csv(
  significativos_fc1,
  "GPL16791_DE_RMS_vs_Normal_FDR005_logFC1_CORREGIDO.csv",
  row.names = FALSE
)

cat("\nArchivos generados:\n")
cat("- GPL16791_DE_RMS_vs_Normal_todos_los_genes_CORREGIDO.csv\n")
cat("- GPL16791_DE_RMS_vs_Normal_FDR005_CORREGIDO.csv\n")
cat("- GPL16791_DE_RMS_vs_Normal_FDR005_logFC1_CORREGIDO.csv\n")

cat("\nRESULTADO FINAL: tabla reconstruida con identidad génica verificada.\n")