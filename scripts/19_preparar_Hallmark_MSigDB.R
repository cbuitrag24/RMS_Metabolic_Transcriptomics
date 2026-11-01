# ============================================================
# 19_preparar_Hallmark_MSigDB.R
#
# OBJETIVO:
# Obtener y auditar la colección Hallmark humana de MSigDB
# para el posterior análisis funcional prerankeado.
#
# En este paso:
# - se obtiene Hallmark para Homo sapiens;
# - se conserva Symbol -> gene set;
# - se controla número de conjuntos y genes;
# - se guarda exactamente la colección utilizada.
#
# NO se realiza todavía ningún enriquecimiento.
# ============================================================


# ------------------------------------------------------------
# 1. Verificar paquete msigdbr
# ------------------------------------------------------------

if (!requireNamespace("msigdbr", quietly = TRUE)) {

  cat(
    "El paquete msigdbr no está instalado.\n"
  )

  cat(
    "Instalando msigdbr desde CRAN...\n\n"
  )

  install.packages(
    "msigdbr",
    repos = "https://cloud.r-project.org"
  )
}


suppressPackageStartupMessages(
  library(msigdbr)
)


cat("Versión de msigdbr:\n")
cat(
  as.character(
    packageVersion("msigdbr")
  ),
  "\n"
)


# ------------------------------------------------------------
# 2. Obtener Hallmark humano
# ------------------------------------------------------------

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)


cat("\nDimensiones originales de Hallmark:\n")
print(
  dim(hallmark)
)


# ------------------------------------------------------------
# 3. Inspeccionar columnas disponibles
# ------------------------------------------------------------

cat("\nColumnas disponibles:\n")
print(
  colnames(hallmark)
)


necesarias <- c(
  "gs_name",
  "gene_symbol"
)

if (!all(necesarias %in% colnames(hallmark))) {

  stop(
    paste(
      "ERROR: faltan columnas necesarias:",
      paste(
        setdiff(
          necesarias,
          colnames(hallmark)
        ),
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 4. Construir tabla mínima reproducible
# ------------------------------------------------------------

hallmark_min <- data.frame(
  gene_set =
    as.character(hallmark$gs_name),

  gene =
    trimws(
      as.character(
        hallmark$gene_symbol
      )
    ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Eliminar únicamente identificadores vacíos
# ------------------------------------------------------------

hallmark_min <- hallmark_min[
  !is.na(hallmark_min$gene) &
    hallmark_min$gene != "" &
    !is.na(hallmark_min$gene_set) &
    hallmark_min$gene_set != "",
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 6. Verificar duplicados exactos
#
# Un mismo gen puede pertenecer legítimamente a distintos
# Hallmarks. Solo nos interesan duplicados exactos:
# mismo gene_set + mismo gene.
# ------------------------------------------------------------

duplicados_exactos <- duplicated(
  hallmark_min[
    ,
    c(
      "gene_set",
      "gene"
    )
  ]
)


cat("\nDuplicados exactos gene_set + gene:\n")
cat(
  sum(duplicados_exactos),
  "\n"
)


if (any(duplicados_exactos)) {

  hallmark_min <- hallmark_min[
    !duplicados_exactos,
    ,
    drop = FALSE
  ]

  cat(
    "Los duplicados exactos fueron eliminados.\n"
  )
}


# ------------------------------------------------------------
# 7. Auditoría general
# ------------------------------------------------------------

n_sets <- length(
  unique(
    hallmark_min$gene_set
  )
)

n_genes <- length(
  unique(
    hallmark_min$gene
  )
)


cat("\n====================================================\n")
cat("AUDITORÍA HALLMARK\n")
cat("====================================================\n")

cat("\nNúmero de Hallmarks:\n")
cat(
  n_sets,
  "\n"
)

cat("\nGenes únicos representados:\n")
cat(
  n_genes,
  "\n"
)

cat("\nAsignaciones gene-set / gene:\n")
cat(
  nrow(hallmark_min),
  "\n"
)


# ------------------------------------------------------------
# 8. Número de genes por Hallmark
# ------------------------------------------------------------

tamano_sets <- aggregate(
  gene ~ gene_set,
  data = hallmark_min,
  FUN = function(x) {
    length(unique(x))
  }
)

colnames(tamano_sets)[2] <-
  "n_genes"


tamano_sets <- tamano_sets[
  order(
    tamano_sets$gene_set
  ),
  ,
  drop = FALSE
]

rownames(tamano_sets) <- NULL


cat("\nResumen del tamaño de los Hallmarks:\n")
print(
  summary(
    tamano_sets$n_genes
  )
)


cat("\nLista de Hallmarks y tamaños:\n")
print(
  tamano_sets,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Verificación del número de Hallmarks
#
# La colección H humana contiene 50 Hallmarks.
# ------------------------------------------------------------

if (n_sets != 50) {

  stop(
    paste(
      "ERROR: se esperaban 50 Hallmarks y se encontraron",
      n_sets,
      ". No continuar con enriquecimiento."
    )
  )
}


# ------------------------------------------------------------
# 10. Verificar nombres
# ------------------------------------------------------------

if (anyDuplicated(tamano_sets$gene_set)) {

  stop(
    "ERROR: existen nombres de Hallmark duplicados."
  )
}


if (!all(
  grepl(
    "^HALLMARK_",
    tamano_sets$gene_set
  )
)) {

  stop(
    "ERROR: existen conjuntos con nombres inesperados."
  )
}


cat(
  "\nVERIFICACIÓN: se recuperaron exactamente",
  "50 Hallmarks humanos.\n"
)


# ------------------------------------------------------------
# 11. Información de versión disponible en msigdbr
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("INFORMACIÓN DE VERSIÓN\n")
cat("====================================================\n")


if ("db_version" %in% colnames(hallmark)) {

  versiones <- unique(
    hallmark$db_version
  )

  cat("\nVersión MSigDB informada por msigdbr:\n")
  print(
    versiones
  )

} else {

  cat(
    "\nLa tabla devuelta no contiene",
    "una columna db_version.\n"
  )
}


if ("db_species" %in% colnames(hallmark)) {

  cat("\nEspecie de la base MSigDB:\n")

  print(
    unique(
      hallmark$db_species
    )
  )
}


# ------------------------------------------------------------
# 12. Guardar colección exacta utilizada
# ------------------------------------------------------------

write.csv(
  hallmark_min,
  "MSigDB_Hallmark_Homo_sapiens.csv",
  row.names = FALSE
)

saveRDS(
  hallmark_min,
  "MSigDB_Hallmark_Homo_sapiens.rds"
)


write.csv(
  tamano_sets,
  "MSigDB_Hallmark_tamanos.csv",
  row.names = FALSE
)


cat("\nArchivos generados:\n")

cat(
  "- MSigDB_Hallmark_Homo_sapiens.csv\n"
)

cat(
  "- MSigDB_Hallmark_Homo_sapiens.rds\n"
)

cat(
  "- MSigDB_Hallmark_tamanos.csv\n"
)


# ------------------------------------------------------------
# 13. Final
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("PREPARACIÓN DE HALLMARK COMPLETADA\n")
cat("====================================================\n")

cat(
  "\nTodavía NO se realizó ningún análisis",
  "de enriquecimiento funcional.\n"
)

cat(
  "\nLa colección Hallmark utilizada quedó",
  "guardada localmente para reproducibilidad.\n"
)