# ============================================================
# 14_preparar_diseno_DE.R
# Preparación y verificación del diseño para expresión diferencial
# GSE108022
# ============================================================

library(edgeR)

# ------------------------------------------------------------
# 1. Leer metadata
# ------------------------------------------------------------

metadata <- read.delim(
  "GSE108022_metadata.txt.gz",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensiones metadata:\n")
print(dim(metadata))

cat("\nColumnas metadata:\n")
print(colnames(metadata))


# ------------------------------------------------------------
# 2. Leer asignación de plataforma obtenida previamente
#    desde el archivo MINiML
# ------------------------------------------------------------

xml_file <- "GSE108022_family.xml"

xml_lines <- readLines(xml_file, warn = FALSE)

# Extraer IDs GSM y plataformas GPL manteniendo el orden
# Identificar cada bloque <Sample> del MINiML
sample_start <- grep("<Sample iid=", xml_lines, fixed = TRUE)
sample_end   <- grep("</Sample>", xml_lines, fixed = TRUE)

if (length(sample_start) != length(sample_end)) {
  stop("ERROR: no coincide el número de aperturas y cierres de bloques Sample.")
}

sample_names <- character(length(sample_start))
gpl <- character(length(sample_start))

# Extraer el Title y la plataforma dentro de cada bloque Sample
for (i in seq_along(sample_start)) {

  block <- xml_lines[sample_start[i]:sample_end[i]]

  title_line <- grep("<Title>", block, value = TRUE, fixed = TRUE)
  platform_line <- grep("<Platform-Ref ref=", block, value = TRUE, fixed = TRUE)

  if (length(title_line) != 1 || length(platform_line) != 1) {
    stop(paste("ERROR al interpretar el bloque de muestra", i))
  }

  sample_names[i] <- sub(
    ".*<Title>(.*)</Title>.*",
    "\\1",
    title_line
  )

  gpl[i] <- sub(
    '.*<Platform-Ref ref="([^"]+)".*',
    "\\1",
    platform_line
  )
}

cat("\nNúmero de muestras encontradas en MINiML:", length(sample_names), "\n")
cat("Número de plataformas encontradas:", length(gpl), "\n")

platform_table <- data.frame(
  samplename = sample_names,
  platform = gpl,
  stringsAsFactors = FALSE
)

cat("\nPrimeras correspondencias muestra-plataforma:\n")
print(head(platform_table, 10))
# ------------------------------------------------------------
# 3. Unir metadata y plataforma
# ------------------------------------------------------------

design_data <- merge(
  metadata,
  platform_table,
  by = "samplename",
  all.x = TRUE,
  sort = FALSE
)

# Recuperar el orden original del metadata
design_data <- design_data[
  match(metadata$samplename, design_data$samplename),
]


# ------------------------------------------------------------
# 4. Definir condición general RMS vs Normal
# ------------------------------------------------------------

design_data$condition <- ifelse(
  design_data$fusion_status == "Normal",
  "Normal",
  "RMS"
)

design_data$condition <- factor(
  design_data$condition,
  levels = c("Normal", "RMS")
)

design_data$platform <- factor(design_data$platform)


# ------------------------------------------------------------
# 5. Controles básicos
# ------------------------------------------------------------

cat("\nNúmero total de muestras:", nrow(design_data), "\n")

cat("\nCondition:\n")
print(table(design_data$condition))

cat("\nPlatform:\n")
print(table(design_data$platform))

cat("\nCondition x Platform:\n")
print(table(design_data$condition, design_data$platform))

cat("\nPlataformas faltantes:", sum(is.na(design_data$platform)), "\n")


# ------------------------------------------------------------
# 6. Construir matriz de diseño
#    Ajuste por plataforma + efecto RMS vs Normal
# ------------------------------------------------------------

design <- model.matrix(
  ~ platform + condition,
  data = design_data
)

cat("\nMatriz de diseño:\n")
print(colnames(design))

cat("\nDimensiones de la matriz de diseño:\n")
print(dim(design))

cat("\nRango de la matriz:", qr(design)$rank, "\n")
cat("Número de columnas:", ncol(design), "\n")

if (qr(design)$rank == ncol(design)) {
  cat("\nRESULTADO: matriz de diseño de rango completo.\n")
} else {
  cat("\nRESULTADO: matriz de diseño NO es de rango completo.\n")
}


# ------------------------------------------------------------
# 7. Guardar tabla de diseño
# ------------------------------------------------------------

write.csv(
  design_data,
  "DE_diseno_muestras.csv",
  row.names = FALSE
)

cat("\nArchivo generado: DE_diseno_muestras.csv\n")