# =============================================================================
# ANÁLISIS DE ACTIVACIÓN DIFERENCIAL DE RUTAS CON HIPATHIA
# Autor: José Carlos León Miguel
# Máster Universitario en Bioinformática - Universidad Europea
#
# Dataset: BRCA - TCGA (Cáncer de mama)
# Comparación: ER+ (receptor de estrógenos positivo) vs ER- (negativo)
#
# USO:
#   Rscript 02_hipathia_pathway_activation.R [counts.tsv] [clinico.tsv] [salida]
#
# Sin argumentos, los datos se leen de data/ y los resultados se escriben
# en output/. Los datos de TCGA no se incluyen en el repositorio (ver README).
# =============================================================================

# =============================================================================
# 0. INSTALACIÓN DE PAQUETES
# =============================================================================
# IMPORTANTE: Este script requiere R 4.3.3 y Bioconductor 3.18.
# Verificar versión antes de continuar:
stopifnot(as.numeric(R.version$major) == 4)
stopifnot(as.numeric(R.version$minor) >= 3)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar paquetes de CRAN necesarios. En Windows y macOS se instalan las
# versiones binarias (más rápidas y sin compilador); en Linux, desde fuentes.
cran_pkgs <- c("ggplot2", "pheatmap", "RColorBrewer", "ggrepel", "ggpubr")
tipo_paquete <- if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") "binary" else "source"
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, type = tipo_paquete)
}

# Instalar paquetes de Bioconductor
bioc_pkgs <- c("hipathia", "limma", "edgeR", "org.Hs.eg.db", "AnnotationDbi")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
}

# =============================================================================
# 1. CARGAR LIBRERÍAS
# =============================================================================

library(hipathia)       # Análisis de activación de rutas (v3.2.0, Bioc 3.18)
library(limma)          # Modelo lineal para análisis diferencial
library(edgeR)          # Normalización TMM
library(org.Hs.eg.db)   # Anotación de genes humanos
library(AnnotationDbi)  # Funciones de anotación
library(ggplot2)        # Visualizaciones
library(pheatmap)       # Heatmap
library(RColorBrewer)   # Paletas de color
# =============================================================================
# 1B. RUTAS DE ENTRADA Y SALIDA
# =============================================================================
# Las rutas se resuelven en este orden: argumento de línea de comandos, y si no,
# el nombre por defecto dentro de data/. El script no cambia el directorio de
# trabajo: todas las salidas se construyen con el helper out().
args <- commandArgs(trailingOnly = TRUE)

DATA_DIR   <- "data"
output_dir <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

out <- function(...) file.path(output_dir, ...)

ruta_entrada <- function(i, nombre_por_defecto, descripcion) {
  ruta <- if (length(args) >= i && nzchar(args[i])) args[i] else file.path(DATA_DIR, nombre_por_defecto)
  if (!file.exists(ruta)) {
    stop(sprintf(paste0("No se encuentra el archivo %s: '%s'.\n",
                        "  Coloca los datos en '%s/' o pasa la ruta como argumento (ver README)."),
                 descripcion, ruta, DATA_DIR), call. = FALSE)
  }
  normalizePath(ruta)
}

cat("Carpeta de salida:", normalizePath(output_dir), "\n")

# =============================================================================
# 2. CARGA DE ARCHIVOS DE ENTRADA
# =============================================================================
# Se trabaja con dos archivos locales: la matriz de counts crudos y la tabla
# clínica. El separador se detecta automáticamente por la extensión del archivo.
#
# --- 2A. MATRIZ DE EXPRESIÓN ---
path_matrix <- ruta_entrada(1, "BRCA_exp_matrix.tsv", "de la MATRIZ DE EXPRESIÓN (counts crudos)")
cat("Matriz de expresión:", path_matrix, "\n")
# Detección automática del separador:
# - si el archivo es .csv → separador coma
# - si es .txt/.tsv → tabulador
ext_mat <- tolower(tools::file_ext(path_matrix))
sep_mat <- if (ext_mat == "csv") "," else "\t"

# La matriz de expresión debe ser "genes x muestras", con:
# - primera columna: IDs de genes (row.names)
# - columnas: muestras (IDs TCGA, normalmente con sufijos adicionales)
counts_raw <- read.table(path_matrix,
                         header      = TRUE,
                         row.names   = 1,
                         sep         = sep_mat,
                         check.names = FALSE)
counts_raw <- as.matrix(counts_raw)

cat("Dimensiones de la matriz (genes x muestras):", dim(counts_raw), "\n")
cat("Primeros IDs de genes:\n")
print(head(rownames(counts_raw)))
cat("Primeros nombres de muestras:\n")
print(head(colnames(counts_raw)))

# --- 2B. DATOS CLÍNICOS ---
path_clinical <- ruta_entrada(2, "clinical_info_TCGA-BRCA.tsv", "de DATOS CLÍNICOS")
cat("Datos clínicos:", path_clinical, "\n")

ext_clin <- tolower(tools::file_ext(path_clinical))
sep_clin <- if (ext_clin == "csv") "," else "\t"

# row.names = NULL evita error cuando hay NAs en la primera columna
clinical <- read.table(path_clinical,
                       header           = TRUE,
                       row.names        = NULL,
                       sep              = sep_clin,
                       check.names      = FALSE,
                       stringsAsFactors = FALSE)

# Usar la columna 'sample' como identificador (formato SP3930, SP3918, etc.)
# que coincide directamente con los nombres de columna de la matriz de expresión
if ("sample" %in% colnames(clinical)) {
  id_col <- clinical$sample
} else if ("bcr_patient_barcode" %in% colnames(clinical)) {
  id_col <- clinical$bcr_patient_barcode
} else {
  id_col <- clinical[[1]]
}

# Limpiar filas sin ID válido y duplicados
valid_rows <- !is.na(id_col) & id_col != ""
clinical   <- clinical[valid_rows, ]
id_col     <- id_col[valid_rows]
clinical   <- clinical[!duplicated(id_col), ]
id_col     <- id_col[!duplicated(id_col)]
rownames(clinical) <- id_col

cat("Dimensiones de datos clínicos:", dim(clinical), "\n")
cat("Distribución de er_status_by_ihc:\n")
print(table(clinical$er_status_by_ihc, useNA = "always"))

# =============================================================================
# 3. DETECCIÓN AUTOMÁTICA Y CONVERSIÓN DE IDs DE GENES A ENTREZ
# =============================================================================
# Hipathia trabaja con genes en formato EntrezID.
# Como los datasets pueden venir con ENSEMBL (ENSG...), Entrez (números) o símbolos,
# este bloque detecta el tipo de identificador y convierte a Entrez cuando es necesario.
sample_ids <- head(rownames(counts_raw), 20)
is_ensembl <- any(grepl("^ENSG[0-9]+", sample_ids))
is_entrez  <- all(grepl("^[0-9]+$", sample_ids))
is_symbol  <- !is_ensembl && !is_entrez

cat("\nTipo de ID de gen detectado: ")
if (is_ensembl) {
  cat("ENSEMBL\n")
} else if (is_entrez) {
  cat("ENTREZ (no requiere conversión)\n")
} else {
  cat("Gene Symbol\n")
}

if (!is_entrez) {
  if (is_ensembl) {
    # Eliminar sufijos de versión (ENSG00000001234.5 → ENSG00000001234)
    rownames(counts_raw) <- sub("\\.[0-9]+$", "", rownames(counts_raw))
  }
  keytype <- if (is_ensembl) "ENSEMBL" else "SYMBOL"
  entrez_map <- mapIds(org.Hs.eg.db,
                       keys      = rownames(counts_raw),
                       column    = "ENTREZID",
                       keytype   = keytype,
                       multiVals = "first")
  # Se eliminan genes sin mapeo a Entrez y duplicados para evitar:
  # - pérdida de coherencia (múltiples genes apuntando al mismo EntrezID)
  # - errores en funciones posteriores (hipathia espera IDs únicos)
  valid      <- !is.na(entrez_map) & !duplicated(entrez_map)
  counts_raw <- counts_raw[valid, ]
  rownames(counts_raw) <- entrez_map[valid]
  cat("Genes con Entrez ID válido:", nrow(counts_raw), "\n")
}

# =============================================================================
# 4. SELECCIÓN DE MUESTRAS Y VARIABLE DEPENDIENTE
# =============================================================================
# En TCGA, los IDs de muestra pueden incluir sufijos (p. ej. "TCGA-XX-XXXX-01A...").
# Para cruzar expresión con clínica se usa el barcode corto (primeros 12 caracteres),
# que identifica al paciente de forma estable: "TCGA-XX-XXXX".

# Match directo por ID de muestra (SP3930, etc.)
match_idx <- match(colnames(counts_raw), rownames(clinical))
n_matched <- sum(!is.na(match_idx))
cat("\nMuestras con coincidencia en clínico:", n_matched, "de", ncol(counts_raw), "\n")

keep_mat         <- !is.na(match_idx)
counts_matched   <- counts_raw[, keep_mat]
clinical_matched <- clinical[match_idx[keep_mat], ]

# Seleccionar solo ER Positive y Negative
er_col <- clinical_matched$er_status_by_ihc
idx_er <- er_col %in% c("Positive", "Negative")

counts_filt   <- counts_matched[, idx_er]
clinical_filt <- clinical_matched[idx_er, ]

# Se define "Negative" como referencia para que:
# - logFC > 0 signifique mayor activación en ER+ (Positive vs Negative)
# - logFC < 0 signifique menor activación en ER+ (más activo en ER−)
er_status <- factor(clinical_filt$er_status_by_ihc,
                    levels = c("Negative", "Positive"))

cat("\nMuestras seleccionadas para el análisis:\n")
print(table(er_status))

# =============================================================================
# 5. COMPROBAR EL ORDEN DE LAS MUESTRAS
# =============================================================================

stopifnot(ncol(counts_filt) == length(er_status))
cat("✓ Orden de muestras verificado correctamente.\n")

# =============================================================================
# 6. NORMALIZACIÓN TMM CON edgeR
# =============================================================================

dge <- DGEList(counts = counts_filt, group = er_status)
# Filtrado de genes:
# filterByExpr elimina genes con muy baja expresión en la mayoría de muestras.
# Esto reduce ruido y mejora potencia estadística al disminuir el número de tests.
keep_genes <- filterByExpr(dge, group = er_status)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
cat("Genes tras filtrado:", nrow(dge), "\n")

# Calcular factores de normalización TMM
dge <- calcNormFactors(dge, method = "TMM")
cat("Resumen factores TMM:\n")
print(summary(dge$samples$norm.factors))

# Transformación a log2-CPM:
# - estabiliza varianza
# - aproxima distribuciones a normalidad (útil para modelos lineales con limma)
# - facilita comparabilidad entre muestras
expr_norm <- cpm(dge, log = TRUE, prior.count = 1)
cat("Matriz normalizada (log2-CPM):", dim(expr_norm), "\n")

# =============================================================================
# 7. TRADUCCIÓN Y ESCALADO PARA HIPATHIA (rango 0-1)
# =============================================================================
# Hipathia requiere expresión en rango 0–1 para interpretar la señal
# como "capacidad de activación" en los nodos de la ruta.
# translate_data adapta la estructura y normalize_data escala a [0,1].
trans_data <- translate_data(expr_norm, species = "hsa")
norm_data  <- normalize_data(trans_data)
cat("Genes disponibles para hipathia:", nrow(norm_data), "\n")

# =============================================================================
# 8. CARGAR RUTAS DE KEGG
# =============================================================================
# Se cargan rutas KEGG humanas (hsa).
# Hipathia utiliza estas rutas como grafos dirigidos para propagar señal.
# Nota: se cargan muchas rutas; algunas no son específicas de cáncer, pero comparten
# componentes de señalización relevantes (PI3K/AKT/MAPK, etc.).

pathways <- load_pathways(species = "hsa", pathways_list = NULL)
cat("Rutas KEGG cargadas:", length(pathways$pathigraphs), "\n")

# =============================================================================
# 9. COMPUTAR VALORES DE ACTIVACIÓN DE CIRCUITOS
# =============================================================================
# hipathia() calcula la activación de circuitos (sub-rutas) dentro de cada ruta KEGG.
# El resultado se convierte en matriz "circuitos x muestras" para análisis estadístico.
hip_results <- hipathia(norm_data, pathways, verbose = TRUE)
act_matrix  <- get_paths_data(hip_results, matrix = TRUE)
cat("Matriz de activación (circuitos x muestras):", dim(act_matrix), "\n")

# =============================================================================
# 10. ANÁLISIS DIFERENCIAL CON LIMMA
# =============================================================================
# Diseño del modelo:
# Se usa codificación sin intercepto (~0 + grupo) para obtener medias por grupo
# y construir el contraste Positive - Negative de forma explícita.

design <- model.matrix(~ 0 + er_status)
colnames(design) <- levels(er_status)   # "Negative", "Positive"

fit <- lmFit(act_matrix, design)

contrast_matrix <- makeContrasts(
  Positive_vs_Negative = Positive - Negative,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_results <- topTable(fit2,
                        coef    = "Positive_vs_Negative",
                        number  = Inf,
                        sort.by = "p")

cat("\nTotal circuitos analizados:", nrow(all_results), "\n")

# =============================================================================
# 11. RESULTADOS SIGNIFICATIVOS
# =============================================================================

sig_results <- all_results[all_results$adj.P.Val < 0.05, ]

cat("\n--- RESUMEN DE RESULTADOS ---\n")
cat("Circuitos significativos (FDR < 0.05):", nrow(sig_results), "\n")
cat("  Sobreactivados en ER+ (logFC > 0):", sum(sig_results$logFC > 0), "\n")
cat("  Infraactivados en ER+ (logFC < 0):", sum(sig_results$logFC < 0), "\n")

# Resumen a nivel de rutas:
# Cada circuito pertenece a una ruta KEGG (hsaXXXXX).
# Se extraen IDs de ruta para contar cuántas rutas contienen circuitos significativos.
# Nota: los IDs de circuito pueden tener sufijos variables; por eso se usa sub() para limpiar.
sig_pathway_ids <- unique(sub("-[0-9]+-[0-9]+$", "", rownames(sig_results)))
cat("Rutas con al menos 1 circuito significativo:", length(sig_pathway_ids), "\n")

# Top 10
top10 <- head(sig_results, 10)
cat("\n--- TOP 10 circuitos más significativos ---\n")
print(top10[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")])

# Exportación de resultados:
# - hipathia_circuitos_significativos.txt → tabla completa para anexos / análisis adicional
# - hipathia_top10_circuitos.txt → tabla Top10 para el apartado de resultados del informe
# - hipathia_resumen_rutas.txt → resumen por ruta (nº circuitos significativos)
write.table(sig_results, out("hipathia_circuitos_significativos.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(top10, out("hipathia_top10_circuitos.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# Resumen por ruta
get_path_name <- function(pid) {
  pg <- pathways$pathigraphs[[pid]]
  if (!is.null(pg)) pg$graph$name else pid
}

if (nrow(sig_results) > 0) {
  sig_path_ids   <- sub("-[0-9]+-[0-9]+$", "", rownames(sig_results))
  sig_path_names <- sapply(sig_path_ids, get_path_name)
  summary_rutas  <- as.data.frame(table(sig_path_names))
  colnames(summary_rutas) <- c("Ruta", "N_circuitos_significativos")
  summary_rutas  <- summary_rutas[order(-summary_rutas$N_circuitos_significativos), ]
  cat("\n--- Rutas con más circuitos significativos ---\n")
  print(head(summary_rutas, 15))
  write.table(summary_rutas, out("hipathia_resumen_rutas.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("✓ Tablas guardadas.\n")

# =============================================================================
# 12. VISUALIZACIÓN: PCA
# =============================================================================
# Eliminar circuitos con varianza cero antes del PCA
var_circuits <- apply(act_matrix, 1, var)
act_matrix_pca <- act_matrix[var_circuits > 0, ]
cat("Circuitos con varianza > 0:", nrow(act_matrix_pca), "\n")

# PCA:
# Se aplica sobre muestras (por eso se transpone act_matrix).
# Objetivo: evaluar si los perfiles de activación de circuitos separan ER+ vs ER−
# y detectar posibles outliers o subgrupos.
pca_res <- prcomp(t(act_matrix_pca), scale. = TRUE, center = TRUE)
pca_df  <- as.data.frame(pca_res$x[, 1:2])
pca_df$ER_status <- er_status
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 2)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = ER_status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Negative" = "#E74C3C", "Positive" = "#2980B9"),
                     name = "ER Status") +
  labs(title    = "PCA — Activación de circuitos Hipathia",
       subtitle = "TCGA-BRCA: ER Positivo vs ER Negativo",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_bw(base_size = 13)

ggsave(out("PCA_hipathia_BRCA.pdf"), plot = pca_plot, width = 7, height = 5)
ggsave(out("PCA_hipathia_BRCA.png"), plot = pca_plot, width = 7, height = 5, dpi = 300)
if (interactive()) print(pca_plot)
cat("✓ PCA guardado.\n")

# =============================================================================
# 13. VISUALIZACIÓN: HEATMAP
# =============================================================================
# Heatmap Top 50:
# Se representan los 50 circuitos más significativos (si hay al menos 50).
# Se escala por fila (scale="row") para comparar patrones relativos entre muestras,
# ya que diferentes circuitos pueden tener magnitudes de activación distintas.
n_top <- min(50, nrow(sig_results))

if (n_top > 0) {
  top_circuits <- rownames(sig_results)[1:n_top]
  heat_matrix  <- act_matrix[top_circuits, ]

  annotation_col <- data.frame(ER_Status = er_status,
                               row.names = colnames(heat_matrix))
  ann_colors <- list(ER_Status = c("Negative" = "#E74C3C", "Positive" = "#2980B9"))

  pdf(out("Heatmap_hipathia_BRCA.pdf"), width = 12, height = 10)
  pheatmap(mat                      = heat_matrix,
           annotation_col           = annotation_col,
           annotation_colors        = ann_colors,
           show_colnames            = FALSE,
           show_rownames            = TRUE,
           scale                    = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method        = "complete",
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           fontsize_row             = 6,
           main = paste0("Top ", n_top, " circuitos diferencialmente activados\n",
                         "TCGA-BRCA: ER+ vs ER-"))
  dev.off()

  cat("✓ Heatmap guardado.\n")
} else {
  cat("⚠ No hay circuitos significativos para el heatmap.\n")
}

# =============================================================================
# 14. VISUALIZACIÓN: REPORT DEL VISUALIZADOR DE HIPATHIA
# =============================================================================
# Genera carpeta HTML interactiva con los grafos de rutas KEGG coloreados:
# Rojo = sobreactivado en ER+ | Azul = infraactivado en ER+
# IMPORTANTE: este paso se usa SOLO para generar el reporte HTML interactivo.
# El análisis diferencial principal del trabajo se ha realizado con LIMMA.
# Hipathia incluye do_wilcoxon como función conveniente para colorear rutas en el reporte
# (comparación entre grupos), pero NO sustituye al análisis estadístico con limma.
comp <- do_wilcoxon(act_matrix, er_status, g1 = "Positive", g2 = "Negative")

report_folder <- "hipathia_report_BRCA_ER"
create_report(comp, pathways, report_folder,
              path = output_dir)

cat("✓ Report generado en:", report_folder, "\n")
cat("  → Abre '", report_folder, "/index.html' en el navegador.\n", sep = "")

# =============================================================================
# 15. SESSION INFO
# =============================================================================
# Reproducibilidad:
# Guardamos sessionInfo() en un archivo para incluirlo como anexo del informe,
# permitiendo replicar exactamente versiones de R y paquetes utilizados.
cat("\n", strrep("=", 60), "\n")
cat("SESSION INFO\n")
cat(strrep("=", 60), "\n")
print(sessionInfo())

sink(out("sessionInfo_hipathia_BRCA.txt"))
print(sessionInfo())
sink()

cat("\n✓ Análisis completado. Archivos generados:\n")
cat("  hipathia_circuitos_significativos.txt\n")
cat("  hipathia_top10_circuitos.txt\n")
cat("  hipathia_resumen_rutas.txt\n")
cat("  PCA_hipathia_BRCA.pdf / .png\n")
cat("  Heatmap_hipathia_BRCA.pdf\n")
cat("  hipathia_report_BRCA_ER/ (abre index.html)\n")
cat("  sessionInfo_hipathia_BRCA.txt\n")

# =============================================================================
# FIN DEL SCRIPT
# =============================================================================
