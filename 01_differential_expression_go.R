############################################################
# TCGA-BRCA - Expresión diferencial + análisis funcional (GO-BP)
# Autor: José Carlos León Miguel
# Máster Universitario en Bioinformática - Universidad Europea
#
# OBJETIVO:
# 1) Cargar matriz de expresión y datos clínicos de TCGA-BRCA.
# 2) Explorar dimensiones (genes x muestras) y nº de variables clínicas.
# 3) Normalizar datos RNA-seq (TMM) y evaluar el efecto (boxplots).
# 4) Análisis de expresión diferencial ER+ vs ER- corrigiendo por edad.
# 5) Resumen: nº genes significativos, up/down, y Top 20 genes.
# 6) Enriquecimiento funcional (GO-BP) con clusterProfiler y dotplot.
#
# USO:
#   Rscript 01_differential_expression_go.R [expresion.tsv] [clinico.tsv] [salida]
#
# Sin argumentos, los datos se leen de data/ y los resultados se escriben
# en output/. Los datos de TCGA no se incluyen en el repositorio (ver README).
############################################################

## =========================
## 0) Paquetes
## =========================
# En este bloque cargo/instalo los paquetes necesarios.
# - CRAN: herramientas de manipulación y gráficos (ggplot2, dplyr, tibble)
# - Bioconductor: análisis de RNA-seq (edgeR, limma) y análisis funcional (clusterProfiler, org.Hs.eg.db)
# Si falta alguno, se instala automáticamente.
required_pkgs_cran <- c("ggplot2", "dplyr", "tibble", "ggrepel")
required_pkgs_bioc <- c("edgeR", "limma", "clusterProfiler", "org.Hs.eg.db")

to_install_cran <- required_pkgs_cran[!sapply(required_pkgs_cran, requireNamespace, quietly = TRUE)]
if (length(to_install_cran) > 0) install.packages(to_install_cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- required_pkgs_bioc[!sapply(required_pkgs_bioc, requireNamespace, quietly = TRUE)]
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE)

library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)

# Semilla para reproducibilidad (por si en algún paso hubiese aleatoriedad)
set.seed(1)

## =========================
## 1) Rutas de entrada
## =========================
# Las rutas se resuelven en este orden: argumento de línea de comandos,
# si no, el nombre por defecto dentro de data/. Así el script es reproducible
# en cualquier sistema operativo y no depende de rutas locales.
args <- commandArgs(trailingOnly = TRUE)

DATA_DIR <- "data"

ruta_entrada <- function(i, nombre_por_defecto, descripcion) {
  ruta <- if (length(args) >= i && nzchar(args[i])) args[i] else file.path(DATA_DIR, nombre_por_defecto)
  if (!file.exists(ruta)) {
    stop(sprintf(paste0("No se encuentra el archivo %s: '%s'.\n",
                        "  Coloca los datos en '%s/' o pasa la ruta como argumento (ver README)."),
                 descripcion, ruta, DATA_DIR), call. = FALSE)
  }
  normalizePath(ruta)
}

expr_file <- ruta_entrada(1, "BRCA_exp_matrix.tsv",          "de EXPRESIÓN")
clin_file <- ruta_entrada(2, "clinical_info_TCGA-BRCA.tsv",  "CLÍNICO")

message("Matriz de expresión: ", expr_file)
message("Datos clínicos     : ", clin_file)

## =========================
## 1b) Carpeta de salida
## =========================
# En esta carpeta se guardan automáticamente:
# - resumen de exploración
# - figuras de normalización (boxplots)
# - resultados completos de expresión diferencial + resumen + top20
# - resultados de GO enrichment + dotplot
outdir <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
message("Salida en: ", normalizePath(outdir))

## =========================
## 2) Carga de datos
## =========================
# Cargo los archivos TSV:
# - expr: matriz de expresión (genes en filas, muestras en columnas)
# - clinical: tabla clínica (muestras en filas y variables clínicas en columnas)
# Importante: check.names = FALSE para no modificar nombres de columnas (IDs de muestra)
expr <- read.delim(expr_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
clinical <- read.delim(clin_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
clinical <- as.data.frame(clinical)

# Compruebo que el clínico tiene una columna "sample" para poder alinear con la matriz
if (!"sample" %in% colnames(clinical)) stop("No encuentro la columna 'sample' en el clínico.")
clinical$sample <- trimws(as.character(clinical$sample))

## =========================
## 3) Alineación por muestra (base R)
## =========================
# En este paso hago lo que pide el enunciado en la parte de "explorar dataset":
# - Alinear expresión y clínico por el mismo conjunto de muestras.
# Uso intersect() para quedarme con las muestras comunes y match() para reordenar.
common_samples <- intersect(colnames(expr), clinical$sample)
if (length(common_samples) < 20) warning("Pocas muestras en común. ¿Elegiste bien los archivos?")

expr <- expr[, common_samples, drop = FALSE]
clinical <- clinical[clinical$sample %in% common_samples, , drop = FALSE]
clinical <- clinical[match(common_samples, clinical$sample), , drop = FALSE]

# Verificación rápida: el orden de clinical debe coincidir con las columnas de expr
stopifnot(all(clinical$sample == colnames(expr)))

## =========================
## 4) ER (Positive vs Negative) + edad
## =========================
# Variable principal del enunciado: estrogen_receptor
# Además, el enunciado pide corregir por edad.
# Aquí:
# - convierto er_status_by_ihc a un factor limpio: Positive / Negative
# - elimino casos indeterminados (Indeterminate u otros valores raros)
# - dejo "Negative" como referencia (ER-), de forma que el coeficiente sea ER+ vs ER-
if (!"er_status_by_ihc" %in% colnames(clinical)) stop("No encuentro 'er_status_by_ihc' en el clínico.")
if (!"age_at_diagnosis" %in% colnames(clinical)) stop("No encuentro 'age_at_diagnosis' en el clínico.")

er_raw <- trimws(as.character(clinical$er_status_by_ihc))
er_raw <- ifelse(tolower(er_raw) %in% c("positive","pos","1","yes"), "Positive", er_raw)
er_raw <- ifelse(tolower(er_raw) %in% c("negative","neg","0","no"), "Negative", er_raw)
clinical$estrogen_receptor <- factor(er_raw)

# Eliminar indeterminate y valores no deseados
keep_er <- clinical$estrogen_receptor %in% c("Positive", "Negative")
clinical <- clinical[keep_er, , drop = FALSE]
expr <- expr[, clinical$sample, drop = FALSE]

# Referencia: Negative
clinical$estrogen_receptor <- factor(clinical$estrogen_receptor)
clinical$estrogen_receptor <- relevel(clinical$estrogen_receptor, ref = "Negative")

# Edad: el clínico está en días, lo paso a años (aprox)
clinical$age_years <- suppressWarnings(as.numeric(clinical$age_at_diagnosis)) / 365.25

# Elimino muestras con edad NA y re-sincronizo la matriz de expresión
keep_age <- complete.cases(clinical$age_years)
clinical <- clinical[keep_age, , drop = FALSE]
expr <- expr[, clinical$sample, drop = FALSE]

## =========================
## 5) Exploración
## =========================
# Enunciado (punto 2): reportar dimensiones y nº de variables clínicas.
# En lugar de ponerlo a mano, lo guardo en un TXT para copiarlo al informe.
writeLines(c(
  "=== EXPLORACIÓN ===",
  paste0("Expresión (genes x muestras): ", nrow(expr), " x ", ncol(expr)),
  paste0("Clínico (muestras x variables): ", nrow(clinical), " x ", ncol(clinical)),
  paste0("Número variables clínicas: ", ncol(clinical)),
  paste0("N ER Positive: ", sum(clinical$estrogen_receptor == "Positive")),
  paste0("N ER Negative: ", sum(clinical$estrogen_receptor == "Negative"))
), con = file.path(outdir, "00_resumen_exploracion.txt"))

## =========================
## 6) Normalización + boxplots
## =========================
# Enunciado (punto 3): normalizar de forma adecuada para RNA-seq.
# Uso edgeR para normalización TMM y paso a logCPM.
# Además, el enunciado pide un boxplot de expresión media por muestra antes/después.
expr_mat <- as.matrix(expr)
storage.mode(expr_mat) <- "numeric"

# Antes normalizar: transformación sencilla log2(count+1) para visualizar
log_before <- log2(expr_mat + 1)

# Normalización TMM + logCPM (edgeR)
y <- DGEList(counts = expr_mat)
y <- calcNormFactors(y, method = "TMM")
logCPM <- cpm(y, log = TRUE, prior.count = 1)

# (A) Boxplot medias por muestra antes/después (punto 4 del enunciado)
mean_before <- colMeans(log_before, na.rm = TRUE)
mean_after  <- colMeans(logCPM, na.rm = TRUE)

df_means <- tibble(
  sample   = rep(colnames(expr_mat), times = 2),
  mean_expr = c(mean_before, mean_after),
  estado   = rep(c("Antes", "Después"), each = length(mean_before))
)

p_means <- ggplot(df_means, aes(x = estado, y = mean_expr)) +
  geom_boxplot() +
  labs(
    title = "Media de expresión por muestra (log2) antes y después de normalizar",
    x = "", y = "Media (log2)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "01_boxplot_medias_antes_vs_despues.png"),
       p_means, width = 8, height = 5, dpi = 200)

# (B) Extra: boxplots por muestra completos (antes/después),
# para ver si las distribuciones se homogeneizan tras normalizar.
png(file.path(outdir, "01b_boxplot_distribucion_por_muestra_ANTES.png"),
    width = 1600, height = 800, res = 160)
boxplot(log_before, las = 2, cex.axis = 0.5,
        main = "Distribución por muestra (log2(count+1)) - ANTES",
        ylab = "log2(count+1)")
dev.off()

png(file.path(outdir, "01c_boxplot_distribucion_por_muestra_DESPUES.png"),
    width = 1600, height = 800, res = 160)
boxplot(logCPM, las = 2, cex.axis = 0.5,
        main = "Distribución por muestra (logCPM) - DESPUÉS",
        ylab = "logCPM")
dev.off()

## =========================
## 7) Expresión diferencial (limma) + plots
## =========================
# Enunciado (punto 5): análisis de expresión diferencial con estrogen_receptor
# corrigiendo por edad.
# Modelo: ~ estrogen_receptor + age_years
# Coeficiente: estrogen_receptorPositive (ER+ frente a ER-, porque ER- es la referencia)
design <- model.matrix(~ estrogen_receptor + age_years, data = clinical)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "estrogen_receptorPositive", number = Inf, sort.by = "P")

# Guardo resultados completos
write.table(res,
            file = file.path(outdir, "02_DE_resultados_completos.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# Enunciado (punto 6): nº significativos y up/down.
# Defino significativos como FDR<0.05 (adj.P.Val).
sig <- res %>% rownames_to_column("gene") %>% filter(adj.P.Val < 0.05)
n_sig <- nrow(sig)
n_up  <- sum(sig$logFC > 0)
n_down <- sum(sig$logFC < 0)

# Resumen en TXT para el informe
writeLines(c(
  "=== RESUMEN DE ===",
  "Coeficiente ER usado: estrogen_receptorPositive",
  paste0("Significativos (FDR<0.05): ", n_sig),
  paste0("Upregulated (logFC>0): ", n_up),
  paste0("Downregulated (logFC<0): ", n_down)
), con = file.path(outdir, "03_resumen_DE.txt"))

# Enunciado (punto 6): los 20 genes más significativos
# (ordenados por adj.P.Val)
top20 <- sig %>% arrange(adj.P.Val) %>% head(20)
write.table(top20,
            file = file.path(outdir, "03_top20_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figuras de apoyo (extras):
# - Histograma de p-values: confirma que hay señal (exceso de p pequeños)
p_pval <- ggplot(as.data.frame(res), aes(x = P.Value)) +
  geom_histogram(bins = 50) +
  labs(title = "Histograma de p-values (ER+ vs ER-, ajustado por edad)",
       x = "P.Value", y = "Frecuencia") +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "03b_histograma_pvalues.png"),
       p_pval, width = 8, height = 5, dpi = 200)

# - Volcano plot: magnitud del cambio (logFC) frente a significación.
#   Se colorea según el criterio estricto (FDR<0.05 y |logFC|>1), que es el que
#   define el conjunto de genes que pasa al análisis de enriquecimiento.
volc_df <- res %>%
  rownames_to_column("gene") %>%
  mutate(
    neglog10p = -log10(P.Value),
    categoria = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Sobreexpresado en ER+",
      adj.P.Val < 0.05 & logFC < -1 ~ "Infraexpresado en ER+",
      TRUE                          ~ "No relevante"
    )
  )

n_estrictos <- sum(volc_df$categoria != "No relevante")

top_genes <- volc_df %>%
  filter(categoria != "No relevante") %>%
  arrange(adj.P.Val) %>%
  head(12)

p_volc <- ggplot(volc_df, aes(x = logFC, y = neglog10p, color = categoria)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(data = top_genes, aes(label = gene),
                           size = 3.2, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(
    "Sobreexpresado en ER+" = "#C0392B",
    "Infraexpresado en ER+" = "#2471A3",
    "No relevante"          = "grey75"
  )) +
  labs(title = "Volcano plot - ER+ vs ER- (ajustado por edad)",
       subtitle = paste0(n_estrictos, " genes con FDR<0.05 y |logFC|>1"),
       x = "log2 Fold Change", y = "-log10(P.Value)", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "03c_volcano_plot.png"),
       p_volc, width = 9, height = 6, dpi = 200)

## =========================
## 8) Enriquecimiento estricto (FDR<0.05 & |logFC|>1)
## =========================
# Enunciado (punto 7): análisis de enriquecimiento con clusterProfiler.
# Para hacerlo más interpretable, uso un conjunto "estricto" de genes:
# - FDR<0.05 y además |logFC|>1 (cambios de expresión más grandes).
sig_strict <- sig %>% filter(abs(logFC) > 1)

# Guardo cuántos genes entran al enriquecimiento (útil para el informe)
writeLines(c(
  paste0("Genes estrictos (FDR<0.05 & |logFC|>1): ", nrow(sig_strict))
), con = file.path(outdir, "04_resumen_enriquecimiento_strict.txt"))

# clusterProfiler trabaja con ENTREZID, así que convierto SYMBOL -> ENTREZID
conv <- bitr(sig_strict$gene,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)

# Enriquecimiento GO (Biological Process)
ego <- enrichGO(
  gene          = unique(conv$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# Resultados en TSV
go_df <- as.data.frame(ego)
write.table(go_df,
            file = file.path(outdir, "05_GO_BP_strict.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Dotplot del enriquecimiento (punto 7: adjuntar gráficas en el informe)
p_go <- dotplot(ego, showCategory = 20) +
  ggtitle("GO BP enrichment (FDR<0.05 & |logFC|>1)")

ggsave(file.path(outdir, "06_GO_BP_dotplot.png"),
       p_go, width = 9, height = 6, dpi = 200)

## =========================
## 9) Fin
## =========================
# Mensaje final para confirmar que el script terminó y dónde está la carpeta
message("Análisis completado correctamente.")
message("Carpeta generada en: ", outdir)
