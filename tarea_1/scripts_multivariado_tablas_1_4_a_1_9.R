# ============================================================
# Scripts completos de análisis descriptivo multivariado
# Tablas 1.4 a 1.9
# ============================================================
# Autor: Alejandro Meneses Portilla.
# Descripción:
#   Este archivo contiene una función independiente para cada
#   conjunto de datos (T1.4 a T1.9), además de funciones auxiliares
#   para cálculo de estadísticas descriptivas, detección simple de
#   outliers multivariados y generación de gráficos.
#
# Requisitos sugeridos:
# install.packages(c(
#   "tidyverse", "GGally", "factoextra", "corrplot",
#   "gridExtra", "reshape2", "cluster"
# ))
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(GGally)
  library(factoextra)
  library(corrplot)
  library(gridExtra)
  library(reshape2)
  library(cluster)
})

# ------------------------------------------------------------
# Funciones auxiliares generales
# ------------------------------------------------------------

crear_carpeta_salida <- function(output_dir = "salidas_multivariado") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  invisible(output_dir)
}

estadisticas_descriptivas <- function(df_num) {
  tibble(
    variable = names(df_num),
    n = sapply(df_num, function(x) sum(!is.na(x))),
    media = sapply(df_num, mean, na.rm = TRUE),
    mediana = sapply(df_num, median, na.rm = TRUE),
    desviacion = sapply(df_num, sd, na.rm = TRUE),
    varianza = sapply(df_num, var, na.rm = TRUE),
    minimo = sapply(df_num, min, na.rm = TRUE),
    q1 = sapply(df_num, quantile, probs = 0.25, na.rm = TRUE),
    q3 = sapply(df_num, quantile, probs = 0.75, na.rm = TRUE),
    maximo = sapply(df_num, max, na.rm = TRUE),
    rango = sapply(df_num, function(x) diff(range(x, na.rm = TRUE))),
    iqr = sapply(df_num, IQR, na.rm = TRUE),
    cv = sapply(df_num, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  )
}

matriz_correlacion <- function(df_num, metodo = "pearson") {
  cor(df_num, use = "pairwise.complete.obs", method = metodo)
}

outliers_mahalanobis <- function(df_num, alpha = 0.975) {
  df_clean <- na.omit(df_num)
  if (nrow(df_clean) < 3) {
    return(list(indices = integer(0), distancias = numeric(0), corte = NA))
  }

  centro <- colMeans(df_clean)
  cov_mat <- cov(df_clean)
  d2 <- mahalanobis(df_clean, center = centro, cov = cov_mat)
  corte <- qchisq(alpha, df = ncol(df_clean))
  idx <- which(d2 > corte)

  list(indices = idx, distancias = d2, corte = corte)
}

preparar_pca <- function(df_num, scale. = TRUE) {
  prcomp(df_num, center = TRUE, scale. = scale.)
}

resumen_pca <- function(pca_obj) {
  var_exp <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  tibble(
    componente = paste0("PC", seq_along(var_exp)),
    varianza_explicada = var_exp,
    varianza_acumulada = cumsum(var_exp)
  )
}

plot_guardar <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(filename = filename, plot = plot_obj, width = width, height = height)
}

# ------------------------------------------------------------
# 1) Tabla 1.4 - The World's 10 Largest Companies
# ------------------------------------------------------------

analizar_t14_companies <- function(
  path = "T1.4.csv",
  output_dir = "salidas_multivariado"
) {
  crear_carpeta_salida(output_dir)

  df <- read.delim("T1.4.txt", header=T)
  rownames(df) <- df$company
  df_num <- df %>% select(-company)

  stats <- estadisticas_descriptivas(df_num)
  corr <- matriz_correlacion(df_num)
  out <- outliers_mahalanobis(df_num)
  pca <- preparar_pca(df_num)
  pca_res <- resumen_pca(pca)

  p_pairs <- ggpairs(
    df,
    columns = 2:4,
    aes()
  ) + ggtitle("Tabla 1.4 - Matriz de dispersión")

  png(file.path(output_dir, "T1_4_pairs.png"), width = 1200, height = 1000, res = 150)
  print(p_pairs)
  dev.off()

  p_pca <- fviz_pca_biplot(
    pca,
    repel = TRUE,
    label = "var",
    habillage = NULL,
    title = "Tabla 1.4 - Biplot PCA"
  )
  plot_guardar(p_pca, file.path(output_dir, "T1_4_pca_biplot.png"))

  png(file.path(output_dir, "T1_4_corrplot.png"), width = 1000, height = 800, res = 150)
  corrplot(corr, method = "color", type = "upper", tl.col = "black", addCoef.col = "black")
  dev.off()

  write.csv(stats, file.path(output_dir, "T1_4_estadisticas.csv"), row.names = FALSE)
  write.csv(pca_res, file.path(output_dir, "T1_4_pca_resumen.csv"), row.names = FALSE)
  write.csv(as.data.frame(corr), file.path(output_dir, "T1_4_correlacion.csv"))

  list(data = df, descriptivos = stats, correlacion = corr, outliers = out, pca = pca, pca_resumen = pca_res)
}

# ------------------------------------------------------------
# 2) Tabla 1.5 - Air-Pollution Data
# ------------------------------------------------------------

analizar_t15_air_pollution <- function(
  path = "T1.5.csv",
  output_dir = "salidas_multivariado"
) {
  crear_carpeta_salida(output_dir)

  df <- read.csv(path, sep = ";")
  df_num <- df

  stats <- estadisticas_descriptivas(df_num)
  corr <- matriz_correlacion(df_num)
  out <- outliers_mahalanobis(df_num)
  pca <- preparar_pca(df_num)
  pca_res <- resumen_pca(pca)

  p_heat <- df_num %>%
    scale() %>%
    as.data.frame() %>%
    mutate(obs = 1:n()) %>%
    pivot_longer(-obs, names_to = "variable", values_to = "z") %>%
    ggplot(aes(x = variable, y = factor(obs), fill = z)) +
    geom_tile() +
    labs(
      title = "Tabla 1.5 - Heatmap de valores estandarizados",
      x = "Variable",
      y = "Observación"
    ) +
    theme_minimal()
  plot_guardar(p_heat, file.path(output_dir, "T1_5_heatmap.png"), width = 9, height = 7)

  p_pca <- fviz_pca_biplot(
    pca,
    repel = TRUE,
    title = "Tabla 1.5 - Biplot PCA"
  )
  plot_guardar(p_pca, file.path(output_dir, "T1_5_pca_biplot.png"))

  png(file.path(output_dir, "T1_5_corrplot.png"), width = 1000, height = 800, res = 150)
  corrplot(corr, method = "color", type = "upper", tl.col = "black", addCoef.col = "black")
  dev.off()

  write.csv(stats, file.path(output_dir, "T1_5_estadisticas.csv"), row.names = FALSE)
  write.csv(pca_res, file.path(output_dir, "T1_5_pca_resumen.csv"), row.names = FALSE)
  write.csv(as.data.frame(corr), file.path(output_dir, "T1_5_correlacion.csv"))

  list(data = df, descriptivos = stats, correlacion = corr, outliers = out, pca = pca, pca_resumen = pca_res)
}

# ------------------------------------------------------------
# 3) Tabla 1.6 - Multiple-Sclerosis Data
# Nota: el archivo viene sin encabezados en el CSV cargado.
# Se asignan nombres genéricos.
# ------------------------------------------------------------

analizar_t16_multiple_sclerosis <- function(
  path = "T1.6.csv",
  output_dir = "salidas_multivariado"
) {
  crear_carpeta_salida(output_dir)

  df <- read.csv(path, sep = ";", header = FALSE)
  names(df) <- c("X1", "X2", "X3", "X4", "X5", "Grupo")
  df$Grupo <- factor(df$Grupo)
  df_num <- df %>% select(-Grupo)

  stats <- estadisticas_descriptivas(df_num)
  corr <- matriz_correlacion(df_num)
  out <- outliers_mahalanobis(df_num)
  pca <- preparar_pca(df_num)
  pca_res <- resumen_pca(pca)

  p_pca <- fviz_pca_ind(
    pca,
    geom = "point",
    habillage = df$Grupo,
    addEllipses = TRUE,
    repel = TRUE,
    title = "Tabla 1.6 - PCA por grupo"
  )
  plot_guardar(p_pca, file.path(output_dir, "T1_6_pca_grupos.png"))

  p_box <- df %>%
    pivot_longer(cols = X1:X5, names_to = "variable", values_to = "valor") %>%
    ggplot(aes(x = variable, y = valor, fill = Grupo)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    labs(title = "Tabla 1.6 - Boxplots por grupo") +
    theme_minimal()
  plot_guardar(p_box, file.path(output_dir, "T1_6_boxplots_grupo.png"), width = 9, height = 6)

  png(file.path(output_dir, "T1_6_corrplot.png"), width = 1000, height = 800, res = 150)
  corrplot(corr, method = "color", type = "upper", tl.col = "black", addCoef.col = "black")
  dev.off()

  write.csv(stats, file.path(output_dir, "T1_6_estadisticas.csv"), row.names = FALSE)
  write.csv(pca_res, file.path(output_dir, "T1_6_pca_resumen.csv"), row.names = FALSE)
  write.csv(as.data.frame(corr), file.path(output_dir, "T1_6_correlacion.csv"))

  list(data = df, descriptivos = stats, correlacion = corr, outliers = out, pca = pca, pca_resumen = pca_res)
}

# ------------------------------------------------------------
# 4) Tabla 1.7 - Radiotherapy Data
# ------------------------------------------------------------

analizar_t17_radiotherapy <- function(
  path = "T1.7.csv",
  output_dir = "salidas_multivariado"
) {
  crear_carpeta_salida(output_dir)

  df <- read.csv(path, sep = ";")
  df$skin_reaction <- factor(df$skin_reaction)
  df_num <- df %>% select(-skin_reaction)

  stats <- estadisticas_descriptivas(df_num)
  corr <- matriz_correlacion(df_num)
  out <- outliers_mahalanobis(df_num)
  pca <- preparar_pca(df_num)
  pca_res <- resumen_pca(pca)

  p_pca <- fviz_pca_ind(
    pca,
    geom = "point",
    habillage = df$skin_reaction,
    addEllipses = TRUE,
    repel = TRUE,
    title = "Tabla 1.7 - PCA por reacción cutánea"
  )
  plot_guardar(p_pca, file.path(output_dir, "T1_7_pca_skin_reaction.png"))

  p_radar_data <- df %>%
    group_by(skin_reaction) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

  # Representación alternativa: perfiles medios por grupo
  p_profile <- p_radar_data %>%
    pivot_longer(-skin_reaction, names_to = "variable", values_to = "media") %>%
    ggplot(aes(x = variable, y = media, group = skin_reaction, color = skin_reaction)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(title = "Tabla 1.7 - Perfil medio por nivel de reacción cutánea") +
    theme_minimal()
  plot_guardar(p_profile, file.path(output_dir, "T1_7_profile_plot.png"), width = 9, height = 6)

  png(file.path(output_dir, "T1_7_corrplot.png"), width = 1000, height = 800, res = 150)
  corrplot(corr, method = "color", type = "upper", tl.col = "black", addCoef.col = "black")
  dev.off()

  write.csv(stats, file.path(output_dir, "T1_7_estadisticas.csv"), row.names = FALSE)
  write.csv(pca_res, file.path(output_dir, "T1_7_pca_resumen.csv"), row.names = FALSE)
  write.csv(as.data.frame(corr), file.path(output_dir, "T1_7_correlacion.csv"))

  list(data = df, descriptivos = stats, correlacion = corr, outliers = out, pca = pca, pca_resumen = pca_res)
}

# ------------------------------------------------------------
# 5) Tabla 1.8 - Mineral Content in Bones
# ------------------------------------------------------------

analizar_t18_bones <- function(
  path = "T1.8.csv",
  output_dir = "salidas_multivariado"
) {
  crear_carpeta_salida(output_dir)

  df <- read.csv(path, sep = ";")
  df_num <- df

  stats <- estadisticas_descriptivas(df_num)
  corr <- matriz_correlacion(df_num)
  out <- outliers_mahalanobis(df_num)
  pca <- preparar_pca(df_num)
  pca_res <- resumen_pca(pca)

  # Diferencias dominante - no dominante
  df_diff <- tibble(
    radius_diff = df$dominant_radius - df$radius,
    humerus_diff = df$dominant_humerus - df$humerus,
    ulna_diff = df$dominant_ulna - df$ulna
  )

  p_diff <- df_diff %>%
    pivot_longer(everything(), names_to = "hueso", values_to = "diferencia") %>%
    ggplot(aes(x = hueso, y = diferencia)) +
    geom_boxplot(fill = "gray80") +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    labs(title = "Tabla 1.8 - Diferencias dominante vs no dominante") +
    theme_minimal()
  plot_guardar(p_diff, file.path(output_dir, "T1_8_diferencias_boxplot.png"))

  p_pca <- fviz_pca_biplot(
    pca,
    repel = TRUE,
    title = "Tabla 1.8 - Biplot PCA"
  )
  plot_guardar(p_pca, file.path(output_dir, "T1_8_pca_biplot.png"))

  png(file.path(output_dir, "T1_8_corrplot.png"), width = 1000, height = 800, res = 150)
  corrplot(corr, method = "color", type = "upper", tl.col = "black", addCoef.col = "black")
  dev.off()

  write.csv(stats, file.path(output_dir, "T1_8_estadisticas.csv"), row.names = FALSE)
  write.csv(pca_res, file.path(output_dir, "T1_8_pca_resumen.csv"), row.names = FALSE)
  write.csv(as.data.frame(corr), file.path(output_dir, "T1_8_correlacion.csv"))
  write.csv(df_diff, file.path(output_dir, "T1_8_diferencias.csv"), row.names = FALSE)

  list(data = df, descriptivos = stats, correlacion = corr, outliers = out, pca = pca, pca_resumen = pca_res, diferencias = df_diff)
}

# ------------------------------------------------------------
# 6) Tabla 1.9 - National Track Records for Women
# ------------------------------------------------------------

analizar_t19_track_records <- function(
  path = "T1.9.csv",
  output_dir = "salidas_multivariado"
) {
  crear_carpeta_salida(output_dir)

  df <- read.csv(path, sep = ";")
  rownames(df) <- df$COUNTRY
  df_num <- df %>% select(-COUNTRY)

  stats <- estadisticas_descriptivas(df_num)
  corr <- matriz_correlacion(df_num)
  out <- outliers_mahalanobis(df_num)
  pca <- preparar_pca(df_num)
  pca_res <- resumen_pca(pca)

  # Clustering jerárquico sobre datos escalados
  dist_mat <- dist(scale(df_num))
  hc <- hclust(dist_mat, method = "ward.D2")

  png(file.path(output_dir, "T1_9_dendrograma.png"), width = 1200, height = 800, res = 150)
  plot(hc, labels = df$COUNTRY, main = "Tabla 1.9 - Dendrograma", xlab = "", sub = "")
  dev.off()

  p_pca <- fviz_pca_biplot(
    pca,
    repel = TRUE,
    label = "var",
    title = "Tabla 1.9 - Biplot PCA"
  )
  plot_guardar(p_pca, file.path(output_dir, "T1_9_pca_biplot.png"), width = 9, height = 7)

  png(file.path(output_dir, "T1_9_corrplot.png"), width = 1000, height = 800, res = 150)
  corrplot(corr, method = "color", type = "upper", tl.col = "black", addCoef.col = "black")
  dev.off()

  write.csv(stats, file.path(output_dir, "T1_9_estadisticas.csv"), row.names = FALSE)
  write.csv(pca_res, file.path(output_dir, "T1_9_pca_resumen.csv"), row.names = FALSE)
  write.csv(as.data.frame(corr), file.path(output_dir, "T1_9_correlacion.csv"))

  list(data = df, descriptivos = stats, correlacion = corr, outliers = out, pca = pca, pca_resumen = pca_res, clustering = hc)
}

# ------------------------------------------------------------
# Función general para ejecutar todos los análisis
# ------------------------------------------------------------

analizar_todo <- function(
  dir_datos = ".",
  output_dir = "salidas_multivariado"
) {
  crear_carpeta_salida(output_dir)

  resultados <- list(
    T1_4 = analizar_t14_companies(file.path(dir_datos, "T1.4.csv"), output_dir),
    T1_5 = analizar_t15_air_pollution(file.path(dir_datos, "T1.5.csv"), output_dir),
    T1_6 = analizar_t16_multiple_sclerosis(file.path(dir_datos, "T1.6.csv"), output_dir),
    T1_7 = analizar_t17_radiotherapy(file.path(dir_datos, "T1.7.csv"), output_dir),
    T1_8 = analizar_t18_bones(file.path(dir_datos, "T1.8.csv"), output_dir),
    T1_9 = analizar_t19_track_records(file.path(dir_datos, "T1.9.csv"), output_dir)
  )

  message("Análisis completado. Revisa la carpeta: ", output_dir)
  invisible(resultados)
}

# ------------------------------------------------------------
# Ejemplos de uso
# ------------------------------------------------------------
# setwd("ruta/donde/estan/los/csv")
# source("scripts_multivariado_tablas_1_4_a_1_9.R")
#
# resultado_14 <- analizar_t14_companies("T1.4.csv")
# resultado_15 <- analizar_t15_air_pollution("T1.5.csv")
# resultado_16 <- analizar_t16_multiple_sclerosis("T1.6.csv")
# resultado_17 <- analizar_t17_radiotherapy("T1.7.csv")
# resultado_18 <- analizar_t18_bones("T1.8.csv")
# resultado_19 <- analizar_t19_track_records("T1.9.csv")
#
# O ejecutar todo:
# resultados <- analizar_todo(dir_datos = ".", output_dir = "salidas_multivariado")
# ------------------------------------------------------------
