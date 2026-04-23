# =========================================================
# Análisis multivariado de las Tablas 1.4 a 1.9
# Versión adaptada a archivos .txt
# Corrige lectura de números para evitar que queden como string
# =========================================================

# -----------------------------
# 0. Paquetes requeridos
# -----------------------------
required_packages <- c(
  "ggplot2", "GGally", "dplyr", "tidyr", "reshape2", "gridExtra"
)

install_if_missing <- function(pkgs) {
  faltantes <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(faltantes) > 0) {
    message("Instalando paquetes faltantes: ", paste(faltantes, collapse = ", "))
    install.packages(faltantes)
  }
}

install_if_missing(required_packages)

library(ggplot2)
library(GGally)
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)

# -----------------------------
# 1. Funciones auxiliares
# -----------------------------

# Convierte columnas tipo texto/factor a numéricas cuando sea posible.
# Soporta:
# - decimales con coma: "108,28"
# - separador de miles con punto
# - espacios sobrantes
# - guiones bajos en nombres de texto (los conserva)
convert_numeric_like <- function(df) {
  out <- df

  for (j in seq_along(out)) {
    x <- out[[j]]

    if (is.factor(x)) x <- as.character(x)

    if (is.character(x)) {
      x_trim <- trimws(x)

      # Intento 1: decimal con punto
      suppressWarnings(num1 <- as.numeric(x_trim))

      # Intento 2: decimal con coma -> reemplazar coma por punto
      x_comma <- gsub("\\.", "", x_trim)      # quita miles 1.234,56 -> 1234,56
      x_comma <- gsub(",", ".", x_comma)      # 1234,56 -> 1234.56
      suppressWarnings(num2 <- as.numeric(x_comma))

      # Elegimos la conversión que produzca más valores numéricos válidos
      valid1 <- sum(!is.na(num1))
      valid2 <- sum(!is.na(num2))

      # Solo convertir si al menos la mitad de la columna parece numérica
      threshold <- ceiling(length(x_trim) * 0.5)

      if (max(valid1, valid2) >= threshold) {
        out[[j]] <- if (valid2 > valid1) num2 else num1
      } else {
        out[[j]] <- x_trim
      }
    }
  }

  out
}

get_numeric_df <- function(df) {
  df %>% dplyr::select(where(is.numeric))
}

resumen_descriptivo <- function(df_num) {
  data.frame(
    Variable = names(df_num),
    Media = sapply(df_num, mean, na.rm = TRUE),
    Mediana = sapply(df_num, median, na.rm = TRUE),
    DesvStd = sapply(df_num, sd, na.rm = TRUE),
    Varianza = sapply(df_num, var, na.rm = TRUE),
    Minimo = sapply(df_num, min, na.rm = TRUE),
    Maximo = sapply(df_num, max, na.rm = TRUE),
    Rango = sapply(df_num, function(x) diff(range(x, na.rm = TRUE))),
    IQR = sapply(df_num, IQR, na.rm = TRUE),
    row.names = NULL
  )
}

correlacion_numerica <- function(df_num) {
  cor(df_num, use = "pairwise.complete.obs")
}

mahalanobis_outliers <- function(df_num, alpha = 0.975) {
  datos_limpios <- na.omit(df_num)
  centro <- colMeans(datos_limpios)
  covarianza <- cov(datos_limpios)
  d2 <- mahalanobis(datos_limpios, center = centro, cov = covarianza)
  cutoff <- qchisq(alpha, df = ncol(datos_limpios))

  data.frame(
    fila = as.integer(rownames(datos_limpios)),
    distancia_mahalanobis = d2,
    outlier = d2 > cutoff
  )
}

ejecutar_pca <- function(df_num) {
  datos_limpios <- na.omit(df_num)
  pca <- prcomp(datos_limpios, scale. = TRUE)

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  data.frame(
    Componente = paste0("PC", seq_along(var_exp)),
    Varianza_Explicada = var_exp,
    Varianza_Acumulada = cumsum(var_exp)
  ) -> tabla_var

  list(modelo = pca, varianza = tabla_var)
}

plot_correlation_heatmap <- function(cor_mat, title = "Matriz de correlación") {
  cor_long <- reshape2::melt(cor_mat)
  ggplot(cor_long, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
    theme_minimal() +
    labs(title = title, x = "", y = "")
}

plot_pca_individuals <- function(pca_obj, title = "PCA - individuos", grupos = NULL) {
  scores <- as.data.frame(pca_obj$x[, 1:2, drop = FALSE])
  if (!is.null(grupos) && length(grupos) == nrow(scores)) {
    scores$Grupo <- as.factor(grupos)
    ggplot(scores, aes(PC1, PC2, color = Grupo)) +
      geom_point(size = 2.5, alpha = 0.85) +
      theme_minimal() +
      labs(title = title)
  } else {
    ggplot(scores, aes(PC1, PC2)) +
      geom_point(size = 2.5, alpha = 0.85) +
      theme_minimal() +
      labs(title = title)
  }
}

plot_boxplots <- function(df_num, title = "Boxplots") {
  df_long <- df_num %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Valor")

  ggplot(df_long, aes(x = Variable, y = Valor)) +
    geom_boxplot(fill = "grey80") +
    theme_minimal() +
    labs(title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_pairs_base <- function(df_num, main = "Matriz de dispersión") {
  pairs(df_num, main = main, pch = 19, cex = 0.7)
}

guardar_resultados_basicos <- function(nombre, resumen, cor_mat, outliers, pca_var, output_dir = "salidas_r") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(resumen, file.path(output_dir, paste0(nombre, "_resumen.csv")), row.names = FALSE)
  write.csv(cor_mat, file.path(output_dir, paste0(nombre, "_correlacion.csv")))
  write.csv(outliers, file.path(output_dir, paste0(nombre, "_outliers_mahalanobis.csv")), row.names = FALSE)
  write.csv(pca_var, file.path(output_dir, paste0(nombre, "_pca_varianza.csv")), row.names = FALSE)
}

ejecutar_analisis_generico <- function(datos, nombre, grupo = NULL, output_dir = "salidas_r",
                                       grafico_extra = NULL) {
  datos <- convert_numeric_like(datos)
  df_num <- get_numeric_df(datos)

  if (ncol(df_num) < 2) stop("No hay suficientes columnas numéricas en ", nombre)

  resumen <- resumen_descriptivo(df_num)
  cor_mat <- correlacion_numerica(df_num)
  outliers <- mahalanobis_outliers(df_num)
  pca_res <- ejecutar_pca(df_num)

  guardar_resultados_basicos(nombre, resumen, cor_mat, outliers, pca_res$varianza, output_dir)

  # Gráficos
  heat <- plot_correlation_heatmap(cor_mat, paste("Correlación -", nombre))
  pca_plot <- plot_pca_individuals(pca_res$modelo, paste("PCA -", nombre), grupos = grupo)
  boxp <- plot_boxplots(df_num, paste("Boxplots -", nombre))

  ggsave(file.path(output_dir, paste0(nombre, "_heatmap_correlacion.png")), heat, width = 7, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(nombre, "_pca.png")), pca_plot, width = 7, height = 5, dpi = 300)
  ggsave(file.path(output_dir, paste0(nombre, "_boxplots.png")), boxp, width = 8, height = 5, dpi = 300)

  if (!is.null(grafico_extra)) {
    grafico_extra(df_num, datos, nombre, output_dir)
  }

  list(
    datos = datos,
    datos_numericos = df_num,
    resumen = resumen,
    correlacion = cor_mat,
    outliers_mahalanobis = outliers,
    pca = pca_res
  )
}

# -----------------------------
# 2. Lectores específicos .txt
# -----------------------------

leer_T14 <- function(archivo = "T1.4.txt") {
  # company  sales  profits  assets
  # citigroup 108,28 17,05 1484,1
  datos <- read.table(
    archivo,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )
  convert_numeric_like(datos)
}

leer_T15 <- function(archivo = "T1.5.txt") {
  # Sin encabezados en el archivo
  datos <- read.table(
    archivo,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(datos) <- c("Wind", "Radiation", "CO", "NO", "NO2", "O3", "SO2")
  convert_numeric_like(datos)
}

leer_T16 <- function(archivo = "T1.6.txt") {
  datos <- read.table(
    archivo,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(datos) <- c("X1", "X2", "X3", "X4", "X5", "Grupo")
  convert_numeric_like(datos)
}

leer_T17 <- function(archivo = "T1.7.txt") {
  datos <- read.table(
    archivo,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(datos) <- c("X1", "X2", "X3", "X4", "X5", "Grupo")
  convert_numeric_like(datos)
}

leer_T18 <- function(archivo = "T1.8.txt") {
  datos <- read.table(
    archivo,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(datos) <- c("X1", "X2", "X3", "X4", "X5", "X6")
  convert_numeric_like(datos)
}

leer_T19 <- function(archivo = "T1.9.txt") {
  datos <- read.table(
    archivo,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )
  convert_numeric_like(datos)
}

# -----------------------------
# 3. Funciones específicas por tabla
# -----------------------------

analizar_T1_4_companies <- function(archivo = "T1.4.txt", output_dir = "salidas_r") {
  datos <- leer_T14(archivo)

  grafico_extra <- function(df_num, datos_originales, nombre, output_dir) {
    g <- GGally::ggpairs(
      df_num,
      upper = list(continuous = GGally::wrap("cor", size = 4, stars = TRUE)),
      lower = list(continuous = GGally::wrap("points", alpha = 0.8, size = 1.8)),
      diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5))
    ) + theme_bw()

    ggsave(file.path(output_dir, paste0(nombre, "_ggpairs.png")), g, width = 8, height = 8, dpi = 300)
  }

  ejecutar_analisis_generico(datos, "T1_4_companies", output_dir = output_dir, grafico_extra = grafico_extra)
}

analizar_T1_5_air_pollution <- function(
  archivo = "T1.5.txt",
  output_dir = "salidas_r",
  variables_grafico = c("Wind", "Radiation", "CO", "NO", "NO2", "O3")
) {
  datos <- leer_T15(archivo)

  grafico_extra <- function(df_num, datos_originales, nombre, output_dir) {
    vars_existentes <- variables_grafico[variables_grafico %in% names(df_num)]
    datos_plot <- df_num[, vars_existentes, drop = FALSE]

    g <- GGally::ggpairs(
      datos_plot,
      upper = list(
        continuous = GGally::wrap("cor", size = 4, stars = TRUE)
      ),
      lower = list(
        continuous = GGally::wrap("points", alpha = 0.9, size = 1.8)
      ),
      diag = list(
        continuous = GGally::wrap("densityDiag", alpha = 0.5)
      )
    ) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "grey85", colour = "grey40"),
        strip.text = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2)
      )

    ggsave(file.path(output_dir, paste0(nombre, "_ggpairs_air_pollution.png")), g,
           width = 10, height = 10, dpi = 300)
  }

  ejecutar_analisis_generico(datos, "T1_5_air_pollution", output_dir = output_dir, grafico_extra = grafico_extra)
}

analizar_T1_6_multiple_sclerosis <- function(archivo = "T1.6.txt", output_dir = "salidas_r") {
  datos <- leer_T16(archivo)
  grupo <- datos$Grupo
  datos_sin_grupo <- datos %>% select(-Grupo)

  grafico_extra <- function(df_num, datos_originales, nombre, output_dir) {
    g <- GGally::ggpairs(
      cbind(df_num, Grupo = as.factor(grupo)),
      columns = 1:ncol(df_num),
      mapping = aes(color = Grupo, alpha = 0.8),
      upper = list(continuous = GGally::wrap("cor", size = 3)),
      lower = list(continuous = GGally::wrap("points", size = 1.6)),
      diag = list(continuous = GGally::wrap("densityDiag"))
    ) + theme_bw()

    ggsave(file.path(output_dir, paste0(nombre, "_ggpairs.png")), g, width = 9, height = 9, dpi = 300)
  }

  ejecutar_analisis_generico(datos_sin_grupo, "T1_6_multiple_sclerosis", grupo = grupo,
                             output_dir = output_dir, grafico_extra = grafico_extra)
}

analizar_T1_7_radiotherapy <- function(archivo = "T1.7.txt", output_dir = "salidas_r") {
  datos <- leer_T17(archivo)
  grupo <- datos$Grupo
  datos_sin_grupo <- datos %>% select(-Grupo)

  grafico_extra <- function(df_num, datos_originales, nombre, output_dir) {
    p <- GGally::ggpairs(
      cbind(df_num, Grupo = as.factor(grupo)),
      columns = 1:ncol(df_num),
      mapping = aes(color = Grupo, alpha = 0.75),
      upper = list(continuous = GGally::wrap("cor", size = 3)),
      lower = list(continuous = GGally::wrap("points", size = 1.5)),
      diag = list(continuous = GGally::wrap("densityDiag"))
    ) + theme_bw()

    ggsave(file.path(output_dir, paste0(nombre, "_ggpairs.png")), p, width = 9, height = 9, dpi = 300)
  }

  ejecutar_analisis_generico(datos_sin_grupo, "T1_7_radiotherapy", grupo = grupo,
                             output_dir = output_dir, grafico_extra = grafico_extra)
}

analizar_T1_8_mineral_bones <- function(archivo = "T1.8.txt", output_dir = "salidas_r") {
  datos <- leer_T18(archivo)

  grafico_extra <- function(df_num, datos_originales, nombre, output_dir) {
    g <- GGally::ggpairs(
      df_num,
      upper = list(continuous = GGally::wrap("cor", size = 3.5)),
      lower = list(continuous = GGally::wrap("smooth_loess", alpha = 0.35, size = 0.3)),
      diag = list(continuous = GGally::wrap("densityDiag"))
    ) + theme_bw()

    ggsave(file.path(output_dir, paste0(nombre, "_ggpairs.png")), g, width = 9, height = 9, dpi = 300)
  }

  ejecutar_analisis_generico(datos, "T1_8_mineral_bones", output_dir = output_dir, grafico_extra = grafico_extra)
}

analizar_T1_9_track_records <- function(archivo = "T1.9.txt", output_dir = "salidas_r") {
  datos <- leer_T19(archivo)
  paises <- datos[[1]]
  datos_num <- datos %>% select(-1)

  grafico_extra <- function(df_num, datos_originales, nombre, output_dir) {
    cor_mat <- cor(df_num, use = "pairwise.complete.obs")
    cor_long <- melt(cor_mat)

    g <- ggplot(cor_long, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 2)), size = 2.7) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
      theme_minimal() +
      labs(title = "Heatmap de correlación - marcas femeninas", x = "", y = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(output_dir, paste0(nombre, "_heatmap_detallado.png")), g,
           width = 9, height = 7, dpi = 300)
  }

  res <- ejecutar_analisis_generico(datos_num, "T1_9_track_records", output_dir = output_dir,
                                    grafico_extra = grafico_extra)
  res$paises <- paises
  res
}

# -----------------------------
# 4. Función para ejecutar todo
# -----------------------------
analizar_todo <- function(
  archivo_t14 = "T1.4.txt",
  archivo_t15 = "T1.5.txt",
  archivo_t16 = "T1.6.txt",
  archivo_t17 = "T1.7.txt",
  archivo_t18 = "T1.8.txt",
  archivo_t19 = "T1.9.txt",
  output_dir = "salidas_r"
) {
  list(
    T1_4 = analizar_T1_4_companies(archivo_t14, output_dir),
    T1_5 = analizar_T1_5_air_pollution(archivo_t15, output_dir),
    T1_6 = analizar_T1_6_multiple_sclerosis(archivo_t16, output_dir),
    T1_7 = analizar_T1_7_radiotherapy(archivo_t17, output_dir),
    T1_8 = analizar_T1_8_mineral_bones(archivo_t18, output_dir),
    T1_9 = analizar_T1_9_track_records(archivo_t19, output_dir)
  )
}

# -----------------------------
# 5. Ejemplos de uso
# -----------------------------
# resultados <- analizar_todo()
# resultado_t15 <- analizar_T1_5_air_pollution("T1.5.txt")
# print(resultado_t15$resumen)
# print(round(resultado_t15$correlacion, 3))
