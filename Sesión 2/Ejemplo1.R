
# Instalar si no lo tienes
# install.packages("keras")
# install.packages("ggplot2")
# install.packages("dplyr")

# Cargar librerías
library(keras)
library(ggplot2)
library(dplyr)

# Cargar el dataset MNIST
# Contiene 60,000 imágenes de entrenamiento y 10,000 de prueba
mnist <- dataset_mnist()

# Extraemos los datos de entrenamiento
x_train <- mnist$train$x
y_train <- mnist$train$y

# Vistazo a las dimensiones
# 60000 imágenes, cada una de 28x28 píxeles
print(dim(x_train))

# Reshape: Convertir cada matriz de 28x28 en un vector de 784
# Dividimos por 255 para escalar los valores de los píxeles entre 0 y 1
x_train_matrix <- array_reshape(x_train, c(nrow(x_train), 784)) / 255

# Aplicamos PCA. 'prcomp' es la función estrella aquí.
# scale. = FALSE porque ya hemos escalado los datos.
# La función centra los datos por defecto (center = TRUE)
pca_result <- prcomp(x_train_matrix, scale. = FALSE)

# Resumen del resultado
# summary(pca_result) # Descomentar para ver la varianza explicada por cada componente

# --- Visualicemos un dígito para confirmar ---
# Función para graficar un dígito desde un vector
plot_digit <- function(vector_784) {
  image_matrix <- matrix(vector_784, nrow = 28, ncol = 28, byrow = TRUE)
  image(t(apply(image_matrix, 2, rev)), col = grey.colors(255))
}

# Graficamos el dígito de la quinta fila (es un '9')
plot_digit(x_train_matrix[5, ])

# Calculamos la varianza explicada por cada componente
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_variance <- cumsum(variance_explained)

# Creamos un data frame para graficar con ggplot2
scree_data <- data.frame(
  component = 1:length(variance_explained),
  variance_explained = variance_explained,
  cumulative_variance = cumulative_variance
)

# Graficamos la varianza acumulada
ggplot(scree_data, aes(x = component, y = cumulative_variance)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "red") +
  labs(
    title = "Varianza Acumulada Explicada por los Componentes Principales",
    x = "Número de Componentes Principales",
    y = "Varianza Acumulada Explicada"
  ) +
  theme_minimal() +
  ylim(0, 1)


# Los componentes están en pca_result$rotation
# Visualicemos los primeros 4 componentes
par(mfrow = c(2, 2)) # Prepara el layout para 4 gráficos
plot_digit(pca_result$rotation[, 1])
title("PC 1")
plot_digit(pca_result$rotation[, 2])
title("PC 2")
plot_digit(pca_result$rotation[, 3])
title("PC 3")
plot_digit(pca_result$rotation[, 4])
title("PC 4")
par(mfrow = c(1, 1)) # Restaura el layout


# Número de componentes a usar para la reconstrucción
k <- 50 

# Proyección y Reconstrucción
# scores (proyección) = datos_centrados %*% rotación
# reconstrucción = (scores %*% t(rotación)) + centro
# prcomp nos da los scores en pca_result$x

# Tomamos una imagen (ej, la número 20, que es un '4')
original_image <- x_train_matrix[20, ]
image_scores <- pca_result$x[20, 1:k]

# Reconstruimos usando solo k componentes
reconstructed_vector <- image_scores %*% t(pca_result$rotation[, 1:k]) + pca_result$center

# Graficamos la original vs. la reconstruida
par(mfrow = c(1, 2))
plot_digit(original_image)
title("Original")
plot_digit(reconstructed_vector)
title(paste("Reconstruida con", k, "Componentes"))
par(mfrow = c(1, 1))
