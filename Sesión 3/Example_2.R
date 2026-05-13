# --- 0. CONFIGURACIÓN: Cargar todas las librerías necesarias ---

# install.packages("psych")        # Para EFA (fa, fa.parallel, KMO)
# install.packages("GPArotation")  # Requerido por 'psych' para rotaciones

library(psych)
library(GPArotation)
library(caret)
library(readr)
data_features_60 <- read_csv("data_features_60.csv")

# --- PASO 2: Análisis Factorial Exploratorio (EFA) ---

cat("--- PASO 2: Ejecutando Análisis Factorial ---\n")

# 2a. Pruebas de adecuación (con N=2000 y p=60, serán perfectas)
kmo_results <- KMO(data_features_60)
cat(paste("Prueba KMO (Adecuación de Muestra):", round(kmo_results$MSA, 3), "\n"))
# cortest.bartlett(data_features_60) # p-value será ~0

# 2b. Determinar número de factores (lo saltamos, ya que sabemos que es 4)
# En un caso real, ejecutaríamos:
# fa.parallel(data_features_60, fm = "minres", fa = "fa")
n_factores_optimo <- 4

# 2c. Ejecutar el modelo EFA
# fm="minres": Mínimos Residuales (robusto)
# rotate="oblimin": Rotación oblicua (asumimos que los temas están correlacionados)
efa_model_nlp <- fa(data_features_60, 
                    nfactors = n_factores_optimo, 
                    rotate = "oblimin", 
                    fm = "minres")

# Imprimir un resumen del modelo.
# Veremos que las 60 variables cargan limpiamente en 4 factores.
cat("Resultados del Modelo EFA (cargas > 0.4):\n")
print(efa_model_nlp, cut = 0.4, digits = 2, sort = TRUE)


# --- PASO 3: Ingeniería de Características (Feature Engineering) ---

cat("--- PASO 3: Extrayendo Puntuaciones Factoriales ---\n")

# efa_model_nlp$scores contiene las 4 nuevas "features" latentes para cada ticket
data_features_4 <- as.data.frame(efa_model_nlp$scores)
colnames(data_features_4) <- c("F1_Latent", "F2_Latent", "F3_Latent", "F4_Latent")

cat(paste("Dataset 'data_features_4' creado. Dimensiones:", 
          nrow(data_features_4), "x", ncol(data_features_4), "\n\n"))
print(head(data_features_4))


