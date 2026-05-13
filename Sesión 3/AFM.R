# ================================================================
# Análisis Factorial (clásico vs robusto) con datos USArrests
# ================================================================

# --- Paquetes necesarios ---
library(psych)        # fa(), KMO, Bartlett
library(GPArotation)  # rotaciones
library(rrcov)        # CovMcd (robusto)
library(robustbase)   # covOGK (robusto rápido)
library(corrplot)     # visualización
library(tidyverse)

# --- 1. Cargar y preparar los datos ---
data("USArrests")
head(USArrests)
X <- scale(USArrests)     # estandarizar

# --- 2. Evaluar adecuación para AF ---
KMO(cor(X))               # > 0.6 recomendable
cortest.bartlett(cor(X), n = nrow(X))

# --- 3. Scree clásico vs robusto ---
R_class <- cor(X)
R_rob   <- cov2cor(covOGK(X)$cov)

ev_c <- eigen(R_class)$values
ev_r <- eigen(R_rob)$values

par(mfrow=c(1,2))
plot(ev_c, type="b", main="Scree clásico", ylab="Autovalores", xlab="Componente")
abline(h=1, col="red", lty=2)
plot(ev_r, type="b", main="Scree robusto (OGK)", ylab="Autovalores", xlab="Componente")
abline(h=1, col="red", lty=2)

# --- 4. Elegir número de factores (2) ---
m <- 2

# --- 5. AF clásico ---
fa_class <- fa(R_class, nfactors = m, rotate = "varimax", fm = "ml", n.obs = nrow(X))
print(fa_class, digits=2, sort=TRUE)
fa.diagram(fa_class)

# --- 6. AF robusto ---
fa_rob <- fa(R_rob, nfactors = m, rotate = "varimax", fm = "ml", n.obs = nrow(X))
print(fa_rob, digits=2, sort=TRUE)
fa.diagram(fa_rob)

# --- 7. Comparación de cargas ---
load_c <- as.data.frame(fa_class$loadings[1:4,])
load_r <- as.data.frame(fa_rob$loadings[1:4,])
comp <- data.frame(
  Variable = rownames(load_c),
  Clasico_F1 = load_c$ML1, Robusto_F1 = load_r$ML1,
  Clasico_F2 = load_c$ML2, Robusto_F2 = load_r$ML2
)
print(comp)

# --- 8. Comunalidades y residuos ---
res_class <- R_class - (fa_class$loadings %*% t(fa_class$loadings) + diag(fa_class$uniquenesses))
res_rob   <- R_rob   - (fa_rob$loadings   %*% t(fa_rob$loadings)   + diag(fa_rob$uniquenesses))
cat("Norma de residuos (clásico):", sqrt(sum(res_class^2)), "\n")
cat("Norma de residuos (robusto):", sqrt(sum(res_rob^2)), "\n")

par(mfrow=c(1,2))
corrplot(res_class, is.corr=FALSE, main="Residuos clásico")
corrplot(res_rob, is.corr=FALSE, main="Residuos robusto")

# --- 9. Interpretación ---
# Factor 1: Violencia (Assault, Murder, UrbanPop)
# Factor 2: Delitos sexuales / dispersión (Rape, Murder)
