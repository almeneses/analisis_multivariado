# 1) Cargar y explorar
data(iris)
df <- iris[, 1:4]                   # solo variables numéricas
summary(df)
pairs(df, col = "grey70", pch = 16)

# 2) Escalado (recomendado para k-means)
X <- scale(df)

# 3) Elección de k: codo + silhouette
set.seed(123)
wss <- sapply(1:10, function(k){
  kmeans(X, centers = k, nstart = 25)$tot.withinss
})
plot(1:10, wss, type="b", pch=19, xlab="k", ylab="WSS (within-cluster sum of squares)",
     main="Método del codo")

# silhouette promedio para k = 2..10
library(cluster)
sil <- sapply(2:10, function(k){
  km <- kmeans(X, centers = k, nstart = 25)
  mean(silhouette(km$cluster, dist(X))[, 3])
})
plot(2:10, sil, type="b", pch=19, xlab="k", ylab="Silhouette promedio",
     main="Calidad de partición vs k")
k_opt <- which.max(sil) + 1; k_opt

# 4) Ajustar k-means con k óptimo (suele salir 3 en iris)
set.seed(123)
km <- kmeans(X, centers = k_opt, nstart = 50)
table(km$cluster)

# 5) Validación e interpretación
# 5.1) Silhouette
sil_km <- silhouette(km$cluster, dist(X))
plot(sil_km, main = sprintf("Silhouette K-means (k=%d)", k_opt))

# 5.2) Centroides (en escala original para interpretar)
centros_std <- km$centers
centros_orig <- sweep(centros_std, 2, attr(X, "scaled:scale"), "*")
centros_orig <- sweep(centros_orig, 2, attr(X, "scaled:center"), "+")
round(centros_orig, 2)

# 5.3) Evaluación externa (solo docente): ¿cómo se alinean con Species?
tab <- table(Cluster = km$cluster, Species = iris$Species)
tab
prop.table(tab, 2)  # proporción por especie

# 5.4) Visualización en 2D vía PCA
pca <- prcomp(X, center = FALSE, scale. = FALSE)
plot(pca$x[,1], pca$x[,2], col = km$cluster, pch = 19,
     xlab = "PC1", ylab = "PC2", main = "K-means proyectado a PCA")
legend("topright", legend = paste("C", 1:k_opt), col = 1:k_opt, pch = 19, bty = "n")

# 6) Alternativa: clustering jerárquico (Ward)
d <- dist(X, method = "euclidean")
hc <- hclust(d, method = "ward.D2")
plot(hc, cex = 0.6, main = "Dendrograma (Ward)")
rect.hclust(hc, k = k_opt, border = 2:5)
grupos_hc <- cutree(hc, k = k_opt)
table(grupos_hc, km$cluster)

# 7) (ideas):
# - Datos, variables y escalado
# - Selección de k (codo + silhouette promedio = …)
# - Resultados K-means: tamaños de clúster, centroides interpretables:
#   * p.ej., un clúster con pétalos muy cortos/ancho pequeño (perfil Setosa),
#     otro intermedio, otro con pétalos largos/anchos (perfil Virginica).
# - Métrica interna (silhouette promedio = mean(sil_km[,3]))
# - Comparación con Species (si aplica en docencia) y confusiones típicas en versicolor/virginica.
# - Visualizaciones (pairs, PCA, dendrograma)
