
library(ggplot2)
library(tidyr)
library(dplyr)

datos <- read.csv("expresion_tumores.csv")
datos <- datos[,-1]

datos_tidy <- gather(data = datos, key = tumor, value = expresion, -1)
ggplot(data = datos_tidy, aes(x = gene, y = expresion, color = tumor)) +
  geom_path(aes(group = tumor)) +
  geom_point() +
  theme_bw() +
  facet_grid(tumor~.) +
  theme(legend.position = "none")


datos_tidy %>% group_by(gene) %>% summarise(varianza = var(expresion)) %>% 
  ggplot(aes(x = gene, y = varianza)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Dado que se quieren estudiar patrones empleando la expresión de los genes,
# cada gen debe de estar en una columna (predictores)

datos_trans <- datos_tidy %>% spread(key = gene, value = expresion)
rownames(datos_trans) <- datos_trans$tumor
datos_trans <- datos_trans[, -1]
head(datos_trans)

pca <- prcomp(datos_trans)

# Cálculo de la varianza explicada acumulada 
prop_varianza <- pca$sdev^2/sum(pca$sdev^2)
prop_varianza_acum <- cumsum(prop_varianza)
ggplot(data = data.frame(prop_varianza_acum, pc = factor(1:9)),
       aes(x = pc, y = prop_varianza_acum, group = 1)) +
  geom_point() +
  geom_line() +
  geom_label(aes(label = round(prop_varianza_acum,2))) +
  theme_bw() +
  labs(x = "Componentes principales", 
       y = "Prop. varianza explicada acumulada")


pca$rotation[, 1:2]

biplot(x = pca, scale = 0, cex = 0.7, col = c("blue4", "brown3"))

