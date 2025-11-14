# MarthaLeivaSanchez-Viñas_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)
install.packages(c("tidyverse"))
library(tidyverse)
datos <- read.csv("datos_biomed.csv", header = TRUE, sep = ",")

# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)
head(datos)
summary(datos)
dim(datos)
str(datos)
ncol(datos)
length(unique(datos$Tratamiento))

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)
ggplot(datos, aes(x = Tratamiento, y = Glucosa, fill = Tratamiento)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots de Glucosa por Tratamiento")
  
# 4. Realiza un violin plot (investiga qué es). (1 pt)
ggplot(datos, aes(x = Tratamiento, y = Glucosa, fill = Tratamiento)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_minimal()

# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)
plot(datos$Glucosa, datos$Presion,
     col = as.factor(datos$Tratamiento), pch = 19,
     xlab = "Glucosa", ylab = "Presión",
     main = "Glucosa vs Presión por tratamiento")
legend("bottomright", legend = levels(factor(datos$Tratamiento)), col = 1:length(levels(factor(datos$Tratamiento))), pch = 19)

# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)
ggplot(datos, aes(x = Colesterol, y = Presion)) +
  geom_point() +
  facet_grid(~ Tratamiento) +
  theme_minimal()

# 7. Realiza un histogramas para cada variable. (0.5 pts)
datos_long <- datos %>%
  pivot_longer(cols = -Tratamiento, names_to = "Variable", values_to = "Valor")
datos_long %>%
  ggplot(aes(x = Valor, fill = Variable)) +
  geom_histogram() +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()

# 8. Crea un factor a partir del tratamiento. Investifa factor(). (1 pt)
datos$Tratamiento <- factor(datos$Tratamiento)
str(datos$Tratamiento)

# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)
aggregate(Glucosa ~ Tratamiento, datos, mean)
aggregate(Glucosa ~ Tratamiento, datos, sd)

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)
placebo  <- subset(datos, Tratamiento == "Placebo")
farmacoA <- subset(datos, Tratamiento == "FarmacoA")
farmacoB <- subset(datos, Tratamiento == "FarmacoB")
#Comprobación de las variables:
head(placebo)
head(farmacoA)
head(farmacoB)

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)
shapiro.test(datos$Glucosa[datos$Tratamiento == "Placebo"])
shapiro.test(datos$Glucosa[datos$Tratamiento == "FarmacoA"])
shapiro.test(datos$Glucosa[datos$Tratamiento == "FarmacoB"])
kruskal.test(Glucosa ~ Tratamiento, data = datos)

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)
anova_glucosa <- aov(Glucosa ~ Tratamiento, data = datos)
summary(anova_glucosa)

TukeyHSD(anova_glucosa)


