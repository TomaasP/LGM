# Cargar paquetes necesarios
library(tidyverse)

# Leer el archivo de cobertura generado por bedtools
coverage_data <- read.table("HB44IC3.coverage_output.txt", header = FALSE, sep = "\t",dec = ",", 
                            col.names = c("chr","spp", "start", "end", "depth","#BasesCubiertas","fraccionVentana", "coverage"))
# Ver los primeros registros para confirmar la estructura
head(coverage_data)

# Calcular la cobertura promedio por cromosoma
coverage_summary <- coverage_data %>%
  group_by(chr,spp) %>%
  summarize(mean_coverage = mean(depth, na.rm = TRUE)) %>%
  arrange(desc(mean_coverage))

# Ver los resultados
print(coverage_summary)


# Crear un gráfico de barras con ggplot2
p1 <- coverage_summary %>% filter(chr != "chrMT") %>%
  ggplot(aes(x = reorder(chr, mean_coverage), y = mean_coverage, fill = spp)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Rotar el gráfico si hay muchos cromosomas
  labs(
    title = "Cobertura promedio por cromosoma",
    x = "Cromosoma",
    y = "Cobertura promedio"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
p1

ggsave("HB44IC3.png", p1, width = 10, height = 8, dpi = 300)

# Seleccionar un cromosoma específico, por ejemplo, "chr1"
chromosome_data <- coverage_data %>%
  filter(chr == "3CBS12357_Chr03_polished:1-308584")
chromosome_data %>%
  ggplot(aes(x = start, y = depth)) +
  geom_point(color = "blue") +
  labs(
    title = "Cobertura a lo largo",
    x = "Posición (bp)",
    y = "Cobertura"
  ) +
  theme_minimal()

# Graficar la cobertura en ventanas de 1 kb para "chr1"
# Crear una nueva columna que combine especie y cromosoma
coverage_data <- coverage_data %>%
  mutate(spp_chr = paste(spp, chr, sep = "_"))

# Graficar con facetas alineadas en el mismo plano
coverage_data %>%
  ggplot(aes(x = start, y = coverage, color = spp)) +  # Colores por especie
  geom_line() +
  labs(
    title = "Cobertura de secuenciación en el híbrido Se x Sc",
    x = "Posición (bp)",
    y = "Cobertura"
  ) +
  ylim(0, 2) +
  facet_wrap(~ spp_chr, scales = "free_x", nrow = 1) +  # Una sola fila de paneles
  theme_minimal() +
  theme(
    strip.text.x = element_text(angle = 90),  # Rotar etiquetas de paneles
    legend.position = "top"                  # Colocar la leyenda arriba
  )
