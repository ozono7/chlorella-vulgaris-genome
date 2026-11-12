# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie: Artículo Edith — genoma de Chlorella vulgaris (LBCVB)
#
# Figuras del artículo para JGEB
# De novo assembly and functional annotation of sequences from a culture
# identified as Chlorella vulgaris: genetic potential relevant to wastewater
# bioremediation
#
# NOTA v3 (2026-11-11): las 6 figuras que genera este script YA ESTÁN
# GENERADAS y se incluyen en el paquete, en la carpeta figuras/ (PNG, 300 dpi).
# Este script es el material de reproducibilidad: permite verificar o
# regenerar esas figuras. NO es una tarea pendiente para nadie.
# Las figuras que faltan (2b, S1, S2, S3–S11) las produce Edith desde los
# archivos originales; ver lista-figuras-para-edith.md.
#
# Cambios v3: todas las etiquetas en inglés; SIN título ni subtítulo dentro
# de la imagen (el título va en el pie del manuscrito); estética sobria
# apta para revista (paleta compatible con visión de color deficiente,
# tamaños de letra mayores, nombres de género y especie en itálica en los
# ejes). Las etiquetas de los CSV también pasaron a inglés (son etiquetas
# categóricas, no mediciones).
#
# Qué hace este script
# Genera, a partir de los valores verificados de ensamblaje y anotación,
# las figuras que se pueden construir en R sin los archivos
# originales de salida de cada herramienta:
#   Figura 1 — Clasificación taxonómica Kraken2 (conteos por taxón)
#   Figura 2 — Identidad ANIm contra el panel de referencia (UTEX 395 sin valor)
#   Figura 3 — Evaluación BUSCO (viridiplantae_odb10, n = 425)
#   Figura 4 — Anotación estructural y funcional (Funannotate + eggNOG-mapper)
#   Figura 5 — Distribución COG (solo categorías con conteo reportado)
#   Figura S0 — Comparación QUAST de los cuatro ensambladores
#
# Datos
# Los CSV de esta carpeta contienen los valores verificados de ensamblaje
# y anotación que se usan para las figuras. Las categorías
# COG sin conteo reportado (O, P, C, J, I) no se grafican y se señalan en
# la Figura 5; el ANI de UTEX 395 no se grafica porque no fue reportado.
#

Sys.setlocale("LC_ALL", "C.UTF-8")

library(ggplot2)
library(dplyr)
library(readr)
library(scales)

# Carpeta del script: los CSV y las figuras viven junto a este archivo.
# Se detecta desde --file= cuando se corre con Rscript, o se usa el
# directorio de trabajo actual (correr el script desde su propia carpeta).
args <- commandArgs(trailingOnly = FALSE)
archivo <- grep("--file=", args, value = TRUE)
if (length(archivo) > 0) {
  carpeta <- dirname(sub("--file=", "", archivo[1]))
} else {
  carpeta <- getwd()
}
setwd(carpeta)

# Tema común: limpio, sin título dentro de la imagen, leyenda abajo,
# títulos de ejes en negrita.
tema_fig <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 11),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

# Paleta compatible con visión de color deficiente (Okabe-Ito).
azul <- "#0072B2"

# Etiquetas de taxones: género y especie en itálica, rangos superiores
# en romana. Devuelve expresiones plotmath, una por etiqueta.
etiqueta_taxon <- function(x) {
  lapply(x, function(et) {
    partes <- strsplit(et, " ", fixed = TRUE)[[1]]
    es_especie <- length(partes) >= 2 &&
      grepl("^[A-Z][a-z]+$", partes[1]) && grepl("^[a-z]+$", partes[2])
    es_abrev <- length(partes) >= 2 &&
      grepl("^[A-Z]\\.$", partes[1]) && grepl("^[a-z]+$", partes[2])
    if (es_especie || es_abrev) {
      nombre <- paste(partes[1], partes[2])
      resto <- paste(partes[-(1:2)], collapse = " ")
      if (nzchar(resto)) {
        parse(text = paste0('italic("', nombre, '")~"', resto, '"'))[[1]]
      } else {
        parse(text = paste0('italic("', nombre, '")'))[[1]]
      }
    } else if (length(partes) == 1 && et %in% c("Pseudomonas")) {
      parse(text = paste0('italic("', et, '")'))[[1]]
    } else {
      parse(text = paste0('"', et, '"'))[[1]]
    }
  })
}

# ---------------------------------------------------------------------------
# Figura 1. Clasificación taxonómica Kraken2
# ---------------------------------------------------------------------------
# Los conteos son de rangos anidados: el gráfico usa escala logarítmica
# porque los valores abarcan de 2,430 a 50,671,528 unidades asignadas.
kraken <- read_csv("datos_kraken.csv", show_col_types = FALSE) %>%
  mutate(taxon = factor(taxon, levels = rev(taxon[order(conteo)])))

colores_rango <- c(
  "Domain" = "#0072B2", "Clade" = "#56B4E9", "Phylum" = "#009E73",
  "Family" = "#E69F00", "Genus" = "#D55E00", "Species" = "#CC79A7",
  "Group" = "#999999"
)

figura1 <- ggplot(kraken, aes(x = conteo, y = taxon, fill = rango)) +
  geom_col() +
  scale_x_log10(
    labels = c("1,000", "10,000", "100,000", "1M", "10M", "50M"),
    breaks = c(1e3, 1e4, 1e5, 1e6, 1e7, 5e7)
  ) +
  scale_y_discrete(labels = etiqueta_taxon) +
  scale_fill_manual(values = colores_rango) +
  labs(
    x = "Assigned units (log10 scale)",
    y = NULL,
    fill = "Taxonomic rank"
  ) +
  tema_fig

ggsave("figura1_kraken.png", figura1, width = 9, height = 6, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 2. Identidad ANIm contra el panel de referencia
# ---------------------------------------------------------------------------
# UTEX 395 no tiene valor reportado y se excluye del gráfico; se señala en
# el pie. La línea en 95% marca el umbral de identidad habitualmente usado
# para delimitar especies bacterianas; aquí se muestra solo como referencia.
ani <- read_csv("datos_ani.csv", show_col_types = FALSE) %>%
  filter(estado == "resuelto") %>%
  mutate(referencia = factor(referencia, levels = referencia[order(ani)]))

colores_grupo <- c(
  "C. vulgaris" = "#0072B2",
  "Other Chlorella" = "#009E73",
  "Other genus" = "#E69F00"
)

figura2 <- ggplot(ani, aes(x = ani, y = referencia, fill = grupo)) +
  geom_col() +
  geom_vline(xintercept = 95, linetype = "dashed", color = "black") +
  annotate("text", x = 95.3, y = 1, label = "95% threshold", size = 3.5, hjust = 0) +
  coord_cartesian(xlim = c(80, 101)) +
  scale_x_continuous(breaks = seq(80, 100, 5)) +
  scale_y_discrete(labels = etiqueta_taxon) +
  scale_fill_manual(values = colores_grupo) +
  labs(
    x = "Average nucleotide identity (ANIm), %",
    y = NULL,
    fill = "Group"
  ) +
  tema_fig

ggsave("figura2_ani.png", figura2, width = 9, height = 6, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 3. Evaluación BUSCO (viridiplantae_odb10, n = 425)
# ---------------------------------------------------------------------------
busco <- read_csv("datos_busco.csv", show_col_types = FALSE) %>%
  mutate(
    categoria = factor(
      categoria,
      levels = c("Missing", "Fragmented", "Complete duplicated", "Complete single-copy")
    )
  )

figura3 <- ggplot(busco, aes(x = "", y = conteo, fill = categoria)) +
  geom_col(width = 0.5) +
  geom_text(
    aes(label = paste0(conteo, " (", porcentaje, "%)")),
    position = position_stack(vjust = 0.5), size = 4
  ) +
  scale_fill_manual(values = c(
    "Complete single-copy" = "#0072B2",
    "Complete duplicated" = "#56B4E9",
    "Fragmented" = "#E69F00",
    "Missing" = "#D55E00"
  )) +
  labs(
    x = NULL, y = "Orthologs (n = 425)",
    fill = "BUSCO category"
  ) +
  tema_fig

ggsave("figura3_busco.png", figura3, width = 8, height = 6, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 4. Anotación estructural y funcional
# ---------------------------------------------------------------------------
anotacion <- read_csv("datos_anotacion.csv", show_col_types = FALSE) %>%
  mutate(etapa = factor(etapa, levels = rev(etapa)))

figura4 <- ggplot(anotacion, aes(x = conteo, y = etapa)) +
  geom_col(fill = azul) +
  geom_text(aes(label = comma(conteo)), hjust = -0.1, size = 4) +
  scale_x_continuous(labels = comma, limits = c(0, 55000)) +
  labs(
    x = "Sequences",
    y = NULL
  ) +
  tema_fig

ggsave("figura4_anotacion.png", figura4, width = 9, height = 5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 5. Distribución COG
# ---------------------------------------------------------------------------
# Solo se grafican las categorías con conteo reportado. Las cinco categorías
# con porcentaje pero sin conteo (O, P, C, J, I) no se grafican: se
# señalan en el pie y en la Tabla 4 del manuscrito.
cog <- read_csv("datos_cog.csv", show_col_types = FALSE) %>%
  filter(estado == "resuelto") %>%
  mutate(categoria = factor(categoria, levels = categoria[order(conteo)]))

figura5 <- ggplot(cog, aes(x = conteo, y = categoria)) +
  geom_col(fill = azul) +
  geom_text(aes(label = paste0(conteo, " (", porcentaje, "%)")), hjust = -0.1, size = 3.5) +
  scale_x_continuous(labels = comma, limits = c(0, 10000)) +
  labs(
    x = "Assigned genes",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

ggsave("figura5_cog.png", figura5, width = 10, height = 7, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura S0 (suplementaria). Comparación QUAST de los cuatro ensambladores
# ---------------------------------------------------------------------------
quast <- read_csv("datos_quast.csv", show_col_types = FALSE) %>%
  mutate(ensamblador = factor(
    ensamblador,
    levels = c("ABySS", "MEGAHIT", "SPAdes", "Unicycler")
  ))

colores <- c(
  "Unicycler" = "#0072B2", "ABySS" = "#E69F00",
  "MEGAHIT" = "#D55E00", "SPAdes" = "#009E73"
)

tema_panel <- theme_minimal(base_size = 11) +
  theme(legend.position = "none")

panel_n50 <- ggplot(quast, aes(x = ensamblador, y = n50, fill = ensamblador)) +
  geom_col() +
  geom_text(aes(label = comma(n50)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = colores) +
  scale_y_continuous(labels = comma, limits = c(0, 28000)) +
  labs(x = NULL, y = "N50 (bp)") +
  tema_panel

panel_contigs <- ggplot(quast, aes(x = ensamblador, y = contigs_500, fill = ensamblador)) +
  geom_col() +
  geom_text(aes(label = comma(contigs_500)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = colores) +
  scale_y_continuous(labels = comma, limits = c(0, 36000)) +
  labs(x = NULL, y = "Contigs >= 500 bp") +
  tema_panel

panel_span <- ggplot(quast, aes(x = ensamblador, y = span_500 / 1e6, fill = ensamblador)) +
  geom_col() +
  geom_text(aes(label = round(span_500 / 1e6, 2)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = colores) +
  scale_y_continuous(limits = c(0, 185)) +
  labs(x = NULL, y = "Span (Mb)") +
  tema_panel

panel_mayor <- ggplot(quast, aes(x = ensamblador, y = contig_mayor, fill = ensamblador)) +
  geom_col() +
  geom_text(aes(label = comma(contig_mayor)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = colores) +
  scale_y_continuous(labels = comma, limits = c(0, 1600000)) +
  labs(x = NULL, y = "Largest contig (bp)") +
  tema_panel

figura_s0 <- gridExtra::grid.arrange(
  panel_n50, panel_contigs, panel_span, panel_mayor,
  ncol = 2
)
ggsave("figuraS0_quast.png", figura_s0, width = 10, height = 8, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Resumen en consola
# ---------------------------------------------------------------------------
cat("Figuras generadas:\n")
cat("  figura1_kraken.png\n")
cat("  figura2_ani.png\n")
cat("  figura3_busco.png\n")
cat("  figura4_anotacion.png\n")
cat("  figura5_cog.png\n")
cat("  figuraS0_quast.png\n")
cat("Todas en PNG, 300 dpi.\n")
