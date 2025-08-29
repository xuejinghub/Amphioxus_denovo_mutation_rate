library(ggplot2)
library(Hmisc)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(tidyr)
library(ggtext)
library(RColorBrewer)
library(readr)


mut_data <- read.csv(
  "species_mutation_rate_removegc.tsv",
  sep = "\t",
  header = TRUE
)

classification_order <- c(
  "Arthropod",
  "Amphioxus",
  "Fish",
  "Reptile",
  "Bird",
  "Mammal"
)
# ========= A. Parental Origin Distribution =========
mut_data <- mut_data %>%
  mutate(OrganismalClassification = factor(OrganismalClassification,
    levels = classification_order
  )
  ) %>%
  arrange(OrganismalClassification, ScientificName) %>%
  mutate(ScientificName = factor(
    ScientificName,
    levels = rev(unique(ScientificName)
    )
  )
  )
plot_mutation_rate <- ggplot(mut_data,
  aes(x = MutationRateMean,
    y = ScientificName,
    color = OrganismalClassification
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = MutationRateLower95CI,
      xmax = MutationRateUpper95CI
    ),
    height = 1,
    linewidth = 0.5
  ) +
  geom_point(size = 2,
             shape = 15) +
  scale_color_manual(
    values = c(
      "Arthropod" = "#8DD3C7",
      "Amphioxus" = "#4d4398",
      "Bird" = "#fdb462",
      "Fish" = "#fb8072",
      "Mammal" = "#b3de69",
      "Reptile" = "#80b1d3"
    )
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x / 1e-8),
    limits = c(NA, 2.4e-8)
  ) +
  labs(x = "Mutation rate per site per-generation (x10<sup>-8</sup>)",
    y = "Species",
  ) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_markdown(),
        legend.position = "none",
        axis.text.y = element_text(size = 7, face = "italic"),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7))

# ========= B. Comparison of mutation rates and Ne =========
df <- read_tsv("Bergeron_ne_mu.tsv")

plot_dot <- ggplot(df,
  aes(x = Ne / 1000000,
    y = m_generation_modeled * 100000000,
    color = Taxonomic
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Amphioxus" = "#4d4398",
      "Bird" = "#fdb462",
      "Fish" = "#fb8072",
      "Mammal" = "#b3de69",
      "Reptile" = "#80b1d3"
    )
  ) +
  labs(
    x = "Effective population size (Ne)",
    y = "Mutation rate per generation"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.spacing.x = unit(0.2, "cm"),
    axis.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

plot_dot_log <- ggplot(df,
  aes(x = Ne_pi,
      y = m_generation_modeled,
      color = Taxonomic)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(
    data = data.frame(
      x = 10^seq(4.2, 6.8, length.out = 100),
      y = 0.001 / (4 * 10^seq(4.2, 6.8, length.out = 100)),
      curve = factor("0.001")
    ) %>%
      rbind(data.frame(
        x = 10^seq(4.2, 6.8, length.out = 100),
        y = 0.01 / (4 * 10^seq(4.2, 6.8, length.out = 100)),
        curve = factor(0.01)
      )) %>%
      rbind(data.frame(
        x = 10^seq(4.2, 6.8, length.out = 100),
        y = 0.04 / (4 * 10^seq(4.2, 6.8, length.out = 100)),
        curve = factor(0.04)
      )),
    aes(x = x, y = y, group = curve, color = curve),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Amphioxus" = "#4d4398",
      "Bird" = "#fdb462",
      "Fish" = "#fb8072",
      "Mammal" = "#b3de69",
      "Reptile" = "#80b1d3",
      "0.01" = "black",
      "0.001" = "black",
      "0.03" = "black",
      "0.04" = "black"
    ),
    name = "Taxonomic/Curve"
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(10^-8.8, 10^-7.8)
  ) +
  labs(
    x = expression("Effective population size (N"[e] * ")"),
    y = expression("Mutation rate per base per generation")
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  ) +
  annotation_logticks(sides = "bl",
    short = unit(0.1, "cm"),
    mid = unit(0.15, "cm"),
    long = unit(0.2, "cm"),
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# ========= Final Plot Composition =========
plot1 <- ggarrange(
  plot_dot_log,
  ncol = 1, nrow = 1,
  labels = c("b", "c"),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL)
)

fig2_out <- "fig2_mutation_rate.pdf"
ggsave(fig2_out, plot1,
  height = 9,
  width = 7,
  dpi = 300,
  units = "cm"
)
print(paste0("Plot saved as ", fig2_out))
fig2_out <- "fig2_species.pdf"
ggsave(fig2_out, plot_mutation_rate,
  height = 172.174,
  width = 104.843,
  dpi = 300,
  units = "mm"
)