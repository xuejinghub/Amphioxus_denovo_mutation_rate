library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)
library(binom)
library(readr)
library(tidyr)
library(patchwork)

# ======= 1. 1mer spectrum ===========

onemer <- read.csv("onemer_statsdata.tsv", header = TRUE, sep = "\t") # nolint: line_length_linter.

onemer_percent <- onemer %>%
  group_by(Species) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup() %>%
  mutate(Species = factor(Species,
    levels = c("human", "mouse", "stickleback", "B. floridae"),
    labels = c(
      expression(italic("H. sapiens") * " (human)"),
      expression(italic("M. musculus") * " (mouse)"),
      expression(italic("G. aculeatus") * " (stickleback)"),
      expression(italic("B. floridae") * " (Amphioxus)")
    ),
  ),
  ordered = TRUE)
plot_1mer <- ggplot(onemer_percent,
                    aes(x = Type, y = Percent, fill = Species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  labs(x = "Mutation Type", y = "Percentage (%)", fill = "Species") +
  scale_x_discrete(
    labels = function(x) {
      case_when(
        x == "A_C" ~ "A>C",
        x == "C_G" ~ "C>G",
        x == "C_T_CpG" ~ "C>T (CpG)",
        x == "C_T_nonCpG" ~ "C>T (nonCpG)",
        x == "C_A" ~ "C>A",
        x == "A_T" ~ "A>T",
        x == "A_G" ~ "A>G",
        TRUE ~ x
      )
    }
  ) +
  theme_classic() +
  facet_grid(Species ~ .,
             scales = "free_x",
             labeller = label_parsed) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ========= Final Plot output =========
output_file2 <- "figS2_evolution.pdf"
ggsave(output_file2, plot_1mer, dpi = 300, width = 9, height = 15, unit = "cm")
cat("Plot saved to:", output_file2, "\n")
