# 加载必要库
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)
library(binom)
library(readr)
library(tidyr)
library(patchwork)
library(purrr)

fisher_results <- tibble(
  analysis = character(),
  feature = character(),
  odds_ratio = numeric(),
  p_value = numeric()
)
# ========= A. Parental Origin Distribution =========

parental <- "parental.tsv"
df_parental <- read.delim(parental, header = FALSE, sep = "\t")

plot_parental <- ggplot(df_parental, aes(x = V1, y = V2, fill = V1)) +
  geom_col(color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Maternal" = "#a72424", "Paternal" = "#0E4B84"),
    guide = "none"
  ) +
  geom_text(
    aes(label = V2),
    vjust = -0.3,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "parental origin",
    y = "Count",
    title = "Parental origin of DNM"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )


# ========= B. Male Bias Comparison =========

male_bias <- "male_bias.tsv"
df_male_bias <- read.delim(male_bias, header = TRUE, sep = "\t")

plot_male_bias <- ggplot(df_male_bias,
  aes(x = Species, y = a, fill = Species)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(
    aes(ymin = lower_CI, ymax = higher_CI),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  geom_hline(yintercept = 1,
             color = "#333333",
             linetype = "dashed",
             size = 0.5) +
  labs(
    title = expression(
      "The male-to-female contribution ratio (" * italic(alpha) * ")"
    ),
    x = "Genomic Feature",
    y = expression("Male bias (" * italic(alpha) * ")"),
    fill = "Parental Origin"
  ) +
  scale_fill_manual(
    values = sapply(unique(df_male_bias$Species), function(x) {
      ifelse(x == "Amphioxus", "#4d4398", "gray50")
    }),
  ) +
  scale_x_discrete(limits = c("Amphioxus", "Stickleback", "Guppy", "Snake",
                              "Chicken", "Cat",
                              "Mouse", "Human")) +
  scale_y_continuous(
    breaks = seq(0, ceiling(max(df_male_bias$a)), by = 1)
  ) +
  theme_classic() +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )
# ========= C. Mutation Fraction by Gene Feature and Parental Origin=========

gene_feature_data <- read_tsv(
  "DNM_list.lo2ref_function_cata.germline.removeGC.removeGC2.tsv"
) %>%
  mutate(feature = case_when(
    feature %in% c("downstream", "upstream") ~ "Intergenic",
    feature == "splicing" ~ "Extronic",
    feature == "CDS" ~ "Extronic",
    feature %in% c("UTR5", "UTR3") ~ "Extronic",
    TRUE ~ feature
  ))

total_parent <- gene_feature_data %>%
  group_by(parental_origin) %>%
  summarise(total = n()) %>%
  tidyr::pivot_wider(names_from = parental_origin, values_from = total)

gene_feature_proportion <- gene_feature_data %>%
  group_by(feature, parental_origin) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = parental_origin,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total = total_parent$M + total_parent$P,
    total_M = total_parent$M,
    total_P = total_parent$P,
    conf_int_M = purrr::map2(M,
      total_M,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    ),
    conf_int_P = purrr::map2(P,
      total_P,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    )
  ) %>%
  tidyr::unnest_wider(c(conf_int_M, conf_int_P),
    names_sep = "_"
  ) %>%
  rename(M_lower = conf_int_M_lower, M_upper = conf_int_M_upper,
         P_lower = conf_int_P_lower, P_upper = conf_int_P_upper)

significance_proportion <- gene_feature_proportion %>%
  group_by(feature) %>%
  summarise(
    M_count = first(M),
    P_count = first(P),
    total_M = first(total_M),
    total_P = first(total_P),
    M_upper = first(M_upper),
    P_upper = first(P_upper),
    test_result = list(
      fisher.test(matrix(c(M_count,
          total_M - M_count,
          P_count,
          total_P - P_count
        ),
        nrow = 2
      ))
    )
  ) %>%
  mutate(
    p.value = map_dbl(test_result, "p.value"),
    odds_ratio = map_dbl(test_result, ~ .x$estimate),
    sig_label = sprintf("italic(p)*'=%.2g'", p.value),
    y.position = pmax(P_upper, M_upper) + 0.05
  )

fisher_results <- significance_proportion %>%
  transmute(
    analysis = "C.Feature Mutation Type",
    feature = feature,
    odds_ratio,
    p_value = p.value
  ) %>%
  bind_rows(fisher_results)

plot_feature_proportion2 <- ggplot(
  gene_feature_proportion %>%
    tidyr::pivot_longer(c(M, P),
      names_to = "parental_origin",
      values_to = "count"
    )
) +
  geom_col(aes(x = feature,
      y = ifelse(parental_origin == "M", count / total_M, count / total_P),
      fill = parental_origin
    ),
    position = position_dodge(0.7),
    width = 0.7
  ) +
  geom_errorbar(aes(x = feature,
                    ymin = ifelse(parental_origin == "M", M_lower, P_lower),
                    ymax = ifelse(parental_origin == "M", M_upper, P_upper),
                    group = parental_origin),
    position = position_dodge(0.7), width = 0.2
  ) +
  geom_text(
    data = significance_proportion,
    aes(x = feature, y = y.position, label = sig_label),
    inherit.aes = FALSE, size = 7 / .pt,
    parse = TRUE,
  ) +
  scale_y_continuous(labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    y = "Fraction of DNMs",
    title = "DNM Proportion by Genomic Feature and Parental Origin"
  ) +
  scale_fill_manual(
    values = c("M" = "#a72424", "P" = "#0E4B84"),
    name = "Parental Origin"
  ) +
  scale_x_discrete(
    labels = function(x) {
      case_when(
        x == "cds" ~ "CDS",
        x == "intergenic" ~ "Intergenic",
        x == "utr" ~ "UTR",
        x == "intronic" ~ "Intronic",
        TRUE ~ x
      )
    }
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

# ========= D. Mutation Rate by Mutation Type and Parental Origin =========

onemer_data <- read_tsv("all.long_format.1mer_spectrum.removegc.tsv")

onemer_avg_rates <- onemer_data %>%
  group_by(Type, Parent) %>%
  summarise(mean_rate = mean(Rate, na.rm = TRUE), .groups = "drop")

onemer_avg_rates <- onemer_data %>%
  group_by(Type, Parent) %>%
  summarise(
    sum_dnm = sum(Count),
    sum_cg = sum(TotalSites),
    mean_rate = sum_dnm / sum_cg,
    .groups = "drop"
  ) %>%
  mutate(
    conf_int = purrr::map2(sum_dnm, sum_cg,
      ~binom::binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    )
  ) %>%
  tidyr::unnest_wider(conf_int)

print(onemer_avg_rates)
significance <- onemer_avg_rates %>%
  select(Type, Parent, sum_dnm, sum_cg, upper) %>%
  tidyr::pivot_wider(
    names_from = Parent,
    values_from = c(sum_dnm, sum_cg, upper)
  ) %>%
  group_by(Type) %>%
  summarise(
    p.value = chisq.test(
      matrix(c(sum_dnm_M,
               sum_cg_M - sum_dnm_M,
               sum_dnm_P,
               sum_cg_P - sum_dnm_P), nrow = 2)
    )$p.value,
    M_upper = first(upper_M),
    P_upper = first(upper_P)
  ) %>%
  mutate(
    sig_label = sprintf("italic(p)*'=%.2g'", p.value),
    max_upper = pmax(M_upper, P_upper)
  )

plot_1mer <- ggplot(onemer_avg_rates,
  aes(x = Type,
    y = mean_rate * 1000000000,
    fill = Parent
  )
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(
    aes(ymin = lower * 1000000000, ymax = upper * 1000000000),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  geom_text(
    data = significance,
    aes(x = Type, y = (max_upper * 1000000000) + 0.3, label = sig_label),
    inherit.aes = FALSE,
    size = 7 / .pt,
    parse = TRUE
  ) +
  labs(title = "Mutation Rate by Mutation Type and Parental Origin",
       x = "Genomic Feature",
       y = expression("Mutation rate (10"^-9 * ")"),
       fill = "Parental Origin") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 7)) +
  scale_fill_manual(
    values = c("M" = "#a72424", "P" = "#0E4B84"),
  ) +
  theme_classic() +
  scale_x_discrete(
    labels = function(x) {
      case_when(
        x == "A_C" ~ "A>C",
        x == "C_G" ~ "C>G",
        x == "C_T_CpG" ~ "C>T\n(CpG)",
        x == "C_T_nonCpG" ~ "C>T\n(nonCpG)",
        x == "C_A" ~ "C>A",
        x == "A_T" ~ "A>T",
        x == "A_G" ~ "A>G",
        TRUE ~ x
      )
    }
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

# ========= E. Mutation Rate by Mutation Type and Parental Origin =========

onemer_data <- read_tsv("all.long_format.1mer_spectrum.removegc.tsv")

total_parent <- onemer_data %>%
  group_by(Parent) %>%
  summarise(total = sum(Count)) %>%
  tidyr::pivot_wider(names_from = Parent, values_from = total)

onemer_avg_rates_proportion <- onemer_data %>%
  group_by(Type, Parent) %>%
  summarise(
    sum_dnm = sum(Count),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Parent,
    values_from = sum_dnm
  ) %>%
  mutate(
    total = total_parent$M + total_parent$P,
    other_M = total_parent$M,
    other_P = total_parent$P,
    conf_int_M =
      purrr::map2(
        M,
        total_parent$M,
        ~binom::binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
      ),
    conf_int_P =
      purrr::map2(P,
        total_parent$P,
        ~binom::binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    ) # nolint
  ) %>%
  tidyr::unnest_wider(c(conf_int_M, conf_int_P), names_sep = "_") %>%
  rename(M_lower = conf_int_M_lower,
         M_upper = conf_int_M_upper,
         P_lower = conf_int_P_lower,
         P_upper = conf_int_P_upper)
significance_proportion <- onemer_avg_rates_proportion %>%
  group_by(Type) %>%
  summarise(
    test_result = list(fisher.test(matrix(c(M,
        other_M - M,
        P,
        other_P - P
      ),
      nrow = 2
    )
    )),
    M_upper = first(M_upper),
    P_upper = first(P_upper)
  ) %>%
  mutate(
    p.value = map_dbl(test_result, "p.value"),
    odds_ratio = map_dbl(test_result, ~ .x$estimate),
    sig_label = sprintf("italic(p)*'=%.2g'", p.value),
    y.position = pmax(P_upper, M_upper) + 0.03
  )

fisher_results <- significance_proportion %>%
  transmute(
    analysis = "E.1mer Mutation Type",
    feature = Type,
    odds_ratio,
    p_value = p.value
  ) %>%
  bind_rows(fisher_results)

plot_1mer_proportion <- ggplot(
  onemer_avg_rates_proportion %>%
    tidyr::pivot_longer(c(M, P),
      names_to = "Parent",
      values_to = "count"
    )
) +
  geom_col(aes(x = Type,
      y = ifelse(Parent == "M",
        count / other_M, count / other_P
      ),
      fill = Parent
    ),
    position = position_dodge(width = 0.7), width = 0.7
  ) +
  geom_errorbar(
    aes(x = Type,
        ymin = ifelse(Parent == "M", M_lower, P_lower),
        ymax = ifelse(Parent == "M", M_upper, P_upper),
        group = Parent),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  geom_text(
    data = significance_proportion,
    aes(x = Type, y = y.position, label = sig_label),
    size = 7 / .pt,
    parse = TRUE
  ) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    y = "Fraction of DNMs",
    title = "DNM Proportion by Mutation Type and Parental Origin"
  ) +
  scale_fill_manual(
    values = c("M" = "#a72424", "P" = "#0E4B84")
  ) +
  scale_x_discrete(
    labels = function(x) {
      case_when(
        x == "A_C" ~ "A>C",
        x == "C_G" ~ "C>G",
        x == "C_T_CpG" ~ "C>T\n(CpG)",
        x == "C_T_nonCpG" ~ "C>T\n(nonCpG)",
        x == "C_A" ~ "C>A",
        x == "A_T" ~ "A>T",
        x == "A_G" ~ "A>G",
        TRUE ~ x
      )
    }
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

# ========= F. Shared denovo Mutation =========

ssdnm <- "DNM_list_removegc.removegc2.tsv"
df_ssdnm <- read.delim(ssdnm, header = TRUE, sep = "\t") %>%
  filter(ssDNM %in% c("shared", "unique"))

total_parent <- df_ssdnm %>%
  count(parental_origin, name = "total") %>%
  tidyr::pivot_wider(names_from = parental_origin, values_from = total)

ssdnm_total_stats <- df_ssdnm %>%
  group_by(ssDNM) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    fraction = count / sum(count)
  )

ssdnm_stats <- df_ssdnm %>%
  group_by(ssDNM, parental_origin) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = parental_origin,
                     values_from = count,
                     values_fill = 0) %>%
  mutate(
    total = total_parent$M + total_parent$P,
    other_M = total_parent$M,
    other_P = total_parent$P,
    conf_int_M = purrr::map2(M, other_M,
      ~ binom.confint(.x, .y, methods = "exact") %>%
        select(lower, upper) %>%
        rename(M_lower = lower, M_upper = upper)
    ),
    conf_int_P = purrr::map2(P, other_P,
      ~ binom.confint(.x, .y, methods = "exact") %>%
        select(lower, upper) %>%
        rename(P_lower = lower, P_upper = upper)
    )
  ) %>%
  tidyr::unnest_wider(c(conf_int_M, conf_int_P))

plot_data <- ssdnm_stats %>%
  tidyr::pivot_longer(c(M, P), names_to = "parent", values_to = "count")

plot_ssdnm <- ggplot(
  plot_data
) +
  geom_col(aes(x = ssDNM,
      y = count,
      fill = parent
    ),
    position = position_dodge(0.7),
    width = 0.7
  ) +
  geom_text(
    aes(x = ssDNM, y = count, label = count, group = parent),
    position = position_dodge(width = 0.7),
    vjust = -0.3,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(
    values = c("M" = "#a72424", "P" = "#0E4B84"),
    name = "Parental Origin"
  ) +
  labs(
    x = "ssDNM Type",
    y = "Counts",
    title = "Shared denovo mutations"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

# ========= G. Individuals shared count ========
shared_count_data <- read_tsv("shared_count.tsv")
shared_count <- ggplot(shared_count_data,
  aes(x = sharedBySamples, y = counts)
) +
  geom_col(fill = "steelblue") +
  labs(x = "Shared Count", y = "Counts") +
  geom_text(
    aes(label = counts),
    vjust = -0.3,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = seq(0, 9, 2),
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, 9)
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

# ========= H. Fraction Mutation Type and Shared Mutation ========
data <- read_tsv("DNM_list_removegc.removegc2.tsv") %>%
  mutate(
    one_mer = case_when(
      one_mer == "T_G" ~ "A_C",
      one_mer == "T_C" ~ "A_G",
      one_mer == "T_A" ~ "A_T",
      one_mer == "G_C" ~ "C_G",
      one_mer == "G_T" ~ "C_A",
      one_mer == "C_T" & CpG == "none" ~ "C_T_nonCpG",
      one_mer == "C_T" & CpG == "CpG_CtoT" ~ "C_T_CpG",
      one_mer == "G_A" & CpG == "none" ~ "C_T_nonCpG",
      one_mer == "G_A" & CpG == "CpG_CtoT" ~ "C_T_CpG",
      TRUE ~ one_mer
    )
  )
total_type <- data %>%
  group_by(ssDNM) %>%
  summarise(total = n()) %>%
  pivot_wider(names_from = ssDNM, values_from = total)

onemer_proportion <- data %>%
  group_by(one_mer, ssDNM) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = ssDNM,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total_unique = total_type$unique,
    total_shared = total_type$shared,
    conf_int_unique = purrr::map2(unique,
      total_unique,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    ),
    conf_int_shared = purrr::map2(shared,
      total_shared,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    )
  ) %>%
  tidyr::unnest_wider(conf_int_unique, names_sep = "_") %>%
  tidyr::unnest_wider(conf_int_shared, names_sep = "_") %>%
  rename(
    conf_int_unique_lower = conf_int_unique_lower,
    conf_int_unique_upper = conf_int_unique_upper,
    conf_int_shared_lower = conf_int_shared_lower,
    conf_int_shared_upper = conf_int_shared_upper
  )

significance_proportion <- onemer_proportion %>%
  group_by(one_mer) %>%
  summarise(
    test_result = list(fisher.test(matrix(c(unique,
        total_unique - unique,
        shared,
        total_shared - shared
      ),
      nrow = 2
    ))
    ),
    y.position = pmax(conf_int_unique_upper, conf_int_shared_upper) + 0.05
  ) %>%
  mutate(
    p.value = map_dbl(test_result, "p.value"),
    odds_ratio = map_dbl(test_result, ~ .x$estimate),
    sig_label = sprintf("italic('p')*'=%.2g'", p.value)
  )

fisher_results <- significance_proportion %>%
  transmute(
    analysis = "H.1mer Mutation Type",
    feature = one_mer,
    odds_ratio,
    p_value = p.value
  ) %>%
  bind_rows(fisher_results)

plot_onemer_proportion_ssDNM <- ggplot( # nolint
  onemer_proportion %>%
    pivot_longer(c(unique, shared),
      names_to = "ssDNM",
      values_to = "count"
    ),
  aes(x = one_mer,
    y = count / ifelse(ssDNM == "unique",
      total_unique,
      total_shared
    )
  )
) +
  geom_col(aes(fill = ssDNM),
           position = position_dodge(0.7), width = 0.7) +
  geom_errorbar(
    aes(
      group = ssDNM,
      ymin = ifelse(ssDNM == "unique",
        conf_int_unique_lower,
        conf_int_shared_lower
      ),
      ymax = ifelse(ssDNM == "unique",
        conf_int_unique_upper,
        conf_int_shared_upper
      )
    ),
    position = position_dodge(0.7), width = 0.2
  ) +
  geom_text(
    data = significance_proportion,
    aes(x = one_mer, y = y.position, label = sig_label),
    size = 7 / .pt, parse = TRUE
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("shared" = "steelblue", "unique" = "grey80")
  ) +
  labs(y = "Fraction of DNMs",
    title = "DNM Fraction by Chromosomes and unique"
  ) +
  scale_x_discrete(
    labels = function(x) {
      case_when(
        x == "A_C" ~ "A>C",
        x == "C_G" ~ "C>G",
        x == "C_T_CpG" ~ "C>T\n(CpG)",
        x == "C_T_nonCpG" ~ "C>T\n(nonCpG)",
        x == "C_A" ~ "C>A",
        x == "A_T" ~ "A>T",
        x == "A_G" ~ "A>G",
        TRUE ~ x 
      )
    }
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.text = element_text(size = 7)
  )

# ========= Final Plot Composition =========
plot1 <- ggarrange(
  ggarrange(
    plot_parental, plot_male_bias, plot_feature_proportion2,
    ncol = 3, nrow = 1,
    labels = c("A", "B", "C"),
    font.label = list(size = 10, color = "black", face = "bold", family = NULL),
    widths = c(1.2, 2, 2)
  ),
  ggarrange(
    plot_1mer, plot_1mer_proportion,
    ncol = 2, nrow = 1,
    labels = c("D", "E"),
    font.label = list(size = 10, color = "black", face = "bold", family = NULL),
    widths = c(1, 1)
  ),
  ggarrange(
    plot_ssdnm, shared_count, plot_onemer_proportion_ssDNM,
    ncol = 3, nrow = 1,
    labels = c("F", "G", "H"),
    font.label = list(size = 10, color = "black", face = "bold", family = NULL),
    widths = c(1, 1, 2)
  ),
  nrow = 3,
  heights = c(1, 1, 1)
)

# ========= Final Plot output =========

output_pdf <- "fig3.pdf"
ggsave(output_pdf, plot1, dpi = 300, width = 18, height = 15, unit = "cm")
cat("Plot saved to:", output_pdf, "\n")

write_tsv(fisher_results, "fisher_test_results.tsv")
print("Fisher's exact test results saved to fisher_test_results.tsv")
