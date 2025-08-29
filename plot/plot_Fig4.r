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


# ========= 1. Parental Origin Distribution =========

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
    title = "Parental origin of\npzDNM"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

# ========= 4. Mutation Fraction by Gene Feature and Parental Origin ========

gene_feature_data <- read_tsv("postzygoticDNM_removeSimilarnonbubblehighDepth.info.result.removegc.mapped.removesimilar.tsv") %>% # nolint: line_length_linter.
  mutate(feature = case_when(
    feature %in% c("downstream", "upstream") ~ "Intergenic",
    feature == "intergenic" ~ "Intergenic",
    feature == "splicing" ~ "Exonic",
    feature == "CDS" ~ "Exonic",
    feature == "intronic" ~ "Intronic",
    feature %in% c("UTR5", "UTR3") ~ "Exonic",
    TRUE ~ feature
  ))

# 计算全局亲本比例
total_parent <- gene_feature_data %>%
  group_by(parental) %>%
  summarise(total = n()) %>%
  tidyr::pivot_wider(names_from = parental, values_from = total)

gene_feature_proportion <- gene_feature_data %>%
  group_by(feature, parental) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = parental,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total = total_parent$M + total_parent$P,
    total_M = total_parent$M,
    total_P = total_parent$P,
    conf_int_M = purrr::map2(
      M,
      total_M,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    ),
    conf_int_P = purrr::map2(
      P,
      total_P,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    )
  ) %>%
  tidyr::unnest_wider(
    c(conf_int_M, conf_int_P),
    names_sep = "_"
  ) %>%  # 添加names_sep参数
  rename(M_lower = conf_int_M_lower, M_upper = conf_int_M_upper,
         P_lower = conf_int_P_lower, P_upper = conf_int_P_upper)

# 显著性检验（卡方检验）
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
      fisher.test(
        matrix(c(M_count, total_M - M_count, P_count, total_P - P_count),
          nrow = 2
        )
      )
    )
  ) %>%
  mutate(
    p.value = map_dbl(test_result, "p.value"),
    odds_ratio = map_dbl(test_result, ~ .x$estimate),
    sig_label = sprintf("italic('p')*'=%.2g'", p.value),
    group1 = "M",
    group2 = "P",
    feature = as.character(feature),
    y.position = pmax(P_upper, M_upper) + 0.02  # 使用显式获取的列
  )

fisher_results <- significance_proportion %>%
  transmute(
    analysis = "B.feature_proportion_sex",
    feature = feature,
    odds_ratio,
    p_value = p.value
  ) %>%
  bind_rows(fisher_results)

plot_feature_proportion_sex <- ggplot(
  gene_feature_proportion %>%
    tidyr::pivot_longer(c(M, P),
                        names_to = "parental",
                        values_to = "count")
) +
  geom_col(aes(x = feature,
               y = ifelse(parental == "M",
                          count / total_M, count / total_P),
               fill = parental),
           position = position_dodge(0.7),
           width = 0.7) +
  geom_errorbar(
    aes(x = feature,
        ymin = ifelse(parental == "M", M_lower, P_lower),
        ymax = ifelse(parental == "M", M_upper, P_upper),
        group = parental),
    position = position_dodge(0.7), width = 0.2
  ) +
  geom_text(
    data = significance_proportion,
    aes(x = feature, y = y.position, label = sig_label),
    inherit.aes = FALSE, size = 7 / .pt,
    parse = TRUE
  ) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    y = "Fraction of PZMs",
    title = "DNM Fraction by Genomic Feature and Parental Origin"
  ) +
  scale_fill_manual(
    values = c("M" = "#a72424", "P" = "#0E4B84"),
    name = "Parental Origin"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

# ========= 6. Mutation Fraction by 1mer and Parental Origin ========

gene_one_mer_data <- read_tsv("postzygoticDNM_removeSimilarnonbubblehighDepth.info.result.removegc.mapped.removesimilar.tsv") %>% # nolint: line_length_linter.
  mutate(
    one_mer = case_when(
      mutation_type == "T_G" ~ "A_C",
      mutation_type == "T_C" ~ "A_G",
      mutation_type == "T_A" ~ "A_T",
      mutation_type == "G_C" ~ "C_G",
      mutation_type == "G_T" ~ "C_A",
      mutation_type == "G_A" & CpGs == "none" ~ "C_T_nonCpG",
      mutation_type == "G_A" & CpGs == "CpG_CtoT" ~ "C_T_CpG",
      mutation_type == "C_T" & CpGs == "none" ~ "C_T_nonCpG",
      mutation_type == "C_T" & CpGs == "CpG_CtoT" ~ "C_T_CpG",
      TRUE ~ mutation_type
    )
  )

total_parent <- gene_one_mer_data %>%
  group_by(parental) %>%
  summarise(total = n()) %>%
  tidyr::pivot_wider(names_from = parental, values_from = total)

gene_one_mer_proportion <- gene_one_mer_data %>%
  group_by(one_mer, parental) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = parental,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total = total_parent$M + total_parent$P,
    total_M = total_parent$M,
    total_P = total_parent$P,
    conf_int_M = purrr::map2(
      M,
      total_M,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    ),
    conf_int_P = purrr::map2(
      P,
      total_P,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    )
  ) %>%
  tidyr::unnest_wider(
    c(conf_int_M, conf_int_P),
    names_sep = "_"
  ) %>%
  rename(M_lower = conf_int_M_lower, M_upper = conf_int_M_upper,
         P_lower = conf_int_P_lower, P_upper = conf_int_P_upper)

significance_proportion <- gene_one_mer_proportion %>%
  group_by(one_mer) %>%
  summarise(
    M_count = first(M),
    P_count = first(P),
    total_M = first(total_M),
    total_P = first(total_P),
    M_upper = first(M_upper),
    P_upper = first(P_upper),
    test_result = list(
      fisher.test(
        matrix(c(M_count, total_M - M_count, P_count, total_P - P_count),
          nrow = 2
        )
      )
    )
  ) %>%
  mutate(
    p.value = map_dbl(test_result, "p.value"),
    odds_ratio = map_dbl(test_result, ~ .x$estimate),
    sig_label = sprintf("italic('p')*'=%.2g'", p.value),
    group1 = "M",
    group2 = "P",
    one_mer = as.character(one_mer),
    y.position = pmax(P_upper, M_upper) + 0.02
  )

fisher_results <- significance_proportion %>%
  transmute(
    analysis = "C.one_mer_proportion_sex",
    feature = one_mer,
    odds_ratio,
    p_value = p.value
  ) %>%
  bind_rows(fisher_results)

plot_one_mer_proportion_sex <- ggplot(
  gene_one_mer_proportion %>%
    tidyr::pivot_longer(c(M, P),
                        names_to = "parental",
                        values_to = "count")
) +
  geom_col(aes(x = one_mer,
               y = ifelse(parental == "M",
                          count / total_M, count / total_P),
               fill = parental),
           position = position_dodge(0.7),
           width = 0.7) +
  geom_errorbar(
    aes(x = one_mer,
        ymin = ifelse(parental == "M", M_lower, P_lower),
        ymax = ifelse(parental == "M", M_upper, P_upper),
        group = parental),
    position = position_dodge(0.7), width = 0.2
  ) +
  geom_text(
    data = significance_proportion,
    aes(x = one_mer, y = y.position, label = sig_label),
    inherit.aes = FALSE, size = 7 / .pt,
    parse = TRUE
  ) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    y = "Fraction of PZMs",
    title = "DNM Fraction by 1mer and Parental Origin"
  ) +
  scale_fill_manual(
    values = c("M" = "#a72424", "P" = "#0E4B84"),
    name = "Parental Origin"
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
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )
# ========= 5. Mutation Fraction by Gene Feature and mutation cata. ========

gene_feature_data <- read_tsv("/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/haplotype_ref/individual/DNM_list.lo2ref_function_cata.removegc.mapped.addsimilar.removeheteroalleles.tsv") %>% # nolint: line_length_linter.
  mutate(feature = case_when(
    feature %in% c("downstream", "upstream") ~ "Intergenic",
    feature == "intergenic" ~ "Intergenic",
    feature == "splicing" ~ "Exonic",
    feature == "CDS" ~ "Exonic",
    feature == "intronic" ~ "Intronic",
    feature %in% c("UTR5", "UTR3", "exonic") ~ "Exonic",
    TRUE ~ feature
  ))

total_cata <- gene_feature_data %>%
  group_by(cata) %>%
  summarise(total = n()) %>%
  tidyr::pivot_wider(names_from = cata, values_from = total)

gene_feature_proportion <- gene_feature_data %>%
  group_by(feature, cata) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = cata,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total = total_cata$germline + total_cata$somatic,
    total_germline = total_cata$germline,
    total_somatic = total_cata$somatic,
    conf_int_germline = purrr::map2(
      germline,
      total_germline,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    ),
    conf_int_somatic = purrr::map2(
      somatic,
      total_somatic,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    )
  ) %>%
  tidyr::unnest_wider(
    c(conf_int_germline, conf_int_somatic),
    names_sep = "_"
  ) %>%
  rename(germline_lower = conf_int_germline_lower,
    germline_upper = conf_int_germline_upper,
    somatic_lower = conf_int_somatic_lower,
    somatic_upper = conf_int_somatic_upper
  )

significance_proportion <- gene_feature_proportion %>%
  group_by(feature) %>%
  summarise(
    germline_count = first(germline),
    somatic_count = first(somatic),
    total_germline = first(total_germline),
    total_somatic = first(total_somatic),
    germline_upper = first(germline_upper),
    somatic_upper = first(somatic_upper),
    test_result = list(
      fisher.test(
        matrix(c(germline_count,
            somatic_count,
            total_germline - germline_count,
            total_somatic - somatic_count
          ),
          nrow = 2
        )
      )
    )
  ) %>%
  mutate(
    p.value = map_dbl(test_result, "p.value"),
    odds_ratio = map_dbl(test_result, ~ .x$estimate),
    sig_label = sprintf("italic('p')*'=%.2g'", p.value),
    feature = as.character(feature),
    y.position = pmax(somatic_upper, germline_upper) + 0.02  # 使用显式获取的列
  )


fisher_results <- significance_proportion %>%
  transmute(
    analysis = "D.feature_proportion_cata",
    feature = as.character(feature),
    odds_ratio,
    p_value = p.value
  ) %>%
  bind_rows(fisher_results)

plot_feature_proportion_cata <- ggplot(
  gene_feature_proportion %>%
    tidyr::pivot_longer(c(germline, somatic),
                        names_to = "cata",
                        values_to = "count")
) +
  geom_col(aes(x = feature,
               y = ifelse(cata == "germline",
                          count / total_germline, count / total_somatic),
               fill = cata),
           position = position_dodge(0.7),
           width = 0.7) +
  geom_errorbar(
    aes(x = feature,
        ymin = ifelse(cata == "germline", germline_lower, somatic_lower),
        ymax = ifelse(cata == "germline", germline_upper, somatic_upper),
        group = cata),
    position = position_dodge(0.7), width = 0.2
  ) +
  geom_text(
    data = significance_proportion,
    aes(x = feature, y = y.position, label = sig_label),
    inherit.aes = FALSE, size = 7 / .pt,
    parse = TRUE
  ) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    y = "Fraction of mutations",
    title = "DNM Fraction by Genomic Feature and Mutation Catagory"
  ) +
  scale_fill_manual(
    values = c("germline" = "steelblue", "somatic" = "grey80"),
    name = "Parental Origin"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.text = element_text(size = 7)
  )
# ========= 7. Mutation Fraction by one_mer and mutation cata. ========

gene_feature_data <- read_tsv("DNM_list.lo2ref_function_cata.removegc.mapped.addsimilar.tsv") %>% # nolint: line_length_linter.
  mutate(
    one_mer = case_when(
      one_mer == "T_G" ~ "A_C",
      one_mer == "T_C" ~ "A_G",
      one_mer == "T_A" ~ "A_T",
      one_mer == "G_C" ~ "C_G",
      one_mer == "G_T" ~ "C_A",
      one_mer == "G_A" & CpGs == "none" ~ "C_T_nonCpG",
      one_mer == "G_A" & CpGs == "CpG_CtoT" ~ "C_T_CpG",
      one_mer == "C_T" & CpGs == "none" ~ "C_T_nonCpG",
      one_mer == "C_T" & CpGs == "CpG_CtoT" ~ "C_T_CpG",
      TRUE ~ one_mer
    )
  )

total_cata <- gene_feature_data %>%
  group_by(cata) %>%
  summarise(total = n()) %>%
  tidyr::pivot_wider(names_from = cata, values_from = total)

gene_feature_proportion <- gene_feature_data %>%
  group_by(one_mer, cata) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = cata,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total = total_cata$germline + total_cata$somatic,
    total_germline = total_cata$germline,
    total_somatic = total_cata$somatic,
    conf_int_germline = purrr::map2(
      germline,
      total_germline,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    ),
    conf_int_somatic = purrr::map2(
      somatic,
      total_somatic,
      ~binom.confint(.x, .y, methods = "exact")[c("lower", "upper")]
    )
  ) %>%
  tidyr::unnest_wider(
    c(conf_int_germline, conf_int_somatic),
    names_sep = "_"
  ) %>%
  rename(germline_lower = conf_int_germline_lower,
    germline_upper = conf_int_germline_upper,
    somatic_lower = conf_int_somatic_lower,
    somatic_upper = conf_int_somatic_upper
  )

significance_proportion <- gene_feature_proportion %>%
  group_by(one_mer) %>%
  summarise(
    germline_count = first(germline),
    somatic_count = first(somatic),
    total_germline = first(total_germline),
    total_somatic = first(total_somatic),
    germline_upper = first(germline_upper),
    somatic_upper = first(somatic_upper),
    test_result = list(
      fisher.test(
        matrix(c(germline_count,
            total_germline - germline_count,
            somatic_count,
            total_somatic - somatic_count
          ),
          nrow = 2
        )
      )
    )
  ) %>%
  mutate(
    p.value = map_dbl(test_result, "p.value"),
    odds_ratio = map_dbl(test_result, ~ .x$estimate),
    sig_label = sprintf("italic('p')*'=%.2g'", p.value),
    one_mer = as.character(one_mer),
    y.position = pmax(somatic_upper, germline_upper) + 0.02
  )

fisher_results <- significance_proportion %>%
  transmute(
    analysis = "E.one_mer_proportion_cata",
    feature = one_mer,
    odds_ratio,
    p_value = p.value
  ) %>%
  bind_rows(fisher_results)

plot_one_mer_cata <- ggplot(
  gene_feature_proportion %>%
    tidyr::pivot_longer(c(germline, somatic),
                        names_to = "cata",
                        values_to = "count")
) +
  geom_col(aes(x = one_mer,
               y = ifelse(cata == "germline",
                          count / total_germline, count / total_somatic),
               fill = cata),
           position = position_dodge(0.7),
           width = 0.7) +
  geom_errorbar(
    aes(x = one_mer,
        ymin = ifelse(cata == "germline", germline_lower, somatic_lower),
        ymax = ifelse(cata == "germline", germline_upper, somatic_upper),
        group = cata),
    position = position_dodge(0.7), width = 0.2
  ) +
  geom_text(
    data = significance_proportion,
    aes(x = one_mer, y = y.position, label = sig_label),
    inherit.aes = FALSE, size = 7 / .pt,
    parse = TRUE
  ) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    y = "Fraction of mutations",
    title = "DNM Fraction by 1mer and Mutation Catagory"
  ) +
  scale_fill_manual(
    values = c("germline" = "steelblue", "somatic" = "grey80"),
    name = "Parental Origin"
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
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.text = element_text(size = 7)
  )

# ========= Final Plot Composition =========
plot1 <- ggarrange(
  ggarrange(
    plot_parental, plot_feature_proportion_sex, plot_one_mer_proportion_sex,
    ncol = 3, nrow = 1,
    labels = c("A", "B", "C"),
    widths = c(1, 1.5, 2.5),
    font.label = list(size = 10, color = "black", face = "bold", family = NULL)
  ),
  ggarrange(
    plot_feature_proportion_cata, plot_one_mer_cata,
    ncol = 2, nrow = 1,
    labels = c("D", "E"),
    font.label = list(size = 10, color = "black", face = "bold", family = NULL),
    widths = c(1, 2)
  ),
  nrow = 2,
  heights = c(1, 1)
)

# ========= Final Plot output =========
output_pdf <- "fig4_pzDNM.pdf"
ggsave(output_pdf, plot1, dpi = 300, width = 18, height = 10, unit = "cm")
cat("Plot saved to:", output_pdf, "\n")
write_tsv(fisher_results, "fisher_test_results.tsv")
print("Fisher's exact test results saved to fisher_test_results.tsv")
