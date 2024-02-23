

# Env -------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(rstatix)
library(ggprism)
library(ggpubr)
library(broom)
library(patchwork)

# Data input ------------------------------------------------------------------------------------------------------

data <- read.delim("../data/tidy/cell_cycle_2d_48h_run_2_all_samples.txt")
metadata <- read_delim("../data/metadata/metadata.txt")

# Merge metadata
full_data <- left_join(data, metadata, by = c("Metadata_Well" = "Metadata_Well")) %>%
  select(contains("Count"), contains("Metadata"))

# Stats analysis --------------------------------------------------------------------------------------------------

group_stats <- full_data %>%
  mutate(across(contains("Count_"),
                ~ round(.x / Count_DAPI * 100, digits = 2),
                .names = "Percent_{.col}"
  ),
  Percent_Count_EdU_PH3 = round(Count_EdU_PH3/Count_PH3, digits = 2),
  Percent_Count_EdU_KI67 = round(Count_EdU_KI67/Count_KI67, digits = 2)) %>%
  group_by(Metadata_CellLine, Metadata_Dose, Metadata_Well, Metadata_Hour, Metadata_Plate) %>%
  filter(if_any(contains("Percent_Count_"), ~ near(., mean(.), tol = sd(.)))) %>%
  summarise(
    EdU = mean(Percent_Count_EdU, na.rm = TRUE),
    PH3 = mean(Percent_Count_PH3, na.rm = TRUE),
    KI67 = mean(Percent_Count_KI67, na.rm = TRUE),
    EdU_KI67 = mean(Percent_Count_EdU_KI67, na.rm = TRUE),
    EdU_PH3 = mean(Count_EdU_PH3, na.rm = TRUE)
  ) %>%
  pivot_longer(EdU:EdU_PH3, names_to = "Gene", values_to = "Count") %>%
  filter(Count < 100) %>%
  mutate(Metadata_Dose = as.factor(Metadata_Dose))

write.table(group_stats, "../data/tidy/cell_cycle_2d_48h_run_2_sumstats.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Cell cycle length -----------------------------------------------------------------------------------------------

# using the Nowakowski equation doi:10.1007/BF01190834

# Tc = length of the cell cycle (24h)
# Tm = length of the M-phase
# Ts = length of S-phase

# Growth fraction (GF) is the proportion of total population which is proliferating.
# It's calculated by determining the proportion of cells labelled at any time at or
# after Tc-Ts: GF = LI0 for t >= Tc-Ts. Estimated as the y-intercept and theoretically
# equal to Ts/Tc.

# The presence of non-proliferating cells is a mathematical dilution of the proliferating
# population. Then, the proportion of cells labelled at t = 0 is reduced from Ts/Tc by
# an amount equal to GF modifying equation 1: LI0 = GF X (Ts/Tc). And the rate of increase
# (i.e. slope) of labelled cells is described by f(t) = (GF/Tc) X t + LI0.
# in summary, the cell cycle length and s-phase length can be calculated with the following formulas:
# TC <- GF / slope

# TS <- intercept * TC / GF


cell_cycle_calculations <- group_stats %>%
  filter(Gene == "EdU", Metadata_CellLine == "H7") %>%
  group_by(Metadata_CellLine, Metadata_Dose, Metadata_Plate) %>%
  do({
    fm <- lm(Count ~ Metadata_Hour, .)
    co <- coef(fm)
    summarise(.,
              intercept = co[1],
              slope = co[2],
              r2 = summary(fm)$r.squared,
              GF = max(Count)
    ) %>%
      summarise(
        TC = GF / slope,
        TS = (intercept * TC) / GF
      )
  }) %>%
  ungroup()

# Dataviz -----------------------------------------------------------------------------------------------

ts_plot <- cell_cycle_calculations %>%
  ggplot() +
  stat_summary(
    aes(as.factor(Metadata_Dose), TS, fill = as.factor(Metadata_Dose)),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(.9)
  ) +
  geom_bar(
    aes(as.factor(Metadata_Dose), TS, fill = as.factor(Metadata_Dose)),
    stat = "summary",
    fun = "mean",
    position = position_dodge(.9)
  ) +
  geom_jitter(aes(as.factor(Metadata_Dose), TS, fill = as.factor(Metadata_Dose)),
              stat = "identity",
              position = position_dodge2(.25),
              size = .75) +
  facet_grid(. ~ Metadata_CellLine) +
  theme(text = element_text(size = 14)) +
  xlab("Dose") +
  ylab("S-phase length (h)") +
  scale_fill_discrete(name = "24(S),25-EC\n dose [µM]") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.title.x = element_blank())


tc_plot <- cell_cycle_calculations %>%
  ggplot() +
  stat_summary(
    aes(as.factor(Metadata_Dose), TC, fill = as.factor(Metadata_Dose)),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(.9)
  ) +
  geom_bar(
    aes(as.factor(Metadata_Dose), TC, fill = as.factor(Metadata_Dose)),
    stat = "summary",
    fun = "mean",
    position = position_dodge(.9)
  ) +
  geom_jitter(aes(as.factor(Metadata_Dose), TC, fill = as.factor(Metadata_Dose)),
              stat = "identity",
              position = position_dodge2(.25),
              size = .75) +
  facet_grid(. ~ Metadata_CellLine) +
  theme(text = element_text(size = 14)) +
  xlab("Dose") +
  ylab("Cell cycle length (h)") +
  scale_fill_discrete(name = "24(S),25-EC\n dose [µM]") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    axis.title.x = element_blank()
  )


ts_plot + tc_plot &
  theme(text = element_text(size = 10)) &
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(filename = "../plots/cell_cycle_2d_48h_triplicates_combined.jpeg",
       device = "jpeg",
       dpi = 600,
       units = "cm",
       width = 11,
       height = 4)
