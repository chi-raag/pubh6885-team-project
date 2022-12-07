library(ggplot2)
library(dplyr)
library(forcats)
library(readr)
library(patchwork)
library(omicsArt)
library(scales)

results <- read.table("output/tweedieverse/significant_results.tsv",
                      sep = "\t",
                      header = T)

# COVID-19 DEGs ----------------------------------------------------------------

up_results <- results |>
  filter(value == "COVID-19") |>
  slice_max(coef, n = 10)

down_results <- results |>
  filter(value == "COVID-19") |>
  slice_min(coef, n = 10)

covid_features <- up_results %>%
  bind_rows(down_results) %>%
  pull(feature)

covid_order <- results %>%
  filter(feature %in% covid_features,
         value == "COVID-19") %>%
  mutate(feature = fct_reorder(feature, coef)) %>%
  pull(feature) %>%
  levels()

c <- results %>%
  filter(feature %in% covid_features,
         value == "COVID-19") %>%
  ggplot(aes(fct_reorder(feature, coef), coef, fill = coef)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "Gene", y = "Effect Size", title = "Covid-19") +
  theme(legend.position = "none")

lc <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "long COVID-19") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "", y = "Effect Size", title = "Long Covid-19") +
  theme(legend.position = "none")

r <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "respiratory failure") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "", y = "Effect Size", title = "Respiratory Failure") +
  theme(legend.position = "none")

c + lc + r +
  plot_annotation("Covid-19 DEGs")

ggsave(
  "output/tweedieverse/viz/c_degs.png",
  width = 3000,
  height = 1000,
  units = "px",
  dpi = "retina"
)

# Long COVID-19 DEGs ----------------------------------------------------------------

up_results <- results |>
  filter(value == "long COVID-19") |>
  slice_max(coef, n = 10)

down_results <- results |>
  filter(value == "long COVID-19") |>
  slice_min(coef, n = 10)

covid_features <- up_results %>%
  bind_rows(down_results) %>%
  pull(feature)

covid_order <- results %>%
  filter(feature %in% covid_features,
         value == "long COVID-19") %>%
  mutate(feature = fct_reorder(feature, coef)) %>%
  pull(feature) %>%
  levels()

c <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "COVID-19") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "Gene", y = "Effect Size", title = "Covid-19") +
  theme(legend.position = "none")

lc <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "long COVID-19") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "", y = "Effect Size", title = "Long Covid-19") +
  theme(legend.position = "none")

r <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "respiratory failure") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "", y = "Effect Size", title = "Respiratory Failure") +
  theme(legend.position = "none")

c + lc + r +
  plot_annotation("Long Covid-19 DEGs")

ggsave(
  "output/tweedieverse/viz/lc_degs.png",
  width = 3000,
  height = 1000,
  units = "px",
  dpi = "retina"
)

# rf DEGs ----------------------------------------------------------------

up_results <- results |>
  filter(value == "respiratory failure") |>
  slice_max(coef, n = 10)

down_results <- results |>
  filter(value == "respiratory failure") |>
  slice_min(coef, n = 10)

covid_features <- up_results %>%
  bind_rows(down_results) %>%
  pull(feature)

covid_order <- results %>%
  filter(feature %in% covid_features,
         value == "respiratory failure") %>%
  mutate(feature = fct_reorder(feature, coef)) %>%
  pull(feature) %>%
  levels()

c <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "COVID-19") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "Gene", y = "Effect Size", title = "Covid-19") +
  theme(legend.position = "none")

lc <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "long COVID-19") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "", y = "Effect Size", title = "Long Covid-19") +
  theme(legend.position = "none")

r <- results %>%
  mutate(feature = factor(feature, levels = covid_order)) %>%
  filter(feature %in% covid_features,
         value == "respiratory failure") %>%
  ggplot(aes(feature, coef, fill = coef)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = F) +
  coord_flip() +
  theme_minimal(base_size = 10)  +
  scale_fill_gradient2(low = muted("red"), high = muted("green")) +
  labs(x = "", y = "Effect Size", title = "Respiratory Failure") +
  theme(legend.position = "none")

c + lc + r +
  plot_annotation("Respiratory Failure DEGs")

ggsave(
  "output/tweedieverse/viz/rf_degs.png",
  width = 3000,
  height = 1000,
  units = "px",
  dpi = "retina"
)
