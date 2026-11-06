library(tidyverse)
library(scales)

final_schools_data <- read_csv("data/processed/schools_data.csv", show_col_types = FALSE)

stats_summary <- final_schools_data %>%
  summarise(
    Mean   = mean(MMR_PR, na.rm = TRUE),
    Median = median(MMR_PR, na.rm = TRUE),
    SD     = sd(MMR_PR, na.rm = TRUE),
    Min    = min(MMR_PR, na.rm = TRUE),
    Q1     = quantile(MMR_PR, 0.25, na.rm = TRUE),
    Q3     = quantile(MMR_PR, 0.75, na.rm = TRUE),
    Max    = max(MMR_PR, na.rm = TRUE)
  )

print(stats_summary)

total_schools  <- nrow(final_schools_data)
below_95_count <- sum(final_schools_data$MMR_PR < 0.95, na.rm = TRUE)
pct_below_95   <- (below_95_count / total_schools) * 100

cat(sprintf("\nTotal Schools: %d\n", total_schools))
cat(sprintf("Schools below 95%% threshold: %d (%.2f%%)\n\n", below_95_count, pct_below_95))

p_hist <- ggplot(final_schools_data, aes(x = MMR_PR)) +
  geom_histogram(binwidth = 0.01, fill = "#5B84B1", color = "white") +
  geom_vline(xintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.95, y = Inf, label = "95% Herd Immunity",
           vjust = 2, hjust = 1.1, color = "#B22222", fontface = "bold", size = 4) +
  labs(
    title = "School-Level Coverage Distribution",
    subtitle = "Distribution of MMR coverage showing coverage gaps",
    x = "Observed Vaccinated Fraction (MMR_PR)",
    y = "Number of Schools"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

print(p_hist)

p_box <- ggplot(final_schools_data, aes(y = MMR_PR, x = "")) +
  geom_boxplot(fill = "#5B84B1", alpha = 0.7, outlier.color = "red", outlier.shape = 16) +
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 1.4, y = 0.955, label = "95% Threshold",
           color = "#B22222", fontface = "bold", size = 3.8) +
  labs(
    title = "Overall School Coverage Distribution",
    x = "",
    y = "Observed Vaccinated Fraction (MMR_PR)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.ticks.x = element_blank()
  )

print(p_box)

lowest_10_schools <- final_schools_data %>%
  arrange(MMR_PR) %>%
  head(10)
p_bar_lowest <- ggplot(lowest_10_schools,
                       aes(x = reorder(school_name, MMR_PR), y = MMR_PR)) +
  geom_col(fill = "#D95F02", alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f%%", MMR_PR * 100)),
            hjust = -0.15, size = 3.5, fontface = "bold") +
  coord_flip(ylim = c(0, 1.05)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Top 10 Lowest Coverage Schools",
    subtitle = "Schools with lowest vaccination rates",
    x = "School Name",
    y = "MMR Coverage Rate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_blank()
  )

print(p_bar_lowest)

p_map <- ggplot(final_schools_data, aes(x = school_long, y = school_lat, color = MMR_PR)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_viridis_c(option = "magma", labels = scales::percent_format(accuracy = 1)) +
  geom_point(data = filter(final_schools_data, MMR_PR < 0.95),
             color = "red", shape = 1, stroke = 1.2, size = 4) +
  labs(
    title = "Geographic Distribution of School Vaccination Rates",
    subtitle = "Red outlined circles highlight schools below the 95% herd immunity threshold",
    x = "Longitude",
    y = "Latitude",
    color = "MMR Rate"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

print(p_map)

final_schools_data <- final_schools_data %>%
  mutate(
    risk_category = case_when(
      MMR_PR >= 0.95 ~ ">=95%",
      MMR_PR >= 0.90 ~ "90–94.9%",
      MMR_PR >= 0.80 ~ "80–89.9%",
      TRUE           ~ "<80%"
    ),
    risk_category = factor(
      risk_category,
      levels = c(">=95%", "90–94.9%", "80–89.9%", "<80%")
    )
  )

p_risk <- ggplot(final_schools_data,
                 aes(x = risk_category, fill = risk_category)) +
  geom_bar(alpha = 0.85) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      ">=95%"    = "#2B83BA",
      "90–94.9%" = "#FDAE61",
      "80–89.9%" = "#D7191C",
      "<80%"     = "#8B0000"
    )
  ) +
  labs(
    title = "Breakdown of Schools Across MMR Vaccination Thresholds",
    x = "Vaccination Coverage Level",
    y = "Number of Schools"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

print(p_risk)

coords <- final_schools_data %>% select(school_long, school_lat)
dist_matrix <- as.matrix(dist(coords))
diag(dist_matrix) <- NA
final_schools_data$dist_to_nearest_school <- apply(dist_matrix, 1, min, na.rm = TRUE)

final_schools_data <- final_schools_data %>%
  mutate(
    school_category = case_when(
      grepl("Elementary|Elem|Primary", school_name, ignore.case = TRUE) ~ "Elementary/Primary",
      grepl("High|Secondary|Prep", school_name, ignore.case = TRUE)      ~ "High School",
      grepl("Charter|Academy", school_name, ignore.case = TRUE)          ~ "Charter/Academy",
      grepl("Christian|Catholic|Saint|St.|Lutheran", school_name, ignore.case = TRUE) ~ "Religious Affiliated",
      TRUE ~ "Other / Unspecified"
    )
  )

p_text_type <- ggplot(final_schools_data, aes(x = school_category, y = MMR_PR, fill = school_category)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "MMR Coverage by School Category",
    subtitle = "Categories extracted from school names",
    x = "School Category",
    y = "MMR Coverage Rate"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), axis.text.x = element_text(angle = 15, hjust = 1))

print(p_text_type)

med_lat  <- median(final_schools_data$school_lat, na.rm = TRUE)
med_long <- median(final_schools_data$school_long, na.rm = TRUE)

final_schools_data <- final_schools_data %>%
  mutate(
    region_quadrant = case_when(
      school_lat >= med_lat  & school_long >= med_long ~ "North-East",
      school_lat >= med_lat  & school_long <  med_long ~ "North-West",
      school_lat <  med_lat  & school_long >= med_long ~ "South-East",
      school_lat <  med_lat  & school_long <  med_long ~ "South-West"
    )
  )

p_quadrant <- ggplot(final_schools_data, aes(x = region_quadrant, y = MMR_PR, fill = region_quadrant)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Regional Coverage Comparison by Quadrant",
    x = "Geographic Quadrant",
    y = "MMR Coverage Rate"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

print(p_quadrant)

ggsave("figures/exploration/coverage_histogram.png", plot = p_hist, width = 8, height = 5, dpi = 300)
ggsave("figures/exploration/coverage_boxplot.png", plot = p_box, width = 5, height = 6, dpi = 300)
ggsave("figures/exploration/lowest_10_schools.png", plot = p_bar_lowest, width = 9, height = 5, dpi = 300)
ggsave("figures/exploration/geographic_coverage_map.png", plot = p_map, width = 8, height = 6, dpi = 300)
ggsave("figures/exploration/risk_stratification.png", plot = p_risk, width = 8, height = 5, dpi = 300)
ggsave("figures/exploration/coverage_by_category.png", plot = p_text_type, width = 8, height = 5, dpi = 300)
ggsave("figures/exploration/coverage_by_quadrant.png", plot = p_quadrant, width = 8, height = 5, dpi = 300)
