library(tidyverse)  # data wrangling and plots
library(scales)  # axis labels

final_schools_data <- read_csv("data/processed/schools_data.csv", show_col_types = FALSE)  # unique schools from both ends of each edge

stats_summary <- final_schools_data %>%  # mean, median, and spread of MMR coverage
  summarise(
    Mean   = mean(MMR_PR, na.rm = TRUE),
    Median = median(MMR_PR, na.rm = TRUE),
    SD     = sd(MMR_PR, na.rm = TRUE),
    Min    = min(MMR_PR, na.rm = TRUE),
    Q1     = quantile(MMR_PR, 0.25, na.rm = TRUE),
    Q3     = quantile(MMR_PR, 0.75, na.rm = TRUE),
    Max    = max(MMR_PR, na.rm = TRUE)
  )

print(stats_summary)  # show the result

total_schools  <- nrow(final_schools_data)
below_95_count <- sum(final_schools_data$MMR_PR < 0.95, na.rm = TRUE)  # schools under 95% MMR coverage
pct_below_95   <- (below_95_count / total_schools) * 100  # percent of schools under 95% coverage

cat(sprintf("\nTotal Schools: %d\n", total_schools))  # print a check
cat(sprintf("Schools below 95%% threshold: %d (%.2f%%)\n\n", below_95_count, pct_below_95))  # print a check

p_hist <- ggplot(final_schools_data, aes(x = MMR_PR)) +  # start the plot
  geom_histogram(binwidth = 0.01, fill = "#5B84B1", color = "white") +  # distribution of coverage
  geom_vline(xintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +  # mark the threshold
  annotate("text", x = 0.95, y = Inf, label = "95% Herd Immunity",  # label the reference line
           vjust = 2, hjust = 1.1, color = "#B22222", fontface = "bold", size = 4) +
  labs(  # titles and axis labels
    title = "School-Level Coverage Distribution",
    subtitle = "Distribution of MMR coverage showing coverage gaps",
    x = "Observed Vaccinated Fraction (MMR_PR)",
    y = "Number of Schools"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

print(p_hist)  # show the result

p_box <- ggplot(final_schools_data, aes(y = MMR_PR, x = "")) +  # start the plot
  geom_boxplot(fill = "#5B84B1", alpha = 0.7, outlier.color = "red", outlier.shape = 16) +  # spread of coverage
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +  # reference line
  annotate("text", x = 1.4, y = 0.955, label = "95% Threshold",  # label the reference line
           color = "#B22222", fontface = "bold", size = 3.8) +
  labs(  # titles and axis labels
    title = "Overall School Coverage Distribution",
    x = "",
    y = "Observed Vaccinated Fraction (MMR_PR)"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.ticks.x = element_blank()
  )

print(p_box)  # show the result

lowest_10_schools <- final_schools_data %>%
  arrange(MMR_PR) %>%
  head(10)
p_bar_lowest <- ggplot(lowest_10_schools,  # start the plot
                       aes(x = reorder(school_name, MMR_PR), y = MMR_PR)) +
  geom_col(fill = "#D95F02", alpha = 0.85, width = 0.7) +  # bar lengths show the values
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +  # reference line
  geom_text(aes(label = sprintf("%.1f%%", MMR_PR * 100)),
            hjust = -0.15, size = 3.5, fontface = "bold") +
  coord_flip(ylim = c(0, 1.05)) +  # horizontal bars
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(  # titles and axis labels
    title = "Top 10 Lowest Coverage Schools",
    subtitle = "Schools with lowest vaccination rates",
    x = "School Name",
    y = "MMR Coverage Rate"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_blank()
  )

print(p_bar_lowest)  # show the result

p_map <- ggplot(final_schools_data, aes(x = school_long, y = school_lat, color = MMR_PR)) +  # start the plot
  geom_point(alpha = 0.8, size = 3) +  # one mark per school
  scale_color_viridis_c(option = "magma", labels = scales::percent_format(accuracy = 1)) +
  geom_point(data = filter(final_schools_data, MMR_PR < 0.95),  # one mark per school
             color = "red", shape = 1, stroke = 1.2, size = 4) +
  labs(  # titles and axis labels
    title = "Geographic Distribution of School Vaccination Rates",
    subtitle = "Red outlined circles highlight schools below the 95% herd immunity threshold",
    x = "Longitude",
    y = "Latitude",
    color = "MMR Rate"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(plot.title = element_text(face = "bold", size = 14))

print(p_map)  # show the result

final_schools_data <- final_schools_data %>%  # unique schools from both ends of each edge
  mutate(
    risk_category = case_when(
      MMR_PR >= 0.80 ~ "High Risk (80-89.9%)",
      MMR_PR >= 0.90 ~ "Moderate Risk (90-94.9%)",
      MMR_PR >= 0.95 ~ "Optimal (>=95%)",
      TRUE           ~ "Critical Risk (<80%)"
    ),
    risk_category = factor(
      risk_category,
      levels = c("Optimal (>=95%)", "Moderate Risk (90-94.9%)", "High Risk (80-89.9%)", "Critical Risk (<80%)")
    )
  )

p_risk <- ggplot(final_schools_data, aes(x = risk_category, fill = risk_category)) +  # start the plot
  geom_bar(alpha = 0.85) +  # bar lengths show the values
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, fontface = "bold") +
  scale_fill_manual(values = c("Optimal (>=95%)" = "#2B83BA",  # set the legend colours or line types
                               "Moderate Risk (90-94.9%)" = "#FDAE61",
                               "High Risk (80-89.9%)" = "#D7191C",
                               "Critical Risk (<80%)" = "#8B0000")) +
  labs(  # titles and axis labels
    title = "School Risk Stratification Profiles",
    x = "Risk Category",
    y = "Number of Schools"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "none")

print(p_risk)  # show the result

coords <- final_schools_data %>% select(school_long, school_lat)
dist_matrix <- as.matrix(dist(coords))
diag(dist_matrix) <- NA
final_schools_data$dist_to_nearest_school <- apply(dist_matrix, 1, min, na.rm = TRUE)

final_schools_data <- final_schools_data %>%  # unique schools from both ends of each edge
  mutate(
    school_category = case_when(
      grepl("Elementary|Elem|Primary", school_name, ignore.case = TRUE) ~ "Elementary/Primary",
      grepl("High|Secondary|Prep", school_name, ignore.case = TRUE)      ~ "High School",
      grepl("Charter|Academy", school_name, ignore.case = TRUE)          ~ "Charter/Academy",
      grepl("Christian|Catholic|Saint|St.|Lutheran", school_name, ignore.case = TRUE) ~ "Religious Affiliated",
      TRUE ~ "Other / Unspecified"
    )
  )

p_text_type <- ggplot(final_schools_data, aes(x = school_category, y = MMR_PR, fill = school_category)) +  # start the plot
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +  # spread of coverage
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +  # reference line
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(  # titles and axis labels
    title = "MMR Coverage by School Category",
    subtitle = "Categories extracted from school names",
    x = "School Category",
    y = "MMR Coverage Rate"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(plot.title = element_text(face = "bold", size = 14), axis.text.x = element_text(angle = 15, hjust = 1))

print(p_text_type)  # show the result

med_lat  <- median(final_schools_data$school_lat, na.rm = TRUE)
med_long <- median(final_schools_data$school_long, na.rm = TRUE)

final_schools_data <- final_schools_data %>%  # unique schools from both ends of each edge
  mutate(
    region_quadrant = case_when(
      school_lat >= med_lat  & school_long >= med_long ~ "North-East",
      school_lat >= med_lat  & school_long <  med_long ~ "North-West",
      school_lat <  med_lat  & school_long >= med_long ~ "South-East",
      school_lat <  med_lat  & school_long <  med_long ~ "South-West"
    )
  )

p_quadrant <- ggplot(final_schools_data, aes(x = region_quadrant, y = MMR_PR, fill = region_quadrant)) +  # start the plot
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +  # spread of coverage
  scale_fill_brewer(palette = "Set3") +
  geom_hline(yintercept = 0.95, color = "#B22222", linetype = "dashed", linewidth = 1) +  # reference line
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(  # titles and axis labels
    title = "Regional Coverage Comparison by Quadrant",
    x = "Geographic Quadrant",
    y = "MMR Coverage Rate"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(plot.title = element_text(face = "bold", size = 14))

print(p_quadrant)  # show the result

ggsave("figures/exploration/coverage_histogram.png", plot = p_hist, width = 8, height = 5, dpi = 300)  # save the figure
ggsave("figures/exploration/coverage_boxplot.png", plot = p_box, width = 5, height = 6, dpi = 300)  # save the figure
ggsave("figures/exploration/lowest_10_schools.png", plot = p_bar_lowest, width = 9, height = 5, dpi = 300)  # save the figure
ggsave("figures/exploration/geographic_coverage_map.png", plot = p_map, width = 8, height = 6, dpi = 300)  # save the figure
ggsave("figures/exploration/risk_stratification.png", plot = p_risk, width = 8, height = 5, dpi = 300)  # save the figure
ggsave("figures/exploration/coverage_by_category.png", plot = p_text_type, width = 8, height = 5, dpi = 300)  # save the figure
ggsave("figures/exploration/coverage_by_quadrant.png", plot = p_quadrant, width = 8, height = 5, dpi = 300)  # save the figure
