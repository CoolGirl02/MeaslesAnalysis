strategy_colours <- c(
  "Uniform Matched" = "#E41A1C",
  "Large Emphasis"  = "#008000",
  "Small Emphasis"  = "#0000FF"
)

final_ts <- final_ts %>%
  mutate(
    Strategy_Category = case_when(
      grepl("^Uniform Matched", Scenario) ~ "Uniform Matched",
      grepl("^Large Emphasis", Scenario)  ~ "Large Emphasis",
      grepl("^Small Emphasis", Scenario)  ~ "Small Emphasis",
      TRUE ~ as.character(Scenario)
    ),

    Strategy_Category = factor(
      Strategy_Category,
      levels = c("Uniform Matched", "Large Emphasis", "Small Emphasis")
    )
  )

p_epidemic_curves <- ggplot(
  final_ts,
  aes(
    x = time,
    y = total_infectious,
    colour = Strategy_Category,
    group = Scenario
  )
) +
  geom_line(
    linetype = "solid",
    linewidth = 1.1
  ) +
  facet_wrap(
    ~Cutoff,
    ncol = 3
  ) +
  scale_colour_manual(
    name = "Vaccination Strategy",
    values = strategy_colours
  ) +
  labs(
    title = "Measles Epidemic Curves Across School-Size Cutoffs",
    subtitle = "Comparison of vaccination allocation strategies with matched aggregate coverage",
    x = "Time (Days)",
    y = "Number of Infectious Students"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 11),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.key.width = grid::unit(1.2, "cm"),

    panel.grid.minor = element_blank()
  )

print(p_epidemic_curves)
