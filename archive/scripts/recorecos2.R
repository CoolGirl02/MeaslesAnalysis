p_2x3_grid <- ggplot() +

  geom_line(
    data = hub_data,
    aes(
      x = time,
      y = total_infectious,
      color = School_Name,
      group = School_Name
    ),
    linetype = "solid",
    linewidth = 1
  ) +
  scale_color_manual(
    name = "Hub Schools",
    values = hub_color_map
  ) +

  new_scale_color() +

  geom_line(
    data = peri_data,
    aes(
      x = time,
      y = total_infectious,
      color = School_Name,
      group = School_Name
    ),
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_color_manual(
    name = "Peripheral Schools",
    values = peri_color_map
  ) +

  facet_grid(
    Seed_Type ~ Scenario,
    scales = "free_y"
  ) +

  labs(
    title = "Epidemic Trajectories Across Coverage Scenarios",
    subtitle = "Top Row: Hubs (Solid) | Bottom Row: Peripherals (Dashed)",
    x = "Time (Days)",
    y = "Number of Infectious Cases"
  ) +

  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 13)
  )

print(p_2x3_grid)
