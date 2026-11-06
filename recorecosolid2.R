# ============================================================
# 13. PLOT 6 PANELS WITH DUAL STACKED LEGENDS
# ============================================================

# --- 1. Update Scenario label from 91.6 to 96.3 in both datasets ---
hub_data$Scenario  <- gsub("91\\.6", "96.3", hub_data$Scenario)
peri_data$Scenario <- gsub("91\\.6", "96.3", peri_data$Scenario)

# If Scenario is a factor, make sure the factor levels preserve the order:
scen_levels <- unique(c(hub_data$Scenario, peri_data$Scenario))
hub_data$Scenario  <- factor(hub_data$Scenario, levels = scen_levels)
peri_data$Scenario <- factor(peri_data$Scenario, levels = scen_levels)

# --- 2. Generate Plot ---
p_2x3_grid <- ggplot() +
  # Top Row: Hub lines (Solid lines, custom color legend for Hubs)
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
  
  # Reset scale to generate a separate legend for Peripherals
  new_scale_color() +
  
  # Bottom Row: Peripheral lines (Changed to SOLID lines)
  geom_line(
    data = peri_data,
    aes(
      x = time,
      y = total_infectious,
      color = School_Name,
      group = School_Name
    ),
    linetype = "solid",      # <-- Changed from "dashed" to "solid"
    linewidth = 1
  ) +
  scale_color_manual(
    name = "Peripheral Schools",
    values = peri_color_map
  ) +
  
  # 2x3 Grid Faceting (Hubs top row, Peripherals bottom row)
  facet_grid(
    Seed_Type ~ Scenario,
    scales = "free_y"
  ) +
  
  labs(
    title = "Epidemic Trajectories Across Coverage Scenarios",
    subtitle = "Top Row: Hubs | Bottom Row: Peripherals",  # <-- Updated subtitle
    x = "Time (Days)",
    y = "Number of Infectious Cases"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",             # Stack the legends onto two separate lines
    legend.box.just = "center",          # Center align both legend lines
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 13)
  )

print(p_2x3_grid)