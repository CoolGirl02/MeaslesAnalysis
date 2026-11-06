# ============================================================
# REQUIRED PACKAGES
# ============================================================
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)
library(ggplot2)
library(dplyr)
library(tibble)

# ============================================================
# 11. DEFINE THE 6 SEED SCHOOLS & 3 COVERAGE SCENARIOS
# ============================================================

# Combine the 3 Hubs and 3 Peripherals into one list of seed schools
seed_schools <- bind_rows(
  hubs_df %>% mutate(Seed_Type = "Hub Schools"),
  periph_df %>% mutate(Seed_Type = "Peripheral Schools")
)

# 3 Vaccination Scenarios (Only modifying the seed school)
vax_scenarios <- tibble(
  Scenario_Name = c("Observed", "74.3%", "96.3%"),
  vax_val       = c(NA, 0.743, 0.963)
)

# Color maps for individual school lines
hub_names  <- unique(hubs_df$school_name)
peri_names <- unique(periph_df$school_name)

hub_color_map  <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A"), hub_names)
peri_color_map <- setNames(c("#984EA3", "#FF7F00", "#A65628"), peri_names)


# ============================================================
# 12. RUN 18 SIMULATIONS (6 SCHOOLS x 3 SCENARIOS)
# ============================================================

all_trajectories <- list()
all_outcomes     <- list()

# Loop through each of the 6 individual seed schools
for (i in seq_len(nrow(seed_schools))) {
  
  s_id   <- seed_schools$school_id[i]
  s_name <- seed_schools$school_name[i]
  s_type <- seed_schools$Seed_Type[i]
  
  # Loop through the 3 vaccination scenarios
  for (j in seq_len(nrow(vax_scenarios))) {
    
    scen_name   <- vax_scenarios$Scenario_Name[j]
    target_vax  <- vax_scenarios$vax_val[j]
    
    # Start with observed coverage for all 89 schools
    vax_current <- df$vax_rate
    
    # Modify ONLY the seed school's coverage (if not "Observed")
    if (!is.na(target_vax)) {
      vax_current[s_id] <- target_vax
    }
    
    # Run SEIR simulation seeded in this specific school
    sim_out <- run_vax_scenario(vax_current, s_id)
    
    # 1. Extract Time Series for Plotting
    df_ts <- get_total_I(sim_out, scen_name, s_type, s_name) %>%
      mutate(
        School_Name = s_name,
        Seed_Type   = s_type,
        Scenario    = scen_name
      )
    all_trajectories[[length(all_trajectories) + 1]] <- df_ts
    
    # 2. Extract Key Epidemic Metrics for the Outcome Table
    metrics <- calculate_outbreak_measures(sim_out)
    
    all_outcomes[[length(all_outcomes) + 1]] <- tibble(
      `Seed School`   = s_name,
      `School Type`   = ifelse(s_type == "Hub Schools", "Hub", "Peripheral"),
      `Scenario`      = scen_name,
      `Seed Coverage` = paste0(round(vax_current[s_id] * 100, 1), "%"),
      `Total Cases`   = round(metrics["Total_Cases"]),
      `Peak Burden`   = round(metrics["Peak_Burden"]),
      `Peak Day`      = round(metrics["Peak_Day"]),
      `Attack Rate`   = paste0(round(metrics["Attack_Rate"], 2), "%")
    )
  }
}

# Combine and lock scenario order
all_ts_data <- bind_rows(all_trajectories) %>%
  mutate(
    Scenario  = factor(Scenario, levels = c("Observed", "74.3%", "96.3%")),
    Seed_Type = factor(Seed_Type, levels = c("Hub Schools", "Peripheral Schools"))
  )

# Separate into Hubs and Peripherals for the dual-legend plot
hub_data  <- all_ts_data %>% filter(Seed_Type == "Hub Schools")
peri_data <- all_ts_data %>% filter(Seed_Type == "Peripheral Schools")

# Combine the outcomes summary table
table_outcomes_summary <- bind_rows(all_outcomes)


# ============================================================
# 13. PLOT 6 PANELS (2x3 GRID) WITH SOLID LINES & DUAL LEGENDS
# ============================================================

p_2x3_grid <- ggplot() +
  # Top Row: Hub lines (Solid lines, custom legend for Hubs)
  geom_line(
    data = hub_data,
    aes(
      x     = time,
      y     = total_infectious,
      color = School_Name,
      group = School_Name
    ),
    linetype  = "solid",
    linewidth = 1
  ) +
  scale_color_manual(
    name   = "Hub Schools",
    values = hub_color_map
  ) +
  
  # Reset scale to generate a separate legend for Peripherals
  new_scale_color() +
  
  # Bottom Row: Peripheral lines (Solid lines, custom legend for Peripherals)
  geom_line(
    data = peri_data,
    aes(
      x     = time,
      y     = total_infectious,
      color = School_Name,
      group = School_Name
    ),
    linetype  = "solid",
    linewidth = 1
  ) +
  scale_color_manual(
    name   = "Peripheral Schools",
    values = peri_color_map
  ) +
  
  # 2x3 Grid: Rows = Seed_Type (Hubs vs Peripherals), Columns = Scenario
  facet_grid(
    Seed_Type ~ Scenario,
    scales = "free_y"
  ) +
  
  labs(
    title    = "Epidemic Trajectories Across Coverage Scenarios",
    subtitle = "Top Row: Hub Schools (Solid) | Bottom Row: Peripheral Schools (Solid)",
    x        = "Time (Days)",
    y        = "Number of Infectious Cases"
  ) +
  
  theme_minimal() +
  theme(
    legend.position   = "bottom",
    legend.box        = "vertical",          # Stack the two legends
    legend.box.just   = "center",
    strip.text        = element_text(face = "bold", size = 11),
    strip.background  = element_rect(fill = "grey95", color = "grey80"),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 10)
  )

# Display the plot
print(p_2x3_grid)


# ============================================================
# 14. DISPLAY & EXPORT OUTCOMES SUMMARY TABLE
# ============================================================

cat("\n\n============================================================\n")
cat("SUMMARY TABLE: HUBS VS PERIPHERALS OUTCOMES (18 SIMULATIONS)\n")
cat("============================================================\n\n")

print(table_outcomes_summary, n = Inf, width = Inf)

# Export to CSV for your thesis records
write.csv(
  table_outcomes_summary,
  "table_4_2_hubs_vs_peripherals_outcomes.csv",
  row.names = FALSE
)

cat("\nOutcomes table exported to 'table_4_2_hubs_vs_peripherals_outcomes.csv'\n")