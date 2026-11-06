# ============================================================
# 11. PAIR THE 3 SELECTED HUBS AND 3 PERIPHERALS
# ============================================================
# Building off 'hubs_df' and 'periph_df' defined in Section 7

# Create 3 distinct Hub-Peripheral pairs 

pair_combinations <- tibble(
  Combo_ID   = 1:3,
  Hub_ID     = hubs_df$school_id,
  Peri_ID    = periph_df$school_id,
  Hub_Name   = hubs_df$school_name,
  Peri_Name  = periph_df$school_name,
  Pair_Label = paste("Pair", 1:3)
)

# Define the 3 vaccination configurations for the grid columns

vax_configs <- tibble(
  Config_Name = c("Observed - Observed", "74.3% - 74.3%", "91.6% - 91.6%"),
  vax_val     = c(NA, 0.743, 0.916)
)


# ============================================================
# 12. RUN SIMULATIONS FOR THE 3x3 MATRIX (9 PLOTS)
# ============================================================

comparison_all_list <- list()
results_table_list  <- list()

# Loop through each Hub-Peripheral pair (Rows)
for (i in seq_len(nrow(pair_combinations))) {
  
  pair_label <- pair_combinations$Pair_Label[i]
  h_id       <- pair_combinations$Hub_ID[i]
  p_id       <- pair_combinations$Peri_ID[i]
  h_name     <- pair_combinations$Hub_Name[i]
  p_name     <- pair_combinations$Peri_Name[i]
  
  # Loop through each Vaccination Configuration (Columns)
  for (j in seq_len(nrow(vax_configs))) {
    
    c_name <- vax_configs$Config_Name[j]
    v_val  <- vax_configs$vax_val[j]
    
    # Base vaccination vector
    vax_current <- df$vax_rate 
    
    # If not "Observed", assign the targeted coverage to BOTH the hub and peripheral
    if (!is.na(v_val)) {
      vax_current[h_id] <- v_val
      vax_current[p_id] <- v_val
    }
    
    # Run ODE Simulations
    out_hub  <- run_vax_scenario(vax_current, h_id)
    out_peri <- run_vax_scenario(vax_current, p_id)
    
    # Extract time series and label them for the facet grid
    ts_combined <- bind_rows(
      get_total_I(out_hub,  c_name, "Network Hub", h_name),
      get_total_I(out_peri, c_name, "Peripheral", p_name)
    ) %>%
      mutate(Pair = pair_label)
    
    comparison_all_list[[length(comparison_all_list) + 1]] <- ts_combined
    
    # Calculate metrics for the summary table
    metrics_hub  <- calculate_outbreak_measures(out_hub)
    metrics_peri <- calculate_outbreak_measures(out_peri)
    
    results_table_list[[length(results_table_list) + 1]] <- tibble(
      Pair = pair_label,
      Scenario = c_name,
      `Seed Location` = c(h_name, p_name),
      `Seed Type` = c("Network Hub", "Peripheral"),
      `Coverage` = c(paste0(round(vax_current[h_id] * 100, 1), "%"), paste0(round(vax_current[p_id] * 100, 1), "%")),
      `Total Cases` = c(round(metrics_hub["Total_Cases"]), round(metrics_peri["Total_Cases"])),
      `Peak Burden` = c(round(metrics_hub["Peak_Burden"]), round(metrics_peri["Peak_Burden"])),
      `Attack Rate` = c(paste0(round(metrics_hub["Attack_Rate"], 2), "%"), paste0(round(metrics_peri["Attack_Rate"], 2), "%"))
    )
  }
}

# Combine all into one dataframe and lock factor levels for proper grid ordering
plot_data_3x3 <- bind_rows(comparison_all_list) %>%
  mutate(
    Scenario = factor(Scenario, levels = vax_configs$Config_Name),
    Pair     = factor(Pair, levels = pair_combinations$Pair_Label)
  )

table_3x3_summary <- bind_rows(results_table_list)


# ============================================================
# 13. PLOT THE 3x3 GRID (9 PANELS)
# ============================================================

p_3x3_grid <- ggplot(
  plot_data_3x3,
  aes(
    x = time,
    y = total_infectious,
    color = `Seed Type`,
    linetype = `Seed Type`
  )
) +
  geom_line(linewidth = 1.1) +
  
  # facet_grid(Rows ~ Columns)
  facet_grid(
    Pair ~ Scenario,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      "Network Hub" = "firebrick",
      "Peripheral"  = "darkblue"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "Network Hub" = "solid",
      "Peripheral"  = "solid"
    )
  ) +
  
  labs(
    title = "Epidemic Curves by Hub-Peripheral Pairs Across Coverage Scenarios",
    subtitle = "Comparing trajectories of the top 3 Hubs vs bottom 3 Peripherals across 3 vaccination strategies",
    x = "Days Since Outbreak Introduction",
    y = "Number of Infectious Cases",
    color = "Seed Category",
    linetype = "Seed Category"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 12)
  )

print(p_3x3_grid)


# ============================================================
# 14. DISPLAY SUMMARY TABLE
# ============================================================

cat("\n\n============================================================\n")
cat("SUMMARY TABLE: 9-PANEL SENSITIVITY ANALYSIS\n")
cat("============================================================\n\n")

print(table_3x3_summary, n = Inf)