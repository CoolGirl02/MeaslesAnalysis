pair_combinations <- tibble(
  Combo_ID   = 1:3,
  Hub_ID     = hubs_df$school_id,
  Peri_ID    = periph_df$school_id,
  Hub_Name   = hubs_df$school_name,
  Peri_Name  = periph_df$school_name,
  Pair_Label = paste("Pair", 1:3)
)

vax_configs <- tibble(
  Config_Name = c("Observed - Observed", "74.3% - 74.3%", "91.6% - 91.6%"),
  vax_val     = c(NA, 0.743, 0.916)
)

comparison_all_list <- list()
results_table_list  <- list()

for (i in seq_len(nrow(pair_combinations))) {  # one school at a time

  pair_label <- pair_combinations$Pair_Label[i]
  h_id       <- pair_combinations$Hub_ID[i]
  p_id       <- pair_combinations$Peri_ID[i]
  h_name     <- pair_combinations$Hub_Name[i]
  p_name     <- pair_combinations$Peri_Name[i]

  for (j in seq_len(nrow(vax_configs))) {  # repeat for each item

    c_name <- vax_configs$Config_Name[j]
    v_val  <- vax_configs$vax_val[j]

    vax_current <- df$vax_rate

    if (!is.na(v_val)) {
      vax_current[h_id] <- v_val
      vax_current[p_id] <- v_val
    }

    out_hub  <- run_vax_scenario(vax_current, h_id)
    out_peri <- run_vax_scenario(vax_current, p_id)

    ts_combined <- bind_rows(  # combine these rows
      get_total_I(out_hub,  c_name, "Network Hub", h_name),
      get_total_I(out_peri, c_name, "Peripheral", p_name)
    ) %>%
      mutate(Pair = pair_label)

    comparison_all_list[[length(comparison_all_list) + 1]] <- ts_combined

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

plot_data_3x3 <- bind_rows(comparison_all_list) %>%  # combine these rows
  mutate(
    Scenario = factor(Scenario, levels = vax_configs$Config_Name),
    Pair     = factor(Pair, levels = pair_combinations$Pair_Label)
  )

table_3x3_summary <- bind_rows(results_table_list)  # combine these rows

p_3x3_grid <- ggplot(  # start the plot
  plot_data_3x3,
  aes(
    x = time,
    y = total_infectious,
    color = `Seed Type`,
    linetype = `Seed Type`
  )
) +
  geom_line(linewidth = 1.1) +  # draw the time series

  facet_grid(  # one panel per group
    Pair ~ Scenario,
    scales = "free_y"
  ) +

  scale_color_manual(  # set the legend colours or line types
    values = c(
      "Network Hub" = "firebrick",
      "Peripheral"  = "darkblue"
    )
  ) +

  scale_linetype_manual(  # set the legend colours or line types
    values = c(
      "Network Hub" = "solid",
      "Peripheral"  = "solid"
    )
  ) +

  labs(  # titles and axis labels
    title = "Epidemic Curves by Hub-Peripheral Pairs Across Coverage Scenarios",
    subtitle = "Comparing trajectories of the top 3 Hubs vs bottom 3 Peripherals across 3 vaccination strategies",
    x = "Days Since Outbreak Introduction",
    y = "Number of Infectious Cases",
    color = "Seed Category",
    linetype = "Seed Category"
  ) +

  theme_minimal() +  # plain plot theme
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 12)
  )

print(p_3x3_grid)  # show the result

cat("\n\n============================================================\n")  # print a check
cat("SUMMARY TABLE: 9-PANEL SENSITIVITY ANALYSIS\n")  # print a check
cat("============================================================\n\n")  # print a check

print(table_3x3_summary, n = Inf)  # show the result
