if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)
library(ggplot2)  # plots
library(dplyr)  # data wrangling
library(tibble)  # simple data frames

seed_schools <- bind_rows(  # combine these rows
  hubs_df %>% mutate(Seed_Type = "Hub Schools"),
  periph_df %>% mutate(Seed_Type = "Peripheral Schools")
)

vax_scenarios <- tibble(
  Scenario_Name = c("Observed", "74.3%", "96.3%"),
  vax_val       = c(NA, 0.743, 0.963)
)

hub_names  <- unique(hubs_df$school_name)
peri_names <- unique(periph_df$school_name)

hub_color_map  <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A"), hub_names)
peri_color_map <- setNames(c("#984EA3", "#FF7F00", "#A65628"), peri_names)

all_trajectories <- list()
all_outcomes     <- list()

for (i in seq_len(nrow(seed_schools))) {  # one school at a time

  s_id   <- seed_schools$school_id[i]
  s_name <- seed_schools$school_name[i]
  s_type <- seed_schools$Seed_Type[i]

  for (j in seq_len(nrow(vax_scenarios))) {  # repeat for each item

    scen_name   <- vax_scenarios$Scenario_Name[j]
    target_vax  <- vax_scenarios$vax_val[j]

    vax_current <- df$vax_rate

    if (!is.na(target_vax)) {
      vax_current[s_id] <- target_vax
    }

    sim_out <- run_vax_scenario(vax_current, s_id)

    df_ts <- get_total_I(sim_out, scen_name, s_type, s_name) %>%
      mutate(
        School_Name = s_name,
        Seed_Type   = s_type,
        Scenario    = scen_name
      )
    all_trajectories[[length(all_trajectories) + 1]] <- df_ts

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

all_ts_data <- bind_rows(all_trajectories) %>%  # combine these rows
  mutate(
    Scenario  = factor(Scenario, levels = c("Observed", "74.3%", "96.3%")),
    Seed_Type = factor(Seed_Type, levels = c("Hub Schools", "Peripheral Schools"))
  )

hub_data  <- all_ts_data %>% filter(Seed_Type == "Hub Schools")
peri_data <- all_ts_data %>% filter(Seed_Type == "Peripheral Schools")

table_outcomes_summary <- bind_rows(all_outcomes)  # combine these rows

p_2x3_grid <- ggplot() +  # start the plot

  geom_line(  # draw the time series
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
  scale_color_manual(  # set the legend colours or line types
    name   = "Hub Schools",
    values = hub_color_map
  ) +

  new_scale_color() +

  geom_line(  # draw the time series
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
  scale_color_manual(  # set the legend colours or line types
    name   = "Peripheral Schools",
    values = peri_color_map
  ) +

  facet_grid(  # one panel per group
    Seed_Type ~ Scenario,
    scales = "free_y"
  ) +

  labs(  # titles and axis labels
    title    = "Epidemic Trajectories Across Coverage Scenarios",
    subtitle = "Top Row: Hub Schools (Solid) | Bottom Row: Peripheral Schools (Solid)",
    x        = "Time (Days)",
    y        = "Number of Infectious Cases"
  ) +

  theme_minimal() +  # plain plot theme
  theme(
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.box.just   = "center",
    strip.text        = element_text(face = "bold", size = 11),
    strip.background  = element_rect(fill = "grey95", color = "grey80"),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 10)
  )

print(p_2x3_grid)  # show the result

cat("\n\n============================================================\n")  # print a check
cat("SUMMARY TABLE: HUBS VS PERIPHERALS OUTCOMES (18 SIMULATIONS)\n")  # print a check
cat("============================================================\n\n")  # print a check

print(table_outcomes_summary, n = Inf, width = Inf)  # show the result

write.csv(  # save the table
  table_outcomes_summary,
  "table_4_2_hubs_vs_peripherals_outcomes.csv",
  row.names = FALSE
)

cat("\nOutcomes table exported to 'table_4_2_hubs_vs_peripherals_outcomes.csv'\n")  # print a check
