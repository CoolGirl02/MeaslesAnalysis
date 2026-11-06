# ============================================================
# 19. POPULATION-WIDE AGGREGATION & PLOTTING (EXCLUDING R)
# ============================================================

# Aggregate counts across all schools for each time point and compartment (excluding R)
pop_summary <- results_long %>%
  filter(Compartment != "R") %>%
  group_by(time, Compartment) %>%
  summarise(Total_Count = sum(Count), .groups = "drop")

# Plot total population dynamics for S, E, and I only
p_pop <- ggplot(
  pop_summary,
  aes(
    x = time,
    y = Total_Count,
    color = Compartment,
    linetype = Compartment
  )
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c(
      "S" = "blue",
      "E" = "orange",
      "I" = "red"
    ),
    breaks = c("S", "E", "I"),
    labels = c(
      "S" = "Susceptible (S)",
      "E" = "Exposed (E)",
      "I" = "Infectious (I)"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "S" = "solid",
      "E" = "dashed",
      "I" = "dotted"
    ),
    breaks = c("S", "E", "I"),
    labels = c(
      "S" = "Susceptible (S)",
      "E" = "Exposed (E)",
      "I" = "Infectious (I)"
    )
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = paste0(
      "Total Population SEI Dynamics\nInfection Introduced at ",
      df$school_name[seed_index]
    ),
    x = "Days",
    y = "Total Number of Students",
    color = "Compartment",
    linetype = "Compartment"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  )

print(p_pop)

# ============================================================
# COUNTY & OVERALL POPULATION VACCINATION SUMMARY
# ============================================================

# 1. Total student population across all schools
total_county_students <- sum(df$total_students, na.rm = TRUE)

# 2. Total vaccinated students (handling raw rates or percentages)
total_vaccinated_students <- sum(df$total_students * df$vax_rate, na.rm = TRUE)

# 3. Overall county vaccination coverage rate
county_vax_coverage <- total_vaccinated_students / total_county_students

# Print overall summary metrics
cat("\n=== COUNTY-WIDE SUMMARY METRICS ===\n")
cat("Total Number of Schools:      ", nrow(df), "\n")
cat("Total Student Population:      ", format(total_county_students, big.mark = ","), "\n")
cat("Total Vaccinated Students:     ", format(round(total_vaccinated_students), big.mark = ","), "\n")
cat("Overall Vaccination Coverage:  ", round(county_vax_coverage * 100, 2), "%\n")
cat("===================================\n\n")