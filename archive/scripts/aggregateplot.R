pop_summary <- results_long %>%
  filter(Compartment != "R") %>%
  group_by(time, Compartment) %>%
  summarise(Total_Count = sum(Count), .groups = "drop")

p_pop <- ggplot(  # start the plot
  pop_summary,
  aes(
    x = time,
    y = Total_Count,
    color = Compartment,
    linetype = Compartment
  )
) +
  geom_line(linewidth = 1.2) +  # draw the time series
  scale_color_manual(  # set the legend colours or line types
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
  scale_linetype_manual(  # set the legend colours or line types
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
  labs(  # titles and axis labels
    title = paste0(
      "Total Population SEI Dynamics\nInfection Introduced at ",
      df$school_name[seed_index]
    ),
    x = "Days",
    y = "Total Number of Students",
    color = "Compartment",
    linetype = "Compartment"
  ) +
  theme_minimal() +  # plain plot theme
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  )

print(p_pop)  # show the result

total_county_students <- sum(df$total_students, na.rm = TRUE)

total_vaccinated_students <- sum(df$total_students * df$vax_rate, na.rm = TRUE)

county_vax_coverage <- total_vaccinated_students / total_county_students

cat("\n=== COUNTY-WIDE SUMMARY METRICS ===\n")  # print a check
cat("Total Number of Schools:      ", nrow(df), "\n")  # print a check
cat("Total Student Population:      ", format(total_county_students, big.mark = ","), "\n")  # print a check
cat("Total Vaccinated Students:     ", format(round(total_vaccinated_students), big.mark = ","), "\n")  # print a check
cat("Overall Vaccination Coverage:  ", round(county_vax_coverage * 100, 2), "%\n")  # print a check
cat("===================================\n\n")  # print a check
