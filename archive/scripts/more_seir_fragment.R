for (group_name in names(analysis_groups)) {

  selected_ids <- analysis_groups[[group_name]]

  if (length(selected_ids) == 0) {
    cat("\nNo schools selected for:", group_name, "\n")
    next
  }

  page_data <- results_long %>%
    filter(
      school_id %in% selected_ids,
      Compartment != "R"
    ) %>%
    mutate(
      plot_order = match(school_id, selected_ids)
    ) %>%
    arrange(plot_order)

  p_plot <- ggplot(
    page_data,
    aes(
      x = time,
      y = Count,
      color = Compartment,
      linetype = Compartment
    )
  ) +

    geom_line(linewidth = 0.8) +

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

    facet_wrap(
      ~school_label,
      ncol = 3,
      scales = "free_y"
    ) +

    scale_y_continuous(
      labels = scales::label_comma(),
      expand = expansion(mult = c(0, 0.05))
    ) +

    labs(
      title = paste0(
        group_name,
        "\nInfection introduced at ",
        df$school_name[seed_index]
      ),
      x = "Days",
      y = "Number of Students",
      color = "Compartment",
      linetype = "Compartment"
    ) +

    theme_minimal() +

    theme(
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  print(p_plot)
}

cat("\nAll five grouped school plots have been generated.\n")
