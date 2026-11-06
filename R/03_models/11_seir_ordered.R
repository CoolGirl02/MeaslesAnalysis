library(tidyverse)  # data wrangling and plots
library(deSolve)  # ODE solver

df <- read.csv("data/processed/school_data_joined.csv")  # read the table

names(df) <- tolower(names(df))

if (!"unit_id" %in% names(df)) {  # add a school id if the column is missing
  df$unit_id <- as.character(seq_len(nrow(df)))
} else {
  df$unit_id <- as.character(df$unit_id)
}

required_columns <- c(
  "unit_id",
  "school_name",
  "total_students",
  "percent_immunized"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(  # stop if this check fails
    paste(
      "Missing required school-data columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

df <- df %>%
  mutate(
    vax_rate = if_else(
      percent_immunized > 1,  # values above 1 are percentages, so divide by 100
      percent_immunized / 100,
      percent_immunized
    ),
    school_id = row_number()  # index the schools
  )

n_schools <- nrow(df)  # number of schools

seed_school_name <- "Global Academy of SC"

school_names_clean <- trimws(tolower(as.character(df$school_name)))

seed_index <- which(
  school_names_clean == trimws(tolower(seed_school_name))
)

if (length(seed_index) == 0) {

  possible_matches <- df$school_name[
    grepl(
      "global academy",
      as.character(df$school_name),
      ignore.case = TRUE
    )
  ]

  if (length(possible_matches) > 0) {
    print(possible_matches)  # show the result
  }

  stop("Seed school not found. Check the school name.")  # stop if this check fails

}

if (length(seed_index) > 1) {
  stop("More than one school matches Global Academy of SC.")  # stop if this check fails
}

cat(  # print a check
  "Infection starts at:",
  df$school_name[seed_index],
  "\n"
)

net_df <- read.csv("data/processed/network_df.csv")  # read the table

names(net_df) <- tolower(names(net_df))

required_network_columns <- c(
  "unit_id_a",
  "unit_id_b",
  "total_contacts"
)

missing_network_columns <- setdiff(
  required_network_columns,
  names(net_df)
)

if (length(missing_network_columns) > 0) {
  stop(  # stop if this check fails
    paste(
      "Missing required network-data columns:",
      paste(missing_network_columns, collapse = ", ")
    )
  )
}

net_df <- net_df %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  )

id_map <- setNames(seq_len(n_schools), df$unit_id)

A_raw <- matrix(  # unscaled school-to-school weights
  0,
  nrow = n_schools,
  ncol = n_schools
)

for (k in seq_len(nrow(net_df))) {  # one edge at a time

  i <- unname(id_map[net_df$unit_id_a[k]])
  j <- unname(id_map[net_df$unit_id_b[k]])

  if (
    length(i) == 1 &&
    length(j) == 1 &&
    !is.na(i) &&
    !is.na(j) &&
    i != j
  ) {

    w <- net_df$total_contacts[k]

    if (is.na(w) || !is.finite(w) || w < 0) {
      next
    }

    A_raw[i, j] <- w  # undirected edge weight
    A_raw[j, i] <- w  # undirected edge weight

  }
}

df$connectivity <- rowSums(A_raw)

cat("\nSchools with the highest weighted connectivity:\n")  # print a check

print(  # show the result
  df %>%
    arrange(desc(connectivity)) %>%
    select(
      school_name,
      percent_immunized,
      total_students,
      connectivity
    ) %>%
    head(10)
)

lowest_vax_ids <- df %>%
  arrange(vax_rate, school_name) %>%
  slice_head(n = min(6, n_schools)) %>%
  pull(school_id)

highest_vax_ids <- df %>%
  arrange(desc(vax_rate), school_name) %>%
  slice_head(n = min(6, n_schools)) %>%
  pull(school_id)

intermediate_ids <- df %>%
  filter(
    vax_rate >= 0.75,
    vax_rate < 0.95,
    !school_id %in% c(lowest_vax_ids, highest_vax_ids)
  ) %>%
  arrange(vax_rate, school_name) %>%
  pull(school_id)

highest_connectivity_ids <- df %>%
  arrange(desc(connectivity), school_name) %>%
  slice_head(n = min(6, n_schools)) %>%
  pull(school_id)

cat("\nSIX SCHOOLS WITH THE LOWEST VACCINATION COVERAGE\n")  # print a check

print(  # show the result
  df %>%
    filter(school_id %in% lowest_vax_ids) %>%
    arrange(vax_rate) %>%
    select(school_name, percent_immunized, total_students)
)

cat("\nSIX SCHOOLS WITH THE HIGHEST VACCINATION COVERAGE\n")  # print a check

print(  # show the result
  df %>%
    filter(school_id %in% highest_vax_ids) %>%
    arrange(desc(vax_rate)) %>%
    select(school_name, percent_immunized, total_students)
)

cat("\nSCHOOLS WITH INTERMEDIATE VACCINATION COVERAGE\n")  # print a check

print(  # show the result
  df %>%
    filter(school_id %in% intermediate_ids) %>%
    arrange(vax_rate) %>%
    select(school_name, percent_immunized, total_students)
)

cat("\nSIX SCHOOLS WITH THE HIGHEST CONNECTIVITY\n")  # print a check

print(  # show the result
  df %>%
    filter(school_id %in% highest_connectivity_ids) %>%
    arrange(desc(connectivity)) %>%
    select(
      school_name,
      percent_immunized,
      total_students,
      connectivity
    )
)

spectral_radius <- max(Mod(eigen(A_raw, only.values = TRUE)$values))  # largest absolute eigenvalue of the network

if (!is.finite(spectral_radius) || spectral_radius <= 0) {  # stop when the spectral radius is not usable
  stop("The network spectral radius is zero or invalid.")  # stop if this check fails
}

A <- A_raw / spectral_radius  # divide every edge by the spectral radius

gamma_rate <- 1 / 7  # recovery rate, 7-day infectious period
sigma_rate <- 1 / 10  # latency rate, 10-day incubation
beta_rate <- 15 * gamma_rate  # within-school transmission, R0 about 15
alpha <- 0.05  # fraction of transmission coming from other schools

times <- seq(0, 200, by = 1)  # days to simulate

N <- df$total_students  # students in each school

R0_count <- round(N * df$vax_rate)  # students already immune

I0_count <- rep(0, n_schools)  # index cases at the start

I0_count[seed_index] <- 1

E0_count <- rep(0, n_schools)  # exposed at the start, usually none

S0_count <- pmax(  # students still susceptible
  0,
  N - I0_count - E0_count - R0_count
)

cat("\nInitial infectious students:", sum(I0_count), "\n")  # print a check

initial_state <- c(  # starting S, E, I, and R
  S0_count,
  E0_count,
  I0_count,
  R0_count
)

names(initial_state) <- c(
  paste0("S_", seq_len(n_schools)),
  paste0("E_", seq_len(n_schools)),
  paste0("I_", seq_len(n_schools)),
  paste0("R_", seq_len(n_schools))
)

params <- list(  # parameters passed into the solver
  beta = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N = N,
  A = A,
  n = n_schools
)

seir_model <- function(time, state, parameters) {  # SEIR equations for one time step

  with(as.list(parameters), {

    S <- state[seq_len(n)]

    E <- state[n + seq_len(n)]

    I <- state[2 * n + seq_len(n)]

    R <- state[3 * n + seq_len(n)]

    I_prop <- I / N  # infectious fraction in each school

    lambda <- beta * I_prop +  # force of infection from this school and its neighbours
      alpha * beta * as.vector(A %*% I_prop)

    dS <- -lambda * S  # infections leave the susceptible class

    dE <- lambda * S - sigma * E  # new exposures minus those becoming infectious

    dI <- sigma * E - gamma * I  # new infectious students minus recoveries

    dR <- gamma * I  # students leaving the infectious class

    return(list(c(dS, dE, dI, dR)))  # hand the derivatives to the solver

  })
}

output <- ode(  # solve the SEIR equations
  y = initial_state,
  times = times,
  func = seir_model,
  parms = params
)

cat("SEIR simulation completed successfully.\n")  # print a check

results_long <- as.data.frame(output) %>%  # one row per compartment and day
  pivot_longer(  # one row per compartment
    cols = -time,
    names_to = "var",
    values_to = "Count"
  ) %>%
  mutate(
    Compartment = str_extract(var, "^[SEIR]"),
    school_id = as.numeric(str_extract(var, "\\d+"))
  ) %>%
  left_join(df, by = "school_id") %>%  # join the tables
  mutate(
    school_label = paste0(
      school_name,
      "\n(Vax: ",
      round(vax_rate * 100, 1),
      "%, Pop: ",
      total_students,
      ")"
    )
  )

analysis_groups <- list(

  "Six Schools with the Lowest Vaccination Coverage" =
    lowest_vax_ids,

  "Six Schools with the Highest Vaccination Coverage" =
    highest_vax_ids,

  "Schools with Intermediate Vaccination Coverage (75% to <95%)" =
    intermediate_ids,

  "Six Schools with the Highest Network Connectivity" =
    highest_connectivity_ids

)

for (group_name in names(analysis_groups)) {  # repeat for each item

  selected_ids <- analysis_groups[[group_name]]

  if (length(selected_ids) == 0) {
    cat("\nNo schools selected for:", group_name, "\n")  # print a check
    next
  }

  school_pages <- split(
    selected_ids,
    ceiling(seq_along(selected_ids) / 6)
  )

  for (p in seq_along(school_pages)) {  # repeat for each item

    page_ids <- school_pages[[p]]

    page_data <- results_long %>%  # schools drawn on this page
      filter(school_id %in% page_ids) %>%
      mutate(
        plot_order = match(school_id, page_ids)
      ) %>%
      arrange(plot_order)

    p_plot <- ggplot(  # start the plot
      page_data,
      aes(
        x = time,
        y = Count,
        color = Compartment,
        linetype = Compartment
      )
    ) +

      geom_line(linewidth = 0.8) +  # draw the time series

      scale_color_manual(  # set the legend colours or line types
        values = c(
          "S" = "blue",
          "E" = "orange",
          "I" = "red",
          "R" = "green"
        ),
        breaks = c("S", "E", "I", "R"),
        labels = c(
          "S" = "Susceptible (S)",
          "E" = "Exposed (E)",
          "I" = "Infectious (I)",
          "R" = "Recovered/Vaccinated (R)"
        )
      ) +

      scale_linetype_manual(  # set the legend colours or line types
        values = c(
          "S" = "solid",
          "E" = "dashed",
          "I" = "dotted",
          "R" = "solid"
        ),
        breaks = c("S", "E", "I", "R"),
        labels = c(
          "S" = "Susceptible (S)",
          "E" = "Exposed (E)",
          "I" = "Infectious (I)",
          "R" = "Recovered/Vaccinated (R)"
        )
      ) +

      facet_wrap(  # one panel per group
        ~school_label,
        ncol = 3,
        scales = "free_y"
      ) +

      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        breaks = scales::breaks_pretty(n = 4),
        labels = scales::label_comma(),
        expand = expansion(mult = c(0, 0.05))
      ) +

      labs(  # titles and axis labels
        title = paste0(
          group_name,
          "\nInfection introduced at ",
          df$school_name[seed_index]
        ),

        x = "Days",
        y = "Number of Students (pseudo-log scale)",
        color = "Compartment",
        linetype = "Compartment"
      ) +

      theme_minimal() +  # plain plot theme

      theme(
        strip.text = element_text(
          size = 8,
          face = "bold"
        ),
        legend.position = "bottom"
      )

    print(p_plot)  # show the result

  }
}

cat("\nAll grouped school plots have been generated.\n")  # print a check
