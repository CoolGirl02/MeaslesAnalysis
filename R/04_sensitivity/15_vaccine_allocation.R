library(deSolve)  # ODE solver
library(dplyr)  # data wrangling
library(readr)  # read and write tables
library(tibble)  # simple data frames
library(ggplot2)  # plots
library(tidyr)  # reshape tables

school_file  <- "data/processed/school_data_joined.csv"  # path to the joined school table
network_file <- "data/processed/network_df.csv"  # path to the edge list

if (!file.exists(school_file) || !file.exists(network_file)) {  # check that the input files exist
  stop("Input CSV files missing from data/processed/")  # stop if this check fails
}

school_data  <- read_csv(school_file, show_col_types = FALSE)  # read the table
network_data <- read_csv(network_file, show_col_types = FALSE)  # read the table

standardize_names <- function(data) {  # define the function
  names(data) <- tolower(trimws(names(data)))
  names(data) <- gsub("[[:space:]-]+", "_", names(data))
  names(data) <- make.unique(names(data), sep = "_duplicate_")
  return(data)
}

school_data  <- standardize_names(school_data)
network_data <- standardize_names(network_data)

rename_if_found <- function(df, target, candidates) {  # define the function
  match <- intersect(candidates, names(df))
  if (length(match) > 0 && !target %in% names(df)) {
    df <- dplyr::rename(df, !!target := all_of(match[1]))
  }
  return(df)
}

school_data <- rename_if_found(school_data, "unit_id", c("unitid", "school_id"))
school_data <- rename_if_found(school_data, "school_name", c("school", "schoolname", "name"))
school_data <- rename_if_found(school_data, "total_students", c("total_student", "enrollment", "student_population"))
school_data <- rename_if_found(school_data, "percent_immunized", c("mmr_pr", "vaccination_coverage", "percent_vaccinated"))

network_data <- rename_if_found(network_data, "unit_id_a", c("unitid_a", "school_id_a"))
network_data <- rename_if_found(network_data, "unit_id_b", c("unitid_b", "school_id_b"))
network_data <- rename_if_found(network_data, "connect_risk", c("connection_risk", "weight", "risk"))

school_data <- school_data %>%
  mutate(
    unit_id = trimws(as.character(unit_id)),
    total_students = as.numeric(total_students),
    percent_immunized = as.numeric(percent_immunized)
  ) %>%
  filter(!is.na(unit_id), total_students > 0, !is.na(percent_immunized)) %>%
  distinct(unit_id, .keep_all = TRUE)

N <- school_data$total_students  # students in each school
n_school <- nrow(school_data)  # number of schools
school_ids <- school_data$unit_id
school_names <- school_data$school_name

small_schools <- N <= 700  # schools with at most 700 students
large_schools <- N > 700  # schools with more than 700 students

network_data <- network_data %>%
  mutate(unit_id_a = trimws(as.character(unit_id_a)),
         unit_id_b = trimws(as.character(unit_id_b)),
         connect_risk = as.numeric(connect_risk)) %>%
  filter(unit_id_a %in% school_ids, unit_id_b %in% school_ids, unit_id_a != unit_id_b)

network_edges <- network_data %>%
  mutate(node_1 = pmin(match(unit_id_a, school_ids), match(unit_id_b, school_ids)),  # smaller school id first
         node_2 = pmax(match(unit_id_a, school_ids), match(unit_id_b, school_ids))) %>%  # larger school id second
  group_by(node_1, node_2) %>%
  summarise(weight = mean(connect_risk, na.rm = TRUE), .groups = "drop")

A_raw <- matrix(0, nrow = n_school, ncol = n_school)  # unscaled school-to-school weights
for (k in seq_len(nrow(network_edges))) {  # one edge at a time
  i <- network_edges$node_1[k]
  j <- network_edges$node_2[k]
  w <- network_edges$weight[k]
  A_raw[i, j] <- w  # undirected edge weight
  A_raw[j, i] <- w  # undirected edge weight
}

spectral_radius <- max(Mod(eigen(A_raw, only.values = TRUE)$values))  # largest absolute eigenvalue of the network
A <- A_raw / spectral_radius  # divide every edge by the spectral radius

gamma_rate <- 1 / 7  # recovery rate, 7-day infectious period
sigma_rate <- 1 / 10  # latency rate, 10-day incubation
beta_rate  <- 15 * gamma_rate  # within-school transmission, R0 about 15
alpha      <- 0.05  # fraction of transmission coming from other schools
simulation_times <- seq(0, 200, by = 1)  # days to simulate
seed_school_index <- 1  # school where the index case starts

network_seir <- function(t, state, parms) {  # SEIR equations for one time step
  with(as.list(parms), {
    S <- state[seq_len(n_school)]
    E <- state[n_school + seq_len(n_school)]
    I <- state[2 * n_school + seq_len(n_school)]
    R <- state[3 * n_school + seq_len(n_school)]

    infectious_fraction <- I / N  # infectious fraction in each school
    lambda <- beta_rate * infectious_fraction + alpha * beta_rate * as.vector(A %*% infectious_fraction)  # force of infection from this school and its neighbours

    dS <- -lambda * S  # infections leave the susceptible class
    dE <- lambda * S - sigma_rate * E  # new exposures minus those becoming infectious
    dI <- sigma_rate * E - gamma_rate * I  # new infectious students minus recoveries
    dR <- gamma_rate * I  # students leaving the infectious class

    return(list(c(dS, dE, dI, dR)))  # hand the derivatives to the solver
  })
}

run_simulation <- function(vax_vector) {  # define the function
  R_init <- N * vax_vector  # students already immune
  E_init <- rep(0, n_school)  # exposed at the start, usually none
  I_init <- rep(0, n_school)  # index cases at the start
  I_init[seed_school_index] <- 1
  S_init <- pmax(0, N - R_init - E_init - I_init)  # students still susceptible

  initial_state <- c(S_init, E_init, I_init, R_init)  # starting S, E, I, and R
  parms <- list(N = N, A = A, beta_rate = beta_rate, sigma_rate = sigma_rate,  # parameters passed into the solver
                gamma_rate = gamma_rate, alpha = alpha, n_school = n_school)

  out <- as.data.frame(ode(y = initial_state, times = simulation_times,
                           func = network_seir, parms = parms, method = "lsoda"))

  I_cols <- paste0("1", seq_len(n_school))
  I_matrix <- as.matrix(out[, (2 * n_school + 2):(3 * n_school + 1)])

  tibble(
    Day = out$time,
    Infectious_Students = rowSums(I_matrix)
  )
}

calc_counterpart <- function(target_agg, fixed_cov, fixed_is_small) {  # define the function
  pop_small <- sum(N[small_schools])
  pop_large <- sum(N[large_schools])
  pop_total <- sum(N)

  if (fixed_is_small) {
    return((target_agg * pop_total - fixed_cov * pop_small) / pop_large)
  } else {
    return((target_agg * pop_total - fixed_cov * pop_large) / pop_small)
  }
}

cutoffs <- c(0.70, 0.80, 0.90)

scenarios_list <- list()

for (cut in cutoffs) {  # repeat for each item
  cut_pct <- paste0(round(cut * 100), "%")

  vax_unif <- rep(cut, n_school)
  sim_unif <- run_simulation(vax_unif)
  sim_unif <- sim_unif %>% mutate(
    Cutoff = factor(cut_pct, levels = c("70%", "80%", "90%")),
    Allocation_Type = "Uniform",
    Scenario_Label = paste("Uniform", cut_pct)
  )
  scenarios_list[[length(scenarios_list) + 1]] <- sim_unif

  high_small <- min(1.0, cut + 0.05)
  adj_large  <- calc_counterpart(cut, high_small, fixed_is_small = TRUE)
  vax_small  <- ifelse(small_schools, high_small, adj_large)
  sim_small  <- run_simulation(vax_small)
  sim_small  <- sim_small %>% mutate(
    Cutoff = factor(cut_pct, levels = c("70%", "80%", "90%")),
    Allocation_Type = "Small School Emphasis",
    Scenario_Label = paste("Small Emphasis", cut_pct)
  )
  scenarios_list[[length(scenarios_list) + 1]] <- sim_small

  high_large <- min(1.0, cut + 0.05)
  adj_small  <- calc_counterpart(cut, high_large, fixed_is_small = FALSE)
  vax_large  <- ifelse(large_schools, high_large, adj_small)
  sim_large  <- run_simulation(vax_large)
  sim_large  <- sim_large %>% mutate(
    Cutoff = factor(cut_pct, levels = c("70%", "80%", "90%")),
    Allocation_Type = "Large School Emphasis",
    Scenario_Label = paste("Large Emphasis", cut_pct)
  )
  scenarios_list[[length(scenarios_list) + 1]] <- sim_large
}

plot_data <- bind_rows(scenarios_list)  # combine these rows

epidemic_plot <- ggplot(
  plot_data,
  aes(
    x = Day,
    y = Infectious_Students,
    color = Allocation_Type,
    linetype = Cutoff
  )
) +
  geom_line(linewidth = 1.0) +  # draw the time series
  facet_wrap(~ Cutoff, nrow = 1) +  # one panel per group
  scale_color_manual(  # set the legend colours or line types
    values = c(
      "Uniform"                = "#2b5c8f",
      "Small School Emphasis" = "#e06d53",
      "Large School Emphasis" = "#2ca02c"
    ),
    name = "Allocation Type"
  ) +
  scale_linetype_manual(  # set the legend colours or line types
    values = c(
      "70%" = "solid",
      "80%" = "dashed",
      "90%" = "dotted"
    ),
    name = "Coverage Cutoff"
  ) +
  labs(  # titles and axis labels
    title = "Epidemic Curves Across Vaccination Cutoffs & Allocation Strategies",
    subtitle = "Panels ordered left-to-right by overall aggregate vaccination cutoff",
    x = "Time (days)",
    y = "Total Infectious Students"
  ) +
  theme_minimal(base_size = 12) +  # plain plot theme
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.text = element_text(face = "bold", size = 11)
  )

print(epidemic_plot)  # show the result

ggsave(  # save the figure
  "figures/sensitivity/vaccine_allocation_sensitivity_facet.png",
  plot = epidemic_plot,
  width = 12,
  height = 5,
  dpi = 300
)
