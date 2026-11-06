# ==============================================================================
# SENSITIVITY ANALYSIS: VACCINE ALLOCATION BY SCHOOL SIZE VS UNIFORM
# ==============================================================================

library(deSolve)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------------------
# 1. LOAD & CLEAN DATA
# ------------------------------------------------------------------------------
school_file  <- "data/processed/school_data_joined.csv"
network_file <- "data/processed/network_df.csv"

if (!file.exists(school_file) || !file.exists(network_file)) {
  stop("Input CSV files missing from data/processed/")
}

school_data  <- read_csv(school_file, show_col_types = FALSE)
network_data <- read_csv(network_file, show_col_types = FALSE)

standardize_names <- function(data) {
  names(data) <- tolower(trimws(names(data)))
  names(data) <- gsub("[[:space:]-]+", "_", names(data))
  names(data) <- make.unique(names(data), sep = "_duplicate_")
  return(data)
}

school_data  <- standardize_names(school_data)
network_data <- standardize_names(network_data)

# Helper for standardizing column names
rename_if_found <- function(df, target, candidates) {
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

N <- school_data$total_students
n_school <- nrow(school_data)
school_ids <- school_data$unit_id
school_names <- school_data$school_name

# Size definitions
small_schools <- N <= 700
large_schools <- N > 700

# ------------------------------------------------------------------------------
# 2. NETWORK ADJACENCY MATRIX
# ------------------------------------------------------------------------------
network_data <- network_data %>%
  mutate(unit_id_a = trimws(as.character(unit_id_a)),
         unit_id_b = trimws(as.character(unit_id_b)),
         connect_risk = as.numeric(connect_risk)) %>%
  filter(unit_id_a %in% school_ids, unit_id_b %in% school_ids, unit_id_a != unit_id_b)

network_edges <- network_data %>%
  mutate(node_1 = pmin(match(unit_id_a, school_ids), match(unit_id_b, school_ids)),
         node_2 = pmax(match(unit_id_a, school_ids), match(unit_id_b, school_ids))) %>%
  group_by(node_1, node_2) %>%
  summarise(weight = mean(connect_risk, na.rm = TRUE), .groups = "drop")

A_raw <- matrix(0, nrow = n_school, ncol = n_school)
for (k in seq_len(nrow(network_edges))) {
  i <- network_edges$node_1[k]
  j <- network_edges$node_2[k]
  w <- network_edges$weight[k]
  A_raw[i, j] <- w
  A_raw[j, i] <- w
}

spectral_radius <- max(Mod(eigen(A_raw, only.values = TRUE)$values))
A <- A_raw / spectral_radius

# ------------------------------------------------------------------------------
# 3. SEIR MODEL & SIMULATION SETUP
# ------------------------------------------------------------------------------
gamma_rate <- 1 / 7
sigma_rate <- 1 / 10
beta_rate  <- 15 * gamma_rate
alpha      <- 0.05
simulation_times <- seq(0, 200, by = 1)
seed_school_index <- 1

network_seir <- function(t, state, parms) {
  with(as.list(parms), {
    S <- state[seq_len(n_school)]
    E <- state[n_school + seq_len(n_school)]
    I <- state[2 * n_school + seq_len(n_school)]
    R <- state[3 * n_school + seq_len(n_school)]
    
    infectious_fraction <- I / N
    lambda <- beta_rate * infectious_fraction + alpha * beta_rate * as.vector(A %*% infectious_fraction)
    
    dS <- -lambda * S
    dE <- lambda * S - sigma_rate * E
    dI <- sigma_rate * E - gamma_rate * I
    dR <- gamma_rate * I
    
    return(list(c(dS, dE, dI, dR)))
  })
}

run_simulation <- function(vax_vector) {
  R_init <- N * vax_vector
  E_init <- rep(0, n_school)
  I_init <- rep(0, n_school)
  I_init[seed_school_index] <- 1
  S_init <- pmax(0, N - R_init - E_init - I_init)
  
  initial_state <- c(S_init, E_init, I_init, R_init)
  parms <- list(N = N, A = A, beta_rate = beta_rate, sigma_rate = sigma_rate, 
                gamma_rate = gamma_rate, alpha = alpha, n_school = n_school)
  
  out <- as.data.frame(ode(y = initial_state, times = simulation_times, 
                           func = network_seir, parms = parms, method = "lsoda"))
  
  # Extract total infectious across all schools per time point
  I_cols <- paste0("1", seq_len(n_school)) # pattern matching for I
  I_matrix <- as.matrix(out[, (2 * n_school + 2):(3 * n_school + 1)])
  
  tibble(
    Day = out$time,
    Infectious_Students = rowSums(I_matrix)
  )
}

# ------------------------------------------------------------------------------
# 4. DEFINE SCENARIOS FOR SENSITIVITY PLOTTING
# ------------------------------------------------------------------------------

# Solve required counterpart coverage matching population aggregate:
calc_counterpart <- function(target_agg, fixed_cov, fixed_is_small) {
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

for (cut in cutoffs) {
  cut_pct <- paste0(round(cut * 100), "%")
  
  # 1. Uniform Scenario
  vax_unif <- rep(cut, n_school)
  sim_unif <- run_simulation(vax_unif)
  sim_unif <- sim_unif %>% mutate(
    Cutoff = factor(cut_pct, levels = c("70%", "80%", "90%")),
    Allocation_Type = "Uniform",
    Scenario_Label = paste("Uniform", cut_pct)
  )
  scenarios_list[[length(scenarios_list) + 1]] <- sim_unif
  
  # 2. Small Emphasis (Small high, Large adjusted)
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
  
  # 3. Large Emphasis (Large high, Small adjusted)
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

plot_data <- bind_rows(scenarios_list)

# ------------------------------------------------------------------------------
# 5. GENERATE & DISPLAY PLOT
# ------------------------------------------------------------------------------

epidemic_plot <- ggplot(
  plot_data,
  aes(
    x = Day,
    y = Infectious_Students,
    color = Allocation_Type,
    linetype = Cutoff
  )
) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ Cutoff, nrow = 1) +  # Displays left to right in order (70%, 80%, 90%)
  scale_color_manual(
    values = c(
      "Uniform"                = "#2b5c8f", # Consistent Blue
      "Small School Emphasis" = "#e06d53", # Consistent Red/Coral
      "Large School Emphasis" = "#2ca02c"  # Consistent Green
    ),
    name = "Allocation Type"
  ) +
  scale_linetype_manual(
    values = c(
      "70%" = "solid",
      "80%" = "dashed",
      "90%" = "dotted"
    ),
    name = "Coverage Cutoff"
  ) +
  labs(
    title = "Epidemic Curves Across Vaccination Cutoffs & Allocation Strategies",
    subtitle = "Panels ordered left-to-right by overall aggregate vaccination cutoff",
    x = "Time (days)",
    y = "Total Infectious Students"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.text = element_text(face = "bold", size = 11)
  )

print(epidemic_plot)

ggsave(
  "figures/sensitivity/vaccine_allocation_sensitivity_facet.png",
  plot = epidemic_plot,
  width = 12,
  height = 5,
  dpi = 300
)