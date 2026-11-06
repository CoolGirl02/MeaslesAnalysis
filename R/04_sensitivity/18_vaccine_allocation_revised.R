required_packages <- c(
  "deSolve",
  "dplyr",
  "readr",
  "tibble",
  "tidyr",
  "ggplot2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(deSolve)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(ggplot2)

school_file  <- "data/processed/school_data_joined.csv"
network_file <- "data/processed/network_df.csv"

if (!file.exists(school_file) || !file.exists(network_file)) {
  stop(
    "Input CSV files not found. Check the working directory and file paths."
  )
}

school_data <- read_csv(
  school_file,
  show_col_types = FALSE
)

network_data <- read_csv(
  network_file,
  show_col_types = FALSE
)

standardize_names <- function(data) {

  names(data) <- tolower(trimws(names(data)))

  names(data) <- gsub(
    "[[:space:]-]+",
    "_",
    names(data)
  )

  names(data) <- make.unique(
    names(data),
    sep = "_duplicate_"
  )

  return(data)
}

school_data  <- standardize_names(school_data)
network_data <- standardize_names(network_data)

standardize_column <- function(data, target_name, possible_names) {

  if (target_name %in% names(data)) {
    return(data)
  }

  matches <- intersect(
    possible_names,
    names(data)
  )

  if (length(matches) == 0) {
    stop(
      "Could not identify column: ",
      target_name,
      ". Available columns are: ",
      paste(names(data), collapse = ", ")
    )
  }

  data <- dplyr::rename(
    data,
    !!target_name := all_of(matches[1])
  )

  return(data)
}

school_data <- standardize_column(
  school_data,
  "unit_id",
  c("unitid", "school_id")
)

school_data <- standardize_column(
  school_data,
  "school_name",
  c("school", "schoolname", "name")
)

school_data <- standardize_column(
  school_data,
  "total_students",
  c("total_student", "enrollment", "school_size")
)

school_data <- standardize_column(
  school_data,
  "percent_immunized",
  c("mmr_pr", "mmr_coverage", "vaccination_coverage")
)

network_data <- standardize_column(
  network_data,
  "unit_id_a",
  c("unitid_a", "school_id_a")
)

network_data <- standardize_column(
  network_data,
  "unit_id_b",
  c("unitid_b", "school_id_b")
)

network_data <- standardize_column(
  network_data,
  "connect_risk",
  c("total_contacts", "weight", "risk")
)

school_data <- school_data %>%
  mutate(
    unit_id = trimws(as.character(unit_id)),
    school_name = as.character(school_name),
    total_students = suppressWarnings(
      as.numeric(total_students)
    ),
    percent_immunized = suppressWarnings(
      as.numeric(percent_immunized)
    )
  ) %>%
  filter(
    !is.na(unit_id),
    unit_id != "",
    !is.na(total_students),
    is.finite(total_students),
    total_students > 0
  ) %>%
  distinct(
    unit_id,
    .keep_all = TRUE
  )

if (nrow(school_data) < 2) {
  stop("Insufficient valid schools in the school dataset.")
}

if (anyDuplicated(school_data$unit_id)) {
  stop("Duplicate school IDs remain after cleaning.")
}

N <- school_data$total_students

n_school <- nrow(school_data)

school_names <- school_data$school_name

school_ids <- school_data$unit_id

total_population <- sum(N)

cat("\n============================================================\n")
cat("SCHOOL DATA SUMMARY\n")
cat("============================================================\n")

cat("Number of schools:", n_school, "\n")

cat("Total school population:", total_population, "\n")

cat(
  "Missing observed vaccination percentages:",
  sum(is.na(school_data$percent_immunized)),
  "\n"
)

network_data <- network_data %>%
  mutate(
    unit_id_a = trimws(as.character(unit_id_a)),
    unit_id_b = trimws(as.character(unit_id_b)),
    connect_risk = suppressWarnings(
      as.numeric(connect_risk)
    )
  ) %>%
  filter(
    unit_id_a %in% school_ids,
    unit_id_b %in% school_ids,
    unit_id_a != unit_id_b,
    is.finite(connect_risk),
    connect_risk > 0
  ) %>%
  mutate(
    index_a = match(unit_id_a, school_ids),
    index_b = match(unit_id_b, school_ids)
  )

network_edges <- network_data %>%
  mutate(
    node_1 = pmin(index_a, index_b),
    node_2 = pmax(index_a, index_b)
  ) %>%
  group_by(node_1, node_2) %>%
  summarise(
    weight = mean(connect_risk),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(weight),
    weight > 0
  )

if (nrow(network_edges) == 0) {
  stop(
    "No positive network edges remain after cleaning. ",
    "Check school IDs and connectivity weights."
  )
}

A_raw <- matrix(
  0,
  nrow = n_school,
  ncol = n_school
)

for (k in seq_len(nrow(network_edges))) {

  i <- network_edges$node_1[k]

  j <- network_edges$node_2[k]

  w <- network_edges$weight[k]

  A_raw[i, j] <- w

  A_raw[j, i] <- w
}

spectral_radius <- max(
  Mod(
    eigen(
      A_raw,
      only.values = TRUE
    )$values
  )
)

if (!is.finite(spectral_radius) || spectral_radius <= 0) {
  stop(
    "The network has a non-positive or invalid spectral radius."
  )
}

A <- A_raw / spectral_radius

cat("\n============================================================\n")
cat("NETWORK SUMMARY\n")
cat("============================================================\n")

cat(
  "Number of positive undirected school connections:",
  nrow(network_edges),
  "\n"
)

cat(
  "Spectral radius before normalization:",
  spectral_radius,
  "\n"
)

seed_school_index <- 1

if (seed_school_index < 1 || seed_school_index > n_school) {
  stop("Invalid seed school index.")
}

gamma_rate <- 1 / 7

sigma_rate <- 1 / 10

beta_rate <- 15 * gamma_rate

alpha <- 0.05

simulation_times <- seq(
  0,
  200,
  by = 1
)

network_seir <- function(t, state, parms) {

  N <- parms$N

  A <- parms$A

  beta_rate <- parms$beta_rate

  sigma_rate <- parms$sigma_rate

  gamma_rate <- parms$gamma_rate

  alpha <- parms$alpha

  n <- length(N)

  S <- state[seq_len(n)]

  E <- state[n + seq_len(n)]

  I <- state[2 * n + seq_len(n)]

  R <- state[3 * n + seq_len(n)]

  infectious_fraction <- I / N

  lambda <- beta_rate * infectious_fraction +
    alpha * beta_rate *
    as.vector(A %*% infectious_fraction)

  dS <- -lambda * S

  dE <- lambda * S - sigma_rate * E

  dI <- sigma_rate * E - gamma_rate * I

  dR <- gamma_rate * I

  return(
    list(
      c(dS, dE, dI, dR)
    )
  )
}

run_scenario <- function(
    scenario_name,
    vax_vector,
    cutoff_val,
    small_schools_vec
) {

  if (length(vax_vector) != n_school) {
    stop(
      "Vaccination vector length does not match the number of schools."
    )
  }

  if (any(!is.finite(vax_vector)) ||
      any(vax_vector < 0 | vax_vector > 1)) {

    stop(
      "Invalid vaccination coverage in scenario: ",
      scenario_name
    )
  }

  R_init <- N * vax_vector

  E_init <- rep(
    0,
    n_school
  )

  I_init <- rep(
    0,
    n_school
  )

  I_init[seed_school_index] <- 1

  S_init <- N - R_init - E_init - I_init

  if (any(S_init < -1e-8)) {

    stop(
      "Negative initial susceptible population in scenario: ",
      scenario_name,
      ". Check vaccination coverage and the seed school."
    )
  }

  S_init[abs(S_init) < 1e-10] <- 0

  initial_population_check <- S_init + E_init + I_init + R_init

  if (!isTRUE(
    all.equal(
      as.numeric(initial_population_check),
      as.numeric(N),
      tolerance = 1e-8
    )
  )) {

    stop(
      "Initial population conservation failed in scenario: ",
      scenario_name
    )
  }

  initial_state <- c(
    S_init,
    E_init,
    I_init,
    R_init
  )

  names(initial_state) <- c(
    paste0("S", seq_len(n_school)),
    paste0("E", seq_len(n_school)),
    paste0("I", seq_len(n_school)),
    paste0("R", seq_len(n_school))
  )

  model_parameters <- list(
    N = N,
    A = A,
    beta_rate = beta_rate,
    sigma_rate = sigma_rate,
    gamma_rate = gamma_rate,
    alpha = alpha
  )

  simulation <- as.data.frame(
    ode(
      y = initial_state,
      times = simulation_times,
      func = network_seir,
      parms = model_parameters,
      method = "lsoda"
    )
  )

  S_matrix <- as.matrix(
    simulation[
      ,
      paste0("S", seq_len(n_school)),
      drop = FALSE
    ]
  )

  I_matrix <- as.matrix(
    simulation[
      ,
      paste0("I", seq_len(n_school)),
      drop = FALSE
    ]
  )

  total_I <- rowSums(I_matrix)

  peak_row <- which.max(total_I)

  peak_burden <- total_I[peak_row]

  peak_day <- simulation$time[peak_row]

  final_S <- S_matrix[nrow(S_matrix), ]

  cumulative_new_cases <- sum(
    S_init - final_S
  )

  county_attack_rate <- (
    cumulative_new_cases / total_population
  ) * 100

  aggregate_coverage <- (
    sum(N * vax_vector) / total_population
  ) * 100

  small_population <- sum(
    N[small_schools_vec]
  )

  large_population <- sum(
    N[!small_schools_vec]
  )

  small_cov <- (
    sum(
      N[small_schools_vec] *
        vax_vector[small_schools_vec]
    ) / small_population
  ) * 100

  large_cov <- (
    sum(
      N[!small_schools_vec] *
        vax_vector[!small_schools_vec]
    ) / large_population
  ) * 100

  summary_row <- tibble(

    `Cutoff Threshold` = cutoff_val,

    Scenario = scenario_name,

    `Small School Vax (%)` = round(
      small_cov,
      2
    ),

    `Large School Vax (%)` = round(
      large_cov,
      2
    ),

    `Aggregate Vax (%)` = round(
      aggregate_coverage,
      2
    ),

    `Total New Cases (C)` = round(
      cumulative_new_cases
    ),

    `Peak Burden (Ipeak)` = round(
      peak_burden
    ),

    `Peak Day` = round(
      peak_day
    ),

    `Attack Rate (AR)` = paste0(
      round(county_attack_rate, 2),
      "%"
    )
  )

  time_series <- tibble(

    time = simulation$time,

    total_infectious = total_I,

    Scenario = scenario_name,

    Cutoff = paste0(
      "Cutoff: ",
      cutoff_val,
      " students"
    )
  )

  return(
    list(
      summary = summary_row,
      time_series = time_series
    )
  )
}

cutoff_thresholds <- c(
  700,
  1000,
  1300
)

summary_results <- list()

time_series_data <- list()

for (cutoff in cutoff_thresholds) {

  cat("\n============================================================\n")

  cat(
    "PROCESSING SCHOOL-SIZE CUTOFF:",
    cutoff,
    "STUDENTS\n"
  )

  cat("============================================================\n")

  small_schools <- N <= cutoff

  large_schools <- N > cutoff

  if (sum(small_schools) == 0 ||
      sum(large_schools) == 0) {

    warning(
      "Skipping cutoff ",
      cutoff,
      " because one school-size group is empty."
    )

    next
  }

  pop_small <- sum(
    N[small_schools]
  )

  pop_large <- sum(
    N[large_schools]
  )

  pop_total <- sum(N)

  cat(
    "Number of small schools:",
    sum(small_schools),
    "\n"
  )

  cat(
    "Number of large schools:",
    sum(large_schools),
    "\n"
  )

  cat(
    "Small-school population:",
    pop_small,
    "\n"
  )

  cat(
    "Large-school population:",
    pop_large,
    "\n"
  )

  vax_uniform_70 <- rep(
    0.70,
    n_school
  )

  vax_uniform_80 <- rep(
    0.80,
    n_school
  )

  vax_uniform_90 <- rep(
    0.90,
    n_school
  )

  small_cov_A <- 0.75

  large_cov_A <- 0.90

  vax_size_A <- ifelse(
    small_schools,
    small_cov_A,
    large_cov_A
  )

  anchor_aggregate <- (
    pop_small * small_cov_A +
      pop_large * large_cov_A
  ) / pop_total

  large_cov_B <- 0.75

  small_cov_B <- (
    anchor_aggregate * pop_total -
      large_cov_B * pop_large
  ) / pop_small

  if (!is.finite(small_cov_B) ||
      small_cov_B < 0 ||
      small_cov_B > 1) {

    stop(
      "Infeasible allocation at cutoff ",
      cutoff,
      ". Calculated small-school vaccination coverage is ",
      round(small_cov_B * 100, 2),
      "%. The allocation strategy must be revised."
    )
  }

  vax_size_B <- ifelse(
    small_schools,
    small_cov_B,
    large_cov_B
  )

  aggregate_A <- sum(
    N * vax_size_A
  ) / pop_total

  aggregate_B <- sum(
    N * vax_size_B
  ) / pop_total

  if (!isTRUE(
    all.equal(
      aggregate_A,
      aggregate_B,
      tolerance = 1e-10
    )
  )) {

    stop(
      "Aggregate vaccination coverage does not match at cutoff ",
      cutoff,
      "."
    )
  }

  cat("\nAllocation scenarios:\n")

  cat(
    "Scenario A - Small schools:",
    round(small_cov_A * 100, 2),
    "%; Large schools:",
    round(large_cov_A * 100, 2),
    "%\n"
  )

  cat(
    "Scenario B - Small schools:",
    round(small_cov_B * 100, 2),
    "%; Large schools:",
    round(large_cov_B * 100, 2),
    "%\n"
  )

  cat(
    "Matched aggregate vaccination coverage:",
    round(aggregate_A * 100, 2),
    "%\n"
  )

  sc_70 <- run_scenario(
    scenario_name = "Uniform 70%",
    vax_vector = vax_uniform_70,
    cutoff_val = cutoff,
    small_schools_vec = small_schools
  )

  sc_80 <- run_scenario(
    scenario_name = "Uniform 80%",
    vax_vector = vax_uniform_80,
    cutoff_val = cutoff,
    small_schools_vec = small_schools
  )

  sc_size_A <- run_scenario(
    scenario_name = "Large Emphasis (Small 75%, Large 90%)",
    vax_vector = vax_size_A,
    cutoff_val = cutoff,
    small_schools_vec = small_schools
  )

  sc_size_B <- run_scenario(
    scenario_name = paste0(
      "Small Emphasis (Small ",
      round(small_cov_B * 100, 2),
      "%, Large 75%)"
    ),
    vax_vector = vax_size_B,
    cutoff_val = cutoff,
    small_schools_vec = small_schools
  )

  sc_90 <- run_scenario(
    scenario_name = "Uniform 90%",
    vax_vector = vax_uniform_90,
    cutoff_val = cutoff,
    small_schools_vec = small_schools
  )

  summary_results[[as.character(cutoff)]] <- bind_rows(

    sc_70$summary,

    sc_80$summary,

    sc_size_A$summary,

    sc_size_B$summary,

    sc_90$summary
  )

  time_series_data[[as.character(cutoff)]] <- bind_rows(

    sc_70$time_series,

    sc_80$time_series,

    sc_size_A$time_series,

    sc_size_B$time_series,

    sc_90$time_series
  )
}

if (length(summary_results) == 0) {
  stop(
    "No simulations were completed. Check the school-size cutoffs and data."
  )
}

final_table <- bind_rows(
  summary_results
)

final_ts <- bind_rows(
  time_series_data
)

small_emphasis_names <- unique(
  final_table$Scenario[
    grepl("^Small Emphasis", final_table$Scenario)
  ]
)

scenario_order <- c(
  "Uniform 70%",
  "Uniform 80%",
  "Large Emphasis (Small 75%, Large 90%)",
  small_emphasis_names,
  "Uniform 90%"
)

final_table <- final_table %>%
  mutate(
    Scenario = factor(
      Scenario,
      levels = scenario_order
    )
  ) %>%
  arrange(
    `Cutoff Threshold`,
    Scenario
  ) %>%
  mutate(
    Scenario = as.character(Scenario)
  )

final_ts <- final_ts %>%
  mutate(
    Scenario = factor(
      Scenario,
      levels = scenario_order
    ),
    Cutoff = factor(
      Cutoff,
      levels = paste0(
        "Cutoff: ",
        cutoff_thresholds,
        " students"
      )
    )
  )

cat("\n\n")

cat("======================================================================\n")

cat("VACCINE ALLOCATION SENSITIVITY ANALYSIS\n")

cat("SCHOOL-SIZE CUTOFFS: 700, 1000, AND 1300 STUDENTS\n")

cat("======================================================================\n\n")

print(
  final_table,
  n = Inf,
  width = Inf
)

summary_output <- "output/tables/vaccine_allocation_threshold_summary_table.csv"

write.csv(
  final_table,
  summary_output,
  row.names = FALSE
)

cat("\nSummary table exported to:\n")

cat(
  normalizePath(
    summary_output,
    mustWork = FALSE
  ),
  "\n"
)

final_ts <- final_ts %>%
  mutate(
    Strategy_Category = case_when(
      grepl("^Uniform 70%", Scenario) ~ "Uniform 70%",
      grepl("^Uniform 80%", Scenario) ~ "Uniform 80%",
      grepl("^Large Emphasis", Scenario) ~ "Large Emphasis",
      grepl("^Small Emphasis", Scenario) ~ "Small Emphasis",
      grepl("^Uniform 90%", Scenario) ~ "Uniform 90%",
      TRUE ~ as.character(Scenario)
    ),

    Line_Style = case_when(
      Strategy_Category == "Small Emphasis" & grepl("700", Cutoff)  ~ "dashed",
      Strategy_Category == "Small Emphasis" & grepl("1000", Cutoff) ~ "dotted",
      Strategy_Category == "Small Emphasis" & grepl("1300", Cutoff) ~ "dotdash",
      TRUE ~ "solid"
    )
  )

strategy_colours <- c(
  "Uniform 70%"    = "#B95F89",
  "Uniform 80%"    = "#E41A1C",
  "Large Emphasis" = "#008000",
  "Small Emphasis" = "#0000FF",
  "Uniform 90%"    = "#5C3317"
)

strategy_linetypes <- c(
  "Uniform 70%"    = "solid",
  "Uniform 80%"    = "solid",
  "Large Emphasis" = "solid",
  "Small Emphasis" = "dashed",
  "Uniform 90%"    = "solid"
)

p_epidemic_curves <- ggplot(
  final_ts,
  aes(
    x = time,
    y = total_infectious,
    colour = Strategy_Category,
    linetype = Line_Style,
    group = Scenario
  )
) +
  geom_line(linewidth = 1.1) +
  facet_wrap(
    ~Cutoff,
    ncol = 3
  ) +
  scale_colour_manual(
    name = "Vaccination Strategy",
    values = strategy_colours
  ) +
  scale_linetype_identity() +
  labs(
    title = "Measles Epidemic Curves Across School-Size Cutoffs",
    subtitle = "Comparison of vaccination allocation strategies with matched aggregate coverage",
    x = "Time (Days)",
    y = "Number of Infectious Students"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 11),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.key.width = grid::unit(1.2, "cm"),

    panel.grid.minor = element_blank()
  )

print(p_epidemic_curves)
