required_packages <- c(
  "deSolve",
  "dplyr",
  "readr",
  "tibble"
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

school_file <- "data/processed/school_data_joined.csv"

network_file <- "data/processed/network_df.csv"

if (!file.exists(school_file)) {

  stop(
    "School data file was not found: ",
    school_file
  )

}

if (!file.exists(network_file)) {

  stop(
    "Network data file was not found: ",
    network_file
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

school_data <- standardize_names(school_data)

network_data <- standardize_names(network_data)

cat("\n====================================================\n")
cat("AVAILABLE DATA COLUMNS\n")
cat("====================================================\n")

cat("\nSchool data columns:\n")

print(names(school_data))

cat("\nNetwork data columns:\n")

print(names(network_data))

standardize_column <- function(
    data,
    target_name,
    possible_names
) {

  if (target_name %in% names(data)) {

    return(data)

  }

  matches <- intersect(
    possible_names,
    names(data)
  )

  if (length(matches) == 0) {

    stop(
      "\nCould not identify the column '",
      target_name,
      "'.\nAvailable columns are: ",
      paste(names(data), collapse = ", "),
      "\nPossible column names checked: ",
      paste(possible_names, collapse = ", ")
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
  c(
    "unitid",
    "unit_id_a",
    "unitid_a",
    "school_id",
    "schoolid",
    "unit_id_school"
  )
)

school_data <- standardize_column(
  school_data,
  "school_name",
  c(
    "school",
    "schoolname",
    "name",
    "school_name_a",
    "schoolname_a"
  )
)

school_data <- standardize_column(
  school_data,
  "total_students",
  c(
    "total_student",
    "total_enrollment",
    "enrollment",
    "student_population",
    "school_size",
    "total_students_a"
  )
)

school_data <- standardize_column(
  school_data,
  "percent_immunized",
  c(
    "mmr_pr",
    "mmr_coverage",
    "vaccination_coverage",
    "percent_vaccinated",
    "percent_mmr",
    "immunization_coverage",
    "mmr_pr_a"
  )
)

network_data <- standardize_column(
  network_data,
  "unit_id_a",
  c(
    "unitid_a",
    "unit_id_a_x",
    "school_id_a",
    "schoolid_a",
    "unit_id_from"
  )
)

network_data <- standardize_column(
  network_data,
  "unit_id_b",
  c(
    "unitid_b",
    "unit_id_b_x",
    "school_id_b",
    "schoolid_b",
    "unit_id_to"
  )
)

network_data <- standardize_column(
  network_data,
  "connect_risk",
  c(
    "connection_risk",
    "connectivity_risk",
    "risk",
    "weight",
    "connection_weight"
  )
)

required_school_columns <- c(
  "unit_id",
  "school_name",
  "total_students",
  "percent_immunized"
)

required_network_columns <- c(
  "unit_id_a",
  "unit_id_b",
  "connect_risk"
)

missing_school_columns <- setdiff(
  required_school_columns,
  names(school_data)
)

missing_network_columns <- setdiff(
  required_network_columns,
  names(network_data)
)

if (length(missing_school_columns) > 0) {

  stop(
    "Missing school columns: ",
    paste(missing_school_columns, collapse = ", ")
  )

}

if (length(missing_network_columns) > 0) {

  stop(
    "Missing network columns: ",
    paste(missing_network_columns, collapse = ", ")
  )

}

cat("\nRequired column names have been identified.\n")

school_data <- school_data %>%

  mutate(

    unit_id = trimws(
      as.character(unit_id)
    ),

    school_name = as.character(
      school_name
    ),

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

    total_students > 0,

    !is.na(percent_immunized),

    is.finite(percent_immunized)

  ) %>%

  distinct(
    unit_id,
    .keep_all = TRUE
  )

if (nrow(school_data) == 0) {

  stop(
    "No valid schools remain after cleaning the school data."
  )

}

if (
  all(
    school_data$percent_immunized >= 0 &
    school_data$percent_immunized <= 1
  )
) {

  school_data$percent_immunized <-
    school_data$percent_immunized * 100

  cat(
    "\nVaccination coverage was converted from proportions ",
    "to percentages.\n",
    sep = ""
  )

}

if (
  any(
    school_data$percent_immunized < 0 |
    school_data$percent_immunized > 100
  )
) {

  stop(
    "Vaccination coverage must be between 0 and 100."
  )

}

school_data$school_name[
  is.na(school_data$school_name) |
    trimws(school_data$school_name) == ""
] <- paste0(
  "School_",
  which(
    is.na(school_data$school_name) |
      trimws(school_data$school_name) == ""
  )
)

N <- school_data$total_students

n_school <- nrow(school_data)

school_names <- school_data$school_name

school_ids <- school_data$unit_id

observed_vax <- school_data$percent_immunized / 100

cat("\n====================================================\n")
cat("SCHOOL DATA SUMMARY\n")
cat("====================================================\n")

cat(
  "Number of schools:",
  n_school,
  "\n"
)

cat(
  "Total school population:",
  sum(N),
  "\n"
)

cat(
  "First school in the dataset:",
  school_names[1],
  "\n"
)

cat(
  "First school's ID:",
  school_ids[1],
  "\n"
)

cat(
  "First school's vaccination coverage:",
  observed_vax[1] * 100,
  "%\n"
)

small_schools <- N <= 700

large_schools <- N > 700

if (!any(small_schools)) {

  stop(
    "No small schools were found using the 700-student cutoff."
  )

}

if (!any(large_schools)) {

  stop(
    "No large schools were found using the 700-student cutoff."
  )

}

cat("\n====================================================\n")
cat("SCHOOL SIZE SUMMARY\n")
cat("====================================================\n")

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
  sum(N[small_schools]),
  "\n"
)

cat(
  "Large-school population:",
  sum(N[large_schools]),
  "\n"
)

network_data <- network_data %>%

  mutate(

    unit_id_a = trimws(
      as.character(unit_id_a)
    ),

    unit_id_b = trimws(
      as.character(unit_id_b)
    ),

    connect_risk = suppressWarnings(
      as.numeric(connect_risk)
    )

  ) %>%

  filter(

    !is.na(unit_id_a),

    !is.na(unit_id_b),

    unit_id_a != "",

    unit_id_b != "",

    !is.na(connect_risk),

    is.finite(connect_risk),

    connect_risk >= 0,

    unit_id_a != unit_id_b

  )

network_data <- network_data %>%

  filter(

    unit_id_a %in% school_ids,

    unit_id_b %in% school_ids

  )

if (nrow(network_data) == 0) {

  stop(
    "No valid network connections remain. ",
    "Check whether the school IDs match between the two files."
  )

}

network_data <- network_data %>%

  mutate(

    index_a = match(
      unit_id_a,
      school_ids
    ),

    index_b = match(
      unit_id_b,
      school_ids
    )

  ) %>%

  filter(

    !is.na(index_a),

    !is.na(index_b),

    index_a != index_b

  )

network_edges <- network_data %>%

  mutate(

    node_1 = pmin(
      index_a,
      index_b
    ),

    node_2 = pmax(
      index_a,
      index_b
    )

  ) %>%

  group_by(
    node_1,
    node_2
  ) %>%

  summarise(

    weight = mean(
      connect_risk,
      na.rm = TRUE
    ),

    .groups = "drop"

  ) %>%

  filter(

    is.finite(weight),

    weight > 0

  )

if (nrow(network_edges) == 0) {

  stop(
    "No positive network weights remain after processing the connections. ",
    "Check the connect_risk values."
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

school_strength <- rowSums(A_raw)

eigenvalues <- eigen(
  A_raw,
  only.values = TRUE
)$values

spectral_radius <- max(
  Mod(eigenvalues)
)

if (
  !is.finite(spectral_radius) ||
  spectral_radius <= 0
) {

  stop(
    "The network spectral radius is zero or invalid."
  )

}

A <- A_raw / spectral_radius

cat("\n====================================================\n")
cat("NETWORK SUMMARY\n")
cat("====================================================\n")

cat(
  "Number of school-to-school connections:",
  nrow(network_edges),
  "\n"
)

cat(
  "Schools with positive network strength:",
  sum(school_strength > 0),
  "\n"
)

cat(
  "Highest school network strength:",
  max(school_strength),
  "\n"
)

seed_school_index <- 1

cat("\n====================================================\n")
cat("INITIAL INFECTION\n")
cat("====================================================\n")

cat(
  "Infection starts in:",
  school_names[seed_school_index],
  "\n"
)

cat(
  "School ID:",
  school_ids[seed_school_index],
  "\n"
)

cat(
  "School size:",
  N[seed_school_index],
  "students\n"
)

cat(
  "Network strength:",
  school_strength[seed_school_index],
  "\n"
)

gamma_rate <- 1 / 7

sigma_rate <- 1 / 10

beta_rate <- 15 * gamma_rate

alpha <- 0.05

simulation_times <- seq(
  0,
  200,
  by = 1
)

cat("\n====================================================\n")
cat("SEIR MODEL PARAMETERS\n")
cat("====================================================\n")

cat(
  "Beta:",
  beta_rate,
  "\n"
)

cat(
  "Sigma:",
  sigma_rate,
  "\n"
)

cat(
  "Gamma:",
  gamma_rate,
  "\n"
)

cat(
  "Alpha:",
  alpha,
  "\n"
)

network_seir <- function(
    t,
    state,
    parms
) {

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
    as.vector(
      A %*% infectious_fraction
    )

  dS <- -lambda * S

  dE <- lambda * S - sigma_rate * E

  dI <- sigma_rate * E - gamma_rate * I

  dR <- gamma_rate * I

  return(
    list(
      c(
        dS,
        dE,
        dI,
        dR
      )
    )
  )

}

run_scenario <- function(
    scenario_name,
    vax_vector
) {

  if (length(vax_vector) != n_school) {

    stop(
      "Vaccination vector has the wrong length in scenario: ",
      scenario_name
    )

  }

  if (
    any(!is.finite(vax_vector)) ||
    any(vax_vector < 0 | vax_vector > 1)
  ) {

    stop(
      "Vaccination coverage must be between 0 and 1 in scenario: ",
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
      ". Check the seed school population and vaccination coverage."
    )

  }

  S_init[abs(S_init) < 1e-10] <- 0

  initial_state <- c(
    S_init,
    E_init,
    I_init,
    R_init
  )

  names(initial_state) <- c(

    paste0(
      "S",
      seq_len(n_school)
    ),

    paste0(
      "E",
      seq_len(n_school)
    ),

    paste0(
      "I",
      seq_len(n_school)
    ),

    paste0(
      "R",
      seq_len(n_school)
    )

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

      method = "lsoda",

      rtol = 1e-8,

      atol = 1e-8

    )

  )

  S_matrix <- as.matrix(

    simulation[
      ,
      paste0(
        "S",
        seq_len(n_school)
      ),
      drop = FALSE
    ]

  )

  E_matrix <- as.matrix(

    simulation[
      ,
      paste0(
        "E",
        seq_len(n_school)
      ),
      drop = FALSE
    ]

  )

  I_matrix <- as.matrix(

    simulation[
      ,
      paste0(
        "I",
        seq_len(n_school)
      ),
      drop = FALSE
    ]

  )

  R_matrix <- as.matrix(

    simulation[
      ,
      paste0(
        "R",
        seq_len(n_school)
      ),
      drop = FALSE
    ]

  )

  minimum_compartment <- min(
    S_matrix,
    E_matrix,
    I_matrix,
    R_matrix
  )

  if (minimum_compartment < -1e-4) {

    warning(
      "A negative compartment value occurred in scenario: ",
      scenario_name,
      ". Minimum value: ",
      signif(minimum_compartment, 5)
    )

  }

  total_infectious <- rowSums(
    I_matrix
  )

  peak_row <- which.max(
    total_infectious
  )

  peak_burden <- total_infectious[
    peak_row
  ]

  peak_day <- simulation$time[
    peak_row
  ]

  final_S <- S_matrix[
    nrow(S_matrix),
  ]

  cumulative_cases_by_school <- pmax(
    0,
    S_init - final_S
  )

  cumulative_new_cases <- sum(
    cumulative_cases_by_school
  )

  county_attack_rate <- (

    cumulative_new_cases / sum(N)

  ) * 100

  aggregate_coverage <- (

    sum(N * vax_vector) / sum(N)

  ) * 100

  small_coverage <- (

    sum(
      N[small_schools] *
        vax_vector[small_schools]
    ) /

      sum(
        N[small_schools]
      )

  ) * 100

  large_coverage <- (

    sum(
      N[large_schools] *
        vax_vector[large_schools]
    ) /

      sum(
        N[large_schools]
      )

  ) * 100

  scenario_summary <- tibble(

    Scenario = scenario_name,

    Small_School_Coverage = small_coverage,

    Large_School_Coverage = large_coverage,

    Aggregate_Coverage = aggregate_coverage,

    Peak_Burden = peak_burden,

    Peak_Day = peak_day,

    Total_New_Cases = cumulative_new_cases,

    Attack_Rate_Percent = county_attack_rate

  )

  school_cases <- tibble(

    unit_id = school_ids,

    school_name = school_names,

    school_size = N,

    school_type = ifelse(
      small_schools,
      "Small",
      "Large"
    ),

    vaccination_coverage = vax_vector * 100,

    network_strength = school_strength,

    initially_infectious =
      seq_len(n_school) == seed_school_index,

    susceptible_initially = S_init,

    susceptible_final = final_S,

    cumulative_new_cases =
      cumulative_cases_by_school

  )

  return(

    list(

      summary = scenario_summary,

      simulation = simulation,

      school_cases = school_cases,

      minimum_compartment = minimum_compartment

    )

  )

}

vax_uniform_70 <- rep(
  0.70,
  n_school
)

vax_uniform_80 <- rep(
  0.80,
  n_school
)

vax_size_A <- ifelse(

  small_schools,

  0.75,

  0.90

)

aggregate_A <- sum(
  N * vax_size_A
) / sum(N)

large_coverage_B <- 0.75

population_small <- sum(
  N[small_schools]
)

population_large <- sum(
  N[large_schools]
)

population_total <- sum(N)

small_coverage_B <- (

  aggregate_A * population_total -

    large_coverage_B * population_large

) / population_small

if (
  !is.finite(small_coverage_B) ||
  small_coverage_B < 0 ||
  small_coverage_B > 1
) {

  stop(

    "The calculated small-school coverage for Scenario 4 is outside ",
    "the valid range of 0% to 100%. Check school sizes and Scenario 3."

  )

}

vax_size_B <- ifelse(

  small_schools,

  small_coverage_B,

  large_coverage_B

)

vax_uniform_90 <- rep(
  0.90,
  n_school
)

calculate_coverage <- function(
    vax_vector,
    school_group
) {

  (

    sum(
      N[school_group] *
        vax_vector[school_group]
    ) /

      sum(
        N[school_group]
      )

  ) * 100

}

scenario_names <- c(

  "Uniform 70%",

  "Uniform 80%",

  "Small 75%, Large 90%",

  paste0(
    "Small ",
    round(small_coverage_B * 100, 2),
    "%, Large 75%"
  ),

  "Uniform 90%"

)

vaccination_vectors <- list(

  vax_uniform_70,

  vax_uniform_80,

  vax_size_A,

  vax_size_B,

  vax_uniform_90

)

scenario_coverage <- bind_rows(

  lapply(

    seq_along(vaccination_vectors),

    function(i) {

      v <- vaccination_vectors[[i]]

      tibble(

        Scenario = scenario_names[i],

        Small_School_Coverage =
          calculate_coverage(
            v,
            small_schools
          ),

        Large_School_Coverage =
          calculate_coverage(
            v,
            large_schools
          ),

        Aggregate_Coverage =
          sum(N * v) / sum(N) * 100

      )

    }

  )

)

cat("\n====================================================\n")
cat("VACCINATION COVERAGE BY SCENARIO\n")
cat("====================================================\n")

print(

  scenario_coverage %>%

    mutate(

      across(

        where(is.numeric),

        ~ round(.x, 2)

      )

    )

)

cat("\n====================================================\n")
cat("RUNNING SEIR SIMULATIONS\n")
cat("====================================================\n")

results_70 <- run_scenario(

  "Uniform 70%",

  vax_uniform_70

)

cat("Completed: Uniform 70%\n")

results_80 <- run_scenario(

  "Uniform 80%",

  vax_uniform_80

)

cat("Completed: Uniform 80%\n")

results_size_A <- run_scenario(

  "Small 75%, Large 90%",

  vax_size_A

)

cat("Completed: Small 75%, Large 90%\n")

results_size_B <- run_scenario(

  scenario_names[4],

  vax_size_B

)

cat("Completed: Adjusted small-school coverage, large schools 75%\n")

results_90 <- run_scenario(

  "Uniform 90%",

  vax_uniform_90

)

cat("Completed: Uniform 90%\n")

results_summary <- bind_rows(

  results_70$summary,

  results_80$summary,

  results_size_A$summary,

  results_size_B$summary,

  results_90$summary

)

cat("\n====================================================\n")
cat("FINAL EPIDEMIC OUTCOMES\n")
cat("====================================================\n")

print(

  results_summary %>%

    mutate(

      across(

        where(is.numeric),

        ~ round(.x, 2)

      )

    )

)

coverage_A <- results_size_A$summary$Aggregate_Coverage

coverage_B <- results_size_B$summary$Aggregate_Coverage

coverage_difference <- abs(
  coverage_A - coverage_B
)

cat("\n====================================================\n")
cat("AGGREGATE VACCINATION COVERAGE CHECK\n")
cat("====================================================\n")

cat(

  "Scenario 3 aggregate coverage:",

  round(
    coverage_A,
    4
  ),

  "%\n"

)

cat(

  "Scenario 4 aggregate coverage:",

  round(
    coverage_B,
    4
  ),

  "%\n"

)

cat(

  "Absolute difference:",

  format(
    coverage_difference,
    scientific = FALSE,
    digits = 10
  ),

  "percentage points\n"

)

if (coverage_difference > 1e-6) {

  warning(

    "The aggregate vaccination coverage of the two size-based scenarios ",
    "does not match exactly."

  )

} else {

  cat(
    "The aggregate vaccination coverage matches.\n"
  )

}

school_comparison <- results_size_A$school_cases %>%

  select(

    unit_id,

    school_name,

    school_size,

    school_type,

    network_strength,

    initially_infectious,

    vaccination_coverage,

    cases_A = cumulative_new_cases

  ) %>%

  left_join(

    results_size_B$school_cases %>%

      select(

        unit_id,

        cases_B = cumulative_new_cases

      ),

    by = "unit_id"

  ) %>%

  mutate(

    Difference_B_minus_A = cases_B - cases_A

  ) %>%

  arrange(
    desc(network_strength)
  )

cat("\n====================================================\n")
cat("SCHOOL-LEVEL CASE COMPARISON\n")
cat("====================================================\n")

print(

  school_comparison %>%

    mutate(

      across(

        where(is.numeric),

        ~ round(.x, 2)

      )

    ),

  n = Inf

)

cases_by_size <- bind_rows(

  results_size_A$school_cases %>%

    group_by(
      school_type
    ) %>%

    summarise(

      Scenario = "Small 75%, Large 90%",

      Schools = n(),

      Population = sum(school_size),

      Total_New_Cases = sum(cumulative_new_cases),

      .groups = "drop"

    ),

  results_size_B$school_cases %>%

    group_by(
      school_type
    ) %>%

    summarise(

      Scenario = scenario_names[4],

      Schools = n(),

      Population = sum(school_size),

      Total_New_Cases = sum(cumulative_new_cases),

      .groups = "drop"

    )

)

cat("\n====================================================\n")
cat("CASES BY SCHOOL SIZE\n")
cat("====================================================\n")

print(

  cases_by_size %>%

    mutate(

      across(

        where(is.numeric),

        ~ round(.x, 2)

      )

    )

)

first_school_results <- bind_rows(

  results_70$school_cases %>%

    slice(1) %>%

    mutate(
      Scenario = "Uniform 70%"
    ),

  results_80$school_cases %>%

    slice(1) %>%

    mutate(
      Scenario = "Uniform 80%"
    ),

  results_size_A$school_cases %>%

    slice(1) %>%

    mutate(
      Scenario = "Small 75%, Large 90%"
    ),

  results_size_B$school_cases %>%

    slice(1) %>%

    mutate(
      Scenario = scenario_names[4]
    ),

  results_90$school_cases %>%

    slice(1) %>%

    mutate(
      Scenario = "Uniform 90%"
    )

) %>%

  select(

    Scenario,

    unit_id,

    school_name,

    school_size,

    vaccination_coverage,

    initially_infectious,

    cumulative_new_cases

  )

cat("\n====================================================\n")
cat("FIRST SCHOOL: SEED SCHOOL CHECK\n")
cat("====================================================\n")

print(

  first_school_results %>%

    mutate(

      across(

        where(is.numeric),

        ~ round(.x, 2)

      )

    )

)

cat("\n====================================================\n")
cat("CUMULATIVE CASE CALCULATION CHECK\n")
cat("====================================================\n")

cat(

  "The initial infectious student is excluded from new cases.\n"

)

cat(

  "Cumulative new cases are calculated from the decrease ",
  "in the susceptible population.\n",
  sep = ""

)

cat(

  "Total population used for the county attack rate:",

  sum(N),

  "\n"

)

cat(

  "The county attack rate is cumulative new cases divided by ",
  "the total school population, multiplied by 100.\n",
  sep = ""

)

write_csv(

  results_summary,

  "output/tables/vaccine_allocation_epidemic_results.csv"

)

write_csv(

  scenario_coverage,

  "output/tables/vaccine_allocation_coverage_summary.csv"

)

write_csv(

  school_comparison,

  "output/tables/school_level_case_comparison.csv"

)

write_csv(

  cases_by_size,

  "output/tables/cases_by_school_size.csv"

)

write_csv(

  first_school_results,

  "output/tables/first_school_seed_check.csv"

)

cat("\n====================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

cat(

  "Number of schools included:",

  n_school,

  "\n"

)

cat(

  "Infection started in:",

  school_names[seed_school_index],

  "\n"

)

cat(

  "Seed school ID:",

  school_ids[seed_school_index],

  "\n"

)

cat(

  "Seed school size:",

  N[seed_school_index],

  "\n"

)

cat(

  "Seed school network strength:",

  school_strength[seed_school_index],

  "\n"

)

cat("\nFiles saved:\n")

cat(
  "1. vaccine_allocation_epidemic_results.csv\n"
)

cat(
  "2. vaccine_allocation_coverage_summary.csv\n"
)

cat(
  "3. school_level_case_comparison.csv\n"
)

cat(
  "4. cases_by_school_size.csv\n"
)

cat(
  "5. first_school_seed_check.csv\n"
)

cat("\nAll five vaccination scenarios have been simulated.\n")

cat(
  "The initial infectious student was not counted as a new case.\n"
)

cat(
  "The infection started in the first school, not the network hub.\n"
)

cat("====================================================\n")

epidemic_plot <- ggplot(
  epidemic_curves,
  aes(
    x = Day,
    y = Infectious_Students,
    color = Scenario,
    linetype = Scenario
  )
) +
  geom_line(linewidth = 0.8) +
  scale_linetype_manual(
    values = c(
      "solid",
      "dashed",
      "dotted",
      "dotdash",
      "longdash"
    )
  ) +
  labs(
    title = "Measles Epidemic Curves Under Different Vaccination Scenarios",
    x = "Time (days)",
    y = "Number of Infectious Students",
    color = "Vaccination Scenario",
    linetype = "Vaccination Scenario"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "bottom"
  )

print(epidemic_plot)

ggsave(
  "figures/simulation/measles_epidemic_curves.png",
  plot = epidemic_plot,
  width = 10,
  height = 7,
  dpi = 300
)

county_attack_results <- results_summary %>%
  select(
    Scenario,
    Total_New_Cases,
    Attack_Rate_Percent
  ) %>%
  mutate(
    Total_New_Cases = round(Total_New_Cases, 2),
    Attack_Rate_Percent = round(Attack_Rate_Percent, 2)
  )

print(county_attack_results)

hub_index <- which.max(school_strength)
peripheral_index <- which.min(school_strength)

hub_peripheral_results <- tibble(
  Position = c("Hub", "Peripheral"),
  School_ID = c(school_ids[hub_index], school_ids[peripheral_index]),
  School_Name = c(school_names[hub_index], school_names[peripheral_index]),
  School_Size = c(N[hub_index], N[peripheral_index]),
  Vaccination_Coverage = c(
    observed_vax[hub_index] * 100,
    observed_vax[peripheral_index] * 100
  ),
  Network_Strength = c(
    school_strength[hub_index],
    school_strength[peripheral_index]
  )
)

print(hub_peripheral_results)

hub_index <- which.max(school_strength)

peripheral_index <- which.min(school_strength)

hub_peripheral_results <- tibble(
  Network_Position = c("Hub", "Peripheral"),
  School_ID = c(
    school_ids[hub_index],
    school_ids[peripheral_index]
  ),
  School_Name = c(
    school_names[hub_index],
    school_names[peripheral_index]
  ),
  School_Size = c(
    N[hub_index],
    N[peripheral_index]
  ),
  Vaccination_Coverage_Percent = c(
    observed_vax[hub_index] * 100,
    observed_vax[peripheral_index] * 100
  ),
  Network_Strength = c(
    school_strength[hub_index],
    school_strength[peripheral_index]
  )
)

cat("\n====================================================\n")
cat("HUB AND PERIPHERAL SCHOOLS\n")
cat("====================================================\n")

print(
  hub_peripheral_results %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 2)
      )
    )
)

write_csv(
  hub_peripheral_results,
  "output/tables/hub_and_peripheral_schools.csv"
)
