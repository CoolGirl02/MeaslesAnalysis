# ==============================================================================
# NETWORK SEIR MEASLES SIMULATION
# AND SENSITIVITY ANALYSIS 5:
# VACCINE ALLOCATION BY SCHOOL SIZE VS. BASELINES
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(igraph)
})

theme_set(
  theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12))
)


# ==============================================================================
# 1. LOAD SCHOOL DATA
# ==============================================================================

df <- read.csv("data/cleanned/school_data_joined.csv") %>%
  rename_with(tolower)

if (!"unit_id" %in% names(df)) {
  df <- df %>%
    mutate(unit_id = as.character(row_number()))
} else {
  df <- df %>%
    mutate(unit_id = as.character(unit_id))
}

df <- df %>%
  mutate(
    vax_rate = if_else(
      percent_immunized > 1,
      percent_immunized / 100,
      percent_immunized
    ),
    school_id = row_number(),
    page = ceiling(school_id / 6)
  )

n_schools <- nrow(df)


# ==============================================================================
# 2. BUILD NETWORK ADJACENCY MATRIX
# ==============================================================================

net_df <- read.csv("data/cleanned/network_df.csv") %>%
  rename_with(tolower)

id_cols <- names(net_df)[
  grep("unit|from|to|source|target|id", names(net_df))
]

if (length(id_cols) >= 2) {
  
  net_df <- net_df %>%
    rename(
      unit_id_a = !!sym(id_cols[1]),
      unit_id_b = !!sym(id_cols[2])
    )
  
} else {
  
  colnames(net_df)[1:2] <- c(
    "unit_id_a",
    "unit_id_b"
  )
}

net_df <- net_df %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  )


# ==============================================================================
# 3. CREATE ADJACENCY MATRIX
# ==============================================================================

id_map <- setNames(
  1:n_schools,
  df$unit_id
)

A <- matrix(
  0,
  nrow = n_schools,
  ncol = n_schools
)

for (k in seq_len(nrow(net_df))) {
  
  i <- id_map[net_df$unit_id_a[k]]
  j <- id_map[net_df$unit_id_b[k]]
  
  if (!is.na(i) && !is.na(j)) {
    
    w <- if ("total_contacts" %in% names(net_df)) {
      net_df$total_contacts[k]
    } else {
      1
    }
    
    A[i, j] <- w
    A[j, i] <- w
  }
}


# ==============================================================================
# 4. SCALE NETWORK
# ==============================================================================

A_raw <- A

spectral_radius <- max(
  Re(eigen(A_raw)$values)
)

A <- A_raw / spectral_radius


# ==============================================================================
# 5. IDENTIFY HIGH-CENTRALITY SCHOOL
# ==============================================================================

school_strength <- rowSums(A_raw)

top_hub_id <- which.max(
  school_strength
)

cat("\nHigh-centrality seed school:", top_hub_id, "\n")
cat(
  "School name:",
  df$school_name[top_hub_id],
  "\n\n"
)


# ==============================================================================
# 6. SEIR PARAMETERS
# ==============================================================================

gamma_rate <- 1 / 7
sigma_rate <- 1 / 10
beta_rate  <- 15 * gamma_rate
alpha      <- 0.05

times <- seq(
  0,
  200,
  by = 1
)


# ==============================================================================
# 7. INITIAL CONDITIONS
# ==============================================================================

N <- df$total_students

R0 <- round(
  N * df$vax_rate
)

I0 <- c(
  1,
  rep(0, n_schools - 1)
)

E0 <- rep(
  0,
  n_schools
)

S0 <- pmax(
  0,
  N - I0 - E0 - R0
)

initial_state <- c(
  S0,
  E0,
  I0,
  R0
)

names(initial_state) <- c(
  paste0("S_", 1:n_schools),
  paste0("E_", 1:n_schools),
  paste0("I_", 1:n_schools),
  paste0("R_", 1:n_schools)
)


# ==============================================================================
# 8. MODEL PARAMETERS
# ==============================================================================

params <- list(
  beta  = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N     = N,
  A     = A,
  n     = n_schools
)


# ==============================================================================
# 9. SEIR MODEL
# ==============================================================================

seir_model <- function(time, state, parameters) {
  
  with(as.list(parameters), {
    
    S <- state[
      1:n
    ]
    
    E <- state[
      (n + 1):(2 * n)
    ]
    
    I <- state[
      (2 * n + 1):(3 * n)
    ]
    
    R <- state[
      (3 * n + 1):(4 * n)
    ]
    
    I_prop <- I / N
    
    lambda <- (
      beta * I_prop
    ) +
      (
        alpha *
          beta *
          as.vector(A %*% I_prop)
      )
    
    dS <- -lambda * S
    
    dE <- lambda * S -
      sigma_rate * E
    
    dI <- sigma_rate * E -
      gamma_rate * I
    
    dR <- gamma_rate * I
    
    list(
      c(
        dS,
        dE,
        dI,
        dR
      )
    )
  })
}


# ==============================================================================
# 10. SENSITIVITY ANALYSIS 5:
#     CLASSIFY SCHOOLS BY SIZE
# ==============================================================================

small_idx <- which(
  df$total_students < 700
)

large_idx <- which(
  df$total_students >= 700
)

N_small <- sum(
  df$total_students[small_idx]
)

N_large <- sum(
  df$total_students[large_idx]
)

N_total <- sum(
  df$total_students
)


# ==============================================================================
# 11. SCENARIO A:
#     LOW COVERAGE IN SMALL SCHOOLS (75%)
#     HIGH COVERAGE IN LARGE SCHOOLS (90%)
# ==============================================================================

vax_small_low <- numeric(
  n_schools
)

vax_small_low[small_idx] <- 0.75

vax_small_low[large_idx] <- 0.90

total_vax_A <- sum(
  df$total_students *
    vax_small_low
)

rate_A <- total_vax_A /
  N_total


# ==============================================================================
# 12. SCENARIO B:
#     LOW COVERAGE IN LARGE SCHOOLS (75%)
#     CALCULATE SMALL-SCHOOL COVERAGE
#     TO MATCH SCENARIO A
# ==============================================================================

vax_small_B_needed <- (
  total_vax_A -
    (0.75 * N_large)
) / N_small

vax_large_low <- numeric(
  n_schools
)

vax_large_low[large_idx] <- 0.75

vax_large_low[small_idx] <-
  vax_small_B_needed

total_vax_B <- sum(
  df$total_students *
    vax_large_low
)

rate_B <- total_vax_B /
  N_total


# ==============================================================================
# 13. UNIFORM BASELINES
# ==============================================================================

vax_70 <- rep(
  0.70,
  n_schools
)

vax_80 <- rep(
  0.80,
  n_schools
)

vax_90 <- rep(
  0.90,
  n_schools
)

total_vax_70 <- sum(
  df$total_students *
    vax_70
)

total_vax_80 <- sum(
  df$total_students *
    vax_80
)

total_vax_90 <- sum(
  df$total_students *
    vax_90
)

rate_70 <- total_vax_70 /
  N_total

rate_80 <- total_vax_80 /
  N_total

rate_90 <- total_vax_90 /
  N_total


# ==============================================================================
# 14. CHECK AGGREGATE COVERAGE
# ==============================================================================

cat("\n============================================\n")
cat("AGGREGATE VACCINATION COVERAGE\n")
cat("============================================\n")

cat(
  "Uniform 70%:",
  round(rate_70 * 100, 4),
  "%\n"
)

cat(
  "Uniform 80%:",
  round(rate_80 * 100, 4),
  "%\n"
)

cat(
  "Low Small / High Large:",
  round(rate_A * 100, 4),
  "%\n"
)

cat(
  "High Small / Low Large:",
  round(rate_B * 100, 4),
  "%\n"
)

cat(
  "Uniform 90%:",
  round(rate_90 * 100, 4),
  "%\n"
)

cat(
  "Difference between the two allocated scenarios:",
  round(abs(rate_A - rate_B) * 100, 8),
  "percentage points\n"
)

cat("============================================\n\n")


# ==============================================================================
# 15. FUNCTION TO RUN EACH SCENARIO
#     AND CALCULATE EPIDEMIC OUTCOMES
# ==============================================================================

run_vax_scenario <- function(
    vax_vector,
    label,
    seed_id = 1) {
  
  # Vaccinated population
  R0_vec <- round(
    N * vax_vector
  )
  
  # Seed infection
  I0_vec <- rep(
    0,
    n_schools
  )
  
  I0_vec[seed_id] <- 1
  
  # Susceptible population
  S0_vec <- pmax(
    0,
    N - I0_vec - R0_vec
  )
  
  # Initial state
  state_vec <- c(
    S0_vec,
    rep(0, n_schools),
    I0_vec,
    R0_vec
  )
  
  names(state_vec) <-
    names(initial_state)
  
  # Run SEIR model
  out_df <- as.data.frame(
    ode(
      y = state_vec,
      times = times,
      func = seir_model,
      parms = params
    )
  )
  
  
  # --------------------------------------------------------------------------
  # TOTAL INFECTIOUS CASES
  # --------------------------------------------------------------------------
  
  total_I <- rowSums(
    out_df[
      ,
      paste0(
        "I_",
        1:n_schools
      )
    ]
  )
  
  
  # --------------------------------------------------------------------------
  # PEAK BURDEN
  # --------------------------------------------------------------------------
  
  peak_burden <- max(
    total_I
  )
  
  
  # --------------------------------------------------------------------------
  # PEAK DAY
  # --------------------------------------------------------------------------
  
  peak_day <- out_df$time[
    which.max(total_I)
  ]
  
  
  # --------------------------------------------------------------------------
  # TOTAL NEW CASES
  # --------------------------------------------------------------------------
  
  initial_S <- sum(
    S0_vec
  )
  
  final_S <- sum(
    out_df[
      nrow(out_df),
      paste0(
        "S_",
        1:n_schools
      )
    ]
  )
  
  total_new_cases <-
    initial_S - final_S
  
  
  # --------------------------------------------------------------------------
  # COUNTY ATTACK RATE
  # --------------------------------------------------------------------------
  
  county_attack_rate <-
    total_new_cases /
    N_total
  
  
  # --------------------------------------------------------------------------
  # RETURN RESULTS
  # --------------------------------------------------------------------------
  
  data.frame(
    Scenario = label,
    time = out_df$time,
    total_infectious = total_I,
    Peak_Burden = peak_burden,
    Peak_Day = peak_day,
    Total_New_Cases = total_new_cases,
    County_Attack_Rate = county_attack_rate
  )
}


# ==============================================================================
# 16. RUN THE FIVE SENSITIVITY-ANALYSIS SCENARIOS
# ==============================================================================

sens5_df <- bind_rows(
  
  run_vax_scenario(
    vax_70,
    "Baseline: Uniform 70%",
    seed_id = top_hub_id
  ),
  
  run_vax_scenario(
    vax_80,
    "Baseline: Uniform 80%",
    seed_id = top_hub_id
  ),
  
  run_vax_scenario(
    vax_small_low,
    "Allocated: Low in Small / High in Large",
    seed_id = top_hub_id
  ),
  
  run_vax_scenario(
    vax_large_low,
    "Allocated: High in Small / Low in Large",
    seed_id = top_hub_id
  ),
  
  run_vax_scenario(
    vax_90,
    "Baseline: Uniform 90%",
    seed_id = top_hub_id
  )
)


# ==============================================================================
# 17. CREATE TABLE 4.2
# ==============================================================================

table4_2 <- sens5_df %>%
  
  group_by(
    Scenario
  ) %>%
  
  summarise(
    
    Peak_Burden =
      first(Peak_Burden),
    
    Peak_Day =
      first(Peak_Day),
    
    Total_New_Cases =
      first(Total_New_Cases),
    
    County_Attack_Rate =
      first(County_Attack_Rate),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    # --------------------------------------------------------------------------
    # SMALL SCHOOL COVERAGE
    # --------------------------------------------------------------------------
    
    Small_School_Coverage = case_when(
      
      Scenario ==
        "Baseline: Uniform 70%" ~ 0.70,
      
      Scenario ==
        "Baseline: Uniform 80%" ~ 0.80,
      
      Scenario ==
        "Allocated: Low in Small / High in Large" ~
        0.75,
      
      Scenario ==
        "Allocated: High in Small / Low in Large" ~
        vax_small_B_needed,
      
      Scenario ==
        "Baseline: Uniform 90%" ~ 0.90
    ),
    
    
    # --------------------------------------------------------------------------
    # LARGE SCHOOL COVERAGE
    # --------------------------------------------------------------------------
    
    Large_School_Coverage = case_when(
      
      Scenario ==
        "Baseline: Uniform 70%" ~ 0.70,
      
      Scenario ==
        "Baseline: Uniform 80%" ~ 0.80,
      
      Scenario ==
        "Allocated: Low in Small / High in Large" ~
        0.90,
      
      Scenario ==
        "Allocated: High in Small / Low in Large" ~
        0.75,
      
      Scenario ==
        "Baseline: Uniform 90%" ~ 0.90
    ),
    
    
    # --------------------------------------------------------------------------
    # AGGREGATE COVERAGE
    # --------------------------------------------------------------------------
    
    Aggregate_Coverage = case_when(
      
      Scenario ==
        "Baseline: Uniform 70%" ~
        rate_70,
      
      Scenario ==
        "Baseline: Uniform 80%" ~
        rate_80,
      
      Scenario ==
        "Allocated: Low in Small / High in Large" ~
        rate_A,
      
      Scenario ==
        "Allocated: High in Small / Low in Large" ~
        rate_B,
      
      Scenario ==
        "Baseline: Uniform 90%" ~
        rate_90
    )
  ) %>%
  
  select(
    Scenario,
    Small_School_Coverage,
    Large_School_Coverage,
    Aggregate_Coverage,
    Peak_Burden,
    Peak_Day,
    Total_New_Cases,
    County_Attack_Rate
  )


# ==============================================================================
# 18. ORDER TABLE 4.2
# ==============================================================================

scenario_order <- c(
  "Baseline: Uniform 70%",
  "Baseline: Uniform 80%",
  "Allocated: High in Small / Low in Large",
  "Allocated: Low in Small / High in Large",
  "Baseline: Uniform 90%"
)

table4_2 <- table4_2 %>%
  
  mutate(
    Scenario = factor(
      Scenario,
      levels = scenario_order
    )
  ) %>%
  
  arrange(
    Scenario
  ) %>%
  
  mutate(
    Scenario = as.character(Scenario)
  )


# ==============================================================================
# 19. DISPLAY TABLE 4.2
# ==============================================================================

table4_2_display <- table4_2 %>%
  
  mutate(
    
    Small_School_Coverage =
      paste0(
        round(
          Small_School_Coverage * 100,
          2
        ),
        "%"
      ),
    
    Large_School_Coverage =
      paste0(
        round(
          Large_School_Coverage * 100,
          2
        ),
        "%"
      ),
    
    Aggregate_Coverage =
      paste0(
        round(
          Aggregate_Coverage * 100,
          2
        ),
        "%"
      ),
    
    Peak_Burden =
      round(
        Peak_Burden,
        2
      ),
    
    Peak_Day =
      round(
        Peak_Day,
        0
      ),
    
    Total_New_Cases =
      round(
        Total_New_Cases,
        2
      ),
    
    County_Attack_Rate =
      paste0(
        round(
          County_Attack_Rate * 100,
          2
        ),
        "%"
      )
  )


cat("\n")
cat("=====================================================================\n")
cat("TABLE 4.2: EPIDEMIC OUTCOMES UNDER SIZE-BASED ALLOCATION\n")
cat("              AND UNIFORM BASELINE STRATEGIES\n")
cat("=====================================================================\n\n")

print(
  table4_2_display,
  row.names = FALSE
)

cat("\n=====================================================================\n")


# ==============================================================================
# 20. COMPARATIVE EPIDEMIC CURVES
# ==============================================================================

ggplot(
  sens5_df,
  aes(
    x = time,
    y = total_infectious,
    color = Scenario,
    linetype = Scenario
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = c(
      
      "Baseline: Uniform 70%" =
        "firebrick",
      
      "Baseline: Uniform 80%" =
        "darkorange",
      
      "Baseline: Uniform 90%" =
        "forestgreen",
      
      "Allocated: Low in Small / High in Large" =
        "purple",
      
      "Allocated: High in Small / Low in Large" =
        "deepskyblue3"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      
      "Baseline: Uniform 70%" =
        "dotted",
      
      "Baseline: Uniform 80%" =
        "dotted",
      
      "Baseline: Uniform 90%" =
        "dotted",
      
      "Allocated: Low in Small / High in Large" =
        "solid",
      
      "Allocated: High in Small / Low in Large" =
        "dashed"
    )
  ) +
  
  labs(
    title =
      "Sensitivity Analysis 1: Size-Based Vaccine Allocation vs. Uniform Baselines",
    
    subtitle =
      "Comparing heterogeneous school coverage against uniform vaccination strategies",
    
    x =
      "Days Since Outbreak Introduction",
    
    y =
      "Total Active Infectious Cases in County"
  ) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  
  guides(
    color =
      guide_legend(
        nrow = 2,
        byrow = TRUE
      )
  )









# ============================================================
# SENSITIVITY ANALYSIS 5:
# VACCINE ALLOCATION BY SCHOOL SIZE VS UNIFORM BASELINES
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(igraph)
})


# ============================================================
# 1. LOAD SCHOOL DATA
# ============================================================

df <- read_csv(
  "data/cleanned/school_data_joined.csv",
  show_col_types = FALSE
) %>%
  rename_with(tolower)


# ============================================================
# 2. CHECK SCHOOL DATA COLUMNS
# ============================================================

required_school_cols <- c(
  "unit_id",
  "school_name",
  "total_students",
  "percent_immunized"
)

missing_school_cols <- setdiff(
  required_school_cols,
  names(df)
)

if (length(missing_school_cols) > 0) {
  
  stop(
    paste(
      "The following school-data columns are missing:",
      paste(
        missing_school_cols,
        collapse = ", "
      ),
      "\nAvailable columns are:",
      paste(
        names(df),
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 3. CREATE VACCINATION RATE
# ============================================================

df <- df %>%
  mutate(
    vax_rate = percent_immunized / 100
  )


# ============================================================
# 4. BASIC SCHOOL INFORMATION
# ============================================================

N <- df$total_students

n_schools <- nrow(df)

N_total <- sum(
  N,
  na.rm = TRUE
)

cat(
  "Number of schools:",
  n_schools,
  "\n"
)

cat(
  "Total students:",
  N_total,
  "\n"
)


# ============================================================
# 5. LOAD NETWORK DATA
# ============================================================

network_df <- read_csv(
  "data/cleanned/network_df.csv",
  show_col_types = FALSE
) %>%
  rename_with(tolower)


# ============================================================
# 6. CHECK NETWORK DATA COLUMNS
# ============================================================

required_network_cols <- c(
  "unit_id_a",
  "unit_id_b",
  "total_contacts",
  "connect_risk"
)

missing_network_cols <- setdiff(
  required_network_cols,
  names(network_df)
)

if (length(missing_network_cols) > 0) {
  
  stop(
    paste(
      "The following network-data columns are missing:",
      paste(
        missing_network_cols,
        collapse = ", "
      ),
      "\nAvailable columns are:",
      paste(
        names(network_df),
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 7. CHECK UNIT ID MATCHING
# ============================================================

network_units <- unique(
  c(
    network_df$unit_id_a,
    network_df$unit_id_b
  )
)

network_units <- network_units[
  !is.na(network_units)
]

school_units <- unique(
  df$unit_id
)

unmatched_network_units <- setdiff(
  network_units,
  school_units
)

cat(
  "Network school IDs:",
  length(network_units),
  "\n"
)

cat(
  "School-data IDs:",
  length(school_units),
  "\n"
)

if (
  length(unmatched_network_units) > 0
) {
  
  warning(
    paste(
      length(unmatched_network_units),
      "network unit IDs were not found in school_data_joined.csv."
    )
  )
}


# ============================================================
# 8. CREATE RAW ADJACENCY MATRIX
# ============================================================

A_raw <- matrix(
  0,
  nrow = n_schools,
  ncol = n_schools
)

rownames(A_raw) <- df$school_name

colnames(A_raw) <- df$school_name


# ============================================================
# 9. BUILD WEIGHTED NETWORK
# ============================================================

for (
  i in seq_len(nrow(network_df))
) {
  
  unit_A <- network_df$unit_id_a[i]
  
  unit_B <- network_df$unit_id_b[i]
  
  risk_weight <- network_df$connect_risk[i]
  
  
  # ----------------------------------------------------------
  # SKIP MISSING VALUES
  # ----------------------------------------------------------
  
  if (
    is.na(unit_A) ||
    is.na(unit_B) ||
    is.na(risk_weight)
  ) {
    
    next
  }
  
  
  # ----------------------------------------------------------
  # MATCH NETWORK UNIT IDs TO SCHOOL DATA
  # ----------------------------------------------------------
  
  a <- match(
    unit_A,
    df$unit_id
  )
  
  b <- match(
    unit_B,
    df$unit_id
  )
  
  
  # ----------------------------------------------------------
  # SKIP UNMATCHED OR SELF CONNECTIONS
  # ----------------------------------------------------------
  
  if (
    is.na(a) ||
    is.na(b) ||
    a == b
  ) {
    
    next
  }
  
  
  # ----------------------------------------------------------
  # CREATE UNDIRECTED WEIGHTED EDGE
  # ----------------------------------------------------------
  
  A_raw[a, b] <- risk_weight
  
  A_raw[b, a] <- risk_weight
}


# Remove self-connections
diag(A_raw) <- 0


# ============================================================
# 10. CHECK NETWORK
# ============================================================

number_of_edges <- sum(
  A_raw > 0
) / 2

cat(
  "Number of network edges:",
  number_of_edges,
  "\n"
)


# ============================================================
# 11. SCALE NETWORK BY SPECTRAL RADIUS
# ============================================================

eig_values <- eigen(
  A_raw,
  only.values = TRUE
)$values

spectral_radius <- max(
  Mod(eig_values),
  na.rm = TRUE
)

cat(
  "Network spectral radius:",
  spectral_radius,
  "\n"
)


if (
  is.finite(spectral_radius) &&
  spectral_radius > 0
) {
  
  A <- A_raw / spectral_radius
  
} else {
  
  A <- A_raw
  
  warning(
    "Spectral radius is zero or invalid. Raw network was used."
  )
}


# ============================================================
# 12. CALCULATE SCHOOL CONNECTIVITY
# ============================================================

school_connectivity <- rowSums(
  A_raw,
  na.rm = TRUE
)


# Highest-connectivity school
top_hub_id <- which.max(
  school_connectivity
)


# Lowest-connectivity school
low_hub_id <- which.min(
  school_connectivity
)


top_hub_name <- df$school_name[
  top_hub_id
]

low_hub_name <- df$school_name[
  low_hub_id
]


cat(
  "\nHighest-connectivity school:",
  top_hub_name,
  "\n"
)

cat(
  "Lowest-connectivity school:",
  low_hub_name,
  "\n"
)


# ============================================================
# 13. SEIR PARAMETERS
# ============================================================

gamma_rate <- 1 / 7

sigma_rate <- 1 / 10

beta_rate <- 15 * gamma_rate

alpha <- 0.05

times <- seq(
  0,
  200,
  by = 1
)


# ============================================================
# 14. SEIR NETWORK MODEL
# ============================================================

seir_network <- function(
    time,
    state,
    parameters
) {
  
  beta <- parameters["beta"]
  
  sigma <- parameters["sigma"]
  
  gamma <- parameters["gamma"]
  
  alpha <- parameters["alpha"]
  
  
  # ----------------------------------------------------------
  # EXTRACT COMPARTMENTS
  # ----------------------------------------------------------
  
  S <- state[
    paste0(
      "S_",
      1:n_schools
    )
  ]
  
  E <- state[
    paste0(
      "E_",
      1:n_schools
    )
  ]
  
  I <- state[
    paste0(
      "I_",
      1:n_schools
    )
  ]
  
  R <- state[
    paste0(
      "R_",
      1:n_schools
    )
  ]
  
  
  # ----------------------------------------------------------
  # INFECTIOUS PROPORTION
  # ----------------------------------------------------------
  
  I_prop <- I / N
  
  
  # ----------------------------------------------------------
  # FORCE OF INFECTION
  # ----------------------------------------------------------
  
  lambda <- (
    beta * I_prop
  ) +
    (
      alpha *
        beta *
        as.vector(
          A %*% I_prop
        )
    )
  
  
  # ----------------------------------------------------------
  # SEIR DIFFERENTIAL EQUATIONS
  # ----------------------------------------------------------
  
  dS <- -lambda * S
  
  dE <- (
    lambda * S
  ) -
    (
      sigma * E
    )
  
  dI <- (
    sigma * E
  ) -
    (
      gamma * I
    )
  
  dR <- gamma * I
  
  
  # ----------------------------------------------------------
  # RETURN DERIVATIVES
  # ----------------------------------------------------------
  
  list(
    c(
      dS,
      dE,
      dI,
      dR
    )
  )
}


# ============================================================
# 15. MODEL PARAMETERS
# ============================================================

parameters <- c(
  beta = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha
)


# ============================================================
# 16. DEFINE SCHOOL SIZE GROUPS
# ============================================================

small_schools <- which(
  N <= 1000
)

large_schools <- which(
  N > 1000
)


N_small <- sum(
  N[small_schools],
  na.rm = TRUE
)

N_large <- sum(
  N[large_schools],
  na.rm = TRUE
)


cat(
  "\nNumber of small schools:",
  length(small_schools),
  "\n"
)

cat(
  "Number of large schools:",
  length(large_schools),
  "\n"
)

cat(
  "Students in small schools:",
  N_small,
  "\n"
)

cat(
  "Students in large schools:",
  N_large,
  "\n"
)


# ============================================================
# 17. SCENARIO A
# HIGHER COVERAGE IN LARGE SCHOOLS
#
# Small schools = 75%
# Large schools = 90%
# ============================================================

vax_A <- rep(
  NA_real_,
  n_schools
)

vax_A[small_schools] <- 0.75

vax_A[large_schools] <- 0.90


# Total vaccinated people
total_vax_A <- sum(
  N * vax_A,
  na.rm = TRUE
)


# Overall vaccination coverage
rate_A <- total_vax_A / N_total


cat(
  "\nScenario A overall coverage:",
  round(
    rate_A * 100,
    2
  ),
  "%\n"
)


# ============================================================
# 18. SCENARIO B
# HIGHER COVERAGE IN SMALL SCHOOLS
#
# Large schools = 75%
#
# Small-school coverage is calculated so that the total
# number vaccinated is the same as Scenario A.
# ============================================================

vax_large_B <- 0.75


vax_small_B_needed <- (
  total_vax_A -
    (
      vax_large_B *
        N_large
    )
) /
  N_small


# Check feasibility
if (
  vax_small_B_needed > 1
) {
  
  warning(
    paste0(
      "Scenario B requires ",
      round(
        vax_small_B_needed * 100,
        2
      ),
      "% vaccination in small schools, ",
      "which is greater than 100%. ",
      "This scenario is NOT epidemiologically feasible."
    )
  )
}


if (
  vax_small_B_needed < 0
) {
  
  warning(
    "Scenario B produces a negative vaccination coverage."
  )
}


# Create Scenario B
vax_B <- rep(
  NA_real_,
  n_schools
)

vax_B[small_schools] <-
  vax_small_B_needed

vax_B[large_schools] <-
  vax_large_B


# Overall coverage
total_vax_B <- sum(
  N * vax_B,
  na.rm = TRUE
)

rate_B <- total_vax_B / N_total


cat(
  "Scenario B small-school coverage:",
  round(
    vax_small_B_needed * 100,
    2
  ),
  "%\n"
)

cat(
  "Scenario B overall coverage:",
  round(
    rate_B * 100,
    2
  ),
  "%\n"
)


# ============================================================
# 19. UNIFORM BASELINE SCENARIOS
# ============================================================

vax_70 <- rep(
  0.70,
  n_schools
)

vax_80 <- rep(
  0.80,
  n_schools
)

vax_90 <- rep(
  0.90,
  n_schools
)


# ============================================================
# 20. FUNCTION TO RUN A VACCINATION SCENARIO
# ============================================================

run_vax_scenario <- function(
    vax_vector,
    seed_id,
    scenario_name
) {
  
  
  # ==========================================================
  # INITIAL VACCINATED / IMMUNE POPULATION
  # ==========================================================
  
  R0_vec <- round(
    N * vax_vector
  )
  
  
  # ==========================================================
  # INITIAL INFECTED
  # ==========================================================
  
  I0_vec <- rep(
    0,
    n_schools
  )
  
  I0_vec[seed_id] <- 1
  
  
  # ==========================================================
  # INITIAL EXPOSED
  # ==========================================================
  
  E0_vec <- rep(
    0,
    n_schools
  )
  
  
  # ==========================================================
  # INITIAL SUSCEPTIBLE
  # ==========================================================
  
  S0_vec <- pmax(
    0,
    N -
      I0_vec -
      E0_vec -
      R0_vec
  )
  
  
  # ==========================================================
  # INITIAL STATE VECTOR
  # ==========================================================
  
  state_vec <- c(
    S0_vec,
    E0_vec,
    I0_vec,
    R0_vec
  )
  
  
  names(state_vec) <- c(
    
    paste0(
      "S_",
      1:n_schools
    ),
    
    paste0(
      "E_",
      1:n_schools
    ),
    
    paste0(
      "I_",
      1:n_schools
    ),
    
    paste0(
      "R_",
      1:n_schools
    )
  )
  
  
  # ==========================================================
  # RUN SEIR MODEL
  # ==========================================================
  
  out <- ode(
    y = state_vec,
    times = times,
    func = seir_network,
    parms = parameters
  )
  
  
  # Convert output to data frame
  out_df <- as.data.frame(out)
  
  
  # ==========================================================
  # PEAK INFECTIOUS BURDEN
  # ==========================================================
  
  total_I <- rowSums(
    out_df[
      ,
      paste0(
        "I_",
        1:n_schools
      )
    ],
    na.rm = TRUE
  )
  
  
  I_peak <- max(
    total_I,
    na.rm = TRUE
  )
  
  
  # ==========================================================
  # PEAK DAY
  # ==========================================================
  
  peak_day <- out_df$time[
    which.max(
      total_I
    )
  ]
  
  
  # ==========================================================
  # TOTAL NEW CASES
  #
  # REQUIRED FORMULA:
  #
  # C = sum_i max(0, R_i(T) - N_i*v_i)
  #
  # Since:
  #
  # R0_vec = round(N_i*v_i)
  #
  # the implementation below is:
  #
  # C = sum_i max(0, R_i(T) - R0_vec)
  # ==========================================================
  
  R_final <- out_df[
    nrow(out_df),
    paste0(
      "R_",
      1:n_schools
    )
  ]
  
  
  # Initial vaccinated population
  vaccinated_initial <- R0_vec
  
  
  # New cases in each school
  new_cases_each_school <- pmax(
    0,
    R_final -
      vaccinated_initial
  )
  
  
  # Total new cases
  total_new_cases <- sum(
    new_cases_each_school,
    na.rm = TRUE
  )
  
  
  # ==========================================================
  # COUNTY ATTACK RATE
  # ==========================================================
  
  county_attack_rate <-
    total_new_cases /
    N_total
  
  
  county_attack_rate_percent <-
    county_attack_rate * 100
  
  
  # ==========================================================
  # RETURN RESULTS
  # ==========================================================
  
  return(
    list(
      
      output = out_df,
      
      scenario = scenario_name,
      
      vaccination = vax_vector,
      
      total_new_cases =
        total_new_cases,
      
      peak_burden =
        I_peak,
      
      peak_day =
        peak_day,
      
      attack_rate =
        county_attack_rate,
      
      attack_rate_percent =
        county_attack_rate_percent,
      
      new_cases_each_school =
        new_cases_each_school
    )
  )
}


# ============================================================
# 21. RUN ALL SCENARIOS
# ============================================================

# ------------------------------------------------------------
# Uniform 70%
# ------------------------------------------------------------

result_70 <- run_vax_scenario(
  vax_vector = vax_70,
  seed_id = top_hub_id,
  scenario_name = "Uniform 70%"
)


# ------------------------------------------------------------
# Uniform 80%
# ------------------------------------------------------------

result_80 <- run_vax_scenario(
  vax_vector = vax_80,
  seed_id = top_hub_id,
  scenario_name = "Uniform 80%"
)


# ------------------------------------------------------------
# Scenario A
# 75% Small / 90% Large
# ------------------------------------------------------------

result_A <- run_vax_scenario(
  vax_vector = vax_A,
  seed_id = top_hub_id,
  scenario_name =
    "75% Small / 90% Large"
)


# ------------------------------------------------------------
# Scenario B
# Calculated Small / 75% Large
# ------------------------------------------------------------

result_B <- run_vax_scenario(
  vax_vector = vax_B,
  seed_id = top_hub_id,
  scenario_name =
    "Calculated Small / 75% Large"
)


# ------------------------------------------------------------
# Uniform 90%
# ------------------------------------------------------------

result_90 <- run_vax_scenario(
  vax_vector = vax_90,
  seed_id = top_hub_id,
  scenario_name = "Uniform 90%"
)


# ============================================================
# 22. CREATE RESULTS TABLE
# ============================================================

results_table <- tibble(
  
  Scenario = c(
    "Baseline",
    "Baseline",
    "Allocated",
    "Allocated",
    "Baseline"
  ),
  
  Allocation = c(
    
    "Uniform 70%",
    
    "Uniform 80%",
    
    "75% Small / 90% Large",
    
    paste0(
      round(
        vax_small_B_needed * 100,
        2
      ),
      "% Small / 75% Large"
    ),
    
    "Uniform 90%"
  ),
  
  Small_School_Coverage = c(
    
    70,
    
    80,
    
    75,
    
    vax_small_B_needed * 100,
    
    90
  ),
  
  Large_School_Coverage = c(
    
    70,
    
    80,
    
    90,
    
    75,
    
    90
  ),
  
  Overall_Coverage = c(
    
    70,
    
    80,
    
    rate_A * 100,
    
    rate_B * 100,
    
    90
  ),
  
  Total_New_Cases_C = c(
    
    result_70$total_new_cases,
    
    result_80$total_new_cases,
    
    result_A$total_new_cases,
    
    result_B$total_new_cases,
    
    result_90$total_new_cases
  ),
  
  Peak_Burden_Ipeak = c(
    
    result_70$peak_burden,
    
    result_80$peak_burden,
    
    result_A$peak_burden,
    
    result_B$peak_burden,
    
    result_90$peak_burden
  ),
  
  Peak_Day = c(
    
    result_70$peak_day,
    
    result_80$peak_day,
    
    result_A$peak_day,
    
    result_B$peak_day,
    
    result_90$peak_day
  ),
  
  Attack_Rate = c(
    
    result_70$attack_rate_percent,
    
    result_80$attack_rate_percent,
    
    result_A$attack_rate_percent,
    
    result_B$attack_rate_percent,
    
    result_90$attack_rate_percent
  )
)


# ============================================================
# 23. ROUND RESULTS FOR DISPLAY
# ============================================================

results_table_display <-
  results_table %>%
  
  mutate(
    
    Small_School_Coverage =
      round(
        Small_School_Coverage,
        2
      ),
    
    Large_School_Coverage =
      round(
        Large_School_Coverage,
        2
      ),
    
    Overall_Coverage =
      round(
        Overall_Coverage,
        2
      ),
    
    Total_New_Cases_C =
      round(
        Total_New_Cases_C,
        0
      ),
    
    Peak_Burden_Ipeak =
      round(
        Peak_Burden_Ipeak,
        0
      ),
    
    Peak_Day =
      round(
        Peak_Day,
        0
      ),
    
    Attack_Rate =
      round(
        Attack_Rate,
        2
      )
  )


# ============================================================
# 24. PRINT FINAL RESULTS TABLE
# ============================================================

cat("\n")
cat("====================================================\n")
cat("FINAL SENSITIVITY ANALYSIS 5 RESULTS\n")
cat("====================================================\n\n")

print(
  results_table_display
)


# ============================================================
# 25. SAVE RESULTS TABLE
# ============================================================

write_csv(
  results_table_display,
  "sensitivity_analysis_5_results.csv"
)


# ============================================================
# 26. PRINT SUMMARY OF EACH SCENARIO
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUMMARY\n")
cat("====================================================\n")


cat("\n--- Uniform 70% ---\n")

cat(
  "Total New Cases (C):",
  round(
    result_70$total_new_cases,
    0
  ),
  "\n"
)

cat(
  "Peak Burden:",
  round(
    result_70$peak_burden,
    0
  ),
  "\n"
)

cat(
  "Peak Day:",
  round(
    result_70$peak_day,
    0
  ),
  "\n"
)

cat(
  "Attack Rate:",
  round(
    result_70$attack_rate_percent,
    2
  ),
  "%\n"
)


cat("\n--- Uniform 80% ---\n")

cat(
  "Total New Cases (C):",
  round(
    result_80$total_new_cases,
    0
  ),
  "\n"
)

cat(
  "Peak Burden:",
  round(
    result_80$peak_burden,
    0
  ),
  "\n"
)

cat(
  "Peak Day:",
  round(
    result_80$peak_day,
    0
  ),
  "\n"
)

cat(
  "Attack Rate:",
  round(
    result_80$attack_rate_percent,
    2
  ),
  "%\n"
)


cat("\n--- Scenario A: 75% Small / 90% Large ---\n")

cat(
  "Total New Cases (C):",
  round(
    result_A$total_new_cases,
    0
  ),
  "\n"
)

cat(
  "Peak Burden:",
  round(
    result_A$peak_burden,
    0
  ),
  "\n"
)

cat(
  "Peak Day:",
  round(
    result_A$peak_day,
    0
  ),
  "\n"
)

cat(
  "Attack Rate:",
  round(
    result_A$attack_rate_percent,
    2
  ),
  "%\n"
)


cat("\n--- Scenario B: Calculated Small / 75% Large ---\n")

cat(
  "Total New Cases (C):",
  round(
    result_B$total_new_cases,
    0
  ),
  "\n"
)

cat(
  "Peak Burden:",
  round(
    result_B$peak_burden,
    0
  ),
  "\n"
)

cat(
  "Peak Day:",
  round(
    result_B$peak_day,
    0
  ),
  "\n"
)

cat(
  "Attack Rate:",
  round(
    result_B$attack_rate_percent,
    2
  ),
  "%\n"
)


cat("\n--- Uniform 90% ---\n")

cat(
  "Total New Cases (C):",
  round(
    result_90$total_new_cases,
    0
  ),
  "\n"
)

cat(
  "Peak Burden:",
  round(
    result_90$peak_burden,
    0
  ),
  "\n"
)

cat(
  "Peak Day:",
  round(
    result_90$peak_day,
    0
  ),
  "\n"
)

cat(
  "Attack Rate:",
  round(
    result_90$attack_rate_percent,
    2
  ),
  "%\n"
)


cat("\n")
cat("====================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("====================================================\n")




























# ============================================================
# CHECK POSSIBLE SCHOOL-SIZE CUTOFFS
# ============================================================

cutoffs <- c(400, 500, 600, 700, 800, 900, 1000)

cutoff_results <- map_dfr(cutoffs, function(cutoff) {
  
  small <- N < cutoff
  large <- N >= cutoff
  
  N_small_test <- sum(N[small], na.rm = TRUE)
  N_large_test <- sum(N[large], na.rm = TRUE)
  
  # Scenario A:
  # Small = 75%
  # Large = 90%
  
  total_vax_A_test <-
    (0.75 * N_small_test) +
    (0.90 * N_large_test)
  
  overall_A_test <-
    total_vax_A_test / N_total
  
  # Scenario B:
  # Small = 90%
  # Calculate large-school coverage
  # to preserve Scenario A overall coverage
  
  required_large <-
    (
      total_vax_A_test -
        (0.90 * N_small_test)
    ) / N_large_test
  
  tibble(
    Cutoff = cutoff,
    Small_Schools = sum(small),
    Large_Schools = sum(large),
    Students_Small = N_small_test,
    Students_Large = N_large_test,
    Percent_Students_Small =
      100 * N_small_test / N_total,
    Percent_Students_Large =
      100 * N_large_test / N_total,
    Scenario_A_Overall =
      100 * overall_A_test,
    Scenario_B_Large_Required =
      100 * required_large,
    Feasible =
      required_large >= 0 &
      required_large <= 1
  )
})

print(cutoff_results)














# ============================================================
# 27. PLOT COMPARATIVE EPIDEMIC CURVES
# ============================================================

# Extract infectious time series from each scenario result
extract_curve <- function(res, scenario_label) {
  out_df <- res$output
  i_cols <- paste0("I_", 1:n_schools)
  
  tibble(
    time = out_df$time,
    total_infectious = rowSums(out_df[, i_cols], na.rm = TRUE),
    Scenario = scenario_label
  )
}

# Combine time-series data for all 5 scenarios
plot_df <- bind_rows(
  extract_curve(result_70, "Baseline: Uniform 70%"),
  extract_curve(result_80, "Baseline: Uniform 80%"),
  extract_curve(result_B,  "Allocated: High in Small / Low in Large"),
  extract_curve(result_A,  "Allocated: Low in Small / High in Large"),
  extract_curve(result_90, "Baseline: Uniform 90%")
) %>%
  mutate(
    Scenario = factor(
      Scenario,
      levels = c(
        "Baseline: Uniform 70%",
        "Baseline: Uniform 80%",
        "Allocated: High in Small / Low in Large",
        "Allocated: Low in Small / High in Large",
        "Baseline: Uniform 90%"
      )
    )
  )

# Build plot object matching your exact image layout
p_vax <- ggplot(
  plot_df,
  aes(
    x = time,
    y = total_infectious,
    color = Scenario,
    linetype = Scenario
  )
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c(
      "Baseline: Uniform 70%"                   = "firebrick",
      "Baseline: Uniform 80%"                   = "darkorange",
      "Allocated: High in Small / Low in Large" = "deepskyblue3",
      "Allocated: Low in Small / High in Large" = "purple",
      "Baseline: Uniform 90%"                   = "forestgreen"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Baseline: Uniform 70%"                   = "dotted",
      "Baseline: Uniform 80%"                   = "dotted",
      "Allocated: High in Small / Low in Large" = "dashed",
      "Allocated: Low in Small / High in Large" = "solid",
      "Baseline: Uniform 90%"                   = "dotted"
    )
  ) +
  labs(
    title    = "Sensitivity Analysis 1: Size-Based Vaccine Allocation vs. Uniform Baselines",
    subtitle = "Comparing heterogeneous school coverage against uniform vaccination strategies",
    x        = "Days Since Outbreak Introduction",
    y        = "Total Active Infectious Cases in County"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.title    = element_blank()
  ) +
  guides(
    color    = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

# Explicitly display graph in RStudio Plots tab
print(p_vax)






