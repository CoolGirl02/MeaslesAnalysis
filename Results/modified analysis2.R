# ============================================================
# NETWORK SEIR MEASLES SIMULATION AND SENSITIVITY ANALYSIS
# HUB VS PERIPHERAL SEEDING
#
# TOTAL NEW CASES:
# C(T) = sum_i [S_i(0) - S_i(T)]
#
# IMPORTANT:
# 1. The baseline infection starts in School 1.
# 2. Hub and peripheral scenarios start in their respective schools.
# 3. The adjacency matrix is normalized by its spectral radius.
# 4. Vaccinated students are included in the initial R compartment.
# 5. Total new cases are calculated using the decrease in S.
# ============================================================


suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(igraph)
})


theme_set(
  theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12)
    )
)


# ============================================================
# 1. LOAD SCHOOL DATA
# ============================================================

df <- read.csv(
  "data/cleanned/school_data_joined.csv"
) %>%
  rename_with(tolower)


# Ensure unit_id exists
if (!"unit_id" %in% names(df)) {
  
  df <- df %>%
    mutate(
      unit_id = as.character(row_number())
    )
  
} else {
  
  df <- df %>%
    mutate(
      unit_id = as.character(unit_id)
    )
  
}


# Prepare school data
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


cat("\nNumber of schools:", n_schools, "\n")


# ============================================================
# 2. BUILD NETWORK ADJACENCY MATRIX
# ============================================================

net_df <- read.csv(
  "data/cleanned/network_df.csv"
) %>%
  rename_with(tolower)


# Identify the two school ID columns
id_cols <- names(net_df)[
  grep(
    "unit|from|to|source|target|id",
    names(net_df)
  )
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


# Ensure school IDs use the same format
net_df <- net_df %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  )


# Map school IDs to row numbers
id_map <- setNames(
  seq_len(n_schools),
  df$unit_id
)


# Initialize adjacency matrix
A <- matrix(
  0,
  nrow = n_schools,
  ncol = n_schools
)


# Fill adjacency matrix using total contacts
for (k in seq_len(nrow(net_df))) {
  
  i <- unname(
    id_map[net_df$unit_id_a[k]]
  )
  
  j <- unname(
    id_map[net_df$unit_id_b[k]]
  )
  
  
  if (
    length(i) == 1 &&
    length(j) == 1 &&
    !is.na(i) &&
    !is.na(j)
  ) {
    
    w <- if (
      "total_contacts" %in% names(net_df)
    ) {
      
      net_df$total_contacts[k]
      
    } else {
      
      1
      
    }
    
    
    if (is.na(w) || !is.finite(w) || w < 0) {
      w <- 0
    }
    
    
    A[i, j] <- w
    A[j, i] <- w
    
  }
  
}


# ============================================================
# 3. NORMALIZE NETWORK BY SPECTRAL RADIUS
# ============================================================

A_raw <- A


spectral_radius <- max(
  Mod(
    eigen(
      A_raw,
      only.values = TRUE
    )$values
  )
)


if (
  !is.finite(spectral_radius) ||
  spectral_radius <= 0
) {
  
  stop(
    "The network spectral radius is zero or invalid. ",
    "Check the school IDs and network contact weights."
  )
  
}


# Normalize the adjacency matrix
A <- A_raw / spectral_radius


cat("\n========================================\n")
cat("NETWORK NORMALIZATION\n")
cat("========================================\n")

cat(
  "Original spectral radius:",
  spectral_radius,
  "\n"
)

cat(
  "Normalized spectral radius:",
  max(
    Mod(
      eigen(
        A,
        only.values = TRUE
      )$values
    )
  ),
  "\n"
)


# ============================================================
# 4. MODEL PARAMETERS
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
# 5. POPULATION AND BASELINE INITIAL CONDITIONS
# ============================================================

N <- df$total_students


# Initial vaccinated/recovered population
R0 <- round(
  N * df$vax_rate
)


# Baseline infection starts in School 1
I0 <- c(
  1,
  rep(0, n_schools - 1)
)


# Initial exposed population
E0 <- rep(
  0,
  n_schools
)


# Initial susceptible population
S0 <- pmax(
  0,
  N - I0 - E0 - R0
)


# Initial state
initial_state <- c(
  S0,
  E0,
  I0,
  R0
)


names(initial_state) <- c(
  paste0("S_", seq_len(n_schools)),
  paste0("E_", seq_len(n_schools)),
  paste0("I_", seq_len(n_schools)),
  paste0("R_", seq_len(n_schools))
)


# ============================================================
# 6. MODEL PARAMETERS LIST
# ============================================================

params <- list(
  beta = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N = N,
  A = A,
  n = n_schools
)


# ============================================================
# 7. NETWORK SEIR MODEL
# ============================================================

seir_model <- function(
    time,
    state,
    parameters
) {
  
  with(
    as.list(parameters),
    {
      
      # Extract compartments
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
      
      
      # Proportion infectious in each school
      I_prop <- I / N
      
      
      # Force of infection
      lambda <- (
        beta * I_prop
      ) + (
        alpha *
          beta *
          as.vector(
            A %*% I_prop
          )
      )
      
      
      # SEIR differential equations
      dS <- -lambda * S
      
      dE <- (
        lambda * S
      ) - (
        sigma * E
      )
      
      dI <- (
        sigma * E
      ) - (
        gamma * I
      )
      
      dR <- gamma * I
      
      
      list(
        c(
          dS,
          dE,
          dI,
          dR
        )
      )
      
    }
  )
  
}


# ============================================================
# 8. FUNCTION TO CALCULATE TOTAL INFECTIOUS CASES OVER TIME
# ============================================================

get_total_I <- function(
    out_df,
    scenario_label,
    seed_type
) {
  
  out_df %>%
    
    pivot_longer(
      cols = -time,
      names_to = "var",
      values_to = "Count"
    ) %>%
    
    filter(
      str_detect(
        var,
        "^I_"
      )
    ) %>%
    
    group_by(time) %>%
    
    summarise(
      total_infectious = sum(
        Count,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    
    mutate(
      Scenario = scenario_label,
      `Seed Type` = seed_type
    )
  
}


# ============================================================
# 9. FUNCTION TO RUN VACCINATION SCENARIO
# ============================================================

run_vax_scenario <- function(
    vax_vector,
    seed_id
) {
  
  # Calculate initial vaccinated/recovered population
  R0_vec <- round(
    N * vax_vector
  )
  
  
  # Initial infected population
  I0_vec <- rep(
    0,
    n_schools
  )
  
  I0_vec[seed_id] <- 1
  
  
  # Initial exposed population
  E0_vec <- rep(
    0,
    n_schools
  )
  
  
  # Initial susceptible population
  S0_vec <- pmax(
    0,
    N - I0_vec - E0_vec - R0_vec
  )
  
  
  # Construct initial state
  state_vec <- c(
    S0_vec,
    E0_vec,
    I0_vec,
    R0_vec
  )
  
  
  names(state_vec) <- names(
    initial_state
  )
  
  
  # Run SEIR model
  out_df <- as.data.frame(
    ode(
      y = state_vec,
      times = times,
      func = seir_model,
      parms = params
    )
  )
  
  
  return(out_df)
  
}


# ============================================================
# 10. BASELINE SEIR SIMULATION
# ============================================================

output <- ode(
  y = initial_state,
  times = times,
  func = seir_model,
  parms = params
)


output <- as.data.frame(output)


# ============================================================
# 11. BASELINE RESULTS BY SCHOOL
# ============================================================

results_long <- output %>%
  
  pivot_longer(
    cols = -time,
    names_to = "var",
    values_to = "Count"
  ) %>%
  
  mutate(
    
    Compartment = str_extract(
      var,
      "^[SEIR]"
    ),
    
    school_id = as.numeric(
      str_extract(
        var,
        "\\d+"
      )
    )
    
  ) %>%
  
  left_join(
    df,
    by = "school_id"
  ) %>%
  
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


# ============================================================
# 12. BASELINE SCHOOL PLOTS
# ============================================================

total_pages <- max(
  df$page
)


for (p in seq_len(total_pages)) {
  
  page_data <- results_long %>%
    filter(
      page == p
    )
  
  
  p_plot <- ggplot(
    page_data,
    aes(
      x = time,
      y = Count,
      color = Compartment
    )
  ) +
    
    geom_line(
      linewidth = 1
    ) +
    
    scale_color_manual(
      values = c(
        "S" = "blue",
        "E" = "orange",
        "I" = "red",
        "R" = "green"
      ),
      
      labels = c(
        "S" = "Susceptible",
        "E" = "Exposed",
        "I" = "Infectious",
        "R" = "Recovered/Vaccinated"
      )
    ) +
    
    facet_wrap(
      ~school_label,
      ncol = 3,
      scales = "free_y"
    ) +
    
    labs(
      title = paste0(
        "Baseline SEIR Measles Simulations (Page ",
        p,
        " of ",
        total_pages,
        ")"
      ),
      
      x = "Days",
      
      y = "Number of Students",
      
      color = "Compartment"
    ) +
    
    theme(
      strip.text = element_text(
        size = 8,
        face = "bold"
      ),
      
      legend.position = "bottom"
    )
  
  
  print(p_plot)
  
}


# ============================================================
# 13. NETWORK CENTRALITY
# ============================================================

net_graph <- graph_from_adjacency_matrix(
  A,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)


# Calculate node strength
df$node_strength <- strength(
  net_graph,
  weights = E(net_graph)$weight
)


# Calculate betweenness centrality
# Unweighted distances are used here
df$betweenness <- betweenness(
  net_graph,
  weights = NA
)


# Hub: school with the highest node strength
top_hub_id <- df$school_id[
  which.max(
    df$node_strength
  )
]


# Peripheral: school with the lowest node strength
low_hub_id <- df$school_id[
  which.min(
    df$node_strength
  )
]


top_hub_name <- df$school_name[
  df$school_id == top_hub_id
]


low_hub_name <- df$school_name[
  df$school_id == low_hub_id
]


cat("\n========================================\n")
cat("NETWORK HUB AND PERIPHERAL SCHOOL\n")
cat("========================================\n")

cat(
  "Network Hub:",
  top_hub_name,
  "\n"
)

cat(
  "Hub node strength:",
  df$node_strength[top_hub_id],
  "\n"
)

cat(
  "Peripheral:",
  low_hub_name,
  "\n"
)

cat(
  "Peripheral node strength:",
  df$node_strength[low_hub_id],
  "\n"
)


# ============================================================
# 14. SCENARIO 1: OBSERVED VACCINATION COVERAGE
# ============================================================

vax_empirical <- df$vax_rate


# Outbreak introduced in hub
out_hub_baseline <- run_vax_scenario(
  vax_empirical,
  top_hub_id
)


# Outbreak introduced in peripheral school
out_peri_baseline <- run_vax_scenario(
  vax_empirical,
  low_hub_id
)


# ============================================================
# 15. SCENARIO 2: BOTH SEED SCHOOLS AT 74% COVERAGE
# ============================================================

vax_74 <- df$vax_rate


vax_74[top_hub_id] <- 0.743

vax_74[low_hub_id] <- 0.743


# Hub seeding
out_hub_74 <- run_vax_scenario(
  vax_74,
  top_hub_id
)


# Peripheral seeding
out_peri_74 <- run_vax_scenario(
  vax_74,
  low_hub_id
)


# ============================================================
# 16. SCENARIO 3: BOTH SEED SCHOOLS AT 91.6% COVERAGE
# ============================================================

vax_916 <- df$vax_rate


vax_916[top_hub_id] <- 0.916

vax_916[low_hub_id] <- 0.916


# Hub seeding
out_hub_916 <- run_vax_scenario(
  vax_916,
  top_hub_id
)


# Peripheral seeding
out_peri_916 <- run_vax_scenario(
  vax_916,
  low_hub_id
)


# ============================================================
# 17. COMBINE ALL SCENARIOS FOR PLOTTING
# ============================================================

comparison_all <- bind_rows(
  
  # Scenario 1: Observed vaccination coverage
  
  get_total_I(
    out_hub_baseline,
    "Scenario 1: Observed Coverage",
    "Network Hub"
  ),
  
  get_total_I(
    out_peri_baseline,
    "Scenario 1: Observed Coverage",
    "Peripheral"
  ),
  
  
  # Scenario 2: Both seed schools at 74%
  
  get_total_I(
    out_hub_74,
    "Scenario 2: Both Seed Schools at 74.3%",
    "Network Hub"
  ),
  
  get_total_I(
    out_peri_74,
    "Scenario 2: Both Seed Schools at 74.3%",
    "Peripheral"
  ),
  
  
  # Scenario 3: Both seed schools at 91.6%
  
  get_total_I(
    out_hub_916,
    "Scenario 3: Both Seed Schools at 91.6%",
    "Network Hub"
  ),
  
  get_total_I(
    out_peri_916,
    "Scenario 3: Both Seed Schools at 91.6%",
    "Peripheral"
  )
  
)


# ============================================================
# 18. PLOT HUB VS PERIPHERAL SEEDING
# ============================================================

p_comparison <- ggplot(
  comparison_all,
  aes(
    x = time,
    y = total_infectious,
    color = `Seed Type`
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  facet_wrap(
    ~Scenario,
    ncol = 3,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      "Network Hub" = "firebrick",
      "Peripheral" = "darkblue"
    )
  ) +
  
  labs(
    title = "Sensitivity Analysis: Hub vs Peripheral Seeding",
    
    subtitle = paste(
      "Comparison across three vaccination coverage scenarios"
    ),
    
    x = "Days Since Outbreak Introduction",
    
    y = "Number of  Infectious Cases",
    
    color = "Seed Location"
  ) +
  
  theme(
    legend.position = "bottom",
    
    strip.text = element_text(
      face = "bold",
      size = 10
    )
  )


print(p_comparison)


# ============================================================
# 19. FUNCTION TO CALCULATE OUTBREAK MEASURES
#
# TOTAL NEW CASES:
#
# C(T) = sum_i [S_i(0) - S_i(T)]
#
# This replaces the previous calculation based on R.
# ============================================================

calculate_outbreak_measures <- function(
    out_df,
    seed_id
) {
  
  # ----------------------------------------------------------
  # 19.1 Identify susceptible and infectious columns
  # ----------------------------------------------------------
  
  S_columns <- grep(
    "^S_",
    names(out_df),
    value = TRUE
  )
  
  
  I_columns <- grep(
    "^I_",
    names(out_df),
    value = TRUE
  )
  
  
  # Check that the expected columns exist
  if (
    length(S_columns) != n_schools ||
    length(I_columns) != n_schools
  ) {
    
    stop(
      "The output does not contain the expected ",
      "susceptible and infectious columns."
    )
    
  }
  
  
  # ----------------------------------------------------------
  # 19.2 Initial susceptible population S_i(0)
  # ----------------------------------------------------------
  
  S_initial <- sapply(
    S_columns,
    function(x) {
      out_df[[x]][1]
    }
  )
  
  
  # ----------------------------------------------------------
  # 19.3 Final susceptible population S_i(T)
  # ----------------------------------------------------------
  
  S_final <- sapply(
    S_columns,
    function(x) {
      out_df[[x]][nrow(out_df)]
    }
  )
  
  
  # ----------------------------------------------------------
  # 19.4 Calculate new cases using the requested formula
  #
  # C(T) = sum_i [S_i(0) - S_i(T)]
  # ----------------------------------------------------------
  
  new_cases_each_school <- S_initial - S_final
  
  
  total_cases <- sum(
    new_cases_each_school,
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # 19.5 Calculate total infectious population over time
  # ----------------------------------------------------------
  
  total_I <- rowSums(
    out_df[
      ,
      I_columns,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # 19.6 Peak infectious burden
  # ----------------------------------------------------------
  
  peak_burden <- max(
    total_I,
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # 19.7 Day of peak infectious burden
  # ----------------------------------------------------------
  
  peak_day <- out_df$time[
    which.max(total_I)
  ]
  
  
  # ----------------------------------------------------------
  # 19.8 Total county population
  # ----------------------------------------------------------
  
  county_population <- sum(
    N,
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # 19.9 County attack rate
  # ----------------------------------------------------------
  
  attack_rate <- (
    total_cases / county_population
  ) * 100
  
  
  # ----------------------------------------------------------
  # 19.10 Return outbreak measures
  # ----------------------------------------------------------
  
  return(
    c(
      Total_Cases = total_cases,
      Peak_Burden = peak_burden,
      Peak_Day = peak_day,
      Attack_Rate = attack_rate
    )
  )
  
}


# ============================================================
# 20. CALCULATE SCENARIO 1 RESULTS
# ============================================================

baseline_hub <- calculate_outbreak_measures(
  out_hub_baseline,
  top_hub_id
)


baseline_peri <- calculate_outbreak_measures(
  out_peri_baseline,
  low_hub_id
)


# ============================================================
# 21. CALCULATE SCENARIO 2 RESULTS
# ============================================================

control_74_hub <- calculate_outbreak_measures(
  out_hub_74,
  top_hub_id
)


control_74_peri <- calculate_outbreak_measures(
  out_peri_74,
  low_hub_id
)


# ============================================================
# 22. CALCULATE SCENARIO 3 RESULTS
# ============================================================

control_916_hub <- calculate_outbreak_measures(
  out_hub_916,
  top_hub_id
)


control_916_peri <- calculate_outbreak_measures(
  out_peri_916,
  low_hub_id
)


# ============================================================
# 23. CREATE TABLE 4.3
# ============================================================

results_4_3 <- tibble(
  
  Scenario = c(
    
    "Scenario 1: Observed Coverage",
    "Scenario 1: Observed Coverage",
    
    "Scenario 2: Both Seed Schools at 74%",
    "Scenario 2: Both Seed Schools at 74%",
    
    "Scenario 3: Both Seed Schools at 91.6%",
    "Scenario 3: Both Seed Schools at 91.6%"
    
  ),
  
  
  `Seed Location` = c(
    
    top_hub_name,
    low_hub_name,
    
    top_hub_name,
    low_hub_name,
    
    top_hub_name,
    low_hub_name
    
  ),
  
  
  `Seed Type` = c(
    
    "Network Hub",
    "Peripheral",
    
    "Network Hub",
    "Peripheral",
    
    "Network Hub",
    "Peripheral"
    
  ),
  
  
  `Vaccination Coverage` = c(
    
    paste0(
      round(
        vax_empirical[top_hub_id] * 100,
        1
      ),
      "%"
    ),
    
    paste0(
      round(
        vax_empirical[low_hub_id] * 100,
        1
      ),
      "%"
    ),
    
    "74.3%",
    "74.3%",
    
    "91.6%",
    "91.6%"
    
  ),
  
  
  `Total New Cases (C)` = c(
    
    baseline_hub["Total_Cases"],
    baseline_peri["Total_Cases"],
    
    control_74_hub["Total_Cases"],
    control_74_peri["Total_Cases"],
    
    control_916_hub["Total_Cases"],
    control_916_peri["Total_Cases"]
    
  ),
  
  
  `Peak Burden (Ipeak)` = c(
    
    baseline_hub["Peak_Burden"],
    baseline_peri["Peak_Burden"],
    
    control_74_hub["Peak_Burden"],
    control_74_peri["Peak_Burden"],
    
    control_916_hub["Peak_Burden"],
    control_916_peri["Peak_Burden"]
    
  ),
  
  
  `Peak Day` = c(
    
    baseline_hub["Peak_Day"],
    baseline_peri["Peak_Day"],
    
    control_74_hub["Peak_Day"],
    control_74_peri["Peak_Day"],
    
    control_916_hub["Peak_Day"],
    control_916_peri["Peak_Day"]
    
  ),
  
  
  `Attack Rate (AR)` = c(
    
    baseline_hub["Attack_Rate"],
    baseline_peri["Attack_Rate"],
    
    control_74_hub["Attack_Rate"],
    control_74_peri["Attack_Rate"],
    
    control_916_hub["Attack_Rate"],
    control_916_peri["Attack_Rate"]
    
  )
  
)


# ============================================================
# 24. ROUND TABLE VALUES
# ============================================================

table_4_3 <- results_4_3 %>%
  
  mutate(
    
    `Total New Cases (C)` = round(
      `Total New Cases (C)`
    ),
    
    `Peak Burden (Ipeak)` = round(
      `Peak Burden (Ipeak)`
    ),
    
    `Peak Day` = round(
      `Peak Day`
    ),
    
    `Attack Rate (AR)` = paste0(
      round(
        `Attack Rate (AR)`,
        2
      ),
      "%"
    )
    
  )


# ============================================================
# 25. PRINT TABLE 4.3
# ============================================================

cat("\n\n")

cat(
  "============================================================\n"
)

cat(
  "TABLE 4.3: COMPARATIVE ANALYSIS OF HUB VS PERIPHERAL\n"
)

cat(
  "SEEDING SCENARIOS\n"
)

cat(
  "============================================================\n\n"
)


print(
  table_4_3,
  n = Inf
)


# ============================================================
# 26. OPTIONAL: DISPLAY TABLE USING KNITR
# ============================================================

if (
  requireNamespace(
    "knitr",
    quietly = TRUE
  )
) {
  
  print(
    knitr::kable(
      table_4_3,
      
      caption = paste(
        "Comparative Analysis of Hub vs Peripheral",
        "Seeding Scenarios"
      ),
      
      align = c(
        "l",
        "l",
        "l",
        "l",
        "r",
        "r",
        "r",
        "r"
      )
    )
  )
  
}


# ============================================================
# 27. SAVE TABLE 4.3
# ============================================================

write.csv(
  table_4_3,
  "table_4_3_hub_vs_peripheral.csv",
  row.names = FALSE
)


cat("\nTable 4.3 saved successfully.\n")


# ============================================================
# 28. DISPLAY TOTAL NEW CASES CALCULATED FROM SUSCEPTIBLE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TOTAL NEW CASES CALCULATED USING SUSCEPTIBLE POPULATION\n")
cat("============================================================\n")

cat(
  "Formula: C(T) = sum_i [S_i(0) - S_i(T)]\n\n"
)


cat(
  "Scenario 1 - Hub:",
  round(baseline_hub["Total_Cases"]),
  "cases\n"
)

cat(
  "Scenario 1 - Peripheral:",
  round(baseline_peri["Total_Cases"]),
  "cases\n"
)


cat(
  "Scenario 2 - Hub:",
  round(control_74_hub["Total_Cases"]),
  "cases\n"
)

cat(
  "Scenario 2 - Peripheral:",
  round(control_74_peri["Total_Cases"]),
  "cases\n"
)


cat(
  "Scenario 3 - Hub:",
  round(control_916_hub["Total_Cases"]),
  "cases\n"
)

cat(
  "Scenario 3 - Peripheral:",
  round(control_916_peri["Total_Cases"]),
  "cases\n"
)


cat("\n")
cat("============================================================\n")
cat("SIMULATION AND SENSITIVITY ANALYSIS COMPLETED\n")
cat("============================================================\n")


