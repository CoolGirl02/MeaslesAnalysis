VACCINE ALLOCATION BY SCHOOL SIZE VS. BASELINES
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



