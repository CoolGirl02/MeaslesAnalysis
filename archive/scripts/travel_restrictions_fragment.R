# ============================================================
# SENSITIVITY ANALYSIS 3: PHASED INTER-SCHOOL TRAVEL RESTRICTIONS
# COMPREHENSIVE METRICS: PEAK BURDEN, PEAK DAY, TOTAL CASES, AR (%)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
})

# ------------------------------------------------------------
# 1. SETUP MODEL & PARAMETERS
# ------------------------------------------------------------

alpha_0 <- 0.05
t_star  <- 14
times   <- seq(0, 200, by = 1)
N       <- df$total_students

# Set seed school index (empirical data preserved across all schools)
seed_school_id <- 1  

R0 <- round(N * df$vax_rate)
I0 <- rep(0, n_schools)
I0[seed_school_id] <- 1
E0 <- rep(0, n_schools)
S0 <- pmax(0, N - I0 - E0 - R0)

initial_state <- c(S0, E0, I0, R0)
names(initial_state) <- c(
  paste0("S_", 1:n_schools),
  paste0("E_", 1:n_schools),
  paste0("I_", 1:n_schools),
  paste0("R_", 1:n_schools)
)

params <- list(
  beta  = (15 / 7),
  sigma = (1 / 10),
  gamma = (1 / 7),
  N     = N,
  A     = A,
  n     = n_schools
)

# ------------------------------------------------------------
# 2. DYNAMIC ODE MODEL (Eq. 13)
# ------------------------------------------------------------

seir_model_phased <- function(time, state, parameters) {
  with(as.list(parameters), {
    # Equation 13: alpha(t) = alpha_0 for t < t*, alpha_0 * r for t >= t*
    current_alpha <- if (time < t_star) alpha_0 else (alpha_0 * r_factor)
    
    S <- state[1:n]
    E <- state[(n + 1):(2 * n)]
    I <- state[(2 * n + 1):(3 * n)]
    R <- state[(3 * n + 1):(4 * n)]
    
    I_prop <- I / N
    lambda <- (beta * I_prop) + (current_alpha * beta * as.vector(A %*% I_prop))
    
    dS <- -lambda * S
    dE <-  lambda * S - (sigma * E)
    dI <- (sigma * E) - (gamma * I)
    dR <- (gamma * I)
    
    list(c(dS, dE, dI, dR))
  })
}

# ------------------------------------------------------------
# 3. RUN SCENARIOS & EXTRACT COMPREHENSIVE EPIDEMIC METRICS
# ------------------------------------------------------------

run_r_scenario_full <- function(r) {
  p_temp <- params
  p_temp$r_factor <- r
  
  out <- ode(
    y      = initial_state,
    times  = times,
    func   = seir_model_phased,
    parms  = p_temp,
    method = "lsoda"
  )
  
  out_df <- as.data.frame(out)
  
  I_cols <- paste0("I_", 1:n_schools)
  R_cols <- paste0("R_", 1:n_schools)
  
  # 1. Peak Burden (Ipeak) & Peak Day (tpeak)
  total_I     <- rowSums(out_df[, I_cols, drop = FALSE], na.rm = TRUE)
  peak_burden <- max(total_I, na.rm = TRUE)
  peak_day    <- out_df$time[which.max(total_I)]
  
  # 2. Total Cumulative Cases (C)
  R_final     <- sapply(R_cols, function(col) out_df[[col]][nrow(out_df)])
  total_cases <- sum(pmax(0, R_final - R0), na.rm = TRUE)
  
  # 3. County Attack Rate (%)
  county_ar   <- (total_cases / sum(N, na.rm = TRUE)) * 100
  
  data.frame(
    r           = r,
    alpha_post  = alpha_0 * r,
    I_peak      = round(peak_burden),
    t_peak      = peak_day,
    C_raw       = total_cases,
    C_formatted = format(round(total_cases), big.mark = ","),
    AR_percent  = county_ar
  )
}

r_vals    <- c(1.0, 0.5, 0.1)
results_df <- bind_rows(lapply(r_vals, run_r_scenario_full))

# ------------------------------------------------------------
# 4. COMPUTE PERCENTAGE REDUCTION & ELASTICITY
# ------------------------------------------------------------

C_control <- results_df$C_raw[1] # Control C(1)

# Eq 14: Percent reduction Delta_C(r)
delta_C <- 100 * (C_control - results_df$C_raw) / C_control

# Eq 15: Elasticity E between adjacent r values
E_1_to_05  <- ((results_df$C_raw[2] - results_df$C_raw[1]) / results_df$C_raw[1]) / ((0.5 - 1.0) / 1.0)
E_05_to_01 <- ((results_df$C_raw[3] - results_df$C_raw[2]) / results_df$C_raw[2]) / ((0.1 - 0.5) / 0.5)

elasticity_col <- c(
  "N/A (Control)",
  sprintf("%.4f", E_1_to_05),
  sprintf("%.4f", E_05_to_01)
)

# ------------------------------------------------------------
# 5. CONSTRUCT & DISPLAY REPORT TABLE
# ------------------------------------------------------------

report_table <- data.frame(
  `r (Factor)`              = results_df$r,
  `Alpha (t >= 14)`         = results_df$alpha_post,
  `Peak Burden (Ipeak)`     = results_df$I_peak,
  `Peak Day (tpeak)`        = results_df$t_peak,
  `Total Cases (C)`         = results_df$C_formatted,
  `County AR (%)`           = sprintf("%.2f%%", results_df$AR_percent),
  `Delta_C(r) (%)`          = sprintf("%.2f%%", delta_C),
  `Elasticity (E)`          = elasticity_col,
  check.names               = FALSE
)

cat("\n========================================================================================\n")
cat("SENSITIVITY ANALYSIS 3: PHASED TRAVEL RESTRICTIONS SUMMARY TABLE\n")
cat("========================================================================================\n")
print(report_table)

cat("\n----------------------------------------------------------------------------------------\n")
cat(sprintf(
  "Report Sentence: \"Reducing inter-school connectivity to r = 0.1 after detection reduced total county cases by %.2f%% relative to no restriction.\"\n",
  delta_C[3]
))
cat("----------------------------------------------------------------------------------------\n")






# Check if values before Day 14 are identical across all scenarios
sens8_df %>%
  filter(time <= 14) %>%
  group_by(time) %>%
  summarise(
    min_I = min(total_infectious),
    max_I = max(total_infectious),
    is_identical = (max_I - min_I) < 1e-6
  ) %>%
  print(n = 15)
# Check if values before Day 14 are identical across all scenarios
sens8_df %>%
  filter(time <= 14) %>%
  group_by(time) %>%
  summarise(
    min_I = min(total_infectious),
    max_I = max(total_infectious),
    is_identical = (max_I - min_I) < 1e-6
  ) %>%
  print(n = 15)
