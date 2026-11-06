suppressPackageStartupMessages({  # load packages quietly
  library(tidyverse)  # data wrangling and plots
  library(deSolve)  # ODE solver
})

alpha_0 <- 0.05  # baseline between-school coupling
t_star  <- 14
times   <- seq(0, 200, by = 1)  # days to simulate
N       <- df$total_students  # students in each school

seed_school_id <- 1  # school where the index case starts

R0 <- round(N * df$vax_rate)  # students already immune
I0 <- rep(0, n_schools)  # index cases at the start
I0[seed_school_id] <- 1
E0 <- rep(0, n_schools)  # exposed at the start, usually none
S0 <- pmax(0, N - I0 - E0 - R0)  # students still susceptible

initial_state <- c(S0, E0, I0, R0)  # starting S, E, I, and R
names(initial_state) <- c(
  paste0("S_", 1:n_schools),
  paste0("E_", 1:n_schools),
  paste0("I_", 1:n_schools),
  paste0("R_", 1:n_schools)
)

params <- list(  # parameters passed into the solver
  beta  = (15 / 7),
  sigma = (1 / 10),
  gamma = (1 / 7),
  N     = N,
  A     = A,
  n     = n_schools
)

seir_model_phased <- function(time, state, parameters) {  # SEIR equations for one time step
  with(as.list(parameters), {

    current_alpha <- if (time < t_star) alpha_0 else (alpha_0 * r_factor)

    S <- state[1:n]
    E <- state[(n + 1):(2 * n)]
    I <- state[(2 * n + 1):(3 * n)]
    R <- state[(3 * n + 1):(4 * n)]

    I_prop <- I / N  # infectious fraction in each school
    lambda <- (beta * I_prop) + (current_alpha * beta * as.vector(A %*% I_prop))  # force of infection from this school and its neighbours

    dS <- -lambda * S  # infections leave the susceptible class
    dE <-  lambda * S - (sigma * E)  # new exposures minus those becoming infectious
    dI <- (sigma * E) - (gamma * I)  # new infectious students minus recoveries
    dR <- (gamma * I)  # students leaving the infectious class

    list(c(dS, dE, dI, dR))
  })
}

run_r_scenario_full <- function(r) {  # define the function
  p_temp <- params
  p_temp$r_factor <- r

  out <- ode(  # solve the SEIR equations
    y      = initial_state,
    times  = times,
    func   = seir_model_phased,
    parms  = p_temp,
    method = "lsoda"
  )

  out_df <- as.data.frame(out)

  I_cols <- paste0("I_", 1:n_schools)
  R_cols <- paste0("R_", 1:n_schools)

  total_I     <- rowSums(out_df[, I_cols, drop = FALSE], na.rm = TRUE)
  peak_burden <- max(total_I, na.rm = TRUE)
  peak_day    <- out_df$time[which.max(total_I)]

  R_final     <- sapply(R_cols, function(col) out_df[[col]][nrow(out_df)])  # define the function
  total_cases <- sum(pmax(0, R_final - R0), na.rm = TRUE)  # larger school id second

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
results_df <- bind_rows(lapply(r_vals, run_r_scenario_full))  # combine these rows

C_control <- results_df$C_raw[1]

delta_C <- 100 * (C_control - results_df$C_raw) / C_control

E_1_to_05  <- ((results_df$C_raw[2] - results_df$C_raw[1]) / results_df$C_raw[1]) / ((0.5 - 1.0) / 1.0)
E_05_to_01 <- ((results_df$C_raw[3] - results_df$C_raw[2]) / results_df$C_raw[2]) / ((0.1 - 0.5) / 0.5)

elasticity_col <- c(
  "N/A (Control)",
  sprintf("%.4f", E_1_to_05),
  sprintf("%.4f", E_05_to_01)
)

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

cat("\n========================================================================================\n")  # print a check
cat("SENSITIVITY ANALYSIS 3: PHASED TRAVEL RESTRICTIONS SUMMARY TABLE\n")  # print a check
cat("========================================================================================\n")  # print a check
print(report_table)  # show the result

cat("\n----------------------------------------------------------------------------------------\n")  # print a check
cat(sprintf(  # print a check
  "Report Sentence: \"Reducing inter-school connectivity to r = 0.1 after detection reduced total county cases by %.2f%% relative to no restriction.\"\n",
  delta_C[3]
))
cat("----------------------------------------------------------------------------------------\n")  # print a check

sens8_df %>%
  filter(time <= 14) %>%
  group_by(time) %>%
  summarise(
    min_I = min(total_infectious),
    max_I = max(total_infectious),
    is_identical = (max_I - min_I) < 1e-6
  ) %>%
  print(n = 15)  # show the result

sens8_df %>%
  filter(time <= 14) %>%
  group_by(time) %>%
  summarise(
    min_I = min(total_infectious),
    max_I = max(total_infectious),
    is_identical = (max_I - min_I) < 1e-6
  ) %>%
  print(n = 15)  # show the result
