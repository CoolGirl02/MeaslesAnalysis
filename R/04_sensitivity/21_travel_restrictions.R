suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(igraph)
})

if (!exists("df")) {
  stop(
    "The school data object 'df' does not exist. ",
    "Please load and prepare your school data first."
  )
}

if (!exists("A")) {
  stop(
    "The adjacency matrix 'A' does not exist. ",
    "Please create the network adjacency matrix first."
  )
}

df <- as.data.frame(df)

required_columns <- c("total_students", "vax_rate")

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in df: ",
    paste(missing_columns, collapse = ", ")
  )
}

df$total_students <- as.numeric(df$total_students)
df$vax_rate <- as.numeric(df$vax_rate)

if (any(is.na(df$total_students)) ||
    any(!is.finite(df$total_students)) ||
    any(df$total_students <= 0)) {
  stop("total_students must contain positive, finite values.")
}

if (any(is.na(df$vax_rate)) ||
    any(!is.finite(df$vax_rate)) ||
    any(df$vax_rate < 0 | df$vax_rate > 1)) {
  stop(
    "vax_rate must contain proportions between 0 and 1. ",
    "For example, use 0.90 for 90% coverage, not 90."
  )
}

n <- nrow(df)

cat("Number of schools:", n, "\n")

A <- as.matrix(A)

storage.mode(A) <- "numeric"

if (nrow(A) != n || ncol(A) != n) {
  stop(
    "The adjacency matrix must have the same number of rows and columns ",
    "as the number of schools in df."
  )
}

if (any(!is.finite(A))) {
  stop("The adjacency matrix contains missing or non-finite values.")
}

if (any(A < 0)) {
  stop("The adjacency matrix cannot contain negative values.")
}

diag(A) <- 0

if (all(A == 0)) {
  stop("The adjacency matrix contains no inter-school connections.")
}

spectral_radius_before <- max(
  Mod(eigen(A, only.values = TRUE)$values)
)

if (!is.finite(spectral_radius_before) ||
    spectral_radius_before <= 0) {
  stop("The spectral radius must be positive.")
}

A <- A / spectral_radius_before

spectral_radius_after <- max(
  Mod(eigen(A, only.values = TRUE)$values)
)

cat("\nNETWORK NORMALIZATION\n")
cat("---------------------\n")
cat(
  "Spectral radius before normalization:",
  spectral_radius_before, "\n"
)
cat(
  "Spectral radius after normalization:",
  spectral_radius_after, "\n"
)

gamma <- 1 / 7
sigma <- 1 / 10

R0_target <- 15
beta <- R0_target * gamma

alpha_0 <- 0.05

t_star <- 14

times <- seq(0, 200, by = 1)

cat("\nMODEL PARAMETERS\n")
cat("----------------\n")
cat("Transmission rate beta:", beta, "\n")
cat("Progression rate sigma:", sigma, "\n")
cat("Recovery rate gamma:", gamma, "\n")
cat("Baseline inter-school multiplier:", alpha_0, "\n")
cat("Restrictions begin on day:", t_star, "\n")

N <- df$total_students

R_vax <- round(N * df$vax_rate)

E0 <- rep(0, n)

I0 <- rep(0, n)

seed_school_id <- 1

if (seed_school_id < 1 || seed_school_id > n) {
  stop("The seed school index is outside the range of available schools.")
}

I0[seed_school_id] <- 1

cat("\nINITIAL INFECTION\n")
cat("-----------------\n")
cat("Infection starts in school number:", seed_school_id, "\n")

S0 <- N - R_vax - I0

if (any(S0 < 0)) {
  stop(
    "Some schools have negative initial susceptible populations. ",
    "Check vaccination coverage, school populations, and the initial infection."
  )
}

R0_initial <- R_vax

S_names <- paste0("S", seq_len(n))
E_names <- paste0("E", seq_len(n))
I_names <- paste0("I", seq_len(n))
R_names <- paste0("R", seq_len(n))

initial_state <- c(
  S0,
  E0,
  I0,
  R0_initial
)

names(initial_state) <- c(
  S_names,
  E_names,
  I_names,
  R_names
)

network_seir_model <- function(time, state, parms) {

  N <- parms$N
  A <- parms$A
  beta <- parms$beta
  sigma <- parms$sigma
  gamma <- parms$gamma
  alpha_0 <- parms$alpha_0
  t_star <- parms$t_star
  r_factor <- parms$r_factor

  n <- length(N)

  S <- as.numeric(state[paste0("S", seq_len(n))])
  E <- as.numeric(state[paste0("E", seq_len(n))])
  I <- as.numeric(state[paste0("I", seq_len(n))])
  R <- as.numeric(state[paste0("R", seq_len(n))])

  if (time < t_star) {
    current_alpha <- alpha_0
  } else {
    current_alpha <- alpha_0 * r_factor
  }

  I_prop <- I / N

  lambda <- (beta * I_prop) +
    (current_alpha * beta * as.vector(A %*% I_prop))

  lambda <- pmax(lambda, 0)

  dS <- -lambda * S
  dE <- lambda * S - sigma * E
  dI <- sigma * E - gamma * I
  dR <- gamma * I

  names(dS) <- paste0("S", seq_len(n))
  names(dE) <- paste0("E", seq_len(n))
  names(dI) <- paste0("I", seq_len(n))
  names(dR) <- paste0("R", seq_len(n))

  list(c(dS, dE, dI, dR))
}

run_r_scenario_full <- function(r_factor) {

  parms <- list(
    N = N,
    A = A,
    beta = beta,
    sigma = sigma,
    gamma = gamma,
    alpha_0 = alpha_0,
    t_star = t_star,
    r_factor = r_factor
  )

  out <- deSolve::ode(
    y = initial_state,
    times = times,
    func = network_seir_model,
    parms = parms,
    method = "lsoda"
  )

  out_df <- as.data.frame(out)

  I_cols <- paste0("I", seq_len(n))

  out_df$total_infectious <- rowSums(
    out_df[, I_cols, drop = FALSE]
  )

  peak_index <- which.max(out_df$total_infectious)

  peak_burden <- out_df$total_infectious[peak_index]
  peak_day <- out_df$time[peak_index]

  S_cols <- paste0("S", seq_len(n))

  S_final <- as.numeric(
    out_df[nrow(out_df), S_cols, drop = TRUE]
  )

  total_cases <- sum(S0 - S_final, na.rm = TRUE)

  county_attack_rate <- (
    total_cases / sum(N)
  ) * 100

  scenario_label <- case_when(
    r_factor == 1.0 ~ "No reduction",
    r_factor == 0.5 ~ "50% reduction",
    r_factor == 0.1 ~ "90% reduction",
    TRUE ~ paste0(
      round((1 - r_factor) * 100),
      "% reduction"
    )
  )

  summary_row <- tibble(
    Scenario = scenario_label,
    Reduction_Factor = r_factor,
    Alpha_After_Day_14 = alpha_0 * r_factor,
    Peak_Infectious_Burden = peak_burden,
    Peak_Day = peak_day,
    Total_New_Cases = total_cases,
    County_Attack_Rate_Percent = county_attack_rate
  )

  curve_df <- tibble(
    time = out_df$time,
    Total_Infectious = out_df$total_infectious,
    Scenario = scenario_label,
    Reduction_Factor = r_factor
  )

  list(
    summary = summary_row,
    curve = curve_df,
    model_output = out_df
  )
}

r_vals <- c(1.0, 0.5, 0.1)

scenario_runs <- lapply(
  r_vals,
  run_r_scenario_full
)

results_df <- bind_rows(
  lapply(scenario_runs, function(x) x$summary)
)

sens8_df <- bind_rows(
  lapply(scenario_runs, function(x) x$curve)
)

baseline_cases <- results_df$Total_New_Cases[
  results_df$Reduction_Factor == 1.0
]

if (length(baseline_cases) != 1 || baseline_cases <= 0) {
  stop(
    "The no-restriction scenario must produce a positive number of cases ",
    "to calculate percentage case reductions."
  )
}

results_df <- results_df %>%
  mutate(
    Cases_Reduced = baseline_cases - Total_New_Cases,

    Percentage_Case_Reduction = (
      (baseline_cases - Total_New_Cases) /
        baseline_cases
    ) * 100
  )

pre_restriction_df <- sens8_df %>%
  filter(time < t_star)

pre_restriction_check <- pre_restriction_df %>%
  group_by(time) %>%
  summarise(
    Minimum_Infectious = min(Total_Infectious),
    Maximum_Infectious = max(Total_Infectious),
    Difference = Maximum_Infectious - Minimum_Infectious,
    .groups = "drop"
  )

maximum_pre_restriction_difference <- max(
  pre_restriction_check$Difference,
  na.rm = TRUE
)

cat("\nPRE-RESTRICTION CHECK\n")
cat("---------------------\n")
cat(
  "Maximum difference in total infectious students before day",
  t_star, ":",
  maximum_pre_restriction_difference, "\n"
)

if (maximum_pre_restriction_difference < 1e-6) {
  cat("CHECK PASSED: Scenarios are identical before day 14.\n")
} else {
  cat(
    "WARNING: Scenarios differ before day 14. ",
    "Check the model implementation.\n"
  )
}

results_table <- results_df %>%
  select(
    Scenario,
    Reduction_Factor,
    Alpha_After_Day_14,
    Peak_Infectious_Burden,
    Peak_Day,
    Total_New_Cases,
    County_Attack_Rate_Percent,
    Cases_Reduced,
    Percentage_Case_Reduction
  ) %>%
  mutate(
    Peak_Infectious_Burden = round(Peak_Infectious_Burden, 2),
    Peak_Day = round(Peak_Day, 0),
    Total_New_Cases = round(Total_New_Cases, 0),
    County_Attack_Rate_Percent = round(
      County_Attack_Rate_Percent, 3
    ),
    Cases_Reduced = round(Cases_Reduced, 0),
    Percentage_Case_Reduction = round(
      Percentage_Case_Reduction, 2
    )
  )

cat("\nSENSITIVITY ANALYSIS RESULTS\n")
cat("============================\n")

print(results_table, n = Inf)

elasticity_df <- tibble(
  From_Factor = c(1.0, 0.5),
  To_Factor = c(0.5, 0.1)
) %>%
  rowwise() %>%
  mutate(
    From_Cases = results_df$Total_New_Cases[
      match(From_Factor, results_df$Reduction_Factor)
    ],

    To_Cases = results_df$Total_New_Cases[
      match(To_Factor, results_df$Reduction_Factor)
    ],

    Percentage_Change_Cases = ifelse(
      From_Cases != 0,
      ((To_Cases - From_Cases) / From_Cases) * 100,
      NA_real_
    ),

    Percentage_Change_Connectivity = (
      (To_Factor - From_Factor) / From_Factor
    ) * 100,

    Elasticity = ifelse(
      Percentage_Change_Connectivity != 0,
      Percentage_Change_Cases /
        Percentage_Change_Connectivity,
      NA_real_
    )
  ) %>%
  ungroup()

cat("\nELASTICITY BETWEEN CONNECTIVITY SCENARIOS\n")
cat("=========================================\n")

print(
  elasticity_df %>%
    mutate(
      Percentage_Change_Cases = round(
        Percentage_Change_Cases, 2
      ),
      Percentage_Change_Connectivity = round(
        Percentage_Change_Connectivity, 2
      ),
      Elasticity = round(Elasticity, 4)
    ),
  n = Inf
)

p_infectious <- ggplot(
  sens8_df,
  aes(
    x = time,
    y = Total_Infectious,
    color = Scenario
  )
) +
  geom_line(linewidth = 1.2) +

  geom_vline(
    xintercept = t_star,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +

  annotate(
    "text",
    x = t_star + 3,
    y = max(sens8_df$Total_Infectious, na.rm = TRUE) * 0.95,
    label = paste0("Restrictions begin (day ", t_star, ")"),
    hjust = 0,
    color = "black"
  ) +

  scale_color_manual(
    values = c(
      "No reduction" = "#0072B2",
      "50% reduction" = "#D55E00",
      "90% reduction" = "#009E73"
    ),
    breaks = c(
      "No reduction",
      "50% reduction",
      "90% reduction"
    )
  ) +

  labs(
    title = "Effect of Inter-School Travel Restrictions on Measles Transmission",
    subtitle = paste0("Restrictions begin on day ", t_star),
    x = "Time (days)",
    y = "Total infectious students",
    color = "Connectivity scenario"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(p_infectious)

ggsave(
  filename = "figures/sensitivity/sensitivity_analysis_3_infectious_curves.png",
  plot = p_infectious,
  width = 10,
  height = 6,
  dpi = 300
)

p_cases <- ggplot(
  results_df,
  aes(
    x = factor(
      Scenario,
      levels = c(
        "No reduction",
        "50% reduction",
        "90% reduction"
      )
    ),
    y = Total_New_Cases
  )
) +
  geom_col() +
  geom_text(
    aes(label = round(Total_New_Cases)),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "Total New Measles Cases Under Different Connectivity Scenarios",
    x = "Inter-school connectivity scenario",
    y = "Total new cases"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

print(p_cases)

p_reduction <- ggplot(
  results_df,
  aes(
    x = factor(
      Scenario,
      levels = c(
        "No reduction",
        "50% reduction",
        "90% reduction"
      )
    ),
    y = Percentage_Case_Reduction
  )
) +
  geom_col() +
  geom_text(
    aes(label = paste0(
      round(Percentage_Case_Reduction, 2),
      "%"
    )),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "Percentage Reduction in Cases Compared with No Restrictions",
    x = "Inter-school connectivity scenario",
    y = "Reduction in total new cases (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

print(p_reduction)

write.csv(
  results_df,
  "output/tables/sensitivity_analysis_3_results.csv",
  row.names = FALSE
)

write.csv(
  sens8_df,
  "output/tables/sensitivity_analysis_3_infectious_curves.csv",
  row.names = FALSE
)

write.csv(
  elasticity_df,
  "output/tables/sensitivity_analysis_3_elasticity.csv",
  row.names = FALSE
)

write.csv(
  pre_restriction_check,
  "output/tables/sensitivity_analysis_3_pre_restriction_check.csv",
  row.names = FALSE
)

ggsave(
  "figures/sensitivity/sensitivity_analysis_3_infectious_curves.png",
  plot = p_infectious,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/sensitivity/sensitivity_analysis_3_total_cases.png",
  plot = p_cases,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/sensitivity/sensitivity_analysis_3_case_reduction.png",
  plot = p_reduction,
  width = 9,
  height = 6,
  dpi = 300
)

cat("\nANALYSIS COMPLETE\n")
cat("=================\n")
cat("The three connectivity scenarios have been simulated.\n")
cat("The infection starts in the first school.\n")
cat("The adjacency matrix was normalized by its spectral radius.\n")
cat("Connectivity restrictions begin on day 14.\n")
cat(
  "Total new cases were calculated using susceptible depletion:\n",
  "C(T) = sum_i [S_i(0) - S_i(T)].\n"
)
cat("Results tables and plots have been saved to the working directory.\n")
