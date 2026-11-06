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

df <- read.csv("data/processed/school_data_joined.csv") %>%
  rename_with(tolower)

if (!"unit_id" %in% names(df)) {
  df <- df %>% mutate(unit_id = as.character(row_number()))
} else {
  df <- df %>% mutate(unit_id = as.character(unit_id))
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
cat("\nNumber of schools:", n_schools, "\n")

net_df <- read.csv("data/processed/network_df.csv") %>%
  rename_with(tolower)

id_cols <- names(net_df)[grep("unit|from|to|source|target|id", names(net_df))]

if (length(id_cols) >= 2) {
  net_df <- net_df %>%
    rename(
      unit_id_a = !!sym(id_cols[1]),
      unit_id_b = !!sym(id_cols[2])
    )
} else {
  colnames(net_df)[1:2] <- c("unit_id_a", "unit_id_b")
}

net_df <- net_df %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  )

id_map <- setNames(seq_len(n_schools), df$unit_id)
A <- matrix(0, nrow = n_schools, ncol = n_schools)

for (k in seq_len(nrow(net_df))) {
  i <- unname(id_map[net_df$unit_id_a[k]])
  j <- unname(id_map[net_df$unit_id_b[k]])

  if (length(i) == 1 && length(j) == 1 && !is.na(i) && !is.na(j)) {
    w <- if ("total_contacts" %in% names(net_df)) net_df$total_contacts[k] else 1
    if (is.na(w) || !is.finite(w) || w < 0) w <- 0
    A[i, j] <- w
    A[j, i] <- w
  }
}

A_raw <- A
spectral_radius <- max(Mod(eigen(A_raw, only.values = TRUE)$values))

if (!is.finite(spectral_radius) || spectral_radius <= 0) {
  stop("The network spectral radius is zero or invalid. Check contact weights.")
}

A <- A_raw / spectral_radius

gamma_rate <- 1 / 7
sigma_rate <- 1 / 10
beta_rate  <- 15 * gamma_rate
alpha      <- 0.05

times <- seq(0, 200, by = 1)
N     <- df$total_students

params <- list(
  beta = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N = N,
  A = A,
  n = n_schools
)

seir_model <- function(time, state, parameters) {
  with(as.list(parameters), {
    S <- state[1:n]
    E <- state[(n + 1):(2 * n)]
    I <- state[(2 * n + 1):(3 * n)]
    R <- state[(3 * n + 1):(4 * n)]

    I_prop <- I / N

    lambda <- (beta * I_prop) + (alpha * beta * as.vector(A %*% I_prop))

    dS <- -lambda * S
    dE <- (lambda * S) - (sigma * E)
    dI <- (sigma * E) - (gamma * I)
    dR <- gamma * I

    list(c(dS, dE, dI, dR))
  })
}

run_vax_scenario <- function(vax_vector, seed_id) {
  R0_vec <- round(N * vax_vector)
  I0_vec <- rep(0, n_schools)
  I0_vec[seed_id] <- 1
  E0_vec <- rep(0, n_schools)
  S0_vec <- pmax(0, N - I0_vec - E0_vec - R0_vec)

  state_vec <- c(S0_vec, E0_vec, I0_vec, R0_vec)
  names(state_vec) <- c(
    paste0("S_", seq_len(n_schools)),
    paste0("E_", seq_len(n_schools)),
    paste0("I_", seq_len(n_schools)),
    paste0("R_", seq_len(n_schools))
  )

  out_df <- as.data.frame(
    ode(y = state_vec, times = times, func = seir_model, parms = params)
  )
  return(out_df)
}

get_total_I <- function(out_df, scenario_label, seed_type, seed_name) {
  out_df %>%
    pivot_longer(cols = -time, names_to = "var", values_to = "Count") %>%
    filter(str_detect(var, "^I_")) %>%
    group_by(time) %>%
    summarise(total_infectious = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Scenario = scenario_label,
      `Seed Type` = seed_type,
      `Seed Name` = seed_name
    )
}

calculate_outbreak_measures <- function(out_df) {
  S_columns <- grep("^S_", names(out_df), value = TRUE)
  I_columns <- grep("^I_", names(out_df), value = TRUE)

  S_initial <- sapply(S_columns, function(x) out_df[[x]][1])
  S_final   <- sapply(S_columns, function(x) out_df[[x]][nrow(out_df)])

  total_cases <- sum(S_initial - S_final, na.rm = TRUE)
  total_I     <- rowSums(out_df[, I_columns, drop = FALSE], na.rm = TRUE)

  peak_burden <- max(total_I, na.rm = TRUE)
  peak_day    <- out_df$time[which.max(total_I)]
  attack_rate <- (total_cases / sum(N, na.rm = TRUE)) * 100

  return(c(
    Total_Cases = total_cases,
    Peak_Burden = peak_burden,
    Peak_Day    = peak_day,
    Attack_Rate = attack_rate
  ))
}

net_graph <- graph_from_adjacency_matrix(
  A, mode = "undirected", weighted = TRUE, diag = FALSE
)

df$node_strength <- strength(net_graph, weights = E(net_graph)$weight)
df$betweenness   <- betweenness(net_graph, weights = NA)

sorted_df <- df %>% arrange(desc(node_strength))

hubs_df <- head(sorted_df, 3) %>% mutate(type = "Network Hub")
periph_df <- tail(sorted_df, 3) %>% arrange(node_strength) %>% mutate(type = "Peripheral")

target_schools <- bind_rows(hubs_df, periph_df)

cat("\n========================================\n")
cat("SELECTED SEEDING SCHOOLS (3 HUBS & 3 PERIPHERALS)\n")
cat("========================================\n")
print(target_schools %>% select(school_id, school_name, node_strength, type))

scenarios_list <- list(
  list(name = "Scenario 1: Observed Coverage", set_vax = NULL),
  list(name = "Scenario 2: Seed Schools at 74%", set_vax = 0.74),
  list(name = "Scenario 3: Seed Schools at 91.6%", set_vax = 0.916)
)

results_table_list <- list()
time_series_list    <- list()

for (sc in scenarios_list) {
  sc_name <- sc$name
  fixed_vax <- sc$set_vax

  for (idx in seq_len(nrow(target_schools))) {
    s_id   <- target_schools$school_id[idx]
    s_name <- target_schools$school_name[idx]
    s_type <- target_schools$type[idx]

    current_vax <- df$vax_rate
    if (!is.null(fixed_vax)) {

      current_vax[target_schools$school_id] <- fixed_vax
    }

    out_df <- run_vax_scenario(current_vax, s_id)

    ts_data <- get_total_I(out_df, sc_name, s_type, s_name)
    time_series_list[[length(time_series_list) + 1]] <- ts_data

    metrics <- calculate_outbreak_measures(out_df)

    results_table_list[[length(results_table_list) + 1]] <- tibble(
      Scenario = sc_name,
      `Seed Location` = s_name,
      `Seed Type` = s_type,
      `Vaccination Coverage` = paste0(round(current_vax[s_id] * 100, 1), "%"),
      `Total New Cases (C)` = round(metrics["Total_Cases"]),
      `Peak Burden (Ipeak)` = round(metrics["Peak_Burden"]),
      `Peak Day` = round(metrics["Peak_Day"]),
      `Attack Rate (AR)` = paste0(round(metrics["Attack_Rate"], 2), "%")
    )
  }
}

comparison_all <- bind_rows(time_series_list)
table_4_3      <- bind_rows(results_table_list)

summary_ts <- comparison_all %>%
  group_by(Scenario, `Seed Type`, time) %>%
  summarise(
    mean_I = mean(total_infectious),
    min_I  = min(total_infectious),
    max_I  = max(total_infectious),
    .groups = "drop"
  )

p_comparison <- ggplot(
  summary_ts,
  aes(x = time, y = mean_I, color = `Seed Type`, fill = `Seed Type`)
) +
  geom_ribbon(aes(ymin = min_I, ymax = max_I), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Scenario, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("Network Hub" = "firebrick", "Peripheral" = "darkblue")) +
  scale_fill_manual(values = c("Network Hub" = "firebrick", "Peripheral" = "darkblue")) +
  labs(
    title = "Sensitivity Analysis: Multi-Hub (n=3) vs Multi-Peripheral (n=3) Seeding",
    subtitle = "Solid lines show mean trajectories across seed choices; shaded regions represent range (min-max)",
    x = "Days Since Outbreak Introduction",
    y = "Number of Infectious Cases",
    color = "Seed Category",
    fill = "Seed Category"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10)
  )

print(p_comparison)

cat("\n\n============================================================\n")
cat("TABLE 4.3: SENSITIVITY ANALYSIS OF 3 HUBS VS 3 PERIPHERALS\n")
cat("============================================================\n\n")

print(table_4_3, n = Inf)

if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(table_4_3, caption = "3 Hubs vs 3 Peripherals Sensitivity Analysis"))
}

write.csv(table_4_3, "output/tables/table_4_3_multi_hub_vs_peripheral.csv", row.names = FALSE)
cat("\nTable saved as 'output/tables/table_4_3_multi_hub_vs_peripheral.csv'.\n")
