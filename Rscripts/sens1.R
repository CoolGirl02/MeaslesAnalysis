library(tidyverse)
library(deSolve)

# ==========================================
# 1. LOAD & PREPARE SCHOOL DATA
# ==========================================
df <- read.csv("data/cleanned/school_data_joined.csv")
names(df) <- tolower(names(df))

if (!"unit_id" %in% names(df)) {
  df$unit_id <- as.character(1:nrow(df))
} else {
  df$unit_id <- as.character(df$unit_id)
}

df <- df %>%
  mutate(
    vax_rate = if_else(percent_immunized > 1, percent_immunized / 100, percent_immunized),
    school_id = row_number(),
    page = ceiling(school_id / 6)
  )

n_schools <- nrow(df)

# ==========================================
# 2. BUILD NETWORK ADJACENCY MATRIX (A)
# ==========================================
net_df <- read.csv("data/cleanned/network_df.csv")
names(net_df) <- tolower(names(net_df))

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

id_map <- setNames(1:n_schools, df$unit_id)
A      <- matrix(0, nrow = n_schools, ncol = n_schools)

for (k in 1:nrow(net_df)) {
  i <- id_map[net_df$unit_id_a[k]]
  j <- id_map[net_df$unit_id_b[k]]
  if (!is.na(i) && !is.na(j)) {
    w <- if("connect_risk" %in% names(net_df)) net_df$connect_risk[k] else 1
    A[i, j] <- w
    A[j, i] <- w
  }
}

row_sums <- rowSums(A)
row_sums[row_sums == 0] <- 1
A <- A / row_sums  

# ==========================================
# 3. BASELINE MODEL PARAMETERS
# ==========================================
gamma_rate <- 1 / 7                    # Infectious period (~7 days)
sigma_rate <- 1 / 10                   # Incubation period (~10 days)
beta_rate  <- 15 * gamma_rate          # Intra-school transmission (R0 = 15)
alpha      <- 0.05                     # Inter-school coupling coefficient (5%)
times      <- seq(0, 200, by = 1)

N <- df$total_students
R0_count <- round(N * df$vax_rate)
I0_count <- c(1, rep(0, n_schools - 1)) # Index case in School 1
E0_count <- rep(0, n_schools)
S0_count <- pmax(0, N - I0_count - E0_count - R0_count)

initial_state <- c(S0_count, E0_count, I0_count, R0_count)
names(initial_state) <- c(
  paste0("S_", 1:n_schools),
  paste0("E_", 1:n_schools),
  paste0("I_", 1:n_schools),
  paste0("R_", 1:n_schools)
)

params <- list(
  beta  = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N     = N,
  A     = A,
  n     = n_schools
)

# ==========================================
# 4. SEIR NETWORK ODE FUNCTION
# ==========================================
seir_model <- function(time, state, parameters) {
  with(as.list(parameters), {
    S <- state[1:n]
    E <- state[(n + 1):(2 * n)]
    I <- state[(2 * n + 1):(3 * n)]
    R <- state[(3 * n + 1):(4 * n)]
    
    I_prop <- I / N
    # Force of infection combining Intra- and Inter-school rates
    lambda <- (beta * I_prop) + (alpha * beta * as.vector(A %*% I_prop))
    
    dS <- - lambda * S
    dE <-   lambda * S - (sigma * E)
    dI <-  (sigma * E) - (gamma * I)
    dR <-  (gamma * I)
    
    return(list(c(dS, dE, dI, dR)))
  })
}

# ==========================================
# 5. BASELINE SIMULATION & PLOTS
# ==========================================
output <- ode(y = initial_state, times = times, func = seir_model, parms = params)

results_long <- as.data.frame(output) %>%
  pivot_longer(cols = -time, names_to = "var", values_to = "Count") %>%
  mutate(
    Compartment = str_extract(var, "^[SEIR]"),
    school_id   = as.numeric(str_extract(var, "\\d+"))
  ) %>%
  left_join(df, by = "school_id") %>%
  mutate(school_label = paste0(school_name, "\n(Vax: ", round(vax_rate * 100, 1), "%, Pop: ", total_students, ")"))

total_pages <- max(df$page)

for (p in 1:total_pages) {
  page_data <- results_long %>% filter(page == p)
  
  p_plot <- ggplot(page_data, aes(x = time, y = Count, color = Compartment)) +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = c("S" = "blue", "E" = "orange", "I" = "red", "R" = "green"),
      labels = c("S" = "Susceptible", "E" = "Exposed", "I" = "Infectious", "R" = "Recovered/Vaccinated")
    ) +
    facet_wrap(~school_label, ncol = 3, scales = "free_y") +
    labs(
      title = paste0("Baseline SEIR Measles Simulations (Page ", p, " of ", total_pages, ")"),
      x = "Days", y = "Number of Students", color = "Compartment"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(size = 8, face = "bold"), legend.position = "bottom")
  
  print(p_plot)
}


# ==============================================================================
# 6. SENSITIVITY ANALYSIS: VACCINATION COVERAGE vs. SCHOOL SIZE
# ==============================================================================
# Here we systematically test different vaccination coverage levels (70% - 98%)
# while keeping full network intra + inter force of infection active.

# A. Categorize schools by student enrollment size
df_sens <- df %>%
  mutate(
    school_size_cat = case_when(
      total_students < 400 ~ "Small (<400)",
      total_students >= 400 & total_students <= 800 ~ "Medium (400-800)",
      total_students > 800 ~ "Large (>800)"
    )
  )

# B. Define sensitivity coverage grid
vax_scenarios <- c(0.70, 0.80, 0.85, 0.90, 0.95, 0.98)
sens_results <- data.frame()

# C. Loop over vaccination scenarios
for (v_cov in vax_scenarios) {
  
  # Set uniform initial vaccination across all schools
  R0_sens <- round(N * v_cov)
  I0_sens <- c(1, rep(0, n_schools - 1)) # Seed in School 1
  E0_sens <- rep(0, n_schools)
  S0_sens <- pmax(0, N - I0_sens - E0_sens - R0_sens)
  
  state_sens <- c(S0_sens, E0_sens, I0_sens, R0_sens)
  names(state_sens) <- names(initial_state)
  
  # Run ODE solver for this scenario
  out_sens <- ode(y = state_sens, times = times, func = seir_model, parms = params)
  
  # Process ODE output to calculate attack rate per school
  res_sens <- as.data.frame(out_sens) %>%
    pivot_longer(cols = -time, names_to = "var", values_to = "Count") %>%
    mutate(
      Compartment = str_extract(var, "^[SEIR]"),
      school_id   = as.numeric(str_extract(var, "\\d+"))
    ) %>%
    filter(time == max(times)) %>% # Take final day to get total recovered/infected
    pivot_wider(names_from = Compartment, values_from = Count) %>%
    left_join(df_sens, by = "school_id") %>%
    mutate(
      total_cases = R - (total_students * v_cov), # Subtract initial vaccinated
      attack_rate = (pmax(0, total_cases) / total_students) * 100,
      scenario_vax = paste0(v_cov * 100, "%")
    )
  
  sens_results <- rbind(sens_results, res_sens)
}

# D. Plot Sensitivity 1: Attack Rate vs. Vaccination Coverage grouped by School Size
sens_summary <- sens_results %>%
  group_by(scenario_vax, school_size_cat) %>%
  summarise(
    avg_attack_rate = mean(attack_rate, na.rm = TRUE),
    total_cases_county = sum(pmax(0, total_cases), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(vax_num = as.numeric(gsub("%", "", scenario_vax)))

ggplot(sens_summary, aes(x = vax_num, y = avg_attack_rate, color = school_size_cat, group = school_size_cat)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "red") +
  annotate("text", x = 95.5, y = max(sens_summary$avg_attack_rate) * 0.9, 
           label = "95% Herd Immunity Target", color = "red", angle = 90) +
  labs(
    title = "Sensitivity Analysis 1: Impact of Vaccination Coverage on Attack Rate",
    subtitle = "Evaluated across School Size Categories with Active Network Transmission",
    x = "County-Wide Vaccination Coverage (%)",
    y = "Average Final Attack Rate (%)",
    color = "School Size Category"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12))
















library(igraph)

# ==============================================================================
# 7. SENSITIVITY ANALYSIS 2: ROLE OF HIGH-CENTRALITY (HUB) SCHOOLS
# ==============================================================================

# A. Build graph and calculate Node Strength centrality
net_graph <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)

df$node_strength <- strength(net_graph)
df$betweenness   <- betweenness(net_graph, weights = NA)

# B. Identify Top Hub and Peripheral Schools
top_hub_id   <- df$school_id[which.max(df$node_strength)]
low_hub_id   <- df$school_id[which.min(df$node_strength)]
top_hub_name <- df$school_name[df$school_id == top_hub_id]
low_hub_name <- df$school_name[df$school_id == low_hub_id]

cat("Top Hub School:", top_hub_name, "| Vax Rate:", round(df$vax_rate[df$school_id == top_hub_id]*100, 1), "%\n")
cat("Peripheral School:", low_hub_name, "| Vax Rate:", round(df$vax_rate[df$school_id == low_hub_id]*100, 1), "%\n")

# C. Define Seeding Scenarios using ACTUAL dataset vaccination coverage (df$vax_rate)
R0_actual <- round(N * df$vax_rate)

# Scenario A: Seed index case (I0 = 1) in Top Centrality Hub School
I0_hub <- rep(0, n_schools)
I0_hub[top_hub_id] <- 1
S0_hub <- pmax(0, N - I0_hub - R0_actual)

state_hub <- c(S0_hub, rep(0, n_schools), I0_hub, R0_actual)
names(state_hub) <- names(initial_state)

# Scenario B: Seed index case (I0 = 1) in Peripheral School
I0_peri <- rep(0, n_schools)
I0_peri[low_hub_id] <- 1
S0_peri <- pmax(0, N - I0_peri - R0_actual)

state_peri <- c(S0_peri, rep(0, n_schools), I0_peri, R0_actual)
names(state_peri) <- names(initial_state)

# D. Run SEIR Network Simulations
out_hub  <- as.data.frame(ode(y = state_hub,  times = times, func = seir_model, parms = params))
out_peri <- as.data.frame(ode(y = state_peri, times = times, func = seir_model, parms = params))

# E. Extract Total Active Infectious Cases County-Wide
get_total_I <- function(out_df, scenario_label) {
  out_df %>%
    pivot_longer(cols = -time, names_to = "var", values_to = "Count") %>%
    filter(str_detect(var, "^I_")) %>%
    group_by(time) %>%
    summarise(total_infectious = sum(Count), .groups = "drop") %>%
    mutate(Scenario = scenario_label)
}

hub_curve  <- get_total_I(out_hub,  paste0("Hub : ", top_hub_name))
peri_curve <- get_total_I(out_peri, paste0("Peripheral : ", low_hub_name))

comparison_df <- rbind(hub_curve, peri_curve)

# F. Plot Epidemic Curves
ggplot(comparison_df, aes(x = time, y = total_infectious, color = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("firebrick", "darkblue")) +
  labs(
    title = "Sensitivity Analysis 2: Outbreak Dynamics by Index Case Location",
    subtitle = "Comparing Index Case in Top Centrality Hub vs. Peripheral School (Actual Vax Rates)",
    x = "Days Since Outbreak Introduction",
    y = "Total Active Infectious Cases in County",
    color = "Initial infectious Location"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )







# Controlled Test: Set BOTH seed schools to 74% coverage to remove the vax rate confounder
fixed_vax_rate <- 0.74 

# Base vaccination array from dataset
R0_controlled <- round(N * df$vax_rate)

# Force BOTH the Hub and Peripheral seed schools to have the exact same 74% coverage
R0_controlled[top_hub_id] <- round(N[top_hub_id] * fixed_vax_rate)
R0_controlled[low_hub_id] <- round(N[low_hub_id] * fixed_vax_rate)

# Scenario A: Seed in Hub (with 74% vax)
I0_hub <- rep(0, n_schools); I0_hub[top_hub_id] <- 1
S0_hub <- pmax(0, N - I0_hub - R0_controlled)
state_hub <- c(S0_hub, rep(0, n_schools), I0_hub, R0_controlled)
names(state_hub) <- names(initial_state)

# Scenario B: Seed in Peripheral (also with 74% vax)
I0_peri <- rep(0, n_schools); I0_peri[low_hub_id] <- 1
S0_peri <- pmax(0, N - I0_peri - R0_controlled)
state_peri <- c(S0_peri, rep(0, n_schools), I0_peri, R0_controlled)
names(state_peri) <- names(initial_state)

# Run Simulations
out_hub_ctrl  <- as.data.frame(ode(y = state_hub,  times = times, func = seir_model, parms = params))
out_peri_ctrl <- as.data.frame(ode(y = state_peri, times = times, func = seir_model, parms = params))

# Extract & Plot
hub_curve_ctrl  <- get_total_I(out_hub_ctrl,  paste0("Hub Seed (74% Vax): ", top_hub_name))
peri_curve_ctrl <- get_total_I(out_peri_ctrl, paste0("Peripheral Seed (74% Vax): ", low_hub_name))

comparison_ctrl_df <- rbind(hub_curve_ctrl, peri_curve_ctrl)

ggplot(comparison_ctrl_df, aes(x = time, y = total_infectious, color = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("firebrick", "darkblue")) +
  labs(
    title = "Controlled Sensitivity Analysis 2: Decoupled Network Effect",
    subtitle = "Both Seed Schools Controlled at 74% Vaccination Coverage",
    x = "Days Since Outbreak Introduction",
    y = "Total Active Infectious Cases in County"
  ) +
  theme_minimal()






cat("================ CONNECTIVITY FIGURES ================\n")
cat("HUB SCHOOL:", top_hub_name, "\n")
cat("  - Node Strength (Risk Score):", df$node_strength[top_hub_id], "\n")
cat("  - Betweenness Centrality:   ", df$betweenness[top_hub_id], "\n")
cat("  - Number of Connected Schools:", sum(A[top_hub_id, ] > 0), "\n\n")

cat("PERIPHERAL SCHOOL:", low_hub_name, "\n")
cat("  - Node Strength (Risk Score):", df$node_strength[low_hub_id], "\n")
cat("  - Betweenness Centrality:   ", df$betweenness[low_hub_id], "\n")
cat("  - Number of Connected Schools:", sum(A[low_hub_id, ] > 0), "\n")
cat("======================================================\n")


# Check St. Paul's non-zero connections
st_paul_edges <- data.frame(
  Neighbor_School = df$school_name,
  Weight = A[low_hub_id, ]
) %>% filter(Weight > 0)

print(st_paul_edges)





# ==============================================================================
# CONNECTIVITY RISK ANALYSIS FOR PERIPHERAL NODES
# ==============================================================================

# 1. Connectivity Details for St. Paul (or current min-strength school)
peripheral_details <- data.frame(
  School_ID      = low_hub_id,
  School_Name    = low_hub_name,
  Node_Strength  = df$node_strength[low_hub_id],
  Betweenness    = df$betweenness[low_hub_id],
  Degree_Count   = sum(A[low_hub_id, ] > 0),
  Max_Edge_Weight= ifelse(sum(A[low_hub_id, ] > 0) > 0, max(A[low_hub_id, ]), 0)
)

print("--- PERIPHERAL NODE SUMMARY ---")
print(peripheral_details)

# 2. View the Bottom 10 Schools Ranked by Connectivity Risk (Node Strength)
bottom_10_connectivity <- df %>%
  select(school_id, school_name, node_strength, vax_rate) %>%
  arrange(node_strength) %>%
  head(10)

print("--- BOTTOM 10 LOWEST CONNECTIVITY SCHOOLS ---")
print(bottom_10_connectivity)

# 3. Check if St. Paul has ANY outgoing or incoming edge connections in raw matrix A
raw_outgoing <- sum(A[low_hub_id, ] > 0)
raw_incoming <- sum(A[, low_hub_id] > 0)

cat("\nOutgoing connections from peripheral node:", raw_outgoing, "\n")
cat("Incoming connections to peripheral node:  ", raw_incoming, "\n")








# ==============================================================================
# 9. SENSITIVITY ANALYSIS 3: REDUCING INTER-SCHOOL CONNECTIVITY
# ==============================================================================

# A. Define Coupling Coefficients (Alpha values representing connectivity levels)
# Baseline alpha is 0.05 (100% mobility)
alpha_baseline <- 0.05
alpha_scenarios <- c(
  "Baseline (0% Reduction)" = alpha_baseline,
  "Moderate Restrictions (50% Reduction)" = alpha_baseline * 0.50,
  "Full School Isolation (100% Reduction)" = 0.0
)

# B. Set up Seeding in High-Centrality Hub School (using actual vax coverage)
R0_base <- round(N * df$vax_rate)
I0_hub  <- rep(0, n_schools)
I0_hub[top_hub_id] <- 1
S0_hub  <- pmax(0, N - I0_hub - R0_base)

state_init <- c(S0_hub, rep(0, n_schools), I0_hub, R0_base)
names(state_init) <- names(initial_state)

# C. Loop through Connectivity Reduction Scenarios
results_list <- list()

for (scen_name in names(alpha_scenarios)) {
  
  # Copy baseline parameters and adjust alpha
  params_scen <- params
  params_scen$alpha <- alpha_scenarios[[scen_name]]
  
  # Run ODE Simulation
  out_scen <- as.data.frame(ode(y = state_init, times = times, func = seir_model, parms = params_scen))
  
  # Calculate Total County-Wide Infectious Cases
  df_I <- get_total_I(out_scen, scen_name)
  results_list[[scen_name]] <- df_I
}

# Combine all scenario results into one dataframe
sens3_df <- bind_rows(results_list)

# D. Plot Comparative Epidemic Curves for Connectivity Reduction
ggplot(sens3_df, aes(x = time, y = total_infectious, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c(
    "Baseline (0% Reduction)" = "firebrick",
    "Moderate Restrictions (50% Reduction)" = "darkorange",
    "Full School Isolation (100% Reduction)" = "forestgreen"
  )) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(
    title = "Sensitivity Analysis 3: Impact of Reducing Inter-School Connectivity",
    subtitle = "Simulating Outbreak Response Strategies Starting from High-Centrality Hub",
    x = "Days Since Outbreak Introduction",
    y = "Total Active Infectious Cases in County",
    color = "Inter-School Mobility Policy",
    linetype = "Inter-School Mobility Policy"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )











# ==============================================================================
# 10. SENSITIVITY ANALYSIS 4: TARGETED VS. RANDOM VACCINATION
# ==============================================================================

set.seed(123) # For reproducible random selection

# Identify Top 5 High-Centrality Hub Schools
top_5_hubs <- df %>%
  arrange(desc(node_strength)) %>%
  slice(1:5) %>%
  pull(school_id)

# Select 5 Random Non-Hub Schools for Comparison
random_5_schools <- df %>%
  filter(!school_id %in% top_5_hubs) %>%
  sample_n(5) %>%
  pull(school_id)

# A. Baseline Setup (Actual Dataset Coverage)
vax_base <- df$vax_rate

# B. Scenario 1: Targeted Vaccination (Boost Top 5 Hubs to 95% Coverage)
vax_targeted <- vax_base
vax_targeted[top_5_hubs] <- pmax(vax_targeted[top_5_hubs], 0.95)

# C. Scenario 2: Random Vaccination (Boost 5 Random Schools to 95% Coverage)
vax_random <- vax_base
vax_random[random_5_schools] <- pmax(vax_random[random_5_schools], 0.95)

# Function to run ODE and get epidemic curve seeded at top hub
run_vax_scenario <- function(vax_vector, label) {
  R0_vec <- round(N * vax_vector)
  I0_vec <- rep(0, n_schools)
  I0_vec[top_hub_id] <- 1
  S0_vec <- pmax(0, N - I0_vec - R0_vec)
  
  state_vec <- c(S0_vec, rep(0, n_schools), I0_vec, R0_vec)
  names(state_vec) <- names(initial_state)
  
  out_df <- as.data.frame(ode(y = state_vec, times = times, func = seir_model, parms = params))
  get_total_I(out_df, label)
}

# D. Run All Three Scenarios
curve_baseline <- run_vax_scenario(vax_base,     "Baseline (Actual Coverage)")
curve_targeted <- run_vax_scenario(vax_targeted, "Targeted (Top 5 Hubs Raised to 95%)")
curve_random   <- run_vax_scenario(vax_random,   "Random (5 Random Schools Raised to 95%)")

sens4_df <- rbind(curve_baseline, curve_targeted, curve_random)

# E. Plot Comparative Epidemic Curves
ggplot(sens4_df, aes(x = time, y = total_infectious, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c(
    "Baseline (Actual Coverage)"              = "firebrick",
    "Random (5 Random Schools Raised to 95%)" = "darkorange",
    "Targeted (Top 5 Hubs Raised to 95%)"    = "forestgreen"
  )) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(
    title = "Sensitivity Analysis 4: Targeted vs. Random Vaccine Allocation",
    subtitle = "Comparing County-Wide Outbreak Impact of Protecting Network Hubs",
    x = "Days Since Outbreak Introduction",
    y = "Total Active Infectious Cases in County",
    color = "Allocation Strategy",
    linetype = "Allocation Strategy"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )














# ==============================================================================
# 11. SENSITIVITY ANALYSIS 6: VACCINATION COVERAGE BY SCHOOL SIZE
# (Maintaining Identical County-Wide Overall Vaccination Coverage)
# ==============================================================================

# A. Calculate Baseline County Average Vaccination Rate
target_county_mean <- mean(df$vax_rate)

# B. Classify Schools into Small (<400) and Large (>=400)
small_idx <- which(df$total_students < 400)
large_idx <- which(df$total_students >= 400)

n_small <- length(small_idx)
n_large <- length(large_idx)

# --- SCENARIO A: Low Vax in Small Schools, High Vax in Large Schools ---
# Set initial targets: 70% for small schools, 95% for large schools
vax_small_low <- numeric(n_schools)
vax_small_low[small_idx] <- 0.70
vax_small_low[large_idx] <- 0.95

# Rescale proportionally to match exact county mean coverage
vax_small_low <- vax_small_low * (target_county_mean / mean(vax_small_low))
vax_small_low <- pmin(0.98, pmax(0.40, vax_small_low))

# --- SCENARIO B: High Vax in Small Schools, Low Vax in Large Schools ---
# Set initial targets: 95% for small schools, 70% for large schools
vax_large_low <- numeric(n_schools)
vax_large_low[small_idx] <- 0.95
vax_large_low[large_idx] <- 0.70

# Rescale proportionally to match exact county mean coverage
vax_large_low <- vax_large_low * (target_county_mean / mean(vax_large_low))
vax_large_low <- pmin(0.98, pmax(0.40, vax_large_low))

# C. Console Verification
cat("=== SENSITIVITY 6 COVERAGE VERIFICATION ===\n")
cat("County Mean Target:                ", round(target_county_mean * 100, 2), "%\n")
cat("Scenario A (Low Small / High Large):", round(mean(vax_small_low) * 100, 2), "%\n")
cat("  - Small Schools Avg Vax:         ", round(mean(vax_small_low[small_idx]) * 100, 2), "%\n")
cat("  - Large Schools Avg Vax:         ", round(mean(vax_small_low[large_idx]) * 100, 2), "%\n")
cat("Scenario B (High Small / Low Large):", round(mean(vax_large_low) * 100, 2), "%\n")
cat("  - Small Schools Avg Vax:         ", round(mean(vax_large_low[small_idx]) * 100, 2), "%\n")
cat("  - Large Schools Avg Vax:         ", round(mean(vax_large_low[large_idx]) * 100, 2), "%\n")
cat("===========================================\n")

# D. Run SEIR Network Simulations (Seeded at Top Hub)
curve_small_low <- run_vax_scenario(vax_small_low, "Low Coverage in Small / High in Large")
curve_large_low <- run_vax_scenario(vax_large_low, "High Coverage in Small / Low in Large")

sens6_df <- rbind(curve_small_low, curve_large_low)

# E. Plot Comparative Epidemic Curves
ggplot(sens6_df, aes(x = time, y = total_infectious, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c(
    "Low Coverage in Small / High in Large" = "darkorange",
    "High Coverage in Small / Low in Large" = "purple"
  )) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    title = "Sensitivity Analysis 6: Impact of Vaccine Distribution by School Size",
    subtitle = "Comparing Low Small/High Large vs. High Small/Low Large (Identical County Average Coverage)",
    x = "Days Since Outbreak Introduction",
    y = "Total Active Infectious Cases in County",
    color = "Vax Distribution Strategy",
    linetype = "Vax Distribution Strategy"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )






























# ==========================================
# HELPER FUNCTIONS FOR INTERVENTION SCENARIOS
# ==========================================

# Standardized function to run events/interventions in deSolve
run_intervention_sim <- function(init_state, parms_list, times_seq, event_df = NULL) {
  if (!is.null(event_df)) {
    out <- ode(y = init_state, times = times_seq, func = seir_model, 
               parms = parms_list, events = list(data = event_df))
  } else {
    out <- ode(y = init_state, times = times_seq, func = seir_model, 
               parms = parms_list)
  }
  return(as.data.frame(out))
}


# ==============================================================================
# SENSITIVITY ANALYSIS 7: REACTIVE SCHOOL CLOSURE (NETWORK EDGE REMOVAL)
# ==============================================================================
# Evaluates dynamic isolation of schools once internal cases breach a threshold.

run_closure_scenario <- function(threshold_cases, label) {
  # Copy parameters
  p_temp <- params
  A_dynamic <- p_temp$A
  
  # Tracking arrays
  times_sim <- seq(0, 200, by = 1)
  curr_state <- initial_state
  
  history_list <- list()
  closed_schools <- c()
  
  for (t_idx in 1:(length(times_sim) - 1)) {
    t_start <- times_sim[t_idx]
    t_end   <- times_sim[t_idx + 1]
    
    # Run 1-day step
    step_out <- as.data.frame(ode(y = curr_state, times = c(t_start, t_end), 
                                  func = seir_model, parms = p_temp))
    
    # Update current state to end of step
    curr_state <- unlist(step_out[2, 2:ncol(step_out)])
    
    # Check infectious counts per school
    I_counts <- curr_state[(2 * n_schools + 1):(3 * n_schools)]
    
    # Identify newly breached schools
    newly_closed <- which(I_counts >= threshold_cases & !(1:n_schools %in% closed_schools))
    
    if (length(newly_closed) > 0 && threshold_cases < Inf) {
      closed_schools <- c(closed_schools, newly_closed)
      for (sch in newly_closed) {
        A_dynamic[sch, ] <- 0
        A_dynamic[, sch] <- 0
      }
      p_temp$A <- A_dynamic
    }
    
    # Store daily snapshot
    history_list[[t_idx]] <- step_out[1, ]
  }
  
  full_out <- bind_rows(history_list)
  get_total_I(full_out, label)
}

curve_no_close   <- run_closure_scenario(Inf, "No Closures (Baseline)")
curve_early_close <- run_closure_scenario(2,   "Early Closure (Trigger: 2 Cases/School)")
curve_late_close  <- run_closure_scenario(10,  "Late Closure (Trigger: 10 Cases/School)")

sens7_df <- rbind(curve_no_close, curve_early_close, curve_late_close)

ggplot(sens7_df, aes(x = time, y = total_infectious, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("firebrick", "forestgreen", "darkorange")) +
  labs(
    title = "Sensitivity Analysis 7: Reactive School Closure Interventions",
    subtitle = "Comparing Dynamic Network Node Isolation at Case Thresholds",
    x = "Days Since Outbreak Introduction", y = "Total Active Infectious Cases",
    color = "Intervention Policy", linetype = "Intervention Policy"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")


# ==============================================================================
# SENSITIVITY ANALYSIS 8: PRIORITY RING VS. MASS REACTION VACCINATION
# ==============================================================================
# Triggers post-exposure vaccination at Day 14 targeting ring neighbors vs large schools.

day_trigger <- 14
vax_boost_amount <- 0.15 # 15% boost to susceptible population

# Base initial conditions seeded at Top Hub
R0_base <- round(N * df$vax_rate)
I0_hub  <- rep(0, n_schools); I0_hub[top_hub_id] <- 1
S0_hub  <- pmax(0, N - I0_hub - R0_base)

# Identify Targets
ring_targets <- as.numeric(neighbors(net_graph, top_hub_id))
size_targets <- df %>% arrange(desc(total_students)) %>% slice(1:length(ring_targets)) %>% pull(school_id)

# Event Generation Function
create_vax_events <- function(target_nodes, trigger_day) {
  events <- data.frame()
  for (node in target_nodes) {
    events <- rbind(events, data.frame(
      var = paste0("S_", node),
      time = trigger_day,
      value = 1 - vax_boost_amount,
      method = "mult"
    ))
    events <- rbind(events, data.frame(
      var = paste0("R_", node),
      time = trigger_day,
      value = vax_boost_amount,
      method = "add"
    ))
  }
  return(events)
}

events_ring <- create_vax_events(ring_targets, day_trigger)
events_size <- create_vax_events(size_targets, day_trigger)

state_init_hub <- c(S0_hub, rep(0, n_schools), I0_hub, R0_base)
names(state_init_hub) <- names(initial_state)

out_base_ring <- ode(y = state_init_hub, times = times, func = seir_model, parms = params)
out_ring      <- ode(y = state_init_hub, times = times, func = seir_model, parms = params, events = list(data = events_ring))
out_size      <- ode(y = state_init_hub, times = times, func = seir_model, parms = params, events = list(data = events_size))

curve_base8 <- get_total_I(as.data.frame(out_base_ring), "Baseline (No Reactive Vax)")
curve_ring  <- get_total_I(as.data.frame(out_ring),      "Ring Strategy (Direct Network Neighbors)")
curve_size  <- get_total_I(as.data.frame(out_size),      "Size Strategy (Largest Schools)")

sens8_df <- rbind(curve_base8, curve_ring, curve_size)

ggplot(sens8_df, aes(x = time, y = total_infectious, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("firebrick", "forestgreen", "darkblue")) +
  geom_vline(xintercept = day_trigger, linetype = "dotted", color = "black") +
  annotate("text", x = day_trigger + 2, y = max(sens8_df$total_infectious) * 0.9, 
           label = "Vaccination Trigger (Day 14)", angle = 90, hjust = 1) +
  labs(
    title = "Sensitivity Analysis 8: Reactive Ring vs. Mass-Size Vaccination",
    subtitle = "Post-Exposure Intervention Triggered at Day 14 (Seeded at Top Hub)",
    x = "Days Since Outbreak Introduction", y = "Total Active Infectious Cases",
    color = "Intervention Policy", linetype = "Intervention Policy"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")


# ==============================================================================
# SENSITIVITY ANALYSIS 9: PHASED INTER-SCHOOL TRAVEL RESTRICTIONS
# ==============================================================================
# Steps down alpha dynamically as the outbreak progresses.

run_phased_alpha <- function(alpha_schedule, label) {
  times_sim <- seq(0, 200, by = 1)
  curr_state <- initial_state
  p_temp <- params
  history_list <- list()
  
  for (t_idx in 1:(length(times_sim) - 1)) {
    t_start <- times_sim[t_idx]
    t_end   <- times_sim[t_idx + 1]
    
    # Dynamically select alpha based on time step
    current_alpha <- alpha_schedule(t_start)
    p_temp$alpha  <- current_alpha
    
    step_out <- as.data.frame(ode(y = curr_state, times = c(t_start, t_end), 
                                  func = seir_model, parms = p_temp))
    
    curr_state <- unlist(step_out[2, 2:ncol(step_out)])
    history_list[[t_idx]] <- step_out[1, ]
  }
  
  full_out <- bind_rows(history_list)
  get_total_I(full_out, label)
}

# Define alpha transition functions
alpha_constant <- function(t) 0.05
alpha_step_mod <- function(t) if (t >= 14) 0.025 else 0.05
alpha_step_str <- function(t) if (t >= 14) 0.005 else 0.05

curve_alpha_const <- run_phased_alpha(alpha_constant, "Constant Mobility (alpha = 0.05)")
curve_alpha_mod   <- run_phased_alpha(alpha_step_mod,   "Moderate Restriction at Day 14 (alpha = 0.025)")
curve_alpha_str   <- run_phased_alpha(alpha_step_str,   "Strict Restriction at Day 14 (alpha = 0.005)")

sens9_df <- rbind(curve_alpha_const, curve_alpha_mod, curve_alpha_str)

ggplot(sens9_df, aes(x = time, y = total_infectious, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("firebrick", "darkorange", "forestgreen")) +
  geom_vline(xintercept = 14, linetype = "dotted", color = "black") +
  labs(
    title = "Sensitivity Analysis 9: Phased Inter-School Travel Restrictions",
    subtitle = "Evaluating Impact of Timed Reductions in Coupling Coefficient (alpha)",
    x = "Days Since Outbreak Introduction", y = "Total Active Infectious Cases",
    color = "Mobility Policy", linetype = "Mobility Policy"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")


