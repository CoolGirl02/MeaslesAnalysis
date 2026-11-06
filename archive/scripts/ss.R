suppressPackageStartupMessages({  # load packages quietly
  library(tidyverse)  # data wrangling and plots
  library(deSolve)  # ODE solver
  library(igraph)  # network tools
})

theme_set(theme_minimal() + theme(plot.title = element_text(face = "bold", size = 12)))  # default plot theme

df <- read.csv("data/processed/school_data_joined.csv") %>%  # read the table
  rename_with(tolower)

if (!"unit_id" %in% names(df)) {  # add a school id if the column is missing
  df <- df %>% mutate(unit_id = as.character(row_number()))  # index the schools
} else {
  df <- df %>% mutate(unit_id = as.character(unit_id))
}

df <- df %>%
  mutate(
    vax_rate  = if_else(percent_immunized > 1, percent_immunized / 100, percent_immunized),  # values above 1 are percentages, so divide by 100
    school_id = row_number(),  # index the schools
    page      = ceiling(school_id / 6)  # six schools per figure
  )

n_schools <- nrow(df)  # number of schools

net_df <- read.csv("data/processed/network_df.csv") %>%  # read the table
  rename_with(tolower)

id_cols <- names(net_df)[grep("unit|from|to|source|target|id", names(net_df))]

if (length(id_cols) >= 2) {  # use the columns that name the two schools
  net_df <- net_df %>%
    rename(unit_id_a = !!sym(id_cols[1]), unit_id_b = !!sym(id_cols[2]))
} else {
  colnames(net_df)[1:2] <- c("unit_id_a", "unit_id_b")
}

net_df <- net_df %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  )

id_map <- setNames(1:n_schools, df$unit_id)
A      <- matrix(0, nrow = n_schools, ncol = n_schools)  # adjacency matrix used in the model

for (k in seq_len(nrow(net_df))) {  # one edge at a time
  i <- id_map[net_df$unit_id_a[k]]
  j <- id_map[net_df$unit_id_b[k]]
  if (!is.na(i) && !is.na(j)) {  # keep the edge only if both schools are in the table
    w <- if ("total_contacts" %in% names(net_df)) net_df$total_contacts[k] else 1
    A[i, j] <- w  # undirected edge weight
    A[j, i] <- w  # undirected edge weight
  }
}

A_raw <- A  # unscaled school-to-school weights

spectral_radius <- max(Re(eigen(A_raw)$values))  # largest absolute eigenvalue of the network

A <- A_raw / spectral_radius  # divide every edge by the spectral radius

gamma_rate <- 1 / 7  # recovery rate, 7-day infectious period
sigma_rate <- 1 / 10  # latency rate, 10-day incubation
beta_rate  <- 15 * gamma_rate  # within-school transmission, R0 about 15
alpha      <- 0.05  # fraction of transmission coming from other schools
times      <- seq(0, 200, by = 1)  # days to simulate

N  <- df$total_students  # students in each school
R0 <- round(N * df$vax_rate)  # students already immune
I0 <- c(1, rep(0, n_schools - 1))  # index cases at the start
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
  beta  = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N     = N,
  A     = A,
  n     = n_schools
)

seir_model <- function(time, state, parameters) {  # SEIR equations for one time step
  with(as.list(parameters), {
    S <- state[1:n]
    E <- state[(n + 1):(2 * n)]
    I <- state[(2 * n + 1):(3 * n)]
    R <- state[(3 * n + 1):(4 * n)]

    I_prop <- I / N  # infectious fraction in each school
    lambda <- (beta * I_prop) + (alpha * beta * as.vector(A %*% I_prop))  # force of infection from this school and its neighbours

    dS <- -lambda * S  # infections leave the susceptible class
    dE <-  lambda * S - (sigma * E)  # new exposures minus those becoming infectious
    dI <- (sigma * E) - (gamma * I)  # new infectious students minus recoveries
    dR <- (gamma * I)  # students leaving the infectious class

    list(c(dS, dE, dI, dR))
  })
}

get_total_I <- function(out_df, scenario_label) {  # define the function
  out_df %>%
    pivot_longer(cols = -time, names_to = "var", values_to = "Count") %>%  # one row per compartment
    filter(str_detect(var, "^I_")) %>%
    group_by(time) %>%
    summarise(total_infectious = sum(Count), .groups = "drop") %>%
    mutate(Scenario = scenario_label)
}

run_vax_scenario <- function(vax_vector, label, seed_id = 1) {  # define the function
  R0_vec <- round(N * vax_vector)
  I0_vec <- rep(0, n_schools)
  I0_vec[seed_id] <- 1
  S0_vec <- pmax(0, N - I0_vec - R0_vec)  # larger school id second

  state_vec <- c(S0_vec, rep(0, n_schools), I0_vec, R0_vec)
  names(state_vec) <- names(initial_state)

  out_df <- as.data.frame(ode(y = state_vec, times = times, func = seir_model, parms = params))
  get_total_I(out_df, label)
}

output <- ode(y = initial_state, times = times, func = seir_model, parms = params)  # solve the SEIR equations

results_long <- as.data.frame(output) %>%  # one row per compartment and day
  pivot_longer(cols = -time, names_to = "var", values_to = "Count") %>%  # one row per compartment
  mutate(
    Compartment = str_extract(var, "^[SEIR]"),
    school_id   = as.numeric(str_extract(var, "\\d+"))
  ) %>%
  left_join(df, by = "school_id") %>%  # join the tables
  mutate(school_label = paste0(school_name, "\n(Vax: ", round(vax_rate * 100, 1), "%, Pop: ", total_students, ")"))

total_pages <- max(df$page)  # number of figure pages

for (p in 1:total_pages) {  # one figure page at a time
  page_data <- results_long %>% filter(page == p)  # schools drawn on this page

  p_plot <- ggplot(page_data, aes(x = time, y = Count, color = Compartment)) +  # start the plot
    geom_line(linewidth = 1) +  # draw the time series
    scale_color_manual(  # set the legend colours or line types
      values = c("S" = "blue", "E" = "orange", "I" = "red", "R" = "green"),
      labels = c("S" = "Susceptible", "E" = "Exposed", "I" = "Infectious", "R" = "Recovered/Vaccinated")
    ) +
    facet_wrap(~school_label, ncol = 3, scales = "free_y") +  # one panel per group
    labs(  # titles and axis labels
      title = paste0("Baseline SEIR Measles Simulations (Page ", p, " of ", total_pages, ")"),
      x = "Days", y = "Number of Students", color = "Compartment"
    ) +
    theme(strip.text = element_text(size = 8, face = "bold"), legend.position = "bottom")

  print(p_plot)  # show the result
}

n_schools <- nrow(df)  # number of schools

small_idx <- which(df$total_students < 400)
large_idx <- which(df$total_students >= 400)

N_small <- sum(df$total_students[small_idx])
N_large <- sum(df$total_students[large_idx])
N_total <- sum(df$total_students)

vax_small_low <- numeric(n_schools)
vax_small_low[small_idx] <- 0.75
vax_small_low[large_idx] <- 0.90

total_vax_A <- sum(df$total_students * vax_small_low)
rate_A      <- total_vax_A / N_total

vax_small_B_needed <- (total_vax_A - (0.75 * N_large)) / N_small

vax_large_low <- numeric(n_schools)
vax_large_low[large_idx] <- 0.75
vax_large_low[small_idx] <- vax_small_B_needed

total_vax_B <- sum(df$total_students * vax_large_low)
rate_B      <- total_vax_B / N_total

vax_70 <- rep(0.70, n_schools)
vax_80 <- rep(0.80, n_schools)
vax_90 <- rep(0.90, n_schools)

net_graph       <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
df$node_strength <- strength(net_graph)  # sum of this school's edge weights
df$betweenness   <- betweenness(net_graph, weights = NA)  # schools that sit on many shortest paths

top_hub_id   <- df$school_id[which.max(df$node_strength)]
low_hub_id   <- df$school_id[which.min(df$node_strength)]
top_hub_name <- df$school_name[df$school_id == top_hub_id]
low_hub_name <- df$school_name[df$school_id == low_hub_id]

R0_actual <- round(N * df$vax_rate)

I0_hub <- rep(0, n_schools); I0_hub[top_hub_id] <- 1
state_hub <- c(pmax(0, N - I0_hub - R0_actual), rep(0, n_schools), I0_hub, R0_actual)  # larger school id second
names(state_hub) <- names(initial_state)

I0_peri <- rep(0, n_schools); I0_peri[low_hub_id] <- 1
state_peri <- c(pmax(0, N - I0_peri - R0_actual), rep(0, n_schools), I0_peri, R0_actual)  # larger school id second
names(state_peri) <- names(initial_state)

out_hub  <- as.data.frame(ode(y = state_hub,  times = times, func = seir_model, parms = params))
out_peri <- as.data.frame(ode(y = state_peri, times = times, func = seir_model, parms = params))

comparison_df <- bind_rows(  # combine these rows
  get_total_I(out_hub,  paste0("Hub: ", top_hub_name)),
  get_total_I(out_peri, paste0("Peripheral: ", low_hub_name))
)

ggplot(comparison_df, aes(x = time, y = total_infectious, color = Scenario)) +  # start the plot
  geom_line(linewidth = 1.3) +  # draw the time series
  scale_color_manual(values = c("firebrick", "darkblue")) +  # set the legend colours or line types
  labs(  # titles and axis labels
    title    = "Sensitivity Analysis 2: Outbreak Dynamics by Index Case Location",
    subtitle = "Comparing Index Case in  Hub vs. Peripheral School",
    x        = "Days Since Outbreak Introduction",
    y        = "Total Active Infectious Cases in County",
    color    = "Initial Location"
  ) +
  theme(legend.position = "bottom")

fixed_vax_rate <- 0.74
R0_controlled  <- round(N * df$vax_rate)
R0_controlled[top_hub_id] <- round(N[top_hub_id] * fixed_vax_rate)
R0_controlled[low_hub_id] <- round(N[low_hub_id] * fixed_vax_rate)

state_hub_ctrl  <- c(pmax(0, N - I0_hub - R0_controlled), rep(0, n_schools), I0_hub, R0_controlled)  # larger school id second
state_peri_ctrl <- c(pmax(0, N - I0_peri - R0_controlled), rep(0, n_schools), I0_peri, R0_controlled)  # larger school id second
names(state_hub_ctrl) <- names(state_peri_ctrl) <- names(initial_state)

out_hub_ctrl  <- as.data.frame(ode(y = state_hub_ctrl,  times = times, func = seir_model, parms = params))
out_peri_ctrl <- as.data.frame(ode(y = state_peri_ctrl, times = times, func = seir_model, parms = params))

comparison_ctrl_df <- bind_rows(  # combine these rows
  get_total_I(out_hub_ctrl,  paste0("Hub Seed (74% Vax): ", top_hub_name)),
  get_total_I(out_peri_ctrl, paste0("Peripheral Seed (74% Vax): ", low_hub_name))
)

ggplot(comparison_ctrl_df, aes(x = time, y = total_infectious, color = Scenario)) +  # start the plot
  geom_line(linewidth = 1.3) +  # draw the time series
  scale_color_manual(values = c("firebrick", "darkblue")) +  # set the legend colours or line types
  labs(  # titles and axis labels
    title    = "Controlled Sensitivity Analysis 2: Decoupled Network Effect",
    subtitle = "Both Seed Schools Controlled at 74% Vaccination Coverage",
    x        = "Days Since Outbreak Introduction",
    y        = "Total Active Infectious Cases in County"

    run_phased_alpha <- function(alpha_schedule, label) {  # define the function
      curr_state   <- initial_state
      p_temp       <- params
      history_list <- list()
      times_sim    <- seq(0, 200, by = 1)

      for (t_idx in 1:(length(times_sim) - 1)) {  # repeat for each item
        t_start      <- times_sim[t_idx]
        t_end        <- times_sim[t_idx + 1]
        p_temp$alpha <- alpha_schedule(t_start)

        step_out   <- as.data.frame(ode(y = curr_state, times = c(t_start, t_end), func = seir_model, parms = p_temp))
        curr_state <- unlist(step_out[2, 2:ncol(step_out)])
        history_list[[t_idx]] <- step_out[1, ]
      }

      full_out <- bind_rows(history_list)  # combine these rows
      get_total_I(full_out, label)
    }

    alpha_constant <- function(t) 0.05  # define the function
    alpha_step_mod <- function(t) if (t >= 14) 0.025 else 0.05  # define the function
    alpha_step_str <- function(t) if (t >= 14) 0.005 else 0.05  # define the function

    sens9_df <- bind_rows(  # combine these rows
      run_phased_alpha(alpha_constant, "Constant Mobility (alpha = 0.05)"),
      run_phased_alpha(alpha_step_mod, "Moderate Restriction at Day 14 (alpha = 0.025)"),
      run_phased_alpha(alpha_step_str, "Strict Restriction at Day 14 (alpha = 0.005)"),

    ggplot(sens9_df, aes(x = time, y = total_infectious, color = Scenario, linetype = Scenario)) +  # start the plot
      geom_line(linewidth = 1.3) +  # draw the time series
      scale_color_manual(values = c("firebrick", "darkorange", "forestgreen")) +  # set the legend colours or line types
      geom_vline(xintercept = 14, linetype = "dotted", color = "black") +  # mark the threshold
      labs(  # titles and axis labels
        title    = "Sensitivity Analysis 8: Phased Inter-School Travel Restrictions",
        subtitle = "Evaluating Post-14 Day Inter-School Mobility Reductions",
        x        = "Days Since Outbreak Introduction",
        y        = "Total Active Infectious Cases in County"
      ) +
      theme(legend.position = "bottom")

  )
