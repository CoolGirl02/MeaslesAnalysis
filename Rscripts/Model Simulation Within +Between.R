

library(tidyverse)
library(deSolve)


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


#  Building Adjacency Matrix (A)

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


# Model Parameters & Initial Conditions

gamma_rate <- 1 / 7                    # Infectious period (~7 days)
sigma_rate <- 1 / 10                   # Incubation period (~10 days)
beta_rate  <- 15 * gamma_rate          # Intra-school transmission (R0 = 15)
alpha      <- 0.05                     # Alpha: Inter-school coupling coefficient (5%)
times      <- seq(0, 200, by = 1)

N <- df$total_students
R0_count <- round(N * df$vax_rate)
I0_count <- c(1, rep(0, n_schools - 1)) # Outbreak starts with 1 case in School 1
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


seir_model <- function(time, state, parameters) {
  with(as.list(parameters), {
    S <- state[1:n]
    E <- state[(n + 1):(2 * n)]
    I <- state[(2 * n + 1):(3 * n)]
    R <- state[(3 * n + 1):(4 * n)]
    
    
    I_prop <- I / N
    lambda <- (beta * I_prop) + (alpha * beta * as.vector(A %*% I_prop))
    
    dS <- - lambda * S
    dE <-   lambda * S - (sigma * E)
    dI <-  (sigma * E) - (gamma * I)
    dR <-  (gamma * I)
    
    return(list(c(dS, dE, dI, dR)))
  })
}


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
      labels = c(
        "S" = "Susceptible (S)",
        "E" = "Exposed (E)",
        "I" = "Infectious (I)",
        "R" = "Recovered/Vaccinated (R)"
      )
    ) +
    facet_wrap(~school_label, ncol = 3, scales = "free_y") +
    labs(
      title = paste0("SEIR Measles Simulations (Page ", p, " of ", total_pages, ")"),
      x = "Days",
      y = "Number of Students",
      color = "Compartment"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p_plot)
}

