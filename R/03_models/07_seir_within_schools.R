library(tidyverse)
library(deSolve)

df <- read.csv("data/processed/school_data_joined.csv") %>%
  mutate(
    vax_rate = if_else(percent_immunized > 1, percent_immunized / 100, percent_immunized),
    school_id = row_number(),
    page = ceiling(school_id / 6)
  )

gamma_rate <- 1 / 7
sigma_rate <- 1 / 10
beta_rate  <- 15 * gamma_rate
times <- seq(0, 200, by = 1)

seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- - (beta * S * I) / N
    dE <-   (beta * S * I) / N - (sigma * E)
    dI <-   (sigma * E) - (gamma * I)
    dR <-   (gamma * I)
    return(list(c(dS, dE, dI, dR)))
  })
}

all_simulations <- map_dfr(1:nrow(df), function(i) {
  school <- df[i, ]

  N <- school$total_students
  vax_rate <- school$vax_rate

  R0_count <- round(N * vax_rate)
  I0_count <- 1
  E0_count <- 0
  S0_count <- max(0, N - I0_count - E0_count - R0_count)

  initial_state <- c(S = S0_count, E = E0_count, I = I0_count, R = R0_count)
  params <- c(beta = beta_rate, sigma = sigma_rate, gamma = gamma_rate, N = N)

  output <- ode(y = initial_state, times = times, func = seir_model, parms = params)

  as.data.frame(output) %>%
    mutate(
      school_name = school$school_name,
      page = school$page,
      school_label = paste0(school$school_name, "\n(Vax: ", round(vax_rate * 100, 1), "%, Pop: ", N, ")")
    )
})

results_long <- all_simulations %>%
  pivot_longer(cols = c("S", "E", "I", "R"), names_to = "Compartment", values_to = "Count")

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
      title = paste0("SEIR Measles Simulations within Schools(Page ", p, " of ", total_pages, ")"),
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
