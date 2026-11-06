library(tidyverse)  # data wrangling and plots
library(deSolve)  # ODE solver

df <- read.csv("data/processed/school_data_joined.csv") %>%  # read the table
  mutate(
    vax_rate = if_else(percent_immunized > 1, percent_immunized / 100, percent_immunized),  # values above 1 are percentages, so divide by 100
    school_id = row_number(),  # index the schools
    page = ceiling(school_id / 6)  # six schools per figure
  )

gamma_rate <- 1 / 7  # recovery rate, 7-day infectious period
sigma_rate <- 1 / 10  # latency rate, 10-day incubation
beta_rate  <- 15 * gamma_rate  # within-school transmission, R0 about 15
times <- seq(0, 200, by = 1)  # days to simulate

seir_model <- function(time, state, parameters) {  # SEIR equations for one time step
  with(as.list(c(state, parameters)), {
    dS <- - (beta * S * I) / N  # infections leave the susceptible class
    dE <-   (beta * S * I) / N - (sigma * E)  # new exposures minus those becoming infectious
    dI <-   (sigma * E) - (gamma * I)  # new infectious students minus recoveries
    dR <-   (gamma * I)  # students leaving the infectious class
    return(list(c(dS, dE, dI, dR)))  # hand the derivatives to the solver
  })
}

all_simulations <- map_dfr(1:nrow(df), function(i) {  # run one school
  school <- df[i, ]

  N <- school$total_students  # students in each school
  vax_rate <- school$vax_rate  # MMR coverage as a proportion

  R0_count <- round(N * vax_rate)  # students already immune
  I0_count <- 1  # index cases at the start
  E0_count <- 0  # exposed at the start, usually none
  S0_count <- max(0, N - I0_count - E0_count - R0_count)  # students still susceptible

  initial_state <- c(S = S0_count, E = E0_count, I = I0_count, R = R0_count)  # starting S, E, I, and R
  params <- c(beta = beta_rate, sigma = sigma_rate, gamma = gamma_rate, N = N)  # parameters passed into the solver

  output <- ode(y = initial_state, times = times, func = seir_model, parms = params)  # solve the SEIR equations

  as.data.frame(output) %>%
    mutate(
      school_name = school$school_name,
      page = school$page,
      school_label = paste0(school$school_name, "\n(Vax: ", round(vax_rate * 100, 1), "%, Pop: ", N, ")")
    )
})

results_long <- all_simulations %>%  # one row per compartment and day
  pivot_longer(cols = c("S", "E", "I", "R"), names_to = "Compartment", values_to = "Count")  # one row per compartment

total_pages <- max(df$page)  # number of figure pages

for (p in 1:total_pages) {  # one figure page at a time
  page_data <- results_long %>% filter(page == p)  # schools drawn on this page

  p_plot <- ggplot(page_data, aes(x = time, y = Count, color = Compartment)) +  # start the plot
    geom_line(linewidth = 1) +  # draw the time series
    scale_color_manual(  # set the legend colours or line types
      values = c("S" = "blue", "E" = "orange", "I" = "red", "R" = "green"),
      labels = c(
        "S" = "Susceptible (S)",
        "E" = "Exposed (E)",
        "I" = "Infectious (I)",
        "R" = "Recovered/Vaccinated (R)"
      )
    ) +
    facet_wrap(~school_label, ncol = 3, scales = "free_y") +  # one panel per group
    labs(  # titles and axis labels
      title = paste0("SEIR Measles Simulations within Schools(Page ", p, " of ", total_pages, ")"),
      x = "Days",
      y = "Number of Students",
      color = "Compartment"
    ) +
    theme_minimal() +  # plain plot theme
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "bottom"
    )

  print(p_plot)  # show the result
}
