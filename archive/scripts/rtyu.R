VACCINE ALLOCATION BY SCHOOL SIZE VS. BASELINES

suppressPackageStartupMessages({  # load packages quietly
  library(tidyverse)  # data wrangling and plots
  library(deSolve)  # ODE solver
  library(igraph)  # network tools
})

theme_set(  # default plot theme
  theme_minimal() +  # plain plot theme
    theme(plot.title = element_text(face = "bold", size = 12))
)

df <- read.csv("data/processed/school_data_joined.csv") %>%  # read the table
  rename_with(tolower)

if (!"unit_id" %in% names(df)) {  # add a school id if the column is missing
  df <- df %>%
    mutate(unit_id = as.character(row_number()))  # index the schools
} else {
  df <- df %>%
    mutate(unit_id = as.character(unit_id))
}

df <- df %>%
  mutate(
    vax_rate = if_else(
      percent_immunized > 1,  # values above 1 are percentages, so divide by 100
      percent_immunized / 100,
      percent_immunized
    ),
    school_id = row_number(),  # index the schools
    page = ceiling(school_id / 6)  # six schools per figure
  )

n_schools <- nrow(df)  # number of schools

net_df <- read.csv("data/processed/network_df.csv") %>%  # read the table
  rename_with(tolower)

id_cols <- names(net_df)[
  grep("unit|from|to|source|target|id", names(net_df))
]

if (length(id_cols) >= 2) {  # use the columns that name the two schools

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

id_map <- setNames(
  1:n_schools,
  df$unit_id
)

A <- matrix(  # adjacency matrix used in the model
  0,
  nrow = n_schools,
  ncol = n_schools
)

for (k in seq_len(nrow(net_df))) {  # one edge at a time

  i <- id_map[net_df$unit_id_a[k]]
  j <- id_map[net_df$unit_id_b[k]]

  if (!is.na(i) && !is.na(j)) {  # keep the edge only if both schools are in the table

    w <- if ("total_contacts" %in% names(net_df)) {
      net_df$total_contacts[k]
    } else {
      1
    }

    A[i, j] <- w  # undirected edge weight
    A[j, i] <- w  # undirected edge weight
  }
}

A_raw <- A  # unscaled school-to-school weights

spectral_radius <- max(  # largest absolute eigenvalue of the network
  Re(eigen(A_raw)$values)
)

A <- A_raw / spectral_radius  # divide every edge by the spectral radius

school_strength <- rowSums(A_raw)

top_hub_id <- which.max(
  school_strength
)

cat("\nHigh-centrality seed school:", top_hub_id, "\n")  # print a check
cat(  # print a check
  "School name:",
  df$school_name[top_hub_id],
  "\n\n"
)

gamma_rate <- 1 / 7  # recovery rate, 7-day infectious period
sigma_rate <- 1 / 10  # latency rate, 10-day incubation
beta_rate  <- 15 * gamma_rate  # within-school transmission, R0 about 15
alpha      <- 0.05  # fraction of transmission coming from other schools

times <- seq(  # days to simulate
  0,
  200,
  by = 1
)

N <- df$total_students  # students in each school

R0 <- round(  # students already immune
  N * df$vax_rate
)

I0 <- c(  # index cases at the start
  1,
  rep(0, n_schools - 1)
)

E0 <- rep(  # exposed at the start, usually none
  0,
  n_schools
)

S0 <- pmax(  # students still susceptible
  0,
  N - I0 - E0 - R0
)

initial_state <- c(  # starting S, E, I, and R
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

    I_prop <- I / N  # infectious fraction in each school

    lambda <- (  # force of infection from this school and its neighbours
      beta * I_prop
    ) +
      (
        alpha *
          beta *
          as.vector(A %*% I_prop)
      )

    dS <- -lambda * S  # infections leave the susceptible class

    dE <- lambda * S -  # new exposures minus those becoming infectious
      sigma_rate * E

    dI <- sigma_rate * E -  # new infectious students minus recoveries
      gamma_rate * I

    dR <- gamma_rate * I  # students leaving the infectious class

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

cat("\n============================================\n")  # print a check
cat("AGGREGATE VACCINATION COVERAGE\n")  # print a check
cat("============================================\n")  # print a check

cat(  # print a check
  "Uniform 70%:",
  round(rate_70 * 100, 4),
  "%\n"
)

cat(  # print a check
  "Uniform 80%:",
  round(rate_80 * 100, 4),
  "%\n"
)

cat(  # print a check
  "Low Small / High Large:",
  round(rate_A * 100, 4),
  "%\n"
)

cat(  # print a check
  "High Small / Low Large:",
  round(rate_B * 100, 4),
  "%\n"
)

cat(  # print a check
  "Uniform 90%:",
  round(rate_90 * 100, 4),
  "%\n"
)

cat(  # print a check
  "Difference between the two allocated scenarios:",
  round(abs(rate_A - rate_B) * 100, 8),
  "percentage points\n"
)

cat("============================================\n\n")  # print a check

run_vax_scenario <- function(  # define the function
    vax_vector,
    label,
    seed_id = 1) {

  R0_vec <- round(
    N * vax_vector
  )

  I0_vec <- rep(
    0,
    n_schools
  )

  I0_vec[seed_id] <- 1

  S0_vec <- pmax(  # larger school id second
    0,
    N - I0_vec - R0_vec
  )

  state_vec <- c(
    S0_vec,
    rep(0, n_schools),
    I0_vec,
    R0_vec
  )

  names(state_vec) <-
    names(initial_state)

  out_df <- as.data.frame(
    ode(  # solve the SEIR equations
      y = state_vec,
      times = times,
      func = seir_model,
      parms = params
    )
  )

  total_I <- rowSums(
    out_df[
      ,
      paste0(
        "I_",
        1:n_schools
      )
    ]
  )

  peak_burden <- max(
    total_I
  )

  peak_day <- out_df$time[
    which.max(total_I)
  ]

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

  county_attack_rate <-
    total_new_cases /
    N_total

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

sens5_df <- bind_rows(  # combine these rows

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

cat("\n")  # print a check
cat("=====================================================================\n")  # print a check
cat("TABLE 4.2: EPIDEMIC OUTCOMES UNDER SIZE-BASED ALLOCATION\n")  # print a check
cat("              AND UNIFORM BASELINE STRATEGIES\n")  # print a check
cat("=====================================================================\n\n")  # print a check

print(  # show the result
  table4_2_display,
  row.names = FALSE
)

cat("\n=====================================================================\n")  # print a check

ggplot(  # start the plot
  sens5_df,
  aes(
    x = time,
    y = total_infectious,
    color = Scenario,
    linetype = Scenario
  )
) +

  geom_line(  # draw the time series
    linewidth = 1.2
  ) +

  scale_color_manual(  # set the legend colours or line types
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

  scale_linetype_manual(  # set the legend colours or line types
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

  labs(  # titles and axis labels
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
