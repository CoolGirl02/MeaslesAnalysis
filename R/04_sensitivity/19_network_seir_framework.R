suppressPackageStartupMessages({  # load packages quietly
  library(tidyverse)  # data wrangling and plots
  library(deSolve)  # ODE solver
  library(igraph)  # network tools
})

theme_set(  # default plot theme
  theme_minimal() +  # plain plot theme
    theme(
      plot.title = element_text(face = "bold", size = 12)
    )
)

df <- read.csv(  # read the table
  "data/processed/school_data_joined.csv"
) %>%
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

net_df <- read.csv(  # read the table
  "data/processed/network_df.csv"
) %>%
  rename_with(tolower)

id_cols <- names(net_df)[
  grep(
    "unit|from|to|source|target|id",
    names(net_df)
  )
]

if (length(id_cols) >= 2) {  # use the columns that name the two schools

  net_df <- net_df %>%
    rename(
      unit_id_a = !!sym(id_cols[1]),
      unit_id_b = !!sym(id_cols[2])
    )

} else {

  colnames(net_df)[1:2] <-
    c("unit_id_a", "unit_id_b")
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

    w <- if (
      "total_contacts" %in% names(net_df)
    ) {
      net_df$total_contacts[k]
    } else {
      1
    }

    A[i, j] <- w  # undirected edge weight
    A[j, i] <- w  # undirected edge weight
  }
}

A_raw <- A  # unscaled school-to-school weights

spectral_radius <-  # largest absolute eigenvalue of the network
  max(
    Re(
      eigen(A_raw)$values
    )
  )

A <- A_raw / spectral_radius  # divide every edge by the spectral radius

gamma_rate <- 1 / 7  # recovery rate, 7-day infectious period

sigma_rate <- 1 / 10  # latency rate, 10-day incubation

beta_rate <- 15 * gamma_rate  # within-school transmission, R0 about 15

alpha <- 0.05  # fraction of transmission coming from other schools

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
  beta = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N = N,
  A = A,
  n = n_schools
)

seir_model <- function(  # SEIR equations for one time step
    time,
    state,
    parameters
) {

  with(
    as.list(parameters),
    {

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

      lambda <-  # force of infection from this school and its neighbours
        (beta * I_prop) +
        (
          alpha *
            beta *
            as.vector(
              A %*% I_prop
            )
        )

      dS <- -lambda * S  # infections leave the susceptible class

      dE <-  # new exposures minus those becoming infectious
        lambda * S -
        (sigma * E)

      dI <-  # new infectious students minus recoveries
        (sigma * E) -
        (gamma * I)

      dR <-  # students leaving the infectious class
        gamma * I

      list(
        c(
          dS,
          dE,
          dI,
          dR
        )
      )
    }
  )
}

get_total_I <- function(  # define the function
    out_df,
    scenario_label,
    seed_type
) {

  out_df %>%

    pivot_longer(  # one row per compartment
      cols = -time,
      names_to = "var",
      values_to = "Count"
    ) %>%

    filter(
      str_detect(
        var,
        "^I_"
      )
    ) %>%

    group_by(time) %>%

    summarise(
      total_infectious =
        sum(
          Count,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) %>%

    mutate(
      Scenario = scenario_label,
      `Seed Type` = seed_type
    )
}

run_vax_scenario <- function(  # define the function
    vax_vector,
    seed_id
) {

  R0_vec <- round(
    N * vax_vector
  )

  I0_vec <- rep(
    0,
    n_schools
  )

  I0_vec[seed_id] <- 1

  E0_vec <- rep(
    0,
    n_schools
  )

  S0_vec <- pmax(  # larger school id second
    0,
    N -
      I0_vec -
      E0_vec -
      R0_vec
  )

  state_vec <- c(
    S0_vec,
    E0_vec,
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

  return(out_df)
}

output <- ode(  # solve the SEIR equations
  y = initial_state,
  times = times,
  func = seir_model,
  parms = params
)

results_long <- as.data.frame(output) %>%  # one row per compartment and day

  pivot_longer(  # one row per compartment
    cols = -time,
    names_to = "var",
    values_to = "Count"
  ) %>%

  mutate(
    Compartment =
      str_extract(
        var,
        "^[SEIR]"
      ),

    school_id =
      as.numeric(
        str_extract(
          var,
          "\\d+"
        )
      )
  ) %>%

  left_join(  # join the tables
    df,
    by = "school_id"
  ) %>%

  mutate(
    school_label =
      paste0(
        school_name,
        "\n(Vax: ",
        round(vax_rate * 100, 1),
        "%, Pop: ",
        total_students,
        ")"
      )
  )

total_pages <- max(  # number of figure pages
  df$page
)

for (p in 1:total_pages) {  # one figure page at a time

  page_data <-  # schools drawn on this page
    results_long %>%
    filter(page == p)

  p_plot <- ggplot(  # start the plot
    page_data,
    aes(
      x = time,
      y = Count,
      color = Compartment
    )
  ) +

    geom_line(  # draw the time series
      linewidth = 1
    ) +

    scale_color_manual(  # set the legend colours or line types
      values = c(
        "S" = "blue",
        "E" = "orange",
        "I" = "red",
        "R" = "green"
      ),

      labels = c(
        "S" = "Susceptible",
        "E" = "Exposed",
        "I" = "Infectious",
        "R" = "Recovered/Vaccinated"
      )
    ) +

    facet_wrap(  # one panel per group
      ~school_label,
      ncol = 3,
      scales = "free_y"
    ) +

    labs(  # titles and axis labels
      title =
        paste0(
          "Baseline SEIR Measles Simulations (Page ",
          p,
          " of ",
          total_pages,
          ")"
        ),

      x = "Days",

      y = "Number of Students",

      color = "Compartment"
    ) +

    theme(
      strip.text =
        element_text(
          size = 8,
          face = "bold"
        ),

      legend.position = "bottom"
    )

  print(p_plot)  # show the result
}

net_graph <-
  graph_from_adjacency_matrix(
    A,
    mode = "undirected",
    weighted = TRUE
  )

df$node_strength <-
  strength(  # sum of this school's edge weights
    net_graph
  )

df$betweenness <-
  betweenness(  # schools that sit on many shortest paths
    net_graph,
    weights = NA
  )

top_hub_id <-
  df$school_id[
    which.max(
      df$node_strength
    )
  ]

low_hub_id <-
  df$school_id[
    which.min(
      df$node_strength
    )
  ]

top_hub_name <-
  df$school_name[
    df$school_id == top_hub_id
  ]

low_hub_name <-
  df$school_name[
    df$school_id == low_hub_id
  ]

cat("\n========================================\n")  # print a check
cat("NETWORK HUB AND PERIPHERAL SCHOOL\n")  # print a check
cat("========================================\n")  # print a check

cat(  # print a check
  "Network Hub:",
  top_hub_name,
  "\n"
)

cat(  # print a check
  "Peripheral:",
  low_hub_name,
  "\n"
)

vax_empirical <-
  df$vax_rate

out_hub_baseline <-
  run_vax_scenario(
    vax_empirical,
    top_hub_id
  )

out_peri_baseline <-
  run_vax_scenario(
    vax_empirical,
    low_hub_id
  )

vax_74 <-
  df$vax_rate

vax_74[top_hub_id] <-
  0.74

vax_74[low_hub_id] <-
  0.74

out_hub_74 <-
  run_vax_scenario(
    vax_74,
    top_hub_id
  )

out_peri_74 <-
  run_vax_scenario(
    vax_74,
    low_hub_id
  )

vax_916 <-
  df$vax_rate

vax_916[top_hub_id] <-
  0.916

vax_916[low_hub_id] <-
  0.916

out_hub_916 <-
  run_vax_scenario(
    vax_916,
    top_hub_id
  )

out_peri_916 <-
  run_vax_scenario(
    vax_916,
    low_hub_id
  )

comparison_all <- bind_rows(  # combine these rows

  get_total_I(
    out_hub_baseline,
    "Scenario 1: Empirical Coverage",
    "Network Hub"
  ),

  get_total_I(
    out_peri_baseline,
    "Scenario 1: Empirical Coverage",
    "Peripheral"
  ),

  get_total_I(
    out_hub_74,
    "Scenario 2: Both Seed Schools at 74%",
    "Network Hub"
  ),

  get_total_I(
    out_peri_74,
    "Scenario 2: Both Seed Schools at 74%",
    "Peripheral"
  ),

  get_total_I(
    out_hub_916,
    "Scenario 3: Both Seed Schools at 91.6%",
    "Network Hub"
  ),

  get_total_I(
    out_peri_916,
    "Scenario 3: Both Seed Schools at 91.6%",
    "Peripheral"
  )
)

ggplot(  # start the plot
  comparison_all,
  aes(
    x = time,
    y = total_infectious,
    color = `Seed Type`
  )
) +

  geom_line(  # draw the time series
    linewidth = 1.2
  ) +

  facet_wrap(  # one panel per group
    ~Scenario,
    ncol = 3,
    scales = "free_y"
  ) +

  scale_color_manual(  # set the legend colours or line types
    values = c(
      "Network Hub" = "firebrick",
      "Peripheral" = "darkblue"
    )
  ) +

  labs(  # titles and axis labels
    title =
      "Sensitivity Analysis: Hub vs Peripheral Seeding",

    subtitle =
      "Comparison across three vaccination coverage scenarios",

    x =
      "Days Since Outbreak Introduction",

    y =
      "Total Active Infectious Cases",

    color =
      "Seed Location"
  ) +

  theme(
    legend.position = "bottom",

    strip.text =
      element_text(
        face = "bold",
        size = 10
      )
  )

calculate_outbreak_measures <- function(  # define the function
    out_df,
    seed_id
) {

  R_columns <-
    grep(
      "^R_",
      names(out_df),
      value = TRUE
    )

  I_columns <-
    grep(
      "^I_",
      names(out_df),
      value = TRUE
    )

  R_final <-
    sapply(
      R_columns,
      function(x) {  # define the function
        out_df[[x]][
          nrow(out_df)
        ]
      }
    )

  R_initial <-
    round(
      sapply(
        R_columns,
        function(x) {  # define the function
          out_df[[x]][1]
        }
      )
    )

  new_cases_each_school <-
    pmax(  # larger school id second
      0,
      R_final - R_initial
    )

  total_cases <-
    sum(
      new_cases_each_school,
      na.rm = TRUE
    )

  total_I <-
    rowSums(
      out_df[
        ,
        I_columns,
        drop = FALSE
      ],
      na.rm = TRUE
    )

  peak_burden <-
    max(
      total_I,
      na.rm = TRUE
    )

  peak_day <-
    out_df$time[
      which.max(total_I)
    ]

  county_population <-
    sum(
      N,
      na.rm = TRUE
    )

  attack_rate <-
    (
      total_cases /
        county_population
    ) * 100

  return(
    c(
      Total_Cases = total_cases,
      Peak_Burden = peak_burden,
      Peak_Day = peak_day,
      Attack_Rate = attack_rate
    )
  )
}

baseline_hub <-
  calculate_outbreak_measures(
    out_hub_baseline,
    top_hub_id
  )

baseline_peri <-
  calculate_outbreak_measures(
    out_peri_baseline,
    low_hub_id
  )

control_74_hub <-
  calculate_outbreak_measures(
    out_hub_74,
    top_hub_id
  )

control_74_peri <-
  calculate_outbreak_measures(
    out_peri_74,
    low_hub_id
  )

control_916_hub <-
  calculate_outbreak_measures(
    out_hub_916,
    top_hub_id
  )

control_916_peri <-
  calculate_outbreak_measures(
    out_peri_916,
    low_hub_id
  )

results_4_3 <- tibble(

  Scenario = c(

    "Scenario 1: Empirical Coverage",
    "Scenario 1: Empirical Coverage",

    "Scenario 2: Both Seed Schools at 74%",
    "Scenario 2: Both Seed Schools at 74%",

    "Scenario 3: Both Seed Schools at 91.6%",
    "Scenario 3: Both Seed Schools at 91.6%"
  ),

  `Seed Location` = c(

    top_hub_name,
    low_hub_name,

    top_hub_name,
    low_hub_name,

    top_hub_name,
    low_hub_name
  ),

  `Seed Type` = c(

    "Network Hub",
    "Peripheral",

    "Network Hub",
    "Peripheral",

    "Network Hub",
    "Peripheral"
  ),

  `Vaccination Coverage` = c(

    paste0(
      round(
        vax_empirical[top_hub_id] * 100,
        1
      ),
      "%"
    ),

    paste0(
      round(
        vax_empirical[low_hub_id] * 100,
        1
      ),
      "%"
    ),

    "74.0%",
    "74.0%",

    "91.6%",
    "91.6%"
  ),

  `Total New Cases (C)` = c(

    baseline_hub["Total_Cases"],
    baseline_peri["Total_Cases"],

    control_74_hub["Total_Cases"],
    control_74_peri["Total_Cases"],

    control_916_hub["Total_Cases"],
    control_916_peri["Total_Cases"]
  ),

  `Peak Burden (Ipeak)` = c(

    baseline_hub["Peak_Burden"],
    baseline_peri["Peak_Burden"],

    control_74_hub["Peak_Burden"],
    control_74_peri["Peak_Burden"],

    control_916_hub["Peak_Burden"],
    control_916_peri["Peak_Burden"]
  ),

  `Peak Day` = c(

    baseline_hub["Peak_Day"],
    baseline_peri["Peak_Day"],

    control_74_hub["Peak_Day"],
    control_74_peri["Peak_Day"],

    control_916_hub["Peak_Day"],
    control_916_peri["Peak_Day"]
  ),

  `Attack Rate (AR)` = c(

    baseline_hub["Attack_Rate"],
    baseline_peri["Attack_Rate"],

    control_74_hub["Attack_Rate"],
    control_74_peri["Attack_Rate"],

    control_916_hub["Attack_Rate"],
    control_916_peri["Attack_Rate"]
  )
)

table_4_3 <- results_4_3 %>%

  mutate(

    `Total New Cases (C)` =
      round(
        `Total New Cases (C)`
      ),

    `Peak Burden (Ipeak)` =
      round(
        `Peak Burden (Ipeak)`
      ),

    `Peak Day` =
      round(
        `Peak Day`
      ),

    `Attack Rate (AR)` =
      paste0(
        round(
          `Attack Rate (AR)`,
          2
        ),
        "%"
      )
  )

cat("\n\n")  # print a check
cat("============================================================\n")  # print a check
cat("TABLE 4.3: COMPARATIVE ANALYSIS OF HUB VS PERIPHERAL\n")  # print a check
cat("SEEDING SCENARIOS\n")  # print a check
cat("============================================================\n\n")  # print a check

print(  # show the result
  table_4_3,
  n = Inf
)

if (
  requireNamespace(
    "knitr",
    quietly = TRUE
  )
) {

  print(  # show the result
    knitr::kable(
      table_4_3,
      caption =
        "Comparative Analysis of Hub vs Peripheral Seeding Scenarios",
      align = c(
        "l",
        "l",
        "l",
        "l",
        "r",
        "r",
        "r",
        "r"
      )
    )
  )
}

write.csv(  # save the table
  table_4_3,
  "output/tables/table_4_3_hub_vs_peripheral.csv",
  row.names = FALSE
)
