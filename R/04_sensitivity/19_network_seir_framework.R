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

df <- read.csv(
  "data/processed/school_data_joined.csv"
) %>%
  rename_with(tolower)

if (!"unit_id" %in% names(df)) {
  df <- df %>%
    mutate(unit_id = as.character(row_number()))
} else {
  df <- df %>%
    mutate(unit_id = as.character(unit_id))
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

net_df <- read.csv(
  "data/processed/network_df.csv"
) %>%
  rename_with(tolower)

id_cols <- names(net_df)[
  grep(
    "unit|from|to|source|target|id",
    names(net_df)
  )
]

if (length(id_cols) >= 2) {

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

A <- matrix(
  0,
  nrow = n_schools,
  ncol = n_schools
)

for (k in seq_len(nrow(net_df))) {

  i <- id_map[net_df$unit_id_a[k]]

  j <- id_map[net_df$unit_id_b[k]]

  if (!is.na(i) && !is.na(j)) {

    w <- if (
      "total_contacts" %in% names(net_df)
    ) {
      net_df$total_contacts[k]
    } else {
      1
    }

    A[i, j] <- w
    A[j, i] <- w
  }
}

A_raw <- A

spectral_radius <-
  max(
    Re(
      eigen(A_raw)$values
    )
  )

A <- A_raw / spectral_radius

gamma_rate <- 1 / 7

sigma_rate <- 1 / 10

beta_rate <- 15 * gamma_rate

alpha <- 0.05

times <- seq(
  0,
  200,
  by = 1
)

N <- df$total_students

R0 <- round(
  N * df$vax_rate
)

I0 <- c(
  1,
  rep(0, n_schools - 1)
)

E0 <- rep(
  0,
  n_schools
)

S0 <- pmax(
  0,
  N - I0 - E0 - R0
)

initial_state <- c(
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

params <- list(
  beta = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N = N,
  A = A,
  n = n_schools
)

seir_model <- function(
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

      I_prop <- I / N

      lambda <-
        (beta * I_prop) +
        (
          alpha *
            beta *
            as.vector(
              A %*% I_prop
            )
        )

      dS <- -lambda * S

      dE <-
        lambda * S -
        (sigma * E)

      dI <-
        (sigma * E) -
        (gamma * I)

      dR <-
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

get_total_I <- function(
    out_df,
    scenario_label,
    seed_type
) {

  out_df %>%

    pivot_longer(
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

run_vax_scenario <- function(
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

  S0_vec <- pmax(
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
    ode(
      y = state_vec,
      times = times,
      func = seir_model,
      parms = params
    )
  )

  return(out_df)
}

output <- ode(
  y = initial_state,
  times = times,
  func = seir_model,
  parms = params
)

results_long <- as.data.frame(output) %>%

  pivot_longer(
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

  left_join(
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

total_pages <- max(
  df$page
)

for (p in 1:total_pages) {

  page_data <-
    results_long %>%
    filter(page == p)

  p_plot <- ggplot(
    page_data,
    aes(
      x = time,
      y = Count,
      color = Compartment
    )
  ) +

    geom_line(
      linewidth = 1
    ) +

    scale_color_manual(
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

    facet_wrap(
      ~school_label,
      ncol = 3,
      scales = "free_y"
    ) +

    labs(
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

  print(p_plot)
}

net_graph <-
  graph_from_adjacency_matrix(
    A,
    mode = "undirected",
    weighted = TRUE
  )

df$node_strength <-
  strength(
    net_graph
  )

df$betweenness <-
  betweenness(
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

cat("\n========================================\n")
cat("NETWORK HUB AND PERIPHERAL SCHOOL\n")
cat("========================================\n")

cat(
  "Network Hub:",
  top_hub_name,
  "\n"
)

cat(
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

comparison_all <- bind_rows(

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

ggplot(
  comparison_all,
  aes(
    x = time,
    y = total_infectious,
    color = `Seed Type`
  )
) +

  geom_line(
    linewidth = 1.2
  ) +

  facet_wrap(
    ~Scenario,
    ncol = 3,
    scales = "free_y"
  ) +

  scale_color_manual(
    values = c(
      "Network Hub" = "firebrick",
      "Peripheral" = "darkblue"
    )
  ) +

  labs(
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

calculate_outbreak_measures <- function(
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
      function(x) {
        out_df[[x]][
          nrow(out_df)
        ]
      }
    )

  R_initial <-
    round(
      sapply(
        R_columns,
        function(x) {
          out_df[[x]][1]
        }
      )
    )

  new_cases_each_school <-
    pmax(
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

cat("\n\n")
cat("============================================================\n")
cat("TABLE 4.3: COMPARATIVE ANALYSIS OF HUB VS PERIPHERAL\n")
cat("SEEDING SCENARIOS\n")
cat("============================================================\n\n")

print(
  table_4_3,
  n = Inf
)

if (
  requireNamespace(
    "knitr",
    quietly = TRUE
  )
) {

  print(
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

write.csv(
  table_4_3,
  "output/tables/table_4_3_hub_vs_peripheral.csv",
  row.names = FALSE
)
