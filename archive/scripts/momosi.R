library(tidyverse)
library(deSolve)

df <- read.csv("data/processed/school_data_joined.csv")

names(df) <- tolower(names(df))

if (!"unit_id" %in% names(df)) {
  df$unit_id <- as.character(seq_len(nrow(df)))
} else {
  df$unit_id <- as.character(df$unit_id)
}

required_columns <- c(
  "unit_id",
  "school_name",
  "total_students",
  "percent_immunized"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required school-data columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

df <- df %>%
  mutate(
    vax_rate = if_else(
      percent_immunized > 1,
      percent_immunized / 100,
      percent_immunized
    ),
    school_id = row_number()
  )

n_schools <- nrow(df)

seed_school_name <- "Global Academy of SC"
school_names_clean <- trimws(tolower(as.character(df$school_name)))

seed_index <- which(
  school_names_clean == trimws(tolower(seed_school_name))
)

if (length(seed_index) == 0) {
  possible_matches <- df$school_name[
    grepl("global academy", as.character(df$school_name), ignore.case = TRUE)
  ]
  if (length(possible_matches) > 0) {
    print(possible_matches)
  }
  stop("Seed school not found. Check the school name in the dataset.")
}

if (length(seed_index) > 1) {
  stop("More than one school matches Global Academy of SC.")
}

cat("Infection starts at:", df$school_name[seed_index], "\n")

net_df <- read.csv("data/processed/network_df.csv")
names(net_df) <- tolower(names(net_df))

required_network_columns <- c("unit_id_a", "unit_id_b", "total_contacts")
missing_network_columns <- setdiff(required_network_columns, names(net_df))

if (length(missing_network_columns) > 0) {
  stop(
    paste(
      "Missing required network-data columns:",
      paste(missing_network_columns, collapse = ", ")
    )
  )
}

net_df <- net_df %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  )

id_map <- setNames(seq_len(n_schools), df$unit_id)
A_raw <- matrix(0, nrow = n_schools, ncol = n_schools)

for (k in seq_len(nrow(net_df))) {
  i <- unname(id_map[net_df$unit_id_a[k]])
  j <- unname(id_map[net_df$unit_id_b[k]])

  if (length(i) == 1 && length(j) == 1 && !is.na(i) && !is.na(j) && i != j) {
    w <- net_df$total_contacts[k]
    if (is.na(w) || !is.finite(w) || w < 0) next
    A_raw[i, j] <- w
    A_raw[j, i] <- w
  }
}

df$connectivity <- rowSums(A_raw)

match_six_schools <- function(patterns, df, fallback_ids) {
  matched <- unlist(lapply(patterns, function(p) {
    idx <- which(grepl(p, df$school_name, ignore.case = TRUE))
    if (length(idx) > 0) idx[1] else NULL
  }))
  matched <- unique(matched)
  if (length(matched) < 6) {
    fill <- setdiff(fallback_ids, matched)
    matched <- c(matched, head(fill, 6 - length(matched)))
  }
  return(head(matched, 6))
}

lowest_vax_ids <- df %>%
  arrange(vax_rate, school_name) %>%
  slice_head(n = 6) %>%
  pull(school_id)

highest_vax_ids <- df %>%
  arrange(desc(vax_rate), school_name) %>%
  slice_head(n = 6) %>%
  pull(school_id)

intermediate_patterns <- c(
  "campobello",
  "cooley springs",
  "fairforest middle",
  "byrnes freshman",
  "sugar ridge",
  "covenant"
)
fallback_intermediate <- df %>%
  filter(vax_rate >= 0.75, vax_rate < 0.95) %>%
  pull(school_id)
intermediate_ids <- match_six_schools(intermediate_patterns, df, fallback_intermediate)

highest_conn_patterns <- c(
  "boiling springs elementary",
  "boiling springs middle",
  "libertas",
  "mccracken",
  "oakland",
  "spartanburg high"
)
fallback_highest_conn <- df %>%
  arrange(desc(connectivity)) %>%
  pull(school_id)
highest_connectivity_ids <- match_six_schools(highest_conn_patterns, df, fallback_highest_conn)

lowest_conn_patterns <- c(
  "abner creek academy",
  "pine st",
  "blind elem",
  "blind high",
  "deaf high",
  "cleveland academy"
)
fallback_lowest_conn <- df %>%
  arrange(connectivity) %>%
  pull(school_id)
lowest_connectivity_ids <- match_six_schools(lowest_conn_patterns, df, fallback_lowest_conn)

spectral_radius <- max(Mod(eigen(A_raw, only.values = TRUE)$values))
if (!is.finite(spectral_radius) || spectral_radius <= 0) {
  stop("The network spectral radius is zero or invalid.")
}
A <- A_raw / spectral_radius

gamma_rate <- 1 / 7
sigma_rate <- 1 / 10
beta_rate  <- 15 * gamma_rate
alpha      <- 0.05

times <- seq(0, 200, by = 1)
N <- df$total_students

R0_count <- round(N * df$vax_rate)

I0_count <- rep(0, n_schools)
I0_count[seed_index] <- 1

E0_count <- rep(0, n_schools)

S0_count <- pmax(0, N - I0_count - E0_count - R0_count)

initial_state <- c(S0_count, E0_count, I0_count, R0_count)

names(initial_state) <- c(
  paste0("S_", seq_len(n_schools)),
  paste0("E_", seq_len(n_schools)),
  paste0("I_", seq_len(n_schools)),
  paste0("R_", seq_len(n_schools))
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
    S <- state[seq_len(n)]
    E <- state[n + seq_len(n)]
    I <- state[2 * n + seq_len(n)]
    R <- state[3 * n + seq_len(n)]

    I_prop <- I / N

    lambda <- beta * I_prop + alpha * beta * as.vector(A %*% I_prop)

    dS <- -lambda * S
    dE <-  lambda * S - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I

    return(list(c(dS, dE, dI, dR)))
  })
}

output <- ode(
  y     = initial_state,
  times = times,
  func  = seir_model,
  parms = params
)

cat("SEIR simulation finished successfully.\n")

results_long <- as.data.frame(output) %>%
  pivot_longer(
    cols      = -time,
    names_to  = "var",
    values_to = "Count"
  ) %>%
  mutate(
    Compartment = str_extract(var, "^[SEIR]"),
    school_id   = as.numeric(str_extract(var, "\\d+"))
  ) %>%

  filter(Compartment %in% c("S", "E", "I")) %>%
  left_join(df, by = "school_id") %>%
  mutate(
    school_label = paste0(
      school_name,
      "\n(Vax: ", round(vax_rate * 100, 1), "%, Pop: ", total_students, ")"
    ),
    Compartment = factor(Compartment, levels = c("S", "E", "I"))
  )

analysis_groups <- list(
  "Six Schools with the Highest Vaccination Coverage" = highest_vax_ids,
  "Six Schools with the Lowest Vaccination Coverage"  = lowest_vax_ids,
  "Six Intermediate-Coverage Schools"                 = intermediate_ids,
  "Six Schools with the Highest Network Connectivity" = highest_connectivity_ids,
  "Six Schools with the Lowest Network Connectivity"  = lowest_connectivity_ids
)

for (group_name in names(analysis_groups)) {

  selected_ids <- analysis_groups[[group_name]]

  plot_data <- results_long %>%
    filter(school_id %in% selected_ids) %>%
    mutate(plot_order = match(school_id, selected_ids)) %>%
    arrange(plot_order)

  p_plot <- ggplot(
    plot_data,
    aes(
      x        = time,
      y        = Count,
      color    = Compartment,
      linetype = Compartment
    )
  ) +
    geom_line(linewidth = 0.85) +

    scale_color_manual(
      values = c(
        "S" = "blue",
        "E" = "orange",
        "I" = "red"
      ),
      labels = c(
        "S" = "Susceptible (S)",
        "E" = "Exposed (E)",
        "I" = "Infectious (I)"
      )
    ) +

    scale_linetype_manual(
      values = c(
        "S" = "solid",
        "E" = "dashed",
        "I" = "dotted"
      ),
      labels = c(
        "S" = "Susceptible (S)",
        "E" = "Exposed (E)",
        "I" = "Infectious (I)"
      )
    ) +

    facet_wrap(
      ~school_label,
      ncol   = 3,
      scales = "free_y"
    ) +

    scale_y_continuous(
      labels = scales::label_comma(),
      expand = expansion(mult = c(0, 0.05))
    ) +

    labs(
      title = paste0(
        group_name,
        "\nInfection Introduced at ", df$school_name[seed_index]
      ),
      x        = "Days",
      y        = "Number of Students",
      color    = "Compartment",
      linetype = "Compartment"
    ) +

    theme_minimal() +
    theme(
      strip.text      = element_text(size = 8.5, face = "bold"),
      legend.position = "bottom",
      plot.title      = element_text(size = 11, face = "bold")
    )

  print(p_plot)
}

cat("\nAll 5 thesis-matching SEI linear plots generated successfully.\n")
