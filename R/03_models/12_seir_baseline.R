# ============================================================
# 1. LOAD LIBRARIES
# ============================================================

library(tidyverse)
library(deSolve)


# ============================================================
# 2. LOAD SCHOOL DATA
# ============================================================

df <- read.csv("data/processed/school_data_joined.csv")

names(df) <- tolower(names(df))

# Ensure that every school has an ID
if (!"unit_id" %in% names(df)) {
  df$unit_id <- as.character(seq_len(nrow(df)))
} else {
  df$unit_id <- as.character(df$unit_id)
}

# Check that the required columns exist
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

# Prepare school data
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


# ============================================================
# 3. FIND GLOBAL ACADEMY OF SC IN THE DATASET
# ============================================================

seed_school_name <- "Global Academy of SC"

# Standardise names to make matching case-insensitive
school_names_clean <- trimws(tolower(as.character(df$school_name)))

seed_index <- which(
  school_names_clean == trimws(tolower(seed_school_name))
)

# If the exact name is not found, display possible matches
if (length(seed_index) == 0) {
  
  possible_matches <- df$school_name[
    grepl(
      "global academy",
      as.character(df$school_name),
      ignore.case = TRUE
    )
  ]
  
  if (length(possible_matches) > 0) {
    cat("Possible matching school names:\n")
    print(possible_matches)
  } else {
    cat("No school containing 'Global Academy' was found.\n")
    cat("Available school names include:\n")
    print(head(df$school_name, 20))
  }
  
  stop("Seed school not found. Check the school name in the dataset.")
}

# Ensure that exactly one school matches
if (length(seed_index) > 1) {
  stop("More than one school matches Global Academy of SC.")
}

cat("Seed school found:", as.character(df$school_name[seed_index]), "\n")


# ============================================================
# 4. LOAD NETWORK DATA
# ============================================================

net_df <- read.csv("data/processed/network_df.csv")

names(net_df) <- tolower(names(net_df))

required_network_columns <- c(
  "unit_id_a",
  "unit_id_b",
  "total_contacts"
)

missing_network_columns <- setdiff(
  required_network_columns,
  names(net_df)
)

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


# ============================================================
# 5. BUILD THE ADJACENCY MATRIX
# ============================================================

id_map <- setNames(seq_len(n_schools), df$unit_id)

A <- matrix(
  0,
  nrow = n_schools,
  ncol = n_schools
)

for (k in seq_len(nrow(net_df))) {
  
  i <- unname(id_map[net_df$unit_id_a[k]])
  j <- unname(id_map[net_df$unit_id_b[k]])
  
  if (
    length(i) == 1 &&
    length(j) == 1 &&
    !is.na(i) &&
    !is.na(j) &&
    i != j
  ) {
    
    w <- net_df$total_contacts[k]
    
    if (is.na(w) || !is.finite(w) || w < 0) {
      next
    }
    
    # Set the connection weight between the two schools
    A[i, j] <- w
    A[j, i] <- w
  }
}

# Row-normalise the adjacency matrix
row_sums <- rowSums(A)

row_sums[row_sums == 0] <- 1

A <- A / row_sums


# ============================================================
# 6. SET MODEL PARAMETERS
# ============================================================

gamma_rate <- 1 / 7
sigma_rate <- 1 / 10
beta_rate <- 15 * gamma_rate
alpha <- 0.05

times <- seq(0, 200, by = 1)

# Total population of each school
N <- df$total_students

# Initially vaccinated students
R0_count <- round(N * df$vax_rate)


# ============================================================
# 7. INITIAL CONDITIONS
#    ONE INFECTIOUS STUDENT AT GLOBAL ACADEMY OF SC
# ============================================================

I0_count <- rep(0, n_schools)

# Only the selected seed school starts with an infectious student
I0_count[seed_index] <- 1

# No initially exposed students
E0_count <- rep(0, n_schools)

# Initially susceptible students
S0_count <- pmax(
  0,
  N - I0_count - E0_count - R0_count
)

cat("Initial infectious students:", sum(I0_count), "\n")
cat(
  "Infection starts at:",
  as.character(df$school_name[seed_index]),
  "\n"
)


# ============================================================
# 8. CREATE INITIAL STATE
# ============================================================

initial_state <- c(
  S0_count,
  E0_count,
  I0_count,
  R0_count
)

names(initial_state) <- c(
  paste0("S_", seq_len(n_schools)),
  paste0("E_", seq_len(n_schools)),
  paste0("I_", seq_len(n_schools)),
  paste0("R_", seq_len(n_schools))
)


# ============================================================
# 9. MODEL PARAMETERS LIST
# ============================================================

params <- list(
  beta = beta_rate,
  sigma = sigma_rate,
  gamma = gamma_rate,
  alpha = alpha,
  N = N,
  A = A,
  n = n_schools
)


# ============================================================
# 10. NETWORK-BASED SEIR MODEL
# ============================================================

seir_model <- function(time, state, parameters) {
  
  with(as.list(parameters), {
    
    S <- state[seq_len(n)]
    
    E <- state[n + seq_len(n)]
    
    I <- state[2 * n + seq_len(n)]
    
    R <- state[3 * n + seq_len(n)]
    
    # Proportion of infectious students in each school
    I_prop <- I / N
    
    # Force of infection:
    # within-school transmission plus between-school transmission
    lambda <- (beta * I_prop) +
      (alpha * beta * as.vector(A %*% I_prop))
    
    # SEIR differential equations
    dS <- -lambda * S
    
    dE <- lambda * S - sigma * E
    
    dI <- sigma * E - gamma * I
    
    dR <- gamma * I
    
    return(list(c(dS, dE, dI, dR)))
  })
}


# ============================================================
# 11. RUN THE SEIR SIMULATION
# ============================================================

output <- ode(
  y = initial_state,
  times = times,
  func = seir_model,
  parms = params
)

cat("SEIR simulation completed successfully.\n")


# ============================================================
# 12. PREPARE RESULTS FOR PLOTTING
# ============================================================

results_long <- as.data.frame(output) %>%
  pivot_longer(
    cols = -time,
    names_to = "var",
    values_to = "Count"
  ) %>%
  mutate(
    Compartment = str_extract(var, "^[SEIR]"),
    school_id = as.numeric(str_extract(var, "\\d+"))
  ) %>%
  left_join(df, by = "school_id") %>%
  mutate(
    school_label = paste0(
      school_name,
      "\n(Vax: ",
      round(vax_rate * 100, 1),
      "%, Pop: ",
      total_students,
      ")"
    )
  )


# ============================================================
# 13. PLOT SEIR DYNAMICS FOR ALL SCHOOLS
# ============================================================

total_pages <- max(df$page)

for (p in seq_len(total_pages)) {
  
  page_data <- results_long %>%
    filter(page == p)
  
  p_plot <- ggplot(
    page_data,
    aes(
      x = time,
      y = Count,
      color = Compartment
    )
  ) +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = c(
        "S" = "blue",
        "E" = "orange",
        "I" = "red",
        "R" = "green"
      ),
      breaks = c("S", "E", "I", "R"),
      labels = c(
        "S" = "Susceptible (S)",
        "E" = "Exposed (E)",
        "I" = "Infectious (I)",
        "R" = "Recovered/Vaccinated (R)"
      )
    ) +
    facet_wrap(
      ~school_label,
      ncol = 3,
      scales = "free_y"
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1)
    ) +
    
  labs(
    title = paste0(
      "SEIR Measles Simulations - Seed: ",
      as.character(df$school_name[seed_index])
    ),
    x = "Days",
    y = "Number of Students",
    color = "Compartment"
  ) +
    theme_minimal() +
    theme(
      strip.text = element_text(
        size = 8,
        face = "bold"
      ),
      legend.position = "bottom"
    )

  
  print(p_plot)
}

cat("All school plots have been generated.\n")