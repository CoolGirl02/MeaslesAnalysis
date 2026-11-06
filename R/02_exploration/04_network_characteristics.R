# Geographic Network Plot for Spartanburg County
# Load packages
library(tidyverse)
library(igraph)
library(ggraph)

# Import data and make all column names lowercase
school_data  <- read_csv("data/processed/school_data_joined.csv", show_col_types = FALSE) %>% 
  rename_with(tolower)
network_data <- read_csv("data/processed/network_df.csv", show_col_types = FALSE) %>% 
  rename_with(tolower)

# Clean school data (nodes)
node_data <- school_data %>%
  mutate(
    unit_id = as.character(unit_id),
    herd_immunity = if_else(percent_immunized >= 0.95, "95%+ MMR", "Below 95% MMR"),
    exposure_halo = if_else(measles_exposure == 1, "Exposed School", "Unexposed School")
  ) %>%
  relocate(unit_id, .before = 1) %>% 
  distinct(unit_id, .keep_all = TRUE)

# Clean network connections (edges)
edge_data <- network_data %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  ) %>%
  filter(unit_id_a %in% node_data$unit_id & unit_id_b %in% node_data$unit_id) %>%
  relocate(unit_id_a, unit_id_b, .before = 1)

# Create network graph object
net <- graph_from_data_frame(d = edge_data, vertices = node_data, directed = FALSE)


# ==============================================================================
# NETWORK CHARACTERISTICS ANALYSIS
# ==============================================================================
# This section computes the structural characteristics of the school contact
# network: basic size, node-level centrality measures, global network
# statistics, and the extraction of the hub and peripheral schools.

# ------------------------------------------------------------------------------
# 1. Basic network size
# ------------------------------------------------------------------------------
n_nodes <- vcount(net)   # number of schools (nodes)
n_edges <- ecount(net)   # number of distinct school-pair connections (edges)

cat("Number of nodes (schools):", n_nodes, "\n")
cat("Number of edges (connections):", n_edges, "\n")

# ------------------------------------------------------------------------------
# 2. Node-level centrality measures
# ------------------------------------------------------------------------------
# Ensure the edge weight column is named "weight" for igraph's weighted
# functions to pick it up automatically. Adjust "total_contacts" if your
# edge_data column is named differently.
if (!"weight" %in% edge_attr_names(net)) {
  E(net)$weight <- E(net)$total_contacts
}

node_data <- node_data %>%
  mutate(
    degree      = degree(net),                                   # number of distinct connections
    strength    = strength(net, weights = E(net)$weight),         # weighted degree (k_i)
    betweenness = betweenness(net, weights = E(net)$weight),       # bridging importance
    closeness   = closeness(net, weights = E(net)$weight),         # average reach / speed to others
    eigen_cent  = eigen_centrality(net, weights = E(net)$weight)$vector  # connected-to-well-connected
  )

# ------------------------------------------------------------------------------
# 3. Global (whole-network) characteristics
# ------------------------------------------------------------------------------
net_density      <- edge_density(net)                                    # proportion of possible edges present
net_components   <- components(net)$no                                   # number of separate connected clusters
net_diameter     <- diameter(net, weights = E(net)$weight)                # longest shortest-path (worst-case reach)
net_avg_pathlen  <- mean_distance(net, weights = E(net)$weight)           # typical shortest-path length
net_assortativity <- assortativity_degree(net)                           # do hubs link to other hubs?
net_transitivity  <- transitivity(net, type = "global")                  # clustering / triangle density

global_summary <- tibble(
  Characteristic = c(
    "Number of nodes (n)",
    "Number of edges (m)",
    "Edge density",
    "Number of connected components",
    "Network diameter",
    "Average path length",
    "Degree assortativity",
    "Global clustering coefficient (transitivity)"
  ),
  Value = c(
    n_nodes,
    n_edges,
    round(net_density, 4),
    net_components,
    round(net_diameter, 2),
    round(net_avg_pathlen, 2),
    round(net_assortativity, 4),
    round(net_transitivity, 4)
  )
)

print(global_summary)

# ------------------------------------------------------------------------------
# 4. Extracting the hub and peripheral schools
# ------------------------------------------------------------------------------
# Hub: highest node strength (k_i) -- the school that exports/receives the
# most contact-weighted connection to the rest of the network.
# Peripheral: lowest node strength AMONG SCHOOLS THAT ARE ACTUALLY CONNECTED
# (strength > 0) -- excludes fully isolated schools, which would otherwise
# trivially "win" the peripheral label with zero network exposure.

isolated_schools <- node_data %>% filter(strength == 0)
cat("Number of fully isolated schools (strength = 0):", nrow(isolated_schools), "\n")

hub_school <- node_data %>%
  filter(strength == max(strength)) %>%
  slice(1)   # break ties by taking the first match

peripheral_school <- node_data %>%
  filter(strength > 0) %>%
  filter(strength == min(strength)) %>%
  slice(1)   # break ties by taking the first match

hub_peripheral_table <- bind_rows(
  hub_school        %>% mutate(role = "Hub (highest strength)"),
  peripheral_school %>% mutate(role = "Peripheral (lowest strength, connected)")
) %>%
  select(role, unit_id, school_name, degree, strength, betweenness,
         closeness, total_students, percent_immunized)

print(hub_peripheral_table)

# ------------------------------------------------------------------------------
# 5. Full node-level characteristics table (for appendix / reference)
# ------------------------------------------------------------------------------
node_characteristics_table <- node_data %>%
  select(unit_id, school_name, total_students, percent_immunized,
         degree, strength, betweenness, closeness, eigen_cent) %>%
  arrange(desc(strength))

# View top 10 most-connected schools
print(head(node_characteristics_table, 10))

# Optional: write full table out to CSV for your appendix
# write_csv(node_characteristics_table, "output/tables/network_characteristics_table.csv")


# ==============================================================================
# NETWORK VISUALIZATION
# ==============================================================================

# Plot network
p_geo_halo <- ggraph(net, layout = layout_geo) +
  
  # Draw connection lines between schools
  geom_edge_link(
    color = "gray75",
    alpha = 0.28,
    width = 0.4
  ) +
  
  # Draw school nodes
  geom_node_point(
    aes(
      fill  = herd_immunity,
      color = exposure_halo
    ),
    shape  = 21,
    size   = 6,
    stroke = 1.8
  ) +
  
  # Colors for immunity status (inside point)
  scale_fill_manual(
    values = c(
      "Below 95% MMR" = "#E41A1C",
      "95%+ MMR"      = "#2CA02C"
    ),
    name = "Herd Immunity Status"
  ) +
  
  # Colors for exposure status (outside border)
  scale_color_manual(
    values = c(
      "Exposed School"   = "#00BFFF",
      "Unexposed School" = "#444444"
    ),
    name = "Exposure Status"
  ) +
  
  # Add school ID numbers on top
  geom_node_text(
    aes(label = name),
    color = "white",
    size = 1.8,
    fontface = "bold"
  ) +
  
  # Add theme and titles
  theme_graph(base_family = "sans") +
  labs(
    title = "Geographic Network Plot (Spartanburg County)"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Display plot
print(p_geo_halo)