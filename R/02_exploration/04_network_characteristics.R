library(tidyverse)  # data wrangling and plots
library(igraph)  # network tools
library(ggraph)  # network plots

school_data  <- read_csv("data/processed/school_data_joined.csv", show_col_types = FALSE) %>%  # read the table
  rename_with(tolower)
network_data <- read_csv("data/processed/network_df.csv", show_col_types = FALSE) %>%  # read the table
  rename_with(tolower)

node_data <- school_data %>%
  mutate(
    unit_id = as.character(unit_id),
    herd_immunity = if_else(percent_immunized >= 0.95, "95%+ MMR", "Below 95% MMR"),
    exposure_halo = if_else(measles_exposure == 1, "Exposed School", "Unexposed School")
  ) %>%
  relocate(unit_id, .before = 1) %>%
  distinct(unit_id, .keep_all = TRUE)

edge_data <- network_data %>%
  mutate(
    unit_id_a = as.character(unit_id_a),
    unit_id_b = as.character(unit_id_b)
  ) %>%
  filter(unit_id_a %in% node_data$unit_id & unit_id_b %in% node_data$unit_id) %>%
  relocate(unit_id_a, unit_id_b, .before = 1)

net <- graph_from_data_frame(d = edge_data, vertices = node_data, directed = FALSE)  # build an undirected school network

n_nodes <- vcount(net)
n_edges <- ecount(net)

cat("Number of nodes (schools):", n_nodes, "\n")  # print a check
cat("Number of edges (connections):", n_edges, "\n")  # print a check

if (!"weight" %in% edge_attr_names(net)) {
  E(net)$weight <- E(net)$total_contacts
}

node_data <- node_data %>%
  mutate(
    degree      = degree(net),  # number of schools linked to this one
    strength    = strength(net, weights = E(net)$weight),  # sum of this school's edge weights
    betweenness = betweenness(net, weights = E(net)$weight),  # schools that sit on many shortest paths
    closeness   = closeness(net, weights = E(net)$weight),  # how near the school is to the rest of the network
    eigen_cent  = eigen_centrality(net, weights = E(net)$weight)$vector  # centrality through well-connected neighbours
  )

net_density      <- edge_density(net)
net_components   <- components(net)$no
net_diameter     <- diameter(net, weights = E(net)$weight)
net_avg_pathlen  <- mean_distance(net, weights = E(net)$weight)
net_assortativity <- assortativity_degree(net)  # whether similar-degree schools connect
net_transitivity  <- transitivity(net, type = "global")

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

print(global_summary)  # show the result

isolated_schools <- node_data %>% filter(strength == 0)
cat("Number of fully isolated schools (strength = 0):", nrow(isolated_schools), "\n")  # print a check

hub_school <- node_data %>%
  filter(strength == max(strength)) %>%
  slice(1)

peripheral_school <- node_data %>%
  filter(strength > 0) %>%
  filter(strength == min(strength)) %>%
  slice(1)

hub_peripheral_table <- bind_rows(  # combine these rows
  hub_school        %>% mutate(role = "Hub (highest strength)"),
  peripheral_school %>% mutate(role = "Peripheral (lowest strength, connected)")
) %>%
  select(role, unit_id, school_name, degree, strength, betweenness,
         closeness, total_students, percent_immunized)

print(hub_peripheral_table)  # show the result

node_characteristics_table <- node_data %>%
  select(unit_id, school_name, total_students, percent_immunized,
         degree, strength, betweenness, closeness, eigen_cent) %>%
  arrange(desc(strength))

print(head(node_characteristics_table, 10))  # show the result

p_geo_halo <- ggraph(net, layout = layout_geo) +  # place the schools for the plot

  geom_edge_link(
    color = "gray75",
    alpha = 0.28,
    width = 0.4
  ) +

  geom_node_point(  # one mark per school
    aes(
      fill  = herd_immunity,
      color = exposure_halo
    ),
    shape  = 21,
    size   = 6,
    stroke = 1.8
  ) +

  scale_fill_manual(  # set the legend colours or line types
    values = c(
      "Below 95% MMR" = "#E41A1C",
      "95%+ MMR"      = "#2CA02C"
    ),
    name = "Herd Immunity Status"
  ) +

  scale_color_manual(  # set the legend colours or line types
    values = c(
      "Exposed School"   = "#00BFFF",
      "Unexposed School" = "#444444"
    ),
    name = "Exposure Status"
  ) +

  geom_node_text(
    aes(label = name),
    color = "white",
    size = 1.8,
    fontface = "bold"
  ) +

  theme_graph(base_family = "sans") +
  labs(  # titles and axis labels
    title = "Geographic Network Plot (Spartanburg County)"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

print(p_geo_halo)  # show the result
