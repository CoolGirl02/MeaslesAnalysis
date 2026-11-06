# Geographic Network Plot for Spartanburg County


# Load packages
library(tidyverse)
library(igraph)
library(ggraph)

# Import data and make all column names lowercase
school_data  <- read_csv("data/cleanned/school_data_joined.csv", show_col_types = FALSE) %>% 
  rename_with(tolower)

network_data <- read_csv("data/cleanned/network_df.csv", show_col_types = FALSE) %>% 
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





