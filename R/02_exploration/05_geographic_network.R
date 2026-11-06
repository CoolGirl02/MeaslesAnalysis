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
