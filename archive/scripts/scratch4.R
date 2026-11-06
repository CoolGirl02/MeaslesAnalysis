df <- read.csv("data/raw/spartanburg_connectivity_data.csv")


library(dplyr)
#extract all the nodes
df_1 <- df %>%
  select(unit_id_A, school_name_A, school_long_A, school_lat_A, MMR_PR_A ) %>%
  distinct()

df_1 <- df_1 %>%
  rename(
    unit_Id = unit_id_A,
    school_name = school_name_A,
    school_long = school_long_A, 
    school_lat = school_lat_A,
    MMR_PR = MMR_PR_A
  )

df_2 <- df %>%
  select(unit_id_B, school_name_B, school_long_B, school_lat_B, MMR_PR_B) %>%
  distinct()

df_2 <- df_2 %>%
  rename(
    unit_Id = unit_id_B,
    school_name = school_name_B,
    school_long = school_long_B, 
    school_lat = school_lat_B,
    MMR_PR = MMR_PR_B
  )

final_schools_data <- bind_rows(df_1, df_2)
final_schools_data <- final_schools_data %>%
  distinct()



#build the network(creating the edges)
network_df <- df %>%
  select(unit_id_A, unit_id_B, total_contacts, connect_risk)%>%
  mutate(
    node1 = pmin(unit_id_A, unit_id_B),
    node2 = pmax(unit_id_A, unit_id_B)
  ) %>%
  # 2. Replace original columns with the standardized ones
  mutate(
    unit_id_A = node1,
    unit_id_B = node2
  ) %>%
  select(-node1, -node2)%>%
  distinct(unit_id_A, unit_id_B, .keep_all = TRUE)



duplicated_rows <- network_df %>%
  group_by(unit_id_A, unit_id_B) %>%
  filter(n() > 1) %>%
  ungroup()


#visualization
library(igraph)
edges <-network_df %>%
  select(unit_id_A, unit_id_B)

net <- graph_from_data_frame(edges, directed = FALSE)

# 1. Expand the plotting canvas (removes huge whitespace margins around the plot)
par(mar = c(0, 0, 0, 0))

# 2. Compute a layout that pushes nodes apart
l <- layout_with_fr(net)

# 3. Plot with aggressively smaller node and label sizes
plot(
  net,
  layout = l,
  
  # --- NODE SIZING ---
  vertex.size = 10,               # Was ~30. Set to 2 or 3 to shrink dots dramatically
  vertex.color = "lightgreen",
  vertex.frame.color = NA,       # Removes the thick black border around nodes
  
  # --- LABEL CONTROL ---
  #vertex.label = NA,             # OPTION A: Set to NA to hide labels entirely (best for huge graphs)
  vertex.label.cex = 0.5,      # OPTION B: Uncomment this if you MUST show text, but make it tiny
  # vertex.label.color = "black",
  # vertex.label.dist = 0.5,
  
  # --- EDGE SIZING ---
  edge.color = adjustcolor("gray60", alpha.f = 0.4), # Transparent thin lines
  edge.width = 0.5
)

#save the files
write.csv(final_schools_data, "data/processed/schools_data.csv", row.names = FALSE)
write.csv(network_df, "data/processed/network_df.csv", row.names = FALSE)



# 1. Ensure unit_Id is character in final_schools_data for exact matching
final_schools_data <- final_schools_data %>%
  mutate(unit_Id = as.character(unit_Id))

# 2. Extract coordinates in the exact vertex order of V(net)$name
layout_geo <- final_schools_data %>%
  slice(match(V(net)$name, unit_Id)) %>%
  select(school_long, school_lat) %>%
  as.matrix()


# 1. Jitter coordinates slightly so overlapping schools push apart visually
layout_geo_jittered <- layout_geo
layout_geo_jittered[, 1] <- jitter(layout_geo[, 1], factor = 3) # Nudge longitude
layout_geo_jittered[, 2] <- jitter(layout_geo[, 2], factor = 3) # Nudge latitude

# 2. Define node colors (Red = Below 95% MMR, Green = 95%+ MMR)
v_ids <- V(net)$name
mmr_vec <- final_schools_data$MMR_PR[match(v_ids, final_schools_data$unit_Id)]
node_colors <- ifelse(mmr_vec < 0.95, "#E41A1C", "#4DAF4A")

# 3. Plot clear geographic network
par(mar = c(1, 1, 2, 1))

plot(
  net,
  layout = layout_geo_jittered,
  rescale = TRUE,
  
  # --- NODE STYLING ---
  vertex.size = 10,                      # Slightly reduced to stop overlap
  vertex.color = node_colors,
  vertex.frame.color = "black",          # Solid black border separates stacked circles
  
  # --- LABELS INSIDE CIRCLES ---
  vertex.label = V(net)$name,
  vertex.label.dist = 0,                 # Centered inside node
  vertex.label.cex = 0.45,               # Smaller text so multi-digit IDs fit
  vertex.label.color = "white",
  vertex.label.font = 2,                 # Bold font
  
  # --- LIGHTEN EDGES ---
  edge.color = adjustcolor("gray70", alpha.f = 0.15), # Faded edges reduce background noise
  edge.width = 0.4,
  
  main = "Geographic Network Plot (Spartanburg County)"
)

# Add Legend
legend(
  "bottomright",
  legend = c("Below 95% MMR", "95%+ MMR"),
  fill = c("#E41A1C", "#4DAF4A"),
  bty = "n",
  title = "Herd Immunity Status"
)











