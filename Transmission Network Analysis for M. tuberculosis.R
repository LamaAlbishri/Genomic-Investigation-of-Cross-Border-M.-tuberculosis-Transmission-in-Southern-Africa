###############################################
# Transmission Network Analysis for M. tuberculosis
# This script processes SNP data, constructs transmission networks,
# and generates visualizations for cross-border transmission analysis.
###############################################

#############################
## Load required libraries ##
#############################

library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)
library(igraph)
library(ggmap)
library(geomtextpath)

#################################
## Data Import and Preparation ##
#################################

# Read SNP matrix and metadata
snp_matrix <- "results/transmission/aligned_pseudogenomes_masked_snps.csv"
metadata <- read_csv("sample_info.csv")

snp_matrix_df <- read_csv(snp_matrix, col_names = FALSE)
snp_matrix_df <- transpose(snp_matrix_df, make.names = 'X1')

# Process SNP matrix into pairwise comparisons
snp_matrix_df_decon <- data.frame(
  t(combn(names(snp_matrix_df), 2)), 
  dist = t(snp_matrix_df)[lower.tri(snp_matrix_df)]
)
colnames(snp_matrix_df_decon) <- c("Taxon1", "Taxon2", "dist")

# Add country information to SNP pairs
metadata_2 <- metadata %>% select(sample, Country)

snp_matrix_df_decon_country <- snp_matrix_df_decon %>% 
  left_join(metadata_2, by = c("Taxon1" = "sample")) %>%
  left_join(metadata_2, by = c("Taxon2" = "sample")) %>% 
  filter(!is.na(Country.x), !is.na(Country.y)) %>%
  mutate(
    transmission = case_when(
      Country.x == "Botswana" & Country.y == "Botswana" ~ "within_country",
      Country.x == "Namibia" & Country.y == "Namibia" ~ "within_country",
      TRUE ~ "cross_border"
    )
  )

##############################################
## Plot histogram of pairwise SNP distances ##
##############################################

snp_histogram <- ggplot(snp_matrix_df_decon, aes(x = dist)) + 
  geom_histogram(bins = 200, alpha = 0.8, position = "identity") + 
  theme_bw() + 
  labs(x = "Pairwise SNP distance", y = "Frequency")

###########################################
## Transmission Network Construction ##
###########################################

# Prepare node metadata
metadata_nodes <- metadata %>%
  mutate(Sample_ID = row_number()) %>% 
  select(Sample_ID, collection_date, sample, Country, District, 
         longitude, latitude, DR_type)

# Set SNP threshold and create edges
threshold <- 12

edges <- snp_matrix_df_decon_country %>%
  filter(dist <= threshold) %>% 
  select(Taxon1, Taxon2, dist, transmission) %>%
  left_join(metadata_nodes[, c(1, 3)], by = c("Taxon1" = "sample")) %>% 
  rename(from.Sample_ID = Sample_ID) %>%
  left_join(metadata_nodes[, c(1, 3)], by = c("Taxon2" = "sample")) %>% 
  select(from.Sample_ID, to.Sample_ID = Sample_ID, dist)

# Create initial network
network <- tbl_graph(nodes = metadata_nodes, edges = edges, directed = TRUE)
network_components <- components(network)

metadata_nodes_networked <- metadata_nodes %>% 
  mutate(network_id = network_components$membership)

# Filter out singletons
metadata_nodes_networked_filtered <- metadata_nodes_networked %>% 
  group_by(network_id) %>%
  filter(n() > 1) %>%
  ungroup() %>% 
  mutate(id = row_number()) %>% 
  select(id, sample, longitude, latitude, District, Country, 
         collection_date, network_id, DR_type)

# Create filtered edges
edges_filtered <- snp_matrix_df_decon_country %>%
  filter(dist <= threshold) %>% 
  select(Taxon1, Taxon2, dist, transmission) %>%
  left_join(metadata_nodes_networked_filtered[, c(1, 2)], by = c("Taxon1" = "sample")) %>% 
  rename(from.id = id) %>%
  left_join(metadata_nodes_networked_filtered[, c(1, 2)], by = c("Taxon2" = "sample")) %>%
  select(from.id, to.id = id, dist, transmission)

# Create filtered network
network_filtered <- tbl_graph(nodes = metadata_nodes_networked_filtered, 
                              edges = edges_filtered, directed = TRUE)

###########################################
## Transmission Network Visualization ##
###########################################

# Plot network
network_plot <- ggraph(network_filtered, layout = "nicely") + 
  geom_edge_link(aes(label = dist, color = transmission),
                 edge_width = 1, 
                 angle_calc = "along", 
                 label_dodge = unit(2.5, 'mm'), 
                 label_size = 3) + 
  geom_node_point(aes(colour = District, shape = Country), size = 3) +
  theme_graph() + 
  labs(colour = "District", shape = "Country") +
  facet_nodes(~network_id, scales = "free")

###############################################
## Transmission Threshold Analysis ##
###############################################

# Function to count transmission types by threshold
get_num_transmission <- function(threshold) {
  edges <- snp_matrix_df_decon_country %>%
    filter(dist <= threshold) %>% 
    select(Taxon1, Taxon2, dist, transmission)
  
  data.frame(
    threshold = threshold,
    cross_border = sum(edges$transmission == 'cross_border'),
    within_country = sum(edges$transmission == 'within_country')
  )
}

# Calculate results for thresholds 1-50
results <- map_dfr(1:50, get_num_transmission) %>%
  mutate(total = cross_border + within_country)

# Plot transmission links by threshold
transmission_plot <- ggplot(results, aes(x = threshold)) +
  geom_line(aes(y = cross_border, color = "Cross-Border")) +
  geom_line(aes(y = within_country, color = "Within-Country")) +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "SNP Threshold", y = "Number of Transmission Links",
       color = "Transmission Type") +
  theme_minimal()

###############################################
## Geospatial Visualization ##
###############################################

# Get map data for Botswana-Namibia region
BotsNam_coords <- c(left = 10, bottom = -29, right = 29.5, top = -16)
register_stadiamaps(key = "your-key-here") # Replace with your actual key

BotsNam_coords_stamenmap <- get_stadiamap(BotsNam_coords, zoom = 6, 
                                          maptype = "stamen_terrain", force = TRUE)

# Create map with sample locations
BotsNam_map_TB <- ggmap(BotsNam_coords_stamenmap) + 
  geom_point(data = metadata, 
             aes(x = longitude, y = latitude, color = DR_type, shape = Country), 
             alpha = 0.8) +
  labs(x = "Longitude", y = "Latitude")

###############################################
## Cross-Border Transmission Focus ##
###############################################

# Focus on the cross-border transmission network (network_id 55 in your example)
cross_border_nodes <- metadata_nodes_networked_filtered %>% 
  filter(network_id == 55) %>%  # Replace with relevant network_id
  mutate(id = row_number())

cross_border_edges <- snp_matrix_df_decon %>%
  filter(dist <= threshold) %>% 
  select(Taxon1, Taxon2, dist) %>%
  left_join(cross_border_nodes[, c(1, 2)], by = c("Taxon1" = "sample")) %>% 
  rename(from.id = id) %>%
  left_join(cross_border_nodes[, c(1, 2)], by = c("Taxon2" = "sample")) %>%
  select(from.id, to.id = id, dist) %>% 
  filter(!is.na(from.id), !is.na(to.id))

cross_border_network <- tbl_graph(nodes = cross_border_nodes, 
                                  edges = cross_border_edges, directed = TRUE)

# Plot focused cross-border transmission
cross_border_plot <- ggraph(cross_border_network, 
                            x = cross_border_nodes$longitude,
                            y = cross_border_nodes$latitude) +
  geom_edge_fan(aes(label = dist), alpha = 0.3) +
  geom_node_point(aes(color = DR_type), size = 3) +
  labs(color = "Drug Resistance Type") +
  theme_bw()

# Print summary statistics
cat("Cross-border transmission links at threshold", threshold, ":", 
    sum(edges_filtered$transmission == 'cross_border'), "\n")
cat("Within-country transmission links:", 
    sum(edges_filtered$transmission == 'within_country'), "\n")