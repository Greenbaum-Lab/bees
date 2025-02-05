library(ape)
library(dplyr)
library(ggplot2)
library(geiger)
library(nlme)
library(phytools)
library(picante)
library(svglite)
library(stringr)
library(tidyr)
library(tibble)
library(purrr)




################################################################################
# Directionality Analysis in 4D Space
# Computation of phenotypic diversification rate for species.
# This script integrates phenotypic distances calculated from four principal
# components (PCs) with phylogenetic branch lengths, replicates steps based on
# branch lengths, and then adjusts step values based on a branching-time split.
#################################################################################

# Load data
lineage_df_all <- read.csv("reconstruced_data.csv")
phylo_tree <- read.newick("tree.nwk")

###############################
## User–Defined Parameters
###############################
plotTree(phylo_tree); nodelabels()
# Define which groups should be classified as "super"
super_groups <- c("api", "bom", "mel")

# Define the root (starting) node for tree paths
root_node <- 81

# Replication multiplier used in the replication function
rep_factor <- 10

# Node used for splitting the groups being compared
split_node <- "103"

###############################
## Data Preparation
###############################

# Function to arrange the data by species and add a combined group indicator.
get_combined_pcs_df <- function(df, super_groups) {
  df %>%
    arrange(species) %>%
    group_by(species) %>%
    mutate(combined_group = ifelse(group %in% super_groups, "super", "simple")) %>%
    ungroup()
}

# 'lineage_df_all' is assumed to be available in the workspace (contains PC1-4, species, group)
combined_pcs_df <- get_combined_pcs_df(lineage_df_all, super_groups)

###############################
## Distance Calculations in 4D
###############################

# (A) Euclidean distance function between two points in 4D space.
calc_distance <- function(point1, point2) {
  sqrt(sum((point1 - point2)^2))
}

# (B) Compute sequential distances between rows.
calc_sequential_distances <- function(species_data) {
  species_data %>%
    mutate(
      Distance = sqrt((PC1 - lead(PC1))^2 +
                        (PC2 - lead(PC2))^2 +
                        (PC3 - lead(PC3))^2 +
                        (PC4 - lead(PC4))^2),
      FromRow = row_number(),
      ToRow = row_number() + 1
    ) %>%
    filter(!is.na(Distance)) %>%
    dplyr::select(FromRow, ToRow, Distance, group, combined_group)
}

# (C) Apply the sequential distance calculation by species and add a step index.
get_distance_results <- function(df) {
  df %>%
    group_by(species) %>%
    do(calc_sequential_distances(.)) %>%
    ungroup() %>%
    group_by(species) %>%
    mutate(step_seq = row_number()) %>%
    ungroup()
}

# Compute the sequential distance results from the combined PCs data.
distance_results <- get_distance_results(combined_pcs_df)

###############################
## Tree Path Calculations
###############################

# (A) Prepare the combined edge matrix from the phylogenetic tree.
get_combined_edge_matrix <- function(tree) {
  edge_matrix <- tree$edge
  edge_lengths <- tree$edge.length
  combined_edge_matrix <- cbind(edge_matrix, edge_lengths)
  colnames(combined_edge_matrix) <- c("Parent", "Child", "Edge_Length")
  combined_edge_matrix
}

combined_edge_matrix <- get_combined_edge_matrix(phylo_tree)

# (B) For each species (tip), compute its evolutionary path (set of edges)
#     from the specified root node to the tip.
calculate_path <- function(species_index, tree, combined_edge_matrix, root_node) {
  # Get the node path from the specified root to the species tip.
  path <- nodepath(tree, from = root_node, to = species_index)
  
  # For each consecutive pair in the path, find the matching row in the edge matrix.
  path_rows <- sapply(1:(length(path) - 1), function(j) {
    which(combined_edge_matrix[, "Parent"] == path[j] & 
            combined_edge_matrix[, "Child"] == path[j + 1])
  })
  
  # Extract and transpose the edge lengths, naming the columns as Edge_1, Edge_2, ...
  path_edge_lengths <- combined_edge_matrix[path_rows, "Edge_Length"]
  setNames(as.data.frame(t(path_edge_lengths)), 
           paste("Edge", seq_along(path_edge_lengths), sep = "_"))
}

# (C) Compute the evolutionary path for each species (tip) and store in a list.
get_species_paths_list <- function(tree, combined_edge_matrix, root_node) {
  lapply(seq_along(tree$tip.label), function(idx) {
    calculate_path(species_index = idx, tree = tree, 
                   combined_edge_matrix = combined_edge_matrix, 
                   root_node = root_node)
  })
}

species_paths_list <- get_species_paths_list(phylo_tree, combined_edge_matrix, root_node)

# (D) Standardize the number of edge columns across species.
standardize_species_paths <- function(species_paths_list) {
  max_path_length <- max(sapply(species_paths_list, ncol))
  col_names <- paste("Edge", 1:max_path_length, sep = "_")
  
  species_paths_list <- lapply(species_paths_list, function(df) {
    num_missing_cols <- max_path_length - ncol(df)
    if (num_missing_cols > 0) {
      missing_cols <- setNames(matrix(NA, nrow = 1, ncol = num_missing_cols),
                               paste("Edge", (ncol(df) + 1):max_path_length, sep = "_"))
      df <- cbind(df, missing_cols)
    }
    # Ensure consistent column names and order.
    missing_cols <- setdiff(col_names, colnames(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- matrix(NA, nrow = nrow(df), ncol = length(missing_cols))
    }
    df[, col_names, drop = FALSE]
  })
  
  list(species_paths_list = species_paths_list, col_names = col_names)
}

standardized_paths <- standardize_species_paths(species_paths_list)
species_paths_list <- standardized_paths$species_paths_list

# (E) Combine the species paths into one data frame and convert to long format.
get_species_edge_df <- function(species_paths_list, tree) {
  all_species_paths_df <- do.call(rbind, species_paths_list)
  rownames(all_species_paths_df) <- tree$tip.label
  
  species_edge_df <- all_species_paths_df %>%
    tibble::rownames_to_column(var = "species") %>%
    pivot_longer(
      cols = starts_with("Edge_"),
      names_to = "edge",
      values_to = "edge_length"
    ) %>%
    filter(!is.na(edge_length)) %>%  # Remove missing edge lengths.
    group_by(species) %>%
    mutate(step_seq = row_number()) %>%
    ungroup()
  
  species_edge_df
}

species_edge_df <- get_species_edge_df(species_paths_list, phylo_tree)

###############################
## Combine Phenotypic & Phylogenetic Data
###############################

# Merge the sequential phenotypic distance results with the corresponding edge data.
combine_distance_edge_data <- function(distance_results, species_edge_df) {
  left_join(distance_results, species_edge_df, by = c("species", "step_seq")) %>%
    mutate(step = Distance / edge_length)
}

combined_results <- combine_distance_edge_data(distance_results, species_edge_df)

# Replicate each row of the data frame based on the product of edge_length and rep_factor.
# This allows a consistent number of data points for each species, 
# regrdless of the amount of edges it has.
replicate_steps_with_group <- function(species_id, df, rep_factor) {
  species_data <- subset(df, species == species_id)
  with(species_data, data.frame(
    species = species_id,
    group = rep(group, times = round(edge_length * rep_factor)),
    combined_group = rep(combined_group, times = round(edge_length * rep_factor)),
    step = rep(step, times = round(edge_length * rep_factor))
  ))
}

replicated_steps_df <- do.call(rbind, lapply(unique(combined_results$species), 
                                             replicate_steps_with_group, 
                                             df = combined_results, rep_factor = rep_factor))

###############################
## Sequence Split & Step Adjustment
###############################

# Adjust the step values based on a branch time split.
adjust_steps_based_on_split <- function(replicated_df, tree, split_node) {
  # Compute the split value from the branching times.
  bt <- branching.times(tree)
  split_value <- (round(100 - bt[split_node], 1)) * 10
  
  # Add a sequence number for each species and remove duplicate step values.
  df <- replicated_df %>%
    group_by(species) %>%
    mutate(sequence = row_number()) %>%
    ungroup() %>%
    group_by(sequence) %>%
    filter(!duplicated(step)) %>%
    ungroup()
  
  # Compute average step size for sequences up to the split value.
  avg_by_sequence <- df %>%
    filter(sequence <= split_value) %>%
    group_by(sequence) %>%
    summarize(average_step = mean(step, na.rm = TRUE), .groups = 'drop')
  
  # Update step values for sequences below or equal to the split value.
  left_join(df, avg_by_sequence, by = "sequence") %>%
    mutate(step = ifelse(sequence <= split_value, average_step, step)) %>%
    dplyr::select(-average_step)
}

result_df <- adjust_steps_based_on_split(replicated_steps_df, phylo_tree, split_node)

###############################
## Summarize & Plot
###############################

# Summarize average step sizes by combined group and sequence, then plot the results.
summarize_and_plot <- function(result_df) {
  avg_step_df <- result_df %>%
    group_by(combined_group, sequence) %>%
    summarize(avg_step = mean(step), .groups = 'drop')
  
  ggplot(avg_step_df, aes(x = sequence, y = avg_step, 
                               group = combined_group, color = combined_group)) +
    geom_line() +
    theme_minimal() +
    labs(x = "Sequence", y = "Average Step Size", color = "Combined Group")
  }

plot_obj <- summarize_and_plot(result_df)
plot_obj


################################################################################
# Calculate Angles of Phenotypic Trajectories
#
# This script computes the Euclidean distances and the angles (in radians) 
# between consecutive points (trajectories) in a 4D phenotypic space for each 
# species. It then joins phylogenetic branch lengths to these steps, computes 
# cumulative “time” along each trajectory, and finally creates a density plot 
# of the adjusted angles.
################################################################################


### Helper Functions ###########################################################
# Calculate the dot product between two numeric vectors.
calc_dot_product <- function(v1, v2) {
  sum(v1 * v2)
}

# Calculate the magnitude (Euclidean norm) of a numeric vector.
calc_magnitude <- function(v) {
  sqrt(sum(v^2))
}

# Calculate the angle (in radians) between two vectors.
calc_angle <- function(v1, v2) {
  dot_product  <- calc_dot_product(v1, v2)
  magnitude_v1 <- calc_magnitude(v1)
  magnitude_v2 <- calc_magnitude(v2)
  # Guard against division by zero
  if (magnitude_v1 == 0 || magnitude_v2 == 0) return(NA_real_)
  cos_theta <- dot_product / (magnitude_v1 * magnitude_v2)
  # Ensure the cosine value is in the valid range [-1, 1] to avoid NaN from acos()
  cos_theta <- min(max(cos_theta, -1), 1)
  angle <- acos(cos_theta)
  return(angle)
}

# Calculate sequential distances and angles for a given species.
# The function loops over consecutive rows (steps) in the species’ phenotypic
# data (columns "PC1", "PC2", "PC3", "PC4") and returns a tibble with:
#  - FromRow: starting step index,
#  - ToRow: ending step index,
#  - Distance: Euclidean distance between the two points,
#  - Angle: Angle (in radians) between the two vectors,
#  - Group: The group label from the data,
#  - CombinedGroup: If available, the combined group label.
calc_sequential_distances_and_angles <- function(species_data) {
  results <- tibble()  # Initialize an empty tibble
  
  if (nrow(species_data) > 1) {
    for (i in 1:(nrow(species_data) - 1)) {
      # Extract phenotypic vectors and ensure they are numeric
      v1 <- as.numeric(species_data[i, c("PC1", "PC2", "PC3", "PC4")])
      v2 <- as.numeric(species_data[i + 1, c("PC1", "PC2", "PC3", "PC4")])
      
      distance <- calc_distance(v1, v2)
      angle    <- calc_angle(v1, v2)
      
      # Create a tibble for the current pair
      new_row <- tibble(
        FromRow       = i,
        ToRow         = i + 1,
        Distance      = distance,
        Angle         = angle,
        Group         = species_data$group[i],
        CombinedGroup = if ("combined_group" %in% names(species_data)) {
          species_data$combined_group[i]
        } else {
          NA_character_
        }
      )
      results <- bind_rows(results, new_row)
    }
  }
  return(results)
}

### Compute Angles and Distances per Species ###################################

# Calculate sequential distances and angles for each species.
# Assumes that 'lineage_df_all' includes columns: PC1, PC2, PC3, PC4, species, and group.
angles <- combined_pcs_df %>%
  group_by(species) %>%
  group_modify(~ calc_sequential_distances_and_angles(.x)) %>%
  ungroup()

# Add a sequential step identifier (step_seq) for each species to enable joining
# with phylogenetic branch lengths. (Assumes that the ordering of the rows in
# the phenotypic trajectory matches that in the phylogenetic data.)
angles <- angles %>%
  group_by(species) %>%
  mutate(step_seq = row_number()) %>%
  ungroup()

### Join Phylogenetic Branch Lengths ###########################################
# Instead of using the old 'all_lengths_vector', we join branch lengths from
# 'species_edge_df' (previously computed in the phylogenetic data processing).
# 'species_edge_df' must contain columns: species, step_seq, and edge_length.
angles <- left_join(angles, species_edge_df %>% dplyr::select(species, step_seq, edge_length),
                    by = c("species", "step_seq"))

### Compute Cumulative “Time” Along Trajectories ##############################
# Here, we use the branch lengths (edge_length) to calculate cumulative time.
compute_cumulative_time <- function(df) {
  df %>%
    group_by(species) %>%
    mutate(time = cumsum(edge_length)) %>%
    ungroup()
}

# Apply the cumulative time calculation.
angles_with_time <- compute_cumulative_time(angles)

# Locate timing of split between super group and remaining species
plotTree(phylo_tree); nodelabels(); axisPhylo()
splitting_time_super = 30

# Create a new grouping variable:
# If the cumulative time is less than 30, label as "all"; otherwise use CombinedGroup.
angles_with_time <- angles_with_time %>%
  mutate(newGroup = ifelse(time < splitting_time_super, "all", CombinedGroup))

### Filter and Adjust Angles ###################################################
# Remove duplicate angle values within each newGroup.
angles_filtered <- angles_with_time %>%
  group_by(newGroup) %>%
  filter(!duplicated(Angle)) %>%
  ungroup()

# Adjust angles by taking the complementary angle (pi - Angle)
angles_filtered <- angles_filtered %>%
  mutate(Angle = pi - Angle)

### Plot Density of Angles #####################################################
# Define tick marks for the x–axis.
ticks <- c(pi / 10, (5 * pi) / 10, pi, (15 * pi) / 10)

# Create a density plot of the angles grouped by newGroup.
p <- ggplot(data = angles_filtered, aes(x = Angle, group = newGroup, fill = newGroup)) +
  geom_density(adjust = 3.5, alpha = 0.2) +
  theme_bw() +
  scale_x_continuous(breaks = ticks, limits = c(pi / 10, (15 * pi) / 10)) +
  labs(x = "Adjusted Angle (radians)", y = "Density", fill = "Group")





################################################################################
# Statistical Tests: Comparing Angles Across Groups
################################################################################

# ------------------------------------------------------------------------------
# Extract angle data for each group from the 'angles_filtered' data frame.
# ------------------------------------------------------------------------------
angles_simple <- angles_filtered %>%
  filter(newGroup == "simple") %>%
  pull(Angle)

angles_super <- angles_filtered %>%
  filter(newGroup == "super") %>%
  pull(Angle)

angles_all <- angles_filtered %>%
  filter(newGroup == "all") %>%
  pull(Angle)

# ------------------------------------------------------------------------------
# Perform a t-test to compare the mean angles between groups.
# ------------------------------------------------------------------------------
t_test_result <- t.test(angles_simple, angles_super,
                        alternative = "two.sided",
                        var.equal = FALSE)
cat("T-test Results (Simple vs. Super):\n")
print(t_test_result)

# Calculate and print the standard deviation for the "super" group.
sd_super <- sd(angles_super)
cat("Standard Deviation (Super group):", sd_super, "\n\n")

# ------------------------------------------------------------------------------
# Perform an F-test to compare the variances between the "simple" and "super" groups.
# ------------------------------------------------------------------------------
var_test_result <- var.test(angles_simple, angles_super)
cat("F-test Results (Simple vs. Super):\n")
print(var_test_result)

# ------------------------------------------------------------------------------
# Calculate variance for each group and print the results.
# ------------------------------------------------------------------------------
variance_simple <- var(angles_simple)
variance_super <- var(angles_super)
variance_all <- var(angles_all)

cat("Variances:\n")
cat("Simple group:", variance_simple, "\n")
cat("Super group:", variance_super, "\n")
cat("All group:", variance_all, "\n")




