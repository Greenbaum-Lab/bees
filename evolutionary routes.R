library(ape)
library(dplyr)
library(ggplot2)
library(geiger)
library(nlme)
library(openxlsx)
library(phytools)
library(picante)
library(svglite)
library(stringr)
library(tidyr)
library(tibble)
library(purrr)


################################################################################
# Ancestral State Reconstruction using Phylogenetic Trees and PCA Scores
#
# This script performs ancestral state reconstruction in one and two dimensions.
# It uses a phylogenetic tree (in Newick format) and PCA scores (from a CSV file)
# to (i) reconstruct the ancestral state for PC1 and plot a phenogram, and (ii)
# create a lineage plot of PC1 vs. PC2 with the evolutionary paths inferred.
################################################################################



################################################################################
# Part 1: Ancestral State Reconstruction (1D) for PC1
################################################################################

# Set outer margins for plotting (bottom, left, top, right)
par(oma = c(4, 4, 0, 0))

# -------------------------- Data Loading --------------------------------------

# Load the phylogenetic tree from a Newick file
phylo_tree <- read.newick("tree.nwk")

# Load PCA scores (assumed to have a column 'X' as rownames)
pca_scores <- read.csv("pca_scores.csv", row.names = "X")

# -------------------------- Node Labeling and Tree Coloring -------------------

# For reference, these comments map node numbers to clade names.
# (Adjust these as necessary for your dataset.)
#   - Node 81: Root
#   - Node 94: Euglossini
#   - Node 103: Apis
#   - Node 107: Melipona
#   - Node 135: Bombus
#   - Node 82: Halictidae
#   - Node 151: Xylocopinae

# Plot the original tree with node labels (useful for verifying node numbers)
plotTree(phylo_tree)
nodelabels()

# Create a copy of the tree that will be colored by clade/state
colored_tree <- phylo_tree

# Paint subtrees for different clades using their node numbers.
# The 'state' parameter here is a label (as a string) used to map to colors.
colored_tree <- paintSubTree(colored_tree, node = 81,  state = "1", stem = FALSE)  # Root
colored_tree <- paintSubTree(colored_tree, node = 94,  state = "2", stem = FALSE)  # Euglossini
colored_tree <- paintSubTree(colored_tree, node = 103, state = "6", stem = TRUE)   # Apis
colored_tree <- paintSubTree(colored_tree, node = 107, state = "5", stem = TRUE)   # Melipona
colored_tree <- paintSubTree(colored_tree, node = 135, state = "4", stem = TRUE)   # Bombus
colored_tree <- paintSubTree(colored_tree, node = 82,  state = "7", stem = TRUE)   # Halictidae
colored_tree <- paintSubTree(colored_tree, node = 151, state = "3", stem = FALSE)  # Xylocopinae

# Define a named vector of colors for each state
state_colors <- c("1" = "#000000",  # Black
                  "2" = "#CC79A7",  # Magenta/pink
                  "3" = "#E7B800",  # Yellow
                  "4" = "#0072B2",  # Blue
                  "5" = "#009E73",  # Greenish
                  "6" = "#00FF00",  # Bright green
                  "7" = "#D55E00")  # Orange

# Plot the colored tree using a simmap representation
plotSimmap(colored_tree, state_colors, lwd = 3, pts = FALSE)

# -------------------------- Ancestral State Reconstruction ----------------------

# Extract PC1 scores from the PCA dataset.
# The names of the vector are set to the tip labels for proper matching.
pca_pc1 <- pca_scores[, 1]
names(pca_pc1) <- rownames(pca_scores)

# Reconstruct ancestral states for PC1 using fastAnc (from phytools)
anc_states_pc1 <- fastAnc(phylo_tree, pca_pc1)

# -------------------------- Plot and Save Phenogram -----------------------------

# Save the phenogram plot to an SVG file
svg("reconstruct_pc1.svg", width = 18, height = 12)

# Plot the phenogram showing the evolution of PC1 over the tree
phenogram(tree = colored_tree,
          x = pca_pc1,
          colors = state_colors,
          spread.labels = TRUE,
          lwd = 4,
          spread.cost = c(1, 0))

# Alternatively, display tree with the confidence intervals
fancyTree(tree = phylo_tree,
          type = 'phenogram95',
          x = pca_pc1,
          colors = state_colors,
          spread.cost = c(1, 0))
dev.off()





################################################################################
# Part 2: Ancestral State Reconstruction in 2D
################################################################################

# -------------------------- Data Loading and Matching -------------------------

# Load additional data with taxa information
data_table <- read.csv("data.csv", row.names = "X")

# (Re)load the phylogenetic tree (if not already in memory)
phylo_tree <- read.newick("tree.nwk")

# Check that tip labels in the tree match the row names in the PCA scores dataset
name.check(phylo_tree, pca_scores)

# -------------------------- Ancestral State Calculation -------------------------

# Reconstruct ancestral states for PC1 and PC2
# for downstream analyses, PC3+4 are also needed to be calculated,
# so adjust the code as necessary
ancestral_pc1 <- as.vector(fastAnc(tree = phylo_tree, x = pca_scores[,"PC1"]))
ancestral_pc2 <- as.vector(fastAnc(tree = phylo_tree, x = pca_scores[,"PC2"]))

# -------------------------- Combine Tip and Ancestral Data ----------------------

# Combine the tip scores with the ancestral reconstructions.
# Assume that the ordering of rows in 'pca_scores' matches the tip order in 'phylo_tree'.
combined_scores <- data.frame(
  PC1 = c(pca_scores$PC1, ancestral_pc1),
  PC2 = c(pca_scores$PC2, ancestral_pc2)
)

# -------------------------- Build Data for Lineage Plot -------------------------

# Extract the node paths (list of node indices from root to each tip)
node_paths <- nodepath(phylo_tree)

# For each tip, create a data frame of its evolutionary path through the tree.
# Each path will include the PC1 and PC2 values along that lineage.
lineage_paths_list <- lapply(seq_along(node_paths), function(tip_index) {
  nodes_along_path <- node_paths[[tip_index]]
  data.frame(
    tip_index = tip_index,
    path_order = seq_along(nodes_along_path),
    PC1 = combined_scores[nodes_along_path, "PC1"],
    PC2 = combined_scores[nodes_along_path, "PC2"]
  )
})

# Combine the list into a single data frame for plotting
lineage_df <- do.call(rbind, lineage_paths_list)

# Add tip labels to the lineage data.
# The first 'n' rows (where n is the number of tips) correspond to the tip labels.
lineage_df <- lineage_df %>%
  mutate(tip_label = phylo_tree$tip.label[tip_index])

# Merge grouping information from the taxa data.
# It is assumed that row names in 'data_table' match the tip labels.
lineage_df$group <- data_table[lineage_df$tip_label, "taxa"]

# -------------------------- Define Colors for Groups ---------------------------

# Define a vector of colors for each group.
group_colors <- c("#00FF00", "#0072B2", "#CC79A7", 
                  "#D55E00", "#000000", "#009E73", "#E7B800")

# -------------------------- Create the Lineage Plot -----------------------------

ggplot() +
  # Plot the evolutionary paths (lineages) with arrows
  geom_path(data = lineage_df,
            aes(x = PC1, y = PC2, colour = group, group = tip_label),
            arrow = arrow(length = unit(0.3, "cm"))) +
  # Add small grey points to indicate internal nodes
  geom_point(data = lineage_df, aes(x = PC1, y = PC2),
             color = "grey57", size = 1) +
  # Overlay the tip (observed) PCA scores as black points
  geom_point(data = pca_scores, aes(x = PC1, y = PC2),
             color = "black", size = 1.5) +
  # Use a minimal theme (no background grid)
  theme_void() +
  labs(x = "PC1", y = "PC2") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "none",
        panel.border = element_rect(fill = NA)) +
  scale_color_manual(values = group_colors)


# Create a data frame that combines four principal components (PC1-PC4),
# along with species and group information.
# In this example, it is assumed that:
#   - Columns 3 and 4 of 'lineage_df' contain PC data (used twice for PC1/PC2 and PC3/PC4).
#   - Columns 5 and 6 of 'lineage_df' contain species and group data.
# lineage_df12 <- lineage_df
# lineage_df34 <- lineage_df
# lineage_df_all <- cbind(
#   lineage_df12[, c(3, 4)],   # Assumed to be PC1 and PC2
#   lineage_df34[, c(3, 4)],   # Assumed to be PC3 and PC4
#   lineage_df12[, c(5, 6)]    # Species and group information
# )
# colnames(lineage_df_all) <- c("PC1", "PC2", "PC3", "PC4", "species", "group")





################################################################################
# Identify Adaptive Shifts using PhylogeneticEM
#
# This script uses the PhylogeneticEM package to detect adaptive shifts on a
# phylogenetic tree. Several sections of the code allow for testing different
# data inputs (e.g., NA sensitivity or trait sensitivity) by uncommenting the
# corresponding lines.
################################################################################

# Load the PhylogeneticEM library
library(PhylogeneticEM)

# ------------------------------------------------------------------------------
# Alternative Data Options:
# Uncomment one of the following sections if you want to adjust the analysis
# for NA sensitivity or trait sensitivity.
#
# For NA sensitivity:
#   Use an alternative tree and/or principal component summary.
# ------------------------------------------------------------------------------
# phy <- phy_test
# dat <- pca_test$S[, 1:2]
# dat <- dataq

# For trait sensitivity:
#   Use trait data from a grouped summary.
# ------------------------------------------------------------------------------
# dat <- grouped_summary[, c(2, 4)]

# ------------------------------------------------------------------------------
# Data Preparation:
# Ensure that any open graphics devices are closed.
dev.off()

# Read in the phylogenetic tree from a Newick file.
phy <- read.newick("tree.nwk")

# Plot the raw tree with node labels for verification.
plotTree(phy)
nodelabels()

# Force the tree to be ultrametric (i.e., all tips equidistant from the root).
phy <- force.ultrametric(phy)

# Load PCA data
pca_scores <- read.csv("pca_scores.csv", row.names = "X")
pca_scores <- pca_scores[, 1:2]  # Choose the number of PCs to use.
pca_scores <- t(pca_scores)



# ------------------------------------------------------------------------------
# Run the PhylogeneticEM Analysis:
res <- PhyloEM(phylo = phy,
               Y_data = pca_scores,
               process = "scOU",
               method.selection = c("LINselect", "Djump"),
               stationary.root = TRUE,
               random.root = TRUE)

# Plot the overall results.
svglite::svglite("shifts.svg")
plot(res)
dev.off()

# Plot the selection criterion (e.g., model selection criterion) for LINselect.
plot_criterion(res, method.selection = "LINselect")

# ------------------------------------------------------------------------------
# Multiple Solutions for 4 Principal Components:
multi_solutions <- params_process(res, K = 4)

# Plot the multi–solutions
svglite::svglite("shifts_4pcs_multi.svg")
multi_solutions <- equivalent_shifts(phy, multi_solutions)
plot(multi_solutions, show_shifts_values = FALSE, shifts_cex = 0.3)
dev.off()


