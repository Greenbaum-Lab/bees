################################################################################
# Script: PCA, UMAP, and Ancestral State Reconstruction Analysis & Plots
# used in the supplementary information
################################################################################


################################################################################
## 2. PCA Plots Using Dummy Data
################################################################################

# ---------------------------
# 2.1 Prepare Data Categories
# ---------------------------

# Extract categories from 'data' and add extra rows (to match unique shapes)
sociality <- data[,"sociality"]
sociality <- c(sociality, rep("hgh", 4), rep("sol", 10))

taxa <- data[,"taxa"]
taxa <- c(taxa, rep("hal", 4), rep("meg", 10))

# Define plotting colors, labels, and shapes
color_values <- c("#CC79A7", "#009E73", "#E7B800", "#0072B2", "#000000", "#D55E00")
sociality_labels <- c("Solitary", "Communal", "Subsocial",
                      "Parasocial", "Primitively eusocial",
                      "Advanced eusocial")
# Note: 'shape_values' has 7 values; adjust if necessary.
shape_values <- c(25, 25, 18, 24, 79, 16, 15)

# ---------------------------
# 2.2 Helper Function: Create PCA Scatter Plot
# ---------------------------
create_pca_plot <- function(pca_obj, axes, sociality, taxa, shape_values,
                            color_values, sociality_labels, point_size = 5) {
  # Generate the PCA individual factor map using FactoMineR/factoextra
  basic_plot <- fviz_pca_ind(pca_obj, axes = axes, label = "all")
  
  # Duplicate specified rows to match extra categories (rows 20:23 and 28:37)
  basic_plot$data <- rbind(basic_plot$data, basic_plot$data[c(20:23, 28:37), ])
  
  # Create the ggplot scatter plot
  p <- ggplot(cbind(basic_plot$data[, 2:6], sociality, taxa),
              aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) +
    geom_point(size = point_size) +
    # Use dynamic axis labels (e.g. "PC1" and "PC2"); if custom labels are needed,
    # define them externally.
    labs(x = paste0("PC", axes[1]),
         y = paste0("PC", axes[2]),
         title = NULL,
         color = "sociality",
         shape = "taxa") +
    theme_light() +
    scale_shape_manual(values = shape_values) +
    scale_color_manual(labels = sociality_labels, values = color_values) +
    scale_fill_manual(labels = sociality_labels, values = color_values) +
    theme_test() +
    guides(color = FALSE, shape = FALSE, fill = FALSE)
  return(p)
}

# ---------------------------
# 2.3 Create PCA Plots for Different Axis Combinations
# ---------------------------
# Define the combinations of principal components to plot
axes_list <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
pca_plots <- list()

# Loop over each combination to create plots
for(i in seq_along(axes_list)) {
  pca_plots[[i]] <- create_pca_plot(pca_obj = pca.dummy,
                                    axes = axes_list[[i]],
                                    sociality = sociality,
                                    taxa = taxa,
                                    shape_values = shape_values,
                                    color_values = color_values,
                                    sociality_labels = sociality_labels,
                                    point_size = 5)
}

# Assign names to the plots (a, b, c, d, e, f)
a <- pca_plots[[1]]
b <- pca_plots[[2]]
c <- pca_plots[[3]]
d <- pca_plots[[4]]
e <- pca_plots[[5]]
f <- pca_plots[[6]]

# ---------------------------
# 2.4 Arrange and Save the Combined PCA Plot
# ---------------------------
combined_pca_plot <- plot_grid(a, b, c, d, e, f, ncol = 2)
ggsave("allpcs.svg", combined_pca_plot, dpi = 300, width = 12, height = 12)

################################################################################
## 3. Uncorrected PCA on Imputed Data
################################################################################

# ---------------------------
# 3.1 Load Imputed Data and Merge with Sociality
# ---------------------------
imputed_data <- read.csv('imputed_data.csv', row.names = 'X')
# NOTE: The object 'data' must already be loaded.
sociality <- data[,"sociality"]
imputed_data <- cbind.data.frame(imputed_data, sociality)

# ---------------------------
# 3.2 Run PCA on the Imputed Data
# ---------------------------
# The column number 18 is assumed to correspond to the qualitative variable 'sociality'
res_pca <- PCA(imputed_data, graph = FALSE, quali.sup = 18, scale.unit = TRUE)

# ---------------------------
# 3.3 Define Custom Axis Labels Based on Explained Variance
# ---------------------------
PC1_label <- paste("PC1 (", round(res_pca[["eig"]][1, 2], digits = 0), "%)", sep = "")
PC2_label <- paste("PC2 (", round(res_pca[["eig"]][2, 2], digits = 0), "%)", sep = "")

# ---------------------------
# 3.4 Prepare Categories (Reusing/Overwriting as Needed)
# ---------------------------
sociality <- data[,"sociality"]
sociality <- c(sociality, rep("hgh", 4), rep("sol", 10))

taxa <- data[,"taxa"]
taxa <- c(taxa, rep("hal", 4), rep("meg", 10))

# (Reuse shape_values, color_values, and sociality_labels defined above.)

# ---------------------------
# 3.5 Create the PCA Individuals Plot for Imputed Data
# ---------------------------
basic_plot <- fviz_pca_ind(res_pca, axes = c(1,2), label = "all")
basic_plot$data <- rbind(basic_plot$data, basic_plot$data[c(20:23, 28:37), ])

pca_uncorrect <- ggplot(cbind(basic_plot$data[, 2:6], sociality, taxa),
                        aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) +
  geom_point(size = 8) +
  labs(x = PC1_label, y = PC2_label, title = NULL,
       color = "sociality", shape = "taxa") +
  theme_light() +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(labels = sociality_labels, values = color_values) +
  scale_fill_manual(labels = sociality_labels, values = color_values) +
  theme_test() +
  guides(color = FALSE, shape = FALSE, fill = FALSE)

# ---------------------------
# 3.6 Plot the PCA Variables
# ---------------------------
pca_var_plot <- plot(res_pca, choix = "var") +
  labs(x = PC1_label, y = PC2_label, title = NULL) +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, size = 12))

# ---------------------------
# 3.7 Arrange and Save the Uncorrected PCA Plot
# ---------------------------
combined_pca_uncorrect <- plot_grid(pca_uncorrect, pca_var_plot)
ggsave("pca_uncorrect.svg", combined_pca_uncorrect, dpi = 300, width = 16, height = 8)

################################################################################
## 4. Uncorrected UMAP Analysis and Plotting
################################################################################

# ---------------------------
# 4.1 Prepare Shape Mapping for UMAP
# ---------------------------
# 'shp' is derived from taxa by replacing short taxon codes with numeric shape codes.
shp <- taxa %>%
  str_replace_all(c("api" = "25", "mel" = "16", 
                    "bom" = "25", "eug" = "18",
                    "xyl" = "15", "hal" = "24",
                    "meg" = "79"))
# Ensure 'sociality' is a factor with specified levels
sociality <- factor(sociality, levels = c("sol", "sub", "com", "par", "prm", "hgh"))

# ---------------------------
# 4.2 Define a Helper Function for UMAP Plotting
# ---------------------------
plot.umap <- function(x, labels, shape,
                      main = "",
                      colors = c("#000000", "#D55E00", "#CC79A7", "#E7B800", "#0072B2", "#009E73"),
                      pad = 0.1, cex = 3.2, pch = 17, add = FALSE, legend.suffix = "",
                      cex.main = 1, cex.legend = 2) {
  
  layout <- x
  if (inherits(x, "umap")) {
    layout <- x$layout
  }
  
  # Determine plotting limits with padding
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2] - xylim[1]) * pad) * c(-0.5, 0.5)
  
  if (!add) {
    par(mar = c(0.2, 0.7, 1.2, 0.7), ps = 10)
    plot(xylim, xylim, type = "n", axes = FALSE, frame = FALSE)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border = "grey", lwd = 1.5)
  }
  
  # Define fillable shapes (e.g., shapes 21 to 25)
  fillable_shapes <- 21:25
  
  # Determine outline and fill colors based on shape and label values
  shape_col <- ifelse(as.integer(shape) %in% fillable_shapes & as.integer(shape) %in% c(24, 25),
                      colors[as.integer(labels)],
                      ifelse(as.integer(shape) %in% fillable_shapes, "black", colors[as.integer(labels)]))
  
  shape_fill <- ifelse(as.integer(shape) %in% fillable_shapes, colors[as.integer(labels)], NA)
  
  # Plot the points
  points(layout[, 1], layout[, 2], col = shape_col, bg = shape_fill,
         cex = cex, pch = as.integer(shape))
  
  mtext(side = 3, main, cex = cex.main)
}

# ---------------------------
# 4.3 Set Custom UMAP Settings and Run UMAP
# ---------------------------
custom.settings <- umap.defaults
custom.settings$n_neighbors <- 20
custom.settings$min_dist <- 0.2

# Reload imputed data (if not already in memory)
imputed_data <- read.csv('imputed_data.csv', row.names = 'X')

# Open an SVG device for saving the UMAP plot
svg("umap_uncorrect.svg", width = 22, height = 14)

# Scale imputed data and run UMAP
data.umap <- umap(scale(imputed_data), preserve.seed = FALSE, config = custom.settings)
# Duplicate specific rows to match extra categories
data.umap$layout <- rbind(data.umap$layout, data.umap$layout[c(20:23, 28:37), ])

# Plot UMAP using the custom function
plot.umap(data.umap, sociality, shape = as.numeric(shp))
dev.off()

################################################################################
## 5. Evolutionary Analyses and Ancestral State Reconstruction (ASR)
################################################################################

# ---------------------------
# 5.1 Load Phylogenetic Tree and PCA Scores
# ---------------------------
# Set outer margins for plots
par(oma = c(4, 4, 2, 2))

# Read the tree (in Newick format) and PCA scores
tree <- read.newick("tree.nwk")      # Ensure the tree file exists
pca_scores <- read.csv("pca_scores.csv", row.names = "X")

# Comments regarding taxon labels (ensure these correspond to your tree)
#   78 - hal (black)
#   79 - apida
#   83 - eug
#   93 - mel
#   120 - bom
#   136 - apis

# ---------------------------
# 5.2 Paint Subtrees with States for ASR
# ---------------------------
tree_col <- tree
tree_col <- paintSubTree(tree_col, node = 81, state = "1", stem = FALSE)  # Root
tree_col <- paintSubTree(tree_col, node = 95, state = "2", stem = TRUE)   # Euglossini
tree_col <- paintSubTree(tree_col, node = 103, state = "6", stem = TRUE)  # Apis
tree_col <- paintSubTree(tree_col, node = 107, state = "5", stem = TRUE)  # Melipona
tree_col <- paintSubTree(tree_col, node = 135, state = "4", stem = TRUE)  # Bombus
tree_col <- paintSubTree(tree_col, node = 82, state = "7", stem = FALSE)  # Halictidae
tree_col <- paintSubTree(tree_col, node = 151, state = "3", stem = FALSE) # Xylocopinae

# ---------------------------
# 5.3 Plot the Colored Phylogenetic Tree
# ---------------------------
plotTree(tree)
nodelabels()
cols <- c("#000000", "#CC79A7", "#E7B800", "#0072B2", "#009E73", "#00FF00", "#D55E00")
names(cols) <- c(1, 2, 3, 4, 5, 6, 7)
plotSimmap(tree_col, cols, lwd = 3, pts = FALSE)

# ---------------------------
# 5.4 Generate and Save Phenograms for PC1 to PC4
# ---------------------------
# PC1 phenogram
svg("tree1S.svg", width = 18, height = 22)
pc1 <- pca_scores[, 1]
names(pc1) <- rownames(pca_scores)
phenogram(tree = tree_col, x = pc1, colors = cols, spread.labels = TRUE,
          lwd = 4, spread.cost = c(1, 0))
fancyTree(tree = tree, type = 'phenogram95', x = pc1, colors = cols, spread.cost = c(1, 0))
dev.off()

# PC2 phenogram
svg("tree2S.svg", width = 18, height = 22)
pc2 <- pca_scores[, 2]
names(pc2) <- rownames(pca_scores)
phenogram(tree = tree_col, x = pc2, colors = cols, spread.labels = TRUE,
          lwd = 4, spread.cost = c(1, 0))
fancyTree(tree = tree, type = 'phenogram95', x = pc2, colors = cols, spread.cost = c(1, 0))
dev.off()

# PC3 phenogram
svg("tree3S.svg", width = 18, height = 22)
pc3 <- pca_scores[, 3]
names(pc3) <- rownames(pca_scores)
phenogram(tree = tree_col, x = pc3, colors = cols, spread.labels = TRUE,
          lwd = 4, spread.cost = c(1, 0))
fancyTree(tree = tree, type = 'phenogram95', x = pc3, colors = cols, spread.cost = c(1, 0))
dev.off()

# PC4 phenogram
svg("tree4S.svg", width = 18, height = 22)
pc4 <- pca_scores[, 4]
names(pc4) <- rownames(pca_scores)
phenogram(tree = tree_col, x = pc4, colors = cols, spread.labels = TRUE,
          lwd = 4, spread.cost = c(1, 0))
fancyTree(tree = tree, type = 'phenogram95', x = pc4, colors = cols, spread.cost = c(1, 0))
dev.off()

################################################################################
## 6. ASR Combined with PCA
################################################################################

# ---------------------------
# 6.1 Load Data and Check Consistency with Phylogeny
# ---------------------------
dat <- read.csv("data.csv", row.names = "X")  # 'dat' is used below (ensure it is defined)
pca_scores <- read.csv("pca_scores.csv", row.names = "X")
phy <- read.newick("tree.nwk")
name.check(phy, pca_scores)

# ---------------------------
# 6.2 Choose PC Axes for ASR (Example uses PC3 and PC4)
# ---------------------------
index <- c(3, 4)

# ---------------------------
# 6.3 Estimate Ancestral States for Selected PCs
# ---------------------------
# Use fastAnc to estimate ancestral states
anc_x <- as.vector(fastAnc(tree = phy, x = pca_scores[, 'PC3']))
anc_y <- as.vector(fastAnc(tree = phy, x = pca_scores[, 'PC4']))

# ---------------------------
# 6.4 Combine Tip and Ancestral Data into a Single Data Frame
# ---------------------------
# NOTE: 'dat' is assumed to contain the original PC scores at the tips.
df <- data.frame(x = c(dat[, index[1]], anc_x),
                 y = c(dat[, index[2]], anc_y))

# ---------------------------
# 6.5 Extract Node Paths and Create a Data Frame for Plotting
# ---------------------------
path <- nodepath(phy)
a_list <- lapply(path, function(i) df[i, 1])
b_list <- lapply(path, function(i) df[i, 2])
demo <- cbind(data.frame(v1 = unlist(a_list),
                         v2 = rep(seq(length(a_list)), lengths(a_list))),
              unlist(b_list))
colnames(demo) <- c("x", "species", "y")
demo <- mutate(demo, x = as.numeric(x), y = as.numeric(y))

# ---------------------------
# 6.6 Merge Taxa Information
# ---------------------------
taxa_info <- as.data.frame(dat$taxa)
rownames(taxa_info) <- rownames(dat)
demo$group <- taxa_info[demo$species, ]

# ---------------------------
# 6.7 Create the ASR + PCA Combined Plot
# ---------------------------
cols_combo <- c("#00FF00", "#0072B2", "#CC79A7", "#D55E00", "#000000", "#009E73", "#E7B800")

p6 <- ggplot() +
  geom_path(data = demo, aes(x = x, y = y, colour = group, group = species),
            arrow = arrow(length = unit(0.3, "cm"))) +
  geom_point(data = demo, aes(x = x, y = y), color = "grey57", size = 0.5) +
  geom_point(data = dat, aes(x = PC3, y = PC4), color = "black", size = 0.5) +
  theme_void() +
  labs(x = "PC1", y = "PC2") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "none",
        panel.border = element_rect(fill = NA)) +
  scale_color_manual(values = cols_combo)

# ---------------------------
# 6.8 Combine Multiple Plots and Save
# ---------------------------
# NOTE: p1, p2, p3, p4, and p5 should be defined earlier in your script.
combined_asr_combo <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
ggsave("SUP_asr_combo.svg", combined_asr_combo, dpi = 300, width = 12, height = 16)
