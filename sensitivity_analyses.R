################################################################################
# Script: Inclusion Threshold Adjustment, PCA Sensitivity, and Taxonomic 
#         Sampling & Phylogenetic Tree Sensitivity Analyses
################################################################################

library(missMDA)
library(FactoMineR)
library(ggplot2)
library(factoextra)
library(PCDimension)
library(dplyr)
library(stringr)
library(umap)
library(phytools)
library(geiger)
library(picante)
library(svglite)
library(ggpubr)
library(gridExtra)
library(abs)

#############################
## 2. Adjusting the Inclusion Threshold & Exploring Missing Species
#############################

# Load phylogenetic tree and species data
tree <- read.newick("tree.nwk")
data <- read.csv("data.csv", row.names = "X")
data_test <- data

# Subset numeric columns from the data for missing data analysis
num_cols <- sapply(data_test, is.numeric)
dataq_test <- data_test[, num_cols]

# Compute overall proportion of missing values (for reporting)
p_na_overall <- sum(is.na(dataq_test)) / (nrow(dataq_test) * ncol(dataq_test))
cat("Overall proportion of missing values:", p_na_overall, "\n")

# Create a vector with the count of NAs per species (row)
p_na <- numeric(nrow(dataq_test))
for (i in 1:nrow(dataq_test)) {
  p_na[i] <- sum(is.na(dataq_test[i, ]))
}

# Identify species with more than a threshold of missing data (e.g., 15% of columns)
discard <- which(p_na > ncol(dataq_test) * 0.15)
cat("Number of species to discard:", length(discard), "\n")

# Optionally, inspect species with >25% missing values:
print(dataq_test[which(p_na > ncol(dataq_test) * 0.25), ])

# Remove species with excessive missing data from both datasets
data_test <- data_test[-discard, ]
dataq_test <- dataq_test[-discard, ]

# Log-transform selected columns (e.g., colony size avg, max, foundation, fecundity)
columns_to_log <- c("CF", "CM", "CA", "QF")
dataq_test[columns_to_log] <- lapply(dataq_test[columns_to_log], log10)

#############################
## 3. Imputation and Estimation of PCA Dimensions
#############################

# Estimate number of components for PCA imputation
n_dim <- estim_ncpPCA(dataq_test, scale = TRUE, method = "Regularized",
                      method.cv = "loo")
# Plot the criterion (e.g., cross-validation error) vs. number of dimensions
plot(n_dim$criterion ~ names(n_dim$criterion), xlab = "Number of dimensions")

# Impute missing data using the estimated number of dimensions
sensitiv <- imputePCA(dataq_test, ncp = n_dim$ncp, scale = TRUE, method = "regularized")
sensitiv <- as.data.frame(sensitiv$completeObs)

# Check for any negative imputed values
print(which(sensitiv < 0))

# Manually adjust specific imputed values based on threshold criteria (here, for 0.15)
sensitiv["M.bicolor", "WS"]       <- 1
sensitiv["E.hyacinthina", "CL"]     <- 0.2
sensitiv["C.australensis", "CL"]    <- 0.2
sensitiv["E.tridentata", "PM"]      <- 0

# Scale imputed data and match with phylogeny
sensitiv <- as.data.frame(scale(sensitiv))
sensitiv <- match.phylo.data(tree, sensitiv)$data
phy_test <- match.phylo.data(tree, sensitiv)$phy
name.check(phy_test, sensitiv)

#############################
## 4. Phylogenetic PCA and Dummy PCA Object for Visualization
#############################

# Perform phylogenetic PCA using the lambda method
pca_test <- phyl.pca(phy_test, sensitiv, mode = "cov", method = "lambda")
# (Optional: Flip PC1 sign if desired; e.g., pca_test$S[,1] <- -pca_test$S[,1])

# Plot a biplot of the phylogenetic PCA results
biplot.phyl.pca(pca_test)

# Create a dummy PCA object (using FactoMineR's PCA) to enable external visualization
pca.dummy_test <- PCA(dataq_test, graph = TRUE)
pca.dummy_test[["ind"]][["coord"]] <- pca_test$S
pca.dummy_test[["var"]][["coord"]] <- pca_test$L
sdv <- as.numeric(unlist(summary(pca_test)[1]))
pca.dummy_test[["svd"]][["vs"]] <- sdv
pca.dummy_test[["var"]][["cos2"]] <- pca_test$L^2

# Create axis labels with percentage variance explained
PC1 <- paste("PC1 (", round(summary(pca_test)[[2]][2] * 100, digits = 0), "%)", sep = "")
PC2 <- paste("PC2 (", round(summary(pca_test)[[2]][5] * 100, digits = 0), "%)", sep = "")
PC3 <- paste("PC3 (", round(summary(pca_test)[[2]][8] * 100, digits = 0), "%)", sep = "")
PC4 <- paste("PC4 (", round(summary(pca_test)[[2]][11] * 100, digits = 0), "%)", sep = "")

# Reorder PCA scores to match the species order in data_test
pca.dummy_test$ind$coord <- pca.dummy_test$ind$coord[match(rownames(data_test),
                                                           rownames(pca.dummy_test$ind$coord)), ]

#############################
## 5. Visualization of Phylogenetic PCA Results
#############################

# Define group indices based on genus for later color coding
melipona <- which(data_test$genus == 'Melipona')
apis     <- which(data_test$genus == 'Apis')

# Create categorical variables for plotting
sociality_test <- data_test[,"sociality"]
sociality_test <- c(sociality_test, rep("hgh", length(apis)), rep("sol", length(melipona)))

taxa_test <- data_test[,"taxa"]
taxa_test <- c(taxa_test, rep("hal", length(apis)), rep("meg", length(melipona)))

# Define plotting aesthetics
col_vals    <- c("#CC79A7", "#009E73", "#E7B800", "#0072B2", "#000000", "#D55E00")
labels_vals <- c("Solitary", "Communal", "Subsocial", "Parasocial", 
                 "Primitively eusocial", "Advanced eusocial")
shape_vals  <- c(25, 25, 18, 24, 79, 16, 15)

# Prepare PCA individual plot using factoextra
basic_plot <- fviz_pca_ind(pca.dummy_test, axes = c(1, 2))
# Append rows corresponding to the additional indices (apis and melipona)
basic_plot$data <- rbind(basic_plot$data, basic_plot$data[c(apis, melipona), ])

# Plot the PCA scores with custom aesthetics
p1 <- ggplot(cbind(basic_plot$data[, 2:6], sociality_test, taxa_test),
             aes(x = x, y = y, col = sociality_test, shape = taxa_test, fill = sociality_test)) +
  geom_point(size = 8) +
  labs(x = PC1, y = PC2, title = NULL, color = "sociality_test", shape = "taxa_test") +
  theme_light() +
  scale_shape_manual(values = shape_vals) +
  scale_color_manual(labels = labels_vals, values = col_vals) +
  scale_fill_manual(labels = labels_vals, values = col_vals) +
  theme_test() +
  guides(color = FALSE, shape = FALSE, fill = FALSE)

# Combine with original PCA scores (assuming pca.phyl object exists with S matrix)
# NOTE: pca.phyl is assumed to be available in your environment.
dat <- pca.phyl$S[, 1:2]
colnames(dat) <- c("x", "y")
taxa_orig <- data[,"taxa"]

social_combo <- c(rep("gr", 80), sociality_test)
taxa_combo   <- c(taxa_orig, taxa_test)
pca_combo    <- rbind(dat, basic_plot$data[, 2:3])
pca_combo    <- as.data.frame(pca_combo)
pca_combo$group <- substring(rownames(pca_combo), 1, 6)

# Define a different color palette for the combined plot
col_combo <- c("#CC79A7", "grey", "#E7B800", "#009E73", "#0072B2", "#000000", "#D55E00")

p2 <- ggplot(cbind(pca_combo, social_combo, taxa_combo),
             aes(x = x, y = y, col = social_combo, shape = taxa_combo, fill = social_combo)) +
  geom_point(size = 8) +
  geom_path(aes(group = group)) +
  labs(x = PC1, y = PC2, title = NULL, color = "social_combo", shape = "taxa_combo") +
  theme_light() +
  scale_shape_manual(values = shape_vals) +
  scale_color_manual(labels = labels_vals, values = col_combo) +
  scale_fill_manual(labels = labels_vals, values = col_combo) +
  theme_test() +
  guides(color = FALSE, shape = FALSE, fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 30))

# Save the combined PCA plot
ggsave("sens_na_6.svg", dpi = 300, width = 16, height = 10)

# Report overall missing proportion again for verification
p_na_overall_new <- sum(is.na(dataq_test)) / (nrow(dataq_test) * ncol(dataq_test))
cat("New overall proportion of missing values:", p_na_overall_new, "\n")

#############################
## 6. Confidence Ellipses for Sensitivity Analysis (Leave-One-Out PCA)
#############################

# Load and scale imputed data for sensitivity analysis
data_imp <- read.csv("imputed_data.csv", row.names = "X")
data_sens <- as.data.frame(scale(data_imp))
tree <- read.newick("tree.nwk")
name.check(tree, data_sens)

# Helper function: Run phylogenetic PCA on a data subset and adjust PC sign
run_pca <- function(data_subset) {
  pca_phyl <- phyl.pca(tree, data_subset, mode = "cov", method = "lambda")
  
  if (pca_phyl[["S"]]["A.mellifera", 1] < 0) {
    pca_phyl$S[, 1] <- pca_phyl$S[, 1] * -1
    print("I flipped PC1")
  }
  
  if (pca_phyl[["S"]]["T.angustula", 2] > 0) {
    pca_phyl$S[, 2] <- pca_phyl$S[, 2] * -1
    print("I flipped PC2")
  }
  
  if (pca_phyl[["S"]]["A.mellifera", 3] > 0 || pca_phyl[["S"]]["M.bicolor", 3] > 4) {
    pca_phyl$S[, 3] <- pca_phyl$S[, 3] * -1
    print("I flipped PC3")
  }
  
  if (pca_phyl[["S"]]["B.terrestris", 4] > 0 && pca_phyl[["S"]]["M.bicolor", 4] < abs(2)) {
    pca_phyl$S[, 4] <- pca_phyl$S[, 4] * -1
    print("I flipped PC4")
  }
  
  return(pca_phyl$S[, 1:4])
}

# Initialize an empty data frame to store leave-one-out PCA results
df_loo <- data.frame()

# Loop through each column (leave-one-out) of the scaled data
for (i in 1:ncol(data_sens)) {
  data_subset <- data_sens[, -i]
  pca_subset <- run_pca(data_subset)
  df_subset <- cbind(pca_subset, group = rownames(pca_subset))
  df_loo <- rbind(df_loo, df_subset)
}

# Ensure PCA columns are numeric
df_loo[, 1:4] <- lapply(df_loo[, 1:4], as.numeric)
plot(df_loo[, 3:4])  # Quick diagnostic plot

# Scatter plot of PC1 vs. PC2 with 95% confidence ellipses using ggpubr
p_conf1 <- ggscatter(df_loo, x = "PC1", y = "PC2", color = "white", alpha = 0.2) +
  theme_bw() + theme(legend.position = "none") +
  stat_conf_ellipse(aes(color = group))
ggsave("sens_pca1.svg", p_conf1, width = 15, height = 10, dpi = 300)

# Scatter plot of PC3 vs. PC4 with 95% confidence ellipses
p_conf2 <- ggscatter(df_loo, x = "PC3", y = "PC4", color = "white", alpha = 0.2) +
  theme_bw() + theme(legend.position = "none") +
  stat_conf_ellipse(aes(color = group))
ggsave("sens_pca1_pc34.svg", p_conf2, width = 15, height = 18, dpi = 300)

#############################
## 7. Testing Taxonomic Sampling Bias
#############################

# Load taxa information from the original data to identify taxa to remove
taxa_data <- read.csv("data.csv", row.names = "X")
taxa_to_remove <- (taxa_data[,"taxa"] == 'api' | taxa_data[,"taxa"] == 'mel' | taxa_data[,"taxa"] == 'bom')
n_remove <- (sum(taxa_to_remove) - sum(!taxa_to_remove))
taxa_remove_indices <- which(taxa_to_remove)

# Load and scale imputed data before PCA
data_imp <- read.csv("imputed_data.csv", row.names = "X")
data_scaled <- as.data.frame(scale(data_imp))
tree <- read.newick("tree.nwk")
name.check(tree, data_scaled)

# Define a function to flip PCA axes based on specified species scores
flip_pca_axes <- function(res_pca) {
  if (res_pca["D.novaeangliae", 1] > 0) {
    res_pca[, 1] <- -res_pca[, 1]
    print("I flipped PC1")
  }
  if (res_pca["X.virginica", 2] < 0) {
    res_pca[, 2] <- -res_pca[, 2]
    print("I flipped PC2")
  }
  if (res_pca["M.genalis", 3] > 0) {
    res_pca[, 3] <- -res_pca[, 3]
    print("I flipped PC3")
  }
  if (res_pca["E.townsendi", 4] > 0) {
    res_pca[, 4] <- -res_pca[, 4]
    print("I flipped PC4")
  }
  return(res_pca)
}

# Loop for equal taxonomic sampling and phylogenetic PCA
final_df_tax <- data.frame()
iterations <- 100
for (i in 1:iterations) {
  rows_remove <- sample(taxa_remove_indices, n_remove, replace = FALSE)
  data_subset <- data_scaled[-rows_remove, ]
  subsets <- match.phylo.data(tree, data_subset)
  pca_phyl <- phyl.pca(subsets$phy, subsets$data, mode = "cov", method = "lambda")
  res_pca <- pca_phyl$S[, 1:4]
  print(i)
  res_pca <- flip_pca_axes(res_pca)
  res_pca <- cbind(res_pca, group = rownames(res_pca))
  final_df_tax <- rbind(final_df_tax, res_pca)
}
final_df_tax[, 1:4] <- lapply(final_df_tax[, 1:4], as.numeric)
plot(final_df_tax[1:56, 3:4])  # Diagnostic plot

# Scatter plot of PC1 vs. PC2 with confidence ellipses for taxonomic bias
p_tax1 <- ggscatter(final_df_tax, x = "PC1", y = "PC2", color = "white", alpha = 0.2) +
  theme_bw() + theme(legend.position = "none") +
  stat_conf_ellipse(aes(color = group))
ggsave("tax_bias12.svg", p_tax1, width = 15, height = 10, dpi = 300)

# Scatter plot of PC3 vs. PC4 for taxonomic bias
p_tax2 <- ggscatter(final_df_tax, x = "PC3", y = "PC4", color = "white", alpha = 0.2) +
  theme_bw() + theme(legend.position = "none") +
  stat_conf_ellipse(aes(color = group))
ggsave("tax_bias34.svg", p_tax2, width = 15, height = 10, dpi = 300)

# Compute group means for PC1 and PC2 for select species and plot them
group_means <- final_df_tax %>%
  group_by(group) %>%
  summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2))

taxa_info <- taxa_data[, c('sociality', 'taxa')]
taxa_info$species <- rownames(taxa_info)
group_means <- group_means %>% left_join(taxa_info, by = c("group" = "species"))

p_means <- ggplot(group_means[c(56:62, 65, 64, 67), ], 
                  aes(x = mean_PC1, y = mean_PC2, color = sociality, shape = taxa, fill = sociality)) +
  geom_point(size = 8) +
  labs(x = PC1, y = PC2, title = NULL, color = "Sociality", shape = "Taxa") +
  theme_light() +
  scale_shape_manual(values = shape_vals) +
  scale_color_manual(labels = labels_vals, values = col_vals) +
  scale_fill_manual(labels = labels_vals, values = col_vals) +
  theme(legend.position = "none")
ggsave("tax_bias12_mean.svg", p_means, width = 15, height = 10, dpi = 300)

#############################
## 8. Phylogenetic Tree Sensitivity Analysis
#############################

# Load multiple trees from file (assumed to be in Newick format)
trees <- read.newick('all_trees.nwk')

# Load and scale imputed data
data_imp <- read.csv("imputed_data.csv", row.names = "X")
data_scaled <- as.data.frame(scale(data_imp))
name.check(tree, data_scaled)

final_df_tree <- data.frame()
iterations <- 100
for (i in 1:iterations) {
  # Sample one tree from the list of trees
  tree_sample <- sample(x = trees, size = 1, replace = FALSE)[[1]]
  pca_phyl <- phyl.pca(tree_sample, data_scaled, mode = "cov", method = "lambda")
  res_pca <- pca_phyl$S[, 1:4]
  print(i)
  res_pca <- flip_pca_axes(res_pca)
  res_pca <- cbind(res_pca, group = rownames(res_pca))
  final_df_tree <- rbind(final_df_tree, res_pca)
}
final_df_tree[, 1:4] <- lapply(final_df_tree[, 1:4], as.numeric)
plot(final_df_tree[1:200, 1:2])  # Diagnostic plot

p_tree <- ggscatter(final_df_tree, x = "PC3", y = "PC4", color = "white",
                    alpha = 0.2, mean.point = TRUE) +
  theme_bw() + theme(legend.position = "none") +
  stat_conf_ellipse(aes(color = group))
ggsave("phyl_bias12.svg", p_tree, width = 15, height = 10, dpi = 300)
