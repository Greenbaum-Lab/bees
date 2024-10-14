library(phytools)
library(tibble)
library(ggplot2)
library(svglite)
library(picante)
library(geiger)
library(ape)
library(corrplot)
library(caper)

# Compute phylogenetic signal
# Prepare data after imputation
data <- read.csv('imputed_data.csv',row.names = 'X')
tree=read.newick('tree.nwk')
name.check(data, tree)

# Initialize a dataframe to store the phylogenetic signal and p-values
phyl.sig <- setNames(data.frame(matrix(ncol = ncol(data), nrow = 4)),
                     colnames(data))
rownames(phyl.sig) <- c("lambda", "lambda.p", "k", "k.p")

# Compute phylogenetic signal and p-value for lambda and k
for (i in 1:ncol(phyl.sig)) {
  phyl.sig[1, i] <- phylosig(tree, data[, i], method = "lambda",)[1]
  phyl.sig[2, i] <- phylosig(tree, data[, i], method = "lambda", test = TRUE)[4]
  phyl.sig[3, i] <- phylosig(tree, data[, i], method = "K")[1]
  phyl.sig[4, i] <- phylosig(tree, data[, i], method = "K", test = TRUE)[2]
}

# Round the values and convert to matrix
phyl.sig <- round(phyl.sig, digits = 3)
phyl.sig <- as.matrix(phyl.sig)
print(phyl.sig)
par(mar=c(4,4,4,4))

# order df
desired_order <- c("WM", "PM", "SX", "OR", "WS", "R", "NC", "QL", "CL", "BC", 
                  "QW", "OG", "QS", "QF", "CF", "CA", "CM")

phyl.sig <- phyl.sig[, desired_order]
# Create a bar plot for phylogenetic signal
lim <- 1.2 * max(phyl.sig)
sig.plot <- barplot(phyl.sig[c(1, 3),], beside = TRUE, cex.lab = 1.5,
                    args.legend = list(x = "right", bty = "n", cex = 1, x.intersp = 0.1,
                                       inset = c(0, 0)),
                    legend.text = TRUE, col = c("lightgrey", "darkgrey"), cex.axis = 1.5, cex.names = 1.2,
                    ylim = c(0, 3), ylab = "Phylogenetic signal") +
  abline(h = 1, lwd = 2, lty = 2)

# Save the plot as an SVG file
svglite("phylo_signal.svg", width = 18, height = 10)
dev.off()

# Check significance after Benjamini–Hochberg correction
print(phyl.sig[4,])
print(round(p.adjust(phyl.sig[4,], "BH"), 6))


##########################################################
####################### Compute correlation matrix ########
tree= read.newick('tree.nwk')
df = read.csv("imputed_data.csv",row.names = 'X')
plotTree(tree)
name.check(tree, df)

# Initialize a list to store the contrasts for each trait
contrasts_list <- list()

# Loop over columns of the data
for (i in 1:ncol(df)) {
  # Calculate independent contrasts for this trait
  contrasts_list[[i]] <- pic(df[, i], tree)
}

# Combine all contrasts into a data frame
contrasts <- do.call(cbind, contrasts_list)
correlation_matrix <- cor(contrasts, method = "spearman")

# Calculate the p-value matrix
p.mat <- cor.mtest(contrasts)$p

# Set column and row names for matrices
colnames(correlation_matrix) <- colnames(df)
rownames(correlation_matrix) <- colnames(df)
colnames(p.mat) <- colnames(df)
rownames(p.mat) <- colnames(df)

correlation_matrix <- correlation_matrix[, desired_order]
correlation_matrix <- correlation_matrix[desired_order, ]
p.mat <- p.mat[, desired_order]
p.mat <- p.mat[desired_order, ]

# Define a function to adjust p-values for multiple comparisons
p.adjust.method <- function(p) p.adjust(p, method = "BH")

# Adjust the p-values
p.mat.adjusted <- matrix(p.adjust.method(p.mat), ncol = ncol(p.mat))
p.mat.adjusted <- as.matrix(p.mat.adjusted)
rownames(p.mat.adjusted) <- colnames(p.mat)
colnames(p.mat.adjusted) <- colnames(p.mat)

# Plot the correlation matrix
svglite("cor.svg",width = 10,height = 10)
corrplot(correlation_matrix, order = "original", method = "square",
         hclust.method = "ward.D",
         p.mat = p.mat.adjusted, insig = "blank",
         col = rev(COL2("RdBu", 200)), cl.pos = "b",
         tl.cex = 2, cl.cex = 1, tl.srt = 45, tl.col = "black")
dev.off()


# make hierarchical clustering
svglite("heatmap.svg",10,10)
heatmap(correlation_matrix)
dev.off()




############# traits barplot
desired_order <- c("WM", "PM", "SX", "OR", "WS", "R", "NC", "QL", "CL", "BC", 
                   "QW", "OG", "QS", "QF", "CF", "CA", "CM")


# Reading the imputed data
data_imputed <- read.csv('imputed_data.csv', row.names = 'X')  # Assuming 'X' is your first column as row names
data_imputed <- data_imputed[, desired_order]
rownames(data_imputed)
data_na <- dataq  # Assuming dataq is your original data with NAs
data_na <- data_na[, desired_order]

# Mask for NA values in the original data
na_mask <- is.na(data_na)

# Create a plot layout
par(mfrow = c(1,ncol(data_imputed)), mar = c(1, 1, 1, 1))  # Minimized margins, stack plots vertically

for (i in seq_along(data_imputed)) {
  # Current column data
  data <- data_imputed[[i]]
  original_na <- na_mask[, i]
  
  # Horizontal bar plot for the current column
  heights <- barplot(data, horiz = TRUE, space = 1, axes = FALSE, border = NA,
                     col = ifelse(original_na, "grey", "black"), main = names(data_imputed)[i])
  
  # Adding a base line on the x-axis to mark the starting point
  abline(v = 0, lwd = 0.5)
 }

