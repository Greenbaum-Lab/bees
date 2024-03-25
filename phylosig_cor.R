# Load necessary libraries
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
dat <- res.comp

# Match the order of the tree and data
dat <- dat[match(tree$tip.label, rownames(dat)),]
name.check(dat, tree)

# Initialize a dataframe to store the phylogenetic signal and p-values
phyl.sig <- setNames(data.frame(matrix(ncol = ncol(res.comp), nrow = 4)),
                     colnames(res.comp))
rownames(phyl.sig) <- c("lambda", "lambda.p", "k", "k.p")

# Compute phylogenetic signal and p-value for lambda and k
for (i in 1:ncol(phyl.sig)) {
  phyl.sig[1, i] <- phylosig(tree, dat[, i], method = "lambda",)[1]
  phyl.sig[2, i] <- phylosig(tree, dat[, i], method = "lambda", test = TRUE)[4]
  phyl.sig[3, i] <- phylosig(tree, dat[, i], method = "K")[1]
  phyl.sig[4, i] <- phylosig(tree, dat[, i], method = "K", test = TRUE)[2]
}

# Round the values and convert to matrix
phyl.sig <- round(phyl.sig, digits = 3)
phyl.sig <- as.matrix(phyl.sig)
print(phyl.sig)

# Create a bar plot for phylogenetic signal
lim <- 1.2 * max(phyl.sig)
sig.plot <- barplot(phyl.sig[c(1, 3),], beside = TRUE, cex.lab = 1.5,
                    args.legend = list(x = "right", bty = "n", cex = 1, x.intersp = 0.1,
                                       inset = c(0, 0)),
                    legend.text = TRUE, col = c("lightgrey", "darkgrey"), cex.axis = 1.5, cex.names = 1.2,
                    ylim = c(0, 3), ylab = "Phylogenetic signal") +
  abline(h = 1, lwd = 2, lty = 2)

# Save the plot as an SVG file
# svglite("phylo_signal.svg", width = 18, height = 10)
# dev.off()

# Check significance after Benjamini–Hochberg correction
print(phyl.sig[4,])
print(round(p.adjust(phyl.sig[4,], "BH"), 6))




####################### Compute covariance matrix
plotTree(tree)

# Match data to tree names
dat <- res.comp
name.check(tree, dat)
dat <- dat[match(tree$tip.label, rownames(dat)),]
dat <- as.data.frame(dat)

# Initialize a list to store the contrasts for each trait
contrasts_list <- list()

# Loop over columns of the data
for (i in 1:ncol(dat)) {
  # Calculate independent contrasts for this trait
  contrasts_list[[i]] <- pic(dat[, i], tree)
}

# Combine all contrasts into a data frame
contrasts <- do.call(cbind, contrasts_list)
correlation_matrix <- cor(contrasts, method = "spearman")

# Calculate the p-value matrix
p.mat <- cor.mtest(contrasts, conf.level = 0.95)$p

# Set column and row names for matrices
colnames(correlation_matrix) <- colnames(dat)
rownames(correlation_matrix) <- colnames(dat)
colnames(p.mat) <- colnames(dat)
rownames(p.mat) <- colnames(dat)

# Define a function to adjust p-values for multiple comparisons
p.adjust.method <- function(p) p.adjust(p, method = "BH")

# Adjust the p-values
p.mat.adjusted <- matrix(p.adjust.method(p.mat), ncol = ncol(p.mat))
p.mat.adjusted <- as.matrix(p.mat.adjusted)
rownames(p.mat.adjusted) <- colnames(dat)
colnames(p.mat.adjusted) <- colnames(dat)

# Plot the correlation matrix
# svglite("cor.svg",width = 10,height = 10)
corrplot(correlation_matrix, order = "hclust", method = "square",
         hclust.method = "ward.D",
         p.mat = p.mat.adjusted, insig = "blank",
         col = rev(COL2("RdBu", 200)), cl.pos = "b",
         tl.cex = 2, cl.cex = 1, tl.srt = 45, tl.col = "black")
# dev.off()

#plot cor matrix
corrplot(correlation_matrix, order="hclust",method = "square",
         hclust.method = "ward.D",
         p.mat = p.mat.adjusted,insig = "blank",
         col=rev(COL2("RdBu",200)),cl.pos = "b" ,
         tl.cex = 2, cl.cex = 1,tl.srt = 45,tl.col = "black")


# make hierarchical clustering
# svglite("heatmap.svg",10,10)
heatmap(correlation_matrix)
# dev.off()

