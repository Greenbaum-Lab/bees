library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(geiger)  
library(missMDA)
library(PCDimension)  
library(phytools)
library(stringr)  
library(svglite)
library(umap)

# Data loading and initial preparation
data <- read.csv("data.csv")
rownames(data) <- data[,"name"]  # Set row names
data <- data[,-4]  # Remove row names column
species <- str_replace_all(paste(data[,"genus"], data[,"species"], sep = "_"), " ", "")

# Data pre-processing for PCA
num_cols <- sapply(data, is.numeric)  # Identify numeric columns more efficiently
dataq <- data[, num_cols]  # Subset numeric columns
dataq <- cbind.data.frame(dataq[,-grep("CF|CM|CA|QF", colnames(dataq))],
                          log10(dataq[,grep("CF|CM|CA|QF", colnames(dataq))]))  # Log-transform selected columns

# Dimension estimation and imputation
n.dim <- estim_ncpPCA(dataq, scale = TRUE, method = "Regularized", method.cv = "loo")
plot(n.dim$criterion ~ names(n.dim$criterion), xlab = "number of dimensions")
res.comp <- imputePCA(dataq, ncp = n.dim$ncp, scale = TRUE, method = "regularized")$completeObs

# Data correction
res.comp <- as.data.frame(res.comp)

# Identify columns with negative values and replace them
which(res.comp<0) #7 entries
res.comp %>% filter_all(any_vars(.< 0))
res.comp$PM[res.comp$PM <0] <- 0
res.comp$CL[res.comp$CL <0] <- 0.1
res.comp$WS[res.comp$WS <0] <- 2
res.comp$SX[res.comp$SX <0] <- 0.01
res.comp$NC[res.comp$NC <0] <- 0

# Manual corrections for non-plausible entries
res.comp["E.hyacinthina","CL"]<- 0.2
res.comp["B.fervidus","CL"]<- 0.3
res.comp["B.vitrea","CL"]<- 0.2
res.comp["B.vitrea","QL"]<- 0.4
res.comp["T.angustula","R"]<- 0.75
res.comp["C.capitata","R"]<- 0.75
res.comp["C.capitata","OR"]<- 1


# PCA and phylogenetic PCA
data.imp <- scale(res.comp)
tree <- read.newick("tree.nwk")
name.check(tree, data.imp)
pca.phyl <- phyl.pca(tree, data.imp, mode = "cov", method = "lambda")
pca.phyl$S <- pca.phyl$S * -1
biplot.phyl.pca(pca.phyl)
write.csv(pca.phyl$S, "pca_scores.csv")

# PCA dimension significance and dummy PCA creation for visualization
# Create PCA object with FactoMineR for easy plotting with factoextra
pca.dummy <- PCA(res.comp, graph = FALSE)
pca.dummy[["ind"]][["coord"]] <- pca.phyl$S  # Insert original scores
pca.dummy[["var"]][["coord"]] <- pca.phyl$L  # Insert original loadings
sdv <- unlist(summary(pca.phyl)[1], use.names = FALSE)
pca.dummy[["svd"]][["vs"]] <- sdv


# Data Visualization with ggplot2
sociality <- data[, "sociality"]
taxa <- data[, "taxa"]
PC1 <- paste("PC1 (", round(summary(pca.phyl)[[2]][2]*100, digits = 0), "%)", sep = "")
PC2 <- paste("PC2 (", round(summary(pca.phyl)[[2]][5]*100, digits = 0), "%)", sep = "")

# Create a color palette and shapes for plotting
sociality=data[,"sociality"]
sociality = c(sociality,rep("hgh",4),rep("sol",10))

taxa=data[,"taxa"]
taxa = c(taxa,rep("hal",4),rep("meg",10))

col <- c("#CC79A7", "#009E73","#E7B800", 
         "#0072B2","#000000", "#D55E00")

labels=c('Solitary','Communal','Subsocial',
         'Parasocial','Primitively eusocial',
         'Advanced eusocial')

shape = c(25,25,18,24,79,16,15)

# arrange data to match unique shapes
basic_plot <- fviz_pca_ind(pca.dummy,axes = c(1,2), label="all")
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(1:4,6,9,17:24),]) 

# plot
ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
       aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 8) +
  labs(x = PC1, y = PC2, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)

# ggsave("pca.svg",dpi=300,width = 22, height = 14)


#plot variables 1+2
plot.PCA(pca.dummy, axes=c(1,2), choix="var",
         label="all",invisible = "quali")+
  labs(x =PC1 , y = PC2,title = NULL)+
  theme(axis.title.x.bottom = element_text(hjust = 0.5,size = 14),
        axis.title.y.left = element_text(hjust = 0.5,size = 14 ))

# ggsave("var12.svg",dpi=300,width = 8, height = 8)

##########################umap 
shp=taxa
shp=shp %>%
  str_replace_all(c("api" = "25", "mel" = "16", 
                    "bom" = "25","eug" = "18",
                    "xyl" = "15","hal" = "24",
                    "meg" = "79"))

sociality<-factor(sociality,
                  levels=c("sol","sub","com", "par",
                           "prm", "hgh"))

plot.umap = function(x, labels, shape,
                     main = "",
                     colors = c("#000000", "#D55E00", "#CC79A7", "#E7B800", "#0072B2", "#009E73"),
                     pad = 0.1, cex = 3.2, pch = 17, add = FALSE, legend.suffix = "",
                     cex.main = 1, cex.legend = 2) {
  
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2] - xylim[1]) * pad) * c(-0.5, 0.5)
  
  if (!add) {
    par(mar = c(0.2, 0.7, 1.2, 0.7), ps = 10)
    plot(xylim, xylim, type = "n", axes = F, frame = F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border = "grey", lwd = 1.5)
  }
  
  fillable_shapes = c(21:25)
  
  # Adjust outline color for shapes 24 and 25
  shape_col = ifelse(as.integer(shape) %in% fillable_shapes & as.integer(shape) %in% c(24, 25), 
                     colors[as.integer(labels)], 
                     ifelse(as.integer(shape) %in% fillable_shapes, "black", colors[as.integer(labels)]))
  
  shape_fill = ifelse(as.integer(shape) %in% fillable_shapes, colors[as.integer(labels)], NA)
  
  points(layout[, 1], layout[, 2], col = shape_col, bg = shape_fill, cex = cex, pch = as.integer(shape))
  
  mtext(side = 3, main, cex = cex.main)
}


custom.settings = umap.defaults
custom.settings$n_neighbors = 20
custom.settings$min_dist = 0.2

# svg("umap.svg", width = 22, height = 14)

#pca umap+add apis circles
data.umap=umap(pca.phyl$S[,1:10],preserve.seed =F,config = custom.settings)
data.umap$layout = rbind(data.umap$layout, data.umap$layout[c(1:4,6,9,17:24),])
plot.umap(data.umap,sociality,shape = as.numeric(shp))

# dev.off()

#######subclusters test to show bombus seperation
# sociality.test=sociality
# sociality.test[c(35,36,38,42,43)]="sol" #bombus
# plot.umap(data.umap,sociality.test,shape = as.numeric(shp))


