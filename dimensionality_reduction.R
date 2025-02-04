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
library(geometry)
library(cowplot)


################### Data pre-processing ################### 
# Data loading and initial preparation
data <- read.csv("data.csv",row.names = 'X')

# Get numeric columns
num_cols <- sapply(data, is.numeric)  
dataq <- data[, num_cols]  # Subset numeric columns

# log10 transformation for highly variable traits 
columns_to_log <- c("CF", "CM", "CA", "QF")
dataq[columns_to_log] <- lapply(dataq[columns_to_log], log10)

# Dimension estimation and imputation using iterative PCA
n.dim <- estim_ncpPCA(dataq, scale = TRUE, method = "Regularized", method.cv = "loo")
plot(n.dim$criterion ~ names(n.dim$criterion), xlab = "number of dimensions")
imputed_data <- imputePCA(dataq, ncp = n.dim$ncp, scale = TRUE, method = "regularized")$completeObs
imputed_data <- as.data.frame(imputed_data)

# Identify columns with negative values and fix them
which(imputed_data<0) #8 entries
imputed_data %>% filter_all(any_vars(.< 0))
imputed_data$PM[imputed_data$PM <0] <- 0
imputed_data$CL[imputed_data$CL <0] <- 0.1
imputed_data$SX[imputed_data$SX <0] <- 0.01
imputed_data$NC[imputed_data$NC <0] <- 0

# Manual corrections for non-plausible entries
imputed_data["E.hyacinthina","CL"]<- 0.2
imputed_data["B.fervidus","CL"]<- 0.3
imputed_data["T.angustula","R"]<- 0.75
imputed_data["C.capitata","R"]<- 0.75
imputed_data["C.capitata","OR"]<- 1

# Write imputed data to csv
write.csv(x = imputed_data,file = "imputed_data.csv")

################### Phylogenetically corrected PCA ################### 

# z score standardization of data
imputed_data_norm <- scale(imputed_data)

# Load phylogenetic tree
tree <- read.newick("tree.nwk")
plot(tree)
name.check(tree, imputed_data_norm)

# perform pPCA using phytools
pca.phyl <- phyl.pca(tree, imputed_data_norm, mode = "cov", method = "lambda")
biplot.phyl.pca(pca.phyl)
write.csv(pca.phyl$S, "pca_scores.csv")

# PCA dimension significance and dummy PCA creation for visualization
# Create PCA object with FactoMineR for easy plotting with factoextra
pca.dummy <- PCA(imputed_data_norm, graph = FALSE)
pca.dummy[["ind"]][["coord"]] <- pca.phyl$S  # Insert original scores
pca.dummy[["var"]][["coord"]] <- pca.phyl$L  # Insert original loadings
sdv <- unlist(summary(pca.phyl)[1], use.names = FALSE)
pca.dummy[["svd"]][["vs"]] <- sdv
#contribution of traits to pcs
pca.dummy[["var"]][["cos2"]]<- pca.phyl$L^2 

# Data Visualization with ggplot2
sociality <- data[, "sociality"]
taxa <- data[, "taxa"]
PC1 <- paste("PC1 (", round(summary(pca.phyl)[[2]][2]*100, digits = 0), "%)", sep = "")
PC2 <- paste("PC2 (", round(summary(pca.phyl)[[2]][5]*100, digits = 0), "%)", sep = "")

# Create a color palette and shapes for plotting
sociality=data[,"sociality"]
# Special shapes for honey bees and melipona
sociality = c(sociality,rep("hgh",4),rep("sol",10))
taxa=data[,"taxa"]
taxa = c(taxa,rep("hal",4),rep("meg",10))

col <- c("#CC79A7", "#009E73","#E7B800", 
         "#0072B2","#000000", "#D55E00")

labels=c('Solitary','Communal','Subsocial',
         'Parasocial','Primitively eusocial',
         'Advanced eusocial')

shape = c(25,25,18,24,79,16,15)

# Arrange data to match unique shapes
basic_plot <- fviz_pca_ind(pca.dummy,axes = c(1,2), label="all")
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 

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

ggsave("pca.svg",dpi=300,width = 22, height = 14)


#plot variables 1+2
plot.PCA(pca.dummy, axes=c(1,2), choix="var",
         label="all",invisible = "quali")+
  labs(x =PC1 , y = PC2,title = NULL)+
  theme(axis.title.x.bottom = element_text(hjust = 0.5,size = 14),
        axis.title.y.left = element_text(hjust = 0.5,size = 14 ))

ggsave("var12.svg",dpi=300,width = 8, height = 8)

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

# Custumize UMAP parameters
custom.settings = umap.defaults
custom.settings$n_neighbors = 10
custom.settings$min_dist = 0.5

#Plot UMAP
svg("umap_0.510.svg", width = 22, height = 14)
data.umap=umap(pca.phyl$S[,1:10],preserve.seed =F,config = custom.settings)
data.umap$layout = rbind(data.umap$layout, data.umap$layout[c(20:23,28:37),])
plot.umap(data.umap,sociality,shape = as.numeric(shp))
dev.off()

#######subclusters test to show bombus seperation
# data.umap=umap(pca.phyl$S[,1:10],preserve.seed =F,config = custom.settings)
# data.umap$layout = rbind(data.umap$layout, data.umap$layout[c(20:23,28:37),])
# sociality.test=sociality
# sociality.test[c(52, 54:57)]="sol" #bombus
# plot.umap(data.umap,sociality.test,shape = as.numeric(shp))


#########################################################
####################### calculate convex hull############

# load data
data <- read.csv("data.csv",row.names = "X")
pca_scores <- read.csv("pca_scores.csv", row.names = "X")

# pca_scores=pca.phyl$S
# verify groups of corbiculate bees and remaining species
group_CD <- data["taxa"]=="api" | data["taxa"]=="mel" | data["taxa"]=="bom"
group_AB <- !group_CD

# subset 4 pcs
group_CD <- as.matrix(pca_scores[group_CD,1:4])
group_AB <- as.matrix(pca_scores[group_AB,1:4])


# Convex hull for unscaled data
hull_AB <- convhulln(group_AB, output.options = TRUE, options = "FA")
hull_CD <- convhulln(group_CD, output.options = TRUE, options = "FA")

# Compute ratios
ratio <- hull_CD$vol / hull_AB$vol
cat("Ratio between superorganisms and remaining species:", ratio, "\n")

#calculate number of significant pca dimensions
summary(pca.phyl)
x=as.numeric(unlist(summary(pca.phyl)[1]))
x = round(x^2/sum(x^2),digits = 2)
bsDimension(x)



################# supplementary information #######################
############### contribution of traits to PCs  ###################
loadings <- pca.phyl$L
contributions <- sapply(1:ncol(loadings), function(i) {
  (loadings[, i]^2) * 100 / sum(loadings[, i]^2)
})

contributions <- as.data.frame(contributions)
pca.dummy[["var"]][["contrib"]]<-contributions

#pc variance
PC1=paste("PC1 (",round(summary(pca.phyl)[[2]][2]*100,digits = 0),"%)",sep = "")
PC2=paste("PC2 (",round(summary(pca.phyl)[[2]][5]*100,digits = 0),"%)",sep = "")
PC3=paste("PC3 (",round(summary(pca.phyl)[[2]][8]*100,digits = 0),"%)",sep = "")
PC4=paste("PC4 (",round(summary(pca.phyl)[[2]][11]*100,digits = 0),"%)",sep = "")

# plot
a = fviz_contrib(pca.dummy, choice = "var", axes = 1, top = 10)+
  labs( y = "Contribution (%)",title = NULL)+
  theme(
    axis.text = element_text(size = 12),
    axis.title.y.left = element_text(hjust = 0.5,size = 14 ),
  )

b = fviz_contrib(pca.dummy, choice = "var", axes = 2, top = 10)+
  labs( y = "Contribution (%)",title = NULL)+
  theme(
    axis.text = element_text(size = 12),
    axis.title.y.left = element_text(hjust = 0.5,size = 14 ),
  )

c = fviz_contrib(pca.dummy, choice = "var", axes = 3, top = 10)+
  labs( y = "Contribution (%)",title = NULL)+
  theme(
    axis.text = element_text(size = 12),
    axis.title.y.left = element_text(hjust = 0.5,size = 14 ),
  )

d = fviz_contrib(pca.dummy, choice = "var", axes = 4, top = 10)+
  labs( y = "Contribution (%)",title = NULL)+
  theme(
    axis.text = element_text(size = 12),
    axis.title.y.left = element_text(hjust = 0.5,size = 14 ),
  )

plot_grid(a, b, c, d, labels = c(PC1, PC2, PC3, PC4),
          label_x=0.4,label_y=1)

ggplot2::ggsave("contribution.svg",dpi=300,width = 8, height = 8)


