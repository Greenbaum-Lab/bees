library(missMDA)
library(FactoMineR)
library(ggplot2)
library(factoextra)
library(PCDimension)
library(dplyr)
library(stringr)
library(umap)
library(phytools)
library(stringr)
library(geiger)
library(PCDimension)
library(geiger)
library(svglite)
library(cowplot)
library(picante)

#load and prepare data
# run dimension reduction to do the pca
#make the data categories
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

# arrange data and to match unique shapes
basic_plot <- fviz_pca_ind(pca.dummy,axes = c(1,2), label="all")
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 


a = ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
           aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 5) +
  labs(x = PC1, y = PC2, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)

# pca all combinations 1-4
basic_plot <- fviz_pca_ind(pca.dummy,axes = c(1,3))
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 

b = ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
          aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 5) +
  labs(x = PC1, y = PC3, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)


basic_plot <- fviz_pca_ind(pca.dummy,axes = c(1,4))
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 

c = ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
           aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 5) +
  labs(x = PC1, y = PC4, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)


basic_plot <- fviz_pca_ind(pca.dummy,axes = c(2,3))
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 

d = ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
           aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 5) +
  labs(x = PC2, y = PC3, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)

basic_plot <- fviz_pca_ind(pca.dummy,axes = c(2,4))
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 

e = ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
           aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 5) +
  labs(x = PC2, y = PC4, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)

basic_plot <- fviz_pca_ind(pca.dummy,axes = c(3,4))
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 

f = ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
           aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 5) +
  labs(x = PC3, y = PC4, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)

plot_grid(a, b, c, d, e, f, ncol = 2)

ggsave("allpcs.svg",dpi=300,width = 12, height = 12)




################uncorrected pca
imputed_data = read.csv('imputed_data.csv',row.names = 'X')
sociality=data[,"sociality"]
imputed_data <- cbind.data.frame(imputed_data,sociality)
#number=column of sociality
res_pca <- PCA(imputed_data, graph=F,quali.sup = 18,scale.unit = T)

#plot by tribe/sociality
PC1=paste("PC1 (",round(res_pca[["eig"]][1,2],digits = 0),"%)",sep = "")
PC2=paste("PC2 (",round(res_pca[["eig"]][2,2],digits = 0),"%)",sep = "")


#make the data categories
sociality=data[,"sociality"]
sociality = c(sociality,rep("hgh",4),rep("sol",10))

taxa=data[,"taxa"]
taxa = c(taxa,rep("hal",4),rep("meg",10))

shape = c(25,25,18,24,79,16,15)

# arrange data and to match unique shapes
basic_plot <- fviz_pca_ind(res_pca,axes = c(1,2), label="all")
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(20:23,28:37),]) 


a =  ggplot(cbind(basic_plot$data[,2:6], sociality, taxa),
            aes(x = x, y = y, col = sociality, shape = taxa, fill = sociality)) + 
  geom_point(size = 8) +
  labs(x = PC1, y = PC2, title = NULL, color = "sociality", shape = "taxa") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)

#plot variables
b = plot(res_pca,choix="var")+
  labs(x =PC1 , y = PC2,title = NULL)+
  theme(
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 12),axis.title.x.bottom = element_text(hjust = 0.5,size = 12),
    axis.title.y.left = element_text(hjust = 0.5,size = 12 ), )


plot_grid(a, b)
ggsave("pca_uncorrect.svg",dpi=300,width = 16, height = 8)


##########################uncorrected umap 
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



umap.defaults
custom.settings = umap.defaults
custom.settings$n_neighbors = 20
custom.settings$min_dist = 0.2

imputed_data = read.csv('imputed_data.csv',row.names = 'X')

svg("umap_uncorrect.svg", width = 22, height = 14)

#ppca umap+add apis circles
data.umap=umap(scale(imputed_data),preserve.seed =F,config = custom.settings)
data.umap$layout = rbind(data.umap$layout, data.umap$layout[c(20:23,28:37),])
plot.umap(data.umap,sociality,shape = as.numeric(shp))
dev.off()
######################################################
#################### Evolution ########################
library(phytools)
library(ctpm)
library(geiger)
library(ape)
library(mvMORPH)
library(openxlsx)
library(nlme)
library(picante)
library(ggplot2)
library(dplyr)
library(svglite)
library(stringr)
library(AFM)


#####################ancestral state reconstruction all pcs
library(phytools)

par(oma = c(4, 4, 2,2))
tree <- read.newick("tree.nwk") # Assuming tree is loaded here for context
pca_scores <- read.csv("pca_scores.csv", row.names = "X")

# 78-hal black
# 79-apida
# 83-eug
# 93-mel
# 120-bom
# 136-apis


# Define sub-trees and their states
tree_col <- tree
tree_col <- paintSubTree(tree_col, node = 81, state = "1", stem = F) # root
tree_col <- paintSubTree(tree_col, node = 95, state = "2", stem = T) # euglossini
tree_col <- paintSubTree(tree_col, node = 103, state = "6", stem = T) # apis
tree_col <- paintSubTree(tree_col, node = 107, state = "5", stem = T) # melipona
tree_col <- paintSubTree(tree_col, node = 135, state = "4", stem = T) # bombus
tree_col <- paintSubTree(tree_col, node = 82, state = "7", stem = F) # halictidae
tree_col <- paintSubTree(tree_col, node = 151, state = "3", stem = F) # xylocopinae

# Define colors for plotting
plotTree(tree);nodelabels()
cols <- c("#000000","#CC79A7","#E7B800",  "#0072B2", "#009E73", "#00FF00",  "#D55E00")
names(cols) <- c(1, 2, 3, 4, 5, 6, 7)
plotSimmap(tree_col, cols, lwd = 3, pts = F)


par(mfrow=c(2,2))

#save plot
svg("tree1S.svg",width = 18,height = 22)

#calculate pc1
pc1=pca_scores[,1]
names(pc1) <- rownames(pca_scores)
phenogram(tree=tree_col,x = pc1, colors = cols, spread.labels = T,
          lwd = 4,spread.cost=c(1,0))
fancyTree(tree = tree,type = 'phenogram95',x=pc1,colors=cols,spread.cost=c(1,0))
dev.off()

svg("tree2S.svg",width = 18,height = 22)
pc2=pca_scores[,2]
names(pc2) <- rownames(pca_scores)
phenogram(tree=tree_col,x = pc2, colors = cols, spread.labels = T,
          lwd = 4,spread.cost=c(1,0))
fancyTree(tree = tree,type = 'phenogram95',x=pc2, colors=cols,spread.cost=c(1,0))
dev.off()

svg("tree3S.svg",width = 18,height = 22)
pc3=pca_scores[,3]
names(pc3) <- rownames(pca_scores)

phenogram(tree=tree_col,x = pc3, colors = cols, spread.labels = T,
          lwd = 4,spread.cost=c(1,0))
fancyTree(tree = tree,type = 'phenogram95',x=pc3,colors=cols,spread.cost=c(1,0))
dev.off()

svg("tree4S.svg",width = 18,height = 22)
pc4=pca_scores[,4]
names(pc4) <- rownames(pca_scores)
phenogram(tree=tree_col,x = pc4, colors = cols, spread.labels = T,
          lwd = 4,spread.cost=c(1,0))
fancyTree(tree = tree,type = 'phenogram95',x=pc4,colors=cols,spread.cost=c(1,0))
dev.off()




##################################   ASR +PCA
###change i and j and index in plot across all 6 combinations
# read data
data <- read.csv("data.csv", row.names = 'X')
pca_scores <- read.csv("pca_scores.csv", row.names = "X")
phy <- read.newick("tree.nwk")
name.check(phy, pca_scores)

index <- c(3,4)

# get ancestral states of pcs
x = as.vector(fastAnc(tree = phy,x=pca_scores[,'PC3']))
y = as.vector(fastAnc(tree = phy,x=pca_scores[,'PC4']))

# Create vectors for axes and combine to dataframe
df <- data.frame(x = c(dat[, index[1]], x), y = c(dat[, index[2]], y))

# Extract paths and create combined dataframe for plotting
path <- nodepath(phy)
a <- lapply(path, function(i) df[i, 1])
b <- lapply(path, function(i) df[i, 2])
demo <- cbind(data.frame(v1 = unlist(a), v2 = rep(seq(length(a)), lengths(a))), unlist(b))
colnames(demo) <- c("x", "species", "y")
demo <- mutate(demo, x = as.numeric(x), y = as.numeric(y))

# Read taxa, match with phylogeny, and create plot
taxa <- as.data.frame(data$taxa)
rownames(taxa) <- rownames(data)
demo$group <- taxa[demo$species,]

# Plot by lineage
cols <- c("#00FF00", "#0072B2", "#CC79A7", "#D55E00", "#000000", "#009E73", "#E7B800")

p6=ggplot() +
  geom_path(data = demo, aes(x = x, y = y, colour = group, group = species),
            arrow = arrow(length = unit(0.3, "cm"))) +
  geom_point(data = demo, aes(x = x, y = y), color = "grey57", size = 0.5) +
  geom_point(data = dat, aes(x = PC3, y = PC4), color = "black", size = 0.5) +
  theme_void() +
  labs(x = "PC1", y = "PC2") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "none",
        panel.border = element_rect(fill = NA)) +
  scale_color_manual(values = cols)


plot_grid(ncol = 2,p1,p2,p3,p4,p5,p6)
ggsave("SUP_asr_combo.svg",dpi=300,width = 12, height = 16)
