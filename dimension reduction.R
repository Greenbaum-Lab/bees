
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

#load and prepare data
data=read.csv("data.csv")
#put row names
rownames(data)=data[,"name"]
#delete row names from columns
data=(data[,-4])
#create sociality vector
species=paste(data[,"genus"],data[,"species"],sep = "_")
species=str_replace_all(string=species, pattern=" ", repl="")

#take numeric columns only
num_cols <- unlist(lapply(data, is.numeric))         # Identify numeric columns
dataq <- data[ , num_cols]                        # Subset numeric columns of data


## log on colony size avg+max+foundation, fecundity)
dataq=cbind.data.frame(dataq[,-c(grep("CF|CM|CA|QF",colnames(dataq)))],
                       log10(dataq[,c(grep("CF|CM|CA|QF",colnames(dataq)))]))

##"loo" for leave- one-out or "Kfold" cross-validation
n.dim <- estim_ncpPCA(dataq,  scale = TRUE, method = "Regularized",
                      method.cv = "loo")

plot(n.dim$criterion~names(n.dim$criterion),xlab="number of dimensions")
n.dim$ncp
#impute using iteration pca the NA, no need to scale before!
#i used the ncp with the least negative values

res.comp <- imputePCA(dataq, ncp = n.dim$ncp,scale = T,method = "regularized")
res.comp=res.comp$completeObs

#which col has negative values
res.comp=as.data.frame(res.comp) 
which(res.comp<0) #7 entries
res.comp %>% filter_all(any_vars(.< 0))
res.comp$PM[res.comp$PM <0] <- 0
res.comp$CL[res.comp$CL <0] <- 0.1
res.comp$WS[res.comp$WS <0] <- 2
res.comp$SX[res.comp$SX <0] <- 0.01
res.comp$NC[res.comp$NC <0] <- 0

#correct manually non plausible entries
res.comp["E.hyacinthina","CL"]<- 0.2
res.comp["B.fervidus","CL"]<- 0.3
res.comp["B.vitrea","CL"]<- 0.2
res.comp["B.vitrea","QL"]<- 0.4
res.comp["T.angustula","R"]<- 0.75
res.comp["C.capitata","R"]<- 0.75
res.comp["C.capitata","OR"]<- 1


#create full data from imputed pca itiration after pre-process
data.imp=scale(res.comp)
tree = read.newick("tree.nwk")
name.check(tree,data.imp)

#phylogenetic pca
pca.phyl=phyl.pca(tree,data.imp, mode= "cov",
                  method = "lambda")
pca.phyl$S=pca.phyl$S*-1
biplot.phyl.pca(pca.phyl)

write.csv(pca.phyl$S,"pca_scores.csv")
#calculate number of significant pca dimensions
summary(pca.phyl)
x=as.numeric(unlist(summary(pca.phyl)[1]))
x = x^2/sum(x^2)
bsDimension(x)

#beacause "phyl.pca" is not compatible for visualization outside the
#phytools package we need to
#create a dummy pca object to put the real data inside it
pca.dummy <-PCA(res.comp, graph =F)
#insert original scores
pca.dummy[["ind"]][["coord"]]<-pca.phyl$S
#insert original loadings
pca.dummy[["var"]][["coord"]]<-pca.phyl$L
#extract and insert sdv of PCs
sdv=as.numeric(unlist(summary(pca.phyl)[1]))
pca.dummy[["svd"]][["vs"]]<-sdv


#pc variance
PC1=paste("PC1 (",round(summary(pca.phyl)[[2]][2]*100,digits = 0),"%)",sep = "")
PC2=paste("PC2 (",round(summary(pca.phyl)[[2]][5]*100,digits = 0),"%)",sep = "")
PC3=paste("PC3 (",round(summary(pca.phyl)[[2]][8]*100,digits = 0),"%)",sep = "")
PC4=paste("PC4 (",round(summary(pca.phyl)[[2]][11]*100,digits = 0),"%)",sep = "")

par(mar=c(4,4,4,4),oma=c(2,2,2,2))
plotTree(tree);axisPhylo()

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

ggsave("C:\\Users\\mukid\\Downloads\\new\\pca.svg",dpi=300,width = 22, height = 14)


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



umap.defaults
custom.settings = umap.defaults
custom.settings$n_neighbors = 20
custom.settings$min_dist = 0.5


svg("umap.svg", width = 22, height = 14)


#ppca umap+add apis circles
data.umap=umap(pca.phyl$S[,1:10],preserve.seed =F,config = custom.settings)
data.umap$layout = rbind(data.umap$layout, data.umap$layout[c(1:4,6,9,17:24),])
plot.umap(data.umap,sociality,shape = as.numeric(shp))

dev.off()

#subclusters test
# sociality.test=sociality
# sociality.test[c(35,36,38,42,43)]="sol" #bombus
# plot.umap(data.umap,sociality.test,shape = as.numeric(shp))

