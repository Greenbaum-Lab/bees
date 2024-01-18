
library(PhylogeneticEM)
library(phytools)


# beta-primary optimum
# alpha-selction strength
# there is no convergence evolution here. each trnasition induces a specific optimum,

### for NA sensitivity
# phy=phy_test
# dat=pca_test$S[,1:2]
# dat=(dataq)
### for trait sensitivity
# dat=grouped_summary[,c(2,4)]

dat=read.csv("pca_scores.csv",row.names = "X")
dat=dat[match(phy$tip.label,rownames(dat)),]
dat=dat[,1:2] # choose number of PCs 
dat=t(dat)


phy=read.newick("tree.nwk")
phy=force.ultrametric(phy)

is.ultrametric(phy)
is.rooted(phy)
is.binary(phy)

plotTree((phy));nodelabels()

phy=rotateNodes(tree = phy,nodes = 113)



res=PhyloEM(phylo = phy,Y_data = dat,process =  method.selection = c("LINselect", "Djump"),
            stationary.root = T,random.root = T)


svglite::svglite("shifts_dfvds.svg")

plot(res) 
dev.off()

## Plot selected solution (DDSE)
par(mar=c(4,4,4,4))
plot(res, method.selection = "LINselect")
plot(res, method.selection = "Djump") 

summary(res)
svglite::svglite("criterion.svg",width = 8,height = 4)

par(mfrow = c(1, 2))

plot_criterion(res, "LINselect",select.col = "green")
plot_criterion(res, "Djump",select.col = "green")


?plot_criterion
par(mar=c(4,4,4,4))
multi_solutions <- params_process(res, K = 10)
svglite::svglite("shifts_4pcs_multi.svg")
multi_solutions <- equivalent_shifts(phy, multi_solutions)

plot(multi_solutions, show_shifts_values = F, shifts_cex = 0.3)

dev.off()










#######calculate space in phenotypic space
library(sp)
library(ggplot2)
library(dplyr)
library(ggforce)

data=read.csv('data.csv')

dat=read.csv("pca_scores.csv",row.names = "X")
dat=cbind(dat[,c(2,3)],data$opt2) # choose number of PCs 
colnames(dat)=c('PC1','PC2','group')


# Function to calculate convex hull and return coordinates
hull_coords <- function(df) {
  ch <- chull(df$PC1, df$PC2)
  df[ch, ]
}

# Calculate convex hull for each group
hulls <- dat %>%
  group_by(group) %>%
  do(hull_coords(.))

# Plot
ggplot(dat, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  geom_polygon(data = hulls, aes(fill = group, group = group), alpha = 0.5) +
  theme_minimal()

# ggsave(plot = p,"pca_polygon_full.svg",width = 10,height = 6)


# Function to calculate the area of the polygon
polygon_area <- function(p) {
  n <- nrow(p)
  area <- 0
  for(i in 1:(n-1)){
    area <- area + p[i, 'PC1']*p[i+1, 'PC2'] - p[i+1, 'PC1']*p[i, 'PC2']
  }
  area <- area + p[n, 'PC1']*p[1, 'PC2'] - p[1, 'PC1']*p[n, 'PC2']
  return(abs(area) / 2)
}


metrics <- dat %>%
  group_by(group) %>%
  do(data.frame(area = polygon_area(hull_coords(.))))

metrics = as.data.frame(metrics)

  # Correctly name the columns
colnames(metrics) <- c("group", "area")

# Print metrics
print(metrics)

a=c(14.026,8.591,9.844,7.915,6.091,7.256)
b=c(12.783,13.094,20.06,18.364,24.813,17.709)

mean(a)
mean(b)

sd(b)
t.test(a,b)
