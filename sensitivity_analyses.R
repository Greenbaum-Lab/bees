
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
library(picante)
library(svglite)
library(ggpubr)
library(gridExtra)
library(abs)

#########################Adjusting the inclusion  threshold 
###explore the missing species
tree=read.newick("tree.nwk")
data=read.csv("data.csv",row.names = "X")

data_test = data

#take numeric columns only
num_cols <- unlist(lapply(data_test, is.numeric))         # Identify numeric columns
dataq_test <- data_test[ , num_cols]                        # Subset numeric columns of data

# create vector with number NA for each species
p.na=sum(is.na(dataq_test))/(nrow(dataq_test)*ncol(dataq_test))
p.na

for (i in 1:nrow(dataq_test)) {
  p.na[i]=sum(is.na(dataq_test[i,]))
  
}

#identify species with more than X NA 
#0.15-6%na; 0.2-8%na; 0.25;10%na
discard=which(p.na>ncol(dataq_test)*0.15)
length(discard)
#check who is removed
dataq_test[which(p.na>ncol(dataq_test)*0.25),]

#remove species from data
data_test=data_test[-c(discard),]
dataq_test=dataq_test[-c(discard),]

## log on colony size avg+max+foundation, fecundity)
columns_to_log <- c("CF", "CM", "CA", "QF")
dataq_test[columns_to_log] <- lapply(dataq_test[columns_to_log], log10)


# impute data
n.dim <- estim_ncpPCA(dataq_test,  scale = TRUE, method = "Regularized",
                      method.cv = "loo")

plot(n.dim$criterion~names(n.dim$criterion),xlab="number of dimensions")
sensitiv <- imputePCA(dataq_test, ncp = n.dim$ncp,scale = T,method = "regularized")
sensitiv = as.data.frame(sensitiv$completeObs) 

which(sensitiv<0) 

# choose based on amount nas
# 0.25
# sensitiv["E.viridissima","CL"]<- 0.2
# sensitiv["L.limao","PM"]<- 0
# sensitiv["T.spinipes","SX"]<- 0.01
# sensitiv["E.hyacinthina","CL"]<- 0.2
# sensitiv["T.angustula","R"]<- 0.75

# 0.2
# sensitiv["E.hyacinthina","CL"]<- 0.2
# sensitiv["C.australensis","CL"]<- 0.2
# sensitiv["B.vitrea","CL"]<- 0.2
# sensitiv["M.bicolor","WS"]<- 1
# sensitiv["B.vitrea","QL"]<- 0.4
# sensitiv["E.tridentata","PM"]<- 0


# 0.15
sensitiv["M.bicolor","WS"]<- 1
sensitiv["E.hyacinthina","CL"]<- 0.2
sensitiv["C.australensis","CL"]<- 0.2
sensitiv["E.tridentata","PM"]<- 0
# 

#scale data and match new data to phylo
sensitiv=as.data.frame(scale(sensitiv)) 
sensitiv = match.phylo.data(tree,sensitiv)$data
phy_test = match.phylo.data(tree,sensitiv)$phy
name.check(phy_test,sensitiv)


#phylogenetic pca
pca_test=phyl.pca(phy_test,sensitiv, mode= "cov",
                  method = "lambda")
# pca_test$S[,1]=pca_test$S[,1]*-1

biplot.phyl.pca(pca_test)

#beacause "phyl.pca" is not compatible for visualization outside the
#phytools package i need to
#create a dummy pca object to sneak the real data inside it
pca.dummy_test <-PCA(dataq_test, graph =T)
#insert original scores
pca.dummy_test[["ind"]][["coord"]]<-pca_test$S
#insert original loadings
pca.dummy_test[["var"]][["coord"]]<-pca_test$L
#extract and insert sdv of PCs
sdv=as.numeric(unlist(summary(pca_test)[1]))
pca.dummy_test[["svd"]][["vs"]]<-sdv
#contribution of traots to pcs
pca.dummy_test[["var"]][["cos2"]]<- pca_test$L^2


#pc variance
PC1=paste("PC1 (",round(summary(pca_test)[[2]][2]*100,digits = 0),"%)",sep = "")
PC2=paste("PC2 (",round(summary(pca_test)[[2]][5]*100,digits = 0),"%)",sep = "")
PC3=paste("PC3 (",round(summary(pca_test)[[2]][8]*100,digits = 0),"%)",sep = "")
PC4=paste("PC4 (",round(summary(pca_test)[[2]][11]*100,digits = 0),"%)",sep = "")

# Reorder datasets to match rows
pca.dummy_test$ind$coord <- 
  pca.dummy_test$ind$coord[match(rownames(data_test),
                                 rownames(pca.dummy_test$ind$coord)), ]


melipona <- which(data_test$genus == 'Melipona')
apis <- which(data_test$genus == 'Apis')

#make the data categories
sociality_test=data_test[,"sociality"]
sociality_test = c(sociality_test,rep("hgh",length(apis)),rep("sol",length(melipona)))

taxa_test=data_test[,"taxa"]
taxa_test = c(taxa_test,rep("hal",length(apis)),rep("meg",length(melipona)))

col <- c("#CC79A7", "#009E73","#E7B800", 
         "#0072B2","#000000", "#D55E00")

labels=c('Solitary','Communal','Subsocial',
         'Parasocial','Primitively eusocial',
         'Advanced eusocial')

shape = c(25,25,18,24,79,16,15)

# arrange data and to match unique shapes
basic_plot <- fviz_pca_ind(pca.dummy_test,axes = c(1,2))
basic_plot$data = rbind(basic_plot$data,basic_plot$data[c(apis,melipona),]) 

# plot
ggplot(cbind(basic_plot$data[,2:6], sociality_test, taxa_test),
       aes(x = x, y = y, col = sociality_test, shape = taxa_test, fill = sociality_test)) + 
  geom_point(size = 8) +
  labs(x = PC1, y = PC2, title = NULL, color = "sociality_test", shape = "taxa_test") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)



# use the original pca to combine and plot both 
dat = pca.phyl$S[,1:2]
colnames(dat) = c("x","y")
taxa = data[,"taxa"]

### change number of species
social_combo = c(rep("gr",80), sociality_test)
taxa_combo = c(taxa, taxa_test)
pca_combo = rbind(dat, basic_plot$data[,2:3])


#add group column to match every species
pca_combo$group <- substring(row.names(pca_combo), 1, 6)


col <- c("#CC79A7","grey","#E7B800", "#009E73", 
         "#0072B2","#000000", "#D55E00")



# plot
ggplot(cbind(pca_combo, social_combo, taxa_combo),
       aes(x = x, y = y, col = social_combo, shape = taxa_combo, fill = social_combo)) + 
  geom_point(size = 8) +
  geom_path(data=pca_combo, aes(group=group ))+
  labs(x = PC1, y = PC2, title = NULL, color = "social_combo", shape = "taxa_combo") + 
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme_test() + 
  guides(color = FALSE, shape = FALSE, fill = FALSE)+ 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        text = element_text(size = 30)) 




ggsave("sens_na_6.svg",dpi=300,width = 16, height = 10)

p.na=sum(is.na(dataq_test))/(nrow(dataq_test)*ncol(dataq_test))
p.na





##########################################################
################## confidence elipses ####################

# Load your data (replace this with your own data)
data = read.csv("imputed_data.csv",row.names = "X")
data_sens <- as.data.frame(scale(data))
tree=read.newick("tree.nwk")
name.check(tree,data_sens)

# Define a function to run PCA on the data and extract the first 4 components
run_pca <- function(data) {
  pca_phyl <- phyl.pca(tree, data, mode = "cov", method = "lambda")

    if (pca_phyl[["S"]]["A.mellifera",1]<0) {
    pca_phyl$S[,1] <- pca_phyl$S[,1]*-1
    print("i flipped pc1")
  }
  
  if (pca_phyl[["S"]]["T.angustula",2]>0) {
    pca_phyl$S[,2] <- pca_phyl$S[,2]*-1
    print("i flipped pc2")
  }
  
  if (pca_phyl[["S"]]["A.mellifera",3]>0 | pca_phyl[["S"]]["M.bicolor",3]>4) {
    pca_phyl$S[,3] <- pca_phyl$S[,3]*-1
    print("i flipped pc3")
  }
  

  if ( pca_phyl[["S"]]["B.terrestris",4]>0 & pca_phyl[["S"]]["M.bicolor",4]< abs(2)) {
      
    pca_phyl$S[,4] <- pca_phyl$S[,4]*-1
    print("i flipped pc4")
  }
  
  pca_phyl$S[,1:4]
  
}

# Define an empty data frame to store the results
df <- data.frame()

# Loop through each column of the original data frame
for (i in 1:ncol(data_sens)) {
  # Remove the i-th column from the data frame
  data_subset <- data_sens[, -i]
  
  # Run PCA on the subset of the data and extract the first 4 components
  pca_subset <- run_pca(data_subset)
  
  # Combine the results with the original data frame
  df_subset <- cbind(pca_subset, group = rownames(pca_subset))
  
  df=rbind(df,df_subset)
}


# View the resulting data frame

df[,1]=as.numeric(df[,1])
df[,2]=as.numeric(df[,2])
df[,3]=as.numeric(df[,3])
df[,4]=as.numeric(df[,4])
plot(df[,3:4])

 # scatter plot of pc12 with confidence ellipses at 95%
ggscatter(df, x = "PC1", y = "PC2",color = "white",alpha=0.2)+
  theme_bw()+theme(legend.position = "none")+
  stat_conf_ellipse(aes(color = group))
ggsave("sens_pca1.svg", width = 15,height = 10,dpi = 300)

# scatter plot of pc34 with confidence ellipses at 95%
 ggscatter(df, x = "PC3", y = "PC4",color = "white",alpha=0.2)+
  theme_bw()+theme(legend.position = "none")+
  stat_conf_ellipse(aes(color = group))


ggsave("sens_pca1.svg", width = 15,height = 18,dpi = 300)






##########################################################
################ test taxonomic sampling bias ############

# Load taxa to remove
taxa=read.csv("data.csv",row.names = "X")
taxa_to_remove = taxa[,'taxa'] == 'api' | taxa[,'taxa'] == 'mel' | taxa[,'taxa'] == 'bom'
n_remove = (sum(taxa_to_remove)-sum(!taxa_to_remove))
taxa_to_remove = which(taxa_to_remove)

# Load data and scale data before pca
data = read.csv("imputed_data.csv",row.names = "X")
data <- as.data.frame(scale(data))
tree=read.newick("tree.nwk")
name.check(tree,data)

# pca_phyl <- phyl.pca(tree, data, mode = "cov", method = "lambda")
# x=pca_phyl$S
# plot(pca_phyl$S[,1:2]*-1)
# median(pca_phyl$S[,2]>0)
# mean(pca_phyl$S[,4]*-1)
# 

flip_pca_axes <- function(res_pca) {
  # Check and flip PC1

  if (res_pca["D.novaeangliae", 1] > 0) {
    res_pca[, 1] <- -res_pca[, 1]
    print("I flipped PC1")
  }
  
  # Check and flip PC2
  if (res_pca['X.virginica', 2] < 0) {
    res_pca[, 2] <- -res_pca[, 2]
    print("I flipped PC2")
  }
  
  # Check and flip PC3
  if (res_pca['M.genalis', 3] > 0) {
    res_pca[, 3] <- -res_pca[, 3]
    print("I flipped PC3")
  }
  
  # Check and flip PC4 based on mean
  if (res_pca['E.townsendi', 4] > 0) {
    res_pca[, 4] <- -res_pca[, 4]
    print("I flipped PC4")
  }
  
  return(res_pca)
}

# make a loop to subset data with equal taxonomic sampling
final_df = data.frame()
iterations = 100

for (i in 1:iterations) {
  
  # sample species to remove based on index
  rows_remove = sample(taxa_to_remove,n_remove,replace = F)
  data_subset = data[-rows_remove,]

  # match data and tree
  subsets = match.phylo.data(tree,data_subset)

  # run pca
  pca_phyl <- phyl.pca(subsets$phy, subsets$data, mode = "cov", method = "lambda")
  res_pca = pca_phyl$S[,1:4]
  print(i)

  res_pca = flip_pca_axes(res_pca)
  # add group colum with species names
  res_pca <- cbind(res_pca, group = rownames(res_pca))
  
  # merge all dataframes
  final_df=rbind(final_df,res_pca)
  
}


for (i in 1:4) {
  final_df[,i]=as.numeric(final_df[,i])
}

plot(final_df[1:56,3:4])


# scatter plot of pc12 with confidence ellipses at 95%
ggscatter(final_df, x = "PC1", y = "PC2",color = "white",alpha=0.2) +
  theme_bw()+theme(legend.position = "none")+
  stat_conf_ellipse(aes(color = group))
ggsave("tax_bias12.svg", width = 15,height = 10,dpi = 300)


# pc3+4 are not very good. need to adjust axis rotation
ggscatter(final_df, x = "PC3", y = "PC4",color = "white",alpha=0.2) +
  theme_bw()+theme(legend.position = "none")+
  stat_conf_ellipse(aes(color = group))
ggsave("tax_bias34.svg", width = 15,height = 10,dpi = 300)



######## plot means
group_means <- final_df %>%
  group_by(group) %>%
  summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2))

taxa=taxa[,c('sociality','taxa')]
taxa$species <- rownames(taxa)

# Join 'group_means' with the 'sociality' and 'taxa' information from 'data'
group_means <- group_means %>%
  left_join(taxa, by = c("group" = "species"))

col <- c("#CC79A7", "#009E73","#E7B800", "#0072B2","#000000", "#D55E00")
labels <- c('Solitary', 'Communal', 'Subsocial', 'Parasocial', 'Primitively eusocial', 'Advanced eusocial')
shape <- c(25, 25, 18, 24, 79, 16, 15)

ggplot(group_means[c(56:62,65,64,67),], aes(x = mean_PC1, y = mean_PC2, color = sociality, shape = taxa, fill = sociality)) +
  geom_point(size = 8) +  # Plot mean points
  labs(x = PC1, y = PC2, title = NULL, color = "Sociality", shape = "Taxa") +
  theme_light() +
  scale_shape_manual(values = shape) +
  scale_color_manual(labels = labels, values = col) +
  scale_fill_manual(labels = labels, values = col) +
  theme(legend.position = "none")  # Adjust legend position as needed
ggsave("tax_bias12_mean.svg", width = 15,height = 10,dpi = 300)

################################################################
########## phylogenetic tree sensitivity#####

trees = read.newick('all_trees.nwk')

# Load data and scale data before pca
data = read.csv("imputed_data.csv",row.names = "X")
data <- as.data.frame(scale(data))

name.check(tree,data)



# make a loop to subset data with equal taxonomic sampling
final_df = data.frame()
iterations = 100

for (i in 1:iterations) {
  
  tree = sample(x=trees,size=1,replace = F)[[1]]

  # match data and tree
  # subsets = match.phylo.data(tree,data)
  
  # run pca
  pca_phyl <- phyl.pca(tree, data, mode = "cov", method = "lambda")
  res_pca = pca_phyl$S[,1:4]
  print(i)
  
  res_pca = flip_pca_axes(res_pca)
  # add group colum with species names
  res_pca <- cbind(res_pca, group = rownames(res_pca))
  
  # merge all dataframes
  final_df=rbind(final_df,res_pca)
  
}


for (i in 1:4) {
  final_df[,i]=as.numeric(final_df[,i])
}

plot(final_df[1:200,1:2])


# scatter plot of pc12 with confidence ellipses at 95%
ggscatter(final_df, x = "PC3", y = "PC4",color = "white",
          alpha=0.2,mean.point = T) +
  theme_bw()+theme(legend.position = "none")+
  stat_conf_ellipse(aes(color = group))
ggsave("phyl_bias12.svg", width = 15,height = 10,dpi = 300)


