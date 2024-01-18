library(phytools)
library(ctpm)
library(ape)
library(openxlsx)
library(nlme)
library(picante)
library(ggplot2)
library(dplyr)
library(svglite)
library(stringr)
library(AFM)
library(geiger)
library(AFM)
library(tidyr)




#####################ancestral state reconstruction PC1
par(oma=c(4,4,0,0))

# 78-hal black
# 79-apida
# 83-eug
# 93-mel
# 120-bom
# 136-apis
x=tree
plotTree(tree);nodelabels()
x<-paintSubTree(x,node = 87,state="1",stem=T,anc.state = "7") #xylocopinae
x<-paintSubTree(x,node = 88,state="2",stem=F) #euglossini
x<-paintSubTree(x,node = 89,state="5",stem=T) #apis
x<-paintSubTree(x,node = 93,state="4",stem=T) #melipona
x<-paintSubTree(x,node = 94,state="3",stem=T) #bombus

x<-paintSubTree(x,node = 79,state="6",stem=F) #halictidae

# x<-paintSubTree(x,node = 92,state="3",stem=F)
# x<-paintSubTree(x,node = 139,state="6",stem=F)
# x<-paintBranches(x,edge = 93,state = "4")

cols<-c("#E7B800","#CC79A7","#0072B2","#009E73","#00FF00","#000000","grey58")
names(cols)<-c(1,2,3,4,5,6,7)
plotSimmap(x,cols,lwd=3,pts=F)

cols<c("#E7B800","#CC79A7","#0072B2","#009E73","#33FF33","#000000","grey58")
#starting value based on m.rotunda
#calculate pc1
pc1=pca.phyl$S[,1]
#take value of m. rotunda
nod = pca.phyl$S["M.rotundata",1]
names(nod) = "78"
reconstruct <-fastAnc(tree ,pc1,anc.states=nod)

par(mar=c(2,2,2,2))
#save plot
svg("reconstruct_pc1.svg", width = 18,height = 12)
phenogram(x,c(pc1,reconstruct), colors=cols,spread.labels=F,lwd=4,fsize = 0)

dev.off()


##################################   ASR +PCA
# read data
data <- read.csv("data.csv", row.names = "name")
dat <- read.csv("pca_scores.csv", row.names = "X")
tree <- read.newick("tree.nwk")

# match data with phylo
demo <- match.phylo.data(phy = tree, data = dat)
dat <- demo$data
phy <- demo$phy
name.check(phy, dat)

# function to create vectors
vec_creator <- function(dat, index, phy, species = "M.rotundata") {
  demo = dat[,index]
  names(demo) = rownames(dat)
  x = as.vector(fastAnc(phy, demo, anc.states = setNames(c(dat[species, index]), c("78"))))
  return(x)
}

index = c(2,3)
# make vectors of axis x and y
x <- vec_creator(dat, index[1], phy)
y <- vec_creator(dat, index[2], phy)

# Combine vectors to df
df <- data.frame(x = c(dat[,index[1]], x), y = c(dat[,index[2]], y))

# Extract the nodes each species go through
path <- nodepath(phy)

# create path for each species for x and y
a <- lapply(path, function(i) df[i, 1])
b <- lapply(path, function(i) df[i, 2])

# combine results
a <- data.frame(v1 = unlist(a), v2 = rep(seq(length(a)), lengths(a)))
b <- data.frame(v1 = unlist(b), v2 = rep(seq(length(b)), lengths(b)))
demo <- cbind(a, b[,1])
colnames(demo) <- c("x", "species", "y")

# make values numeric
demo <- mutate(demo, x = as.numeric(x), y = as.numeric(y))

# read taxa
taxa <- as.data.frame(data$taxa)
rownames(taxa) <- rownames(data)
taxa <- match.phylo.data(phy = tree, data = taxa)$data
taxa <- as.vector(taxa$`data[res$phy$tip.label, ]`)

# create a new 'species' column based on 'taxa'
demo=demo
demo$group <- taxa[demo$species]

# plot by lineage
cols <- c("#00FF00", "#0072B2", "#CC79A7", "#000000", "grey58", "#009E73", "#E7B800")




ggplot()+
  geom_path(data = demo,
            aes(x = x, y = y, colour = group,group=species),
            arrow = arrow(length = unit(0.3, "cm"))) +
  
  geom_point(data = demo, aes(x = x, y = y), color = "grey57", size = 0.5) +  # first scatter plot
  geom_point(data = dat, aes(x = PC1, y = PC2), color = "black", size = 0.5) + # second scatter plot
  theme_void() +
  labs(x = "PC1" , y = "PC2", title = NULL) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none",
        panel.border = element_rect(fill = NA),
        legend.text = element_text(size = 0),
        legend.title = element_text(size = 0),
        axis.text = element_text(size = 0),
        axis.title.x.bottom = element_text(hjust = 0.5, size = 14),
        axis.title.y.left = element_text(hjust = 0.5, angle = 90, size = 14)) +
  scale_color_manual(values = cols)





############size of steps
# Calculate vector lengths for each species

df_with_vectors <- demo %>%
  arrange(species) %>%
  group_by(species) %>%
  mutate(dx = lead(x) - x, dy = lead(y) - y,
         vector_length = sqrt(dx^2 + dy^2)) %>%
  slice(-n()) %>%
  mutate(combined_group = ifelse(group %in% c("api", "bom", "mel"), "super", "simple"))


# Extract the edge matrix and edge lengths
edge_matrix <- tree$edge
edge_lengths <- tree$edge.length

# Combine them into a new data frame
combined_edge_matrix <- cbind(edge_matrix, edge_lengths)
colnames(combined_edge_matrix) <- c("Parent", "Child", "Edge_Length")

# Initialize an empty list for storing species path data frames
species_paths_df_list <- vector("list", length(tree$tip.label))

# Function to calculate the path for each species
calculate_path <- function(species_index, tree, combined_edge_matrix) {
  # Find the node path from the root to the current species
  path <- nodepath(tree, from = 78, to = species_index)  # Adjust 'from' as needed
  
  # Corresponding rows in the combined_edge_matrix for each segment in the path
  path_rows <- sapply(1:(length(path) - 1), function(j) {
    which(combined_edge_matrix[, "Parent"] == path[j] & combined_edge_matrix[, "Child"] == path[j + 1])
  })
  
  # Extract edge lengths for this path
  path_edge_lengths <- combined_edge_matrix[path_rows, "Edge_Length"]
  
  # Create a data frame for this species path
  species_path_df <- setNames(as.data.frame(t(path_edge_lengths)), paste("Edge", seq_along(path_edge_lengths), sep = "_"))
  
  return(species_path_df)
}

# Calculate paths for each species and find the maximum path length
species_paths_df_list <- lapply(seq_along(tree$tip.label), calculate_path, tree, combined_edge_matrix)
max_path_length <- max(sapply(species_paths_df_list, ncol))

# Adjust all data frames to have the same number of columns
species_paths_df_list <- lapply(species_paths_df_list, function(df) {
  num_missing_cols <- max_path_length - ncol(df)
  if (num_missing_cols > 0) {
    missing_cols <- setNames(matrix(NA, nrow = 1, ncol = num_missing_cols), paste("Edge", (ncol(df) + 1):max_path_length, sep = "_"))
    return(cbind(df, missing_cols))
  } else {
    return(df)
  }
})



# Ensure all data frames in the list have the same column names
max_cols <- max(sapply(species_paths_df_list, ncol))  # Maximum number of columns in any data frame
col_names <- paste("Edge", 1:max_cols, sep = "_")  # Create column names

species_paths_df_list <- lapply(species_paths_df_list, function(df) {
  missing_cols <- setdiff(col_names, colnames(df))
  if (length(missing_cols) > 0) {
    # Add missing columns as NA
    df[missing_cols] <- matrix(NA, nrow = nrow(df), ncol = length(missing_cols))
  }
  df <- df[, col_names]  # Ensure consistent column order
  return(df)
})

# Now combine all species paths into a single data frame
all_species_paths_df <- do.call(rbind, species_paths_df_list)
rownames(all_species_paths_df) <- tree$tip.label  # Set row names



# Vectorized approach to create all_lengths_vector
all_lengths_vector <- unlist(lapply(seq_len(nrow(all_species_paths_df)), function(i) {
  na.omit(unlist(all_species_paths_df[i, ]))
}))

# Adding the vector to the existing data frame
new <- cbind(df_with_vectors, all_lengths_vector)

# Convert species to character and rounding edge_length
names(new)[names(new) == "...9"] <- "edge_length"
new$species <- as.character(new$species)
new$edge_length <- round(new$edge_length,digits = 1)

# Renaming and creating new columns

new <- transform(new, step = vector_length / edge_length)
names(new)[names(new) == "V10"] <- "step"


# Replace edge_length values of 0 with 1
# new$edge_length[new$edge_length == 0] <- 1

# Optimized replicate_steps_with_group function
replicate_steps_with_group <- function(species_id, df) {
  species_data <- subset(df, species == species_id)
  return(with(species_data, data.frame(
    species = species_id,
    group = rep(group, edge_length*10),
    combined_group = rep(combined_group, edge_length*10),
    step = rep(step, edge_length*10)
  )))
}

# Create the new DataFrame using lapply and do.call
new_df <- do.call(rbind, lapply(unique(new$species), function(species_id) replicate_steps_with_group(species_id, new)))


# Calculate sequence and average step
new_df <- new_df %>%
  group_by(species) %>%
  mutate(sequence = row_number()) %>%
  ungroup()

avg_step_df <- new_df %>%
  group_by(combined_group, sequence) %>%
  summarize(avg_step = mean(step), .groups = 'drop')

split=(round(100-branching.times(tree)["88"][1],1))*10

average_step_per_sequence <- avg_step_df %>%
  group_by(sequence) %>%
  summarize(aveg_step = mean(avg_step)) %>%
  ungroup()

# Joining the data frames on the sequence column
updated_avg_step_df <- avg_step_df %>%
  left_join(average_step_per_sequence, by = "sequence")

# Conditional replacement of avg_step values for sequence 1 to 31
updated_avg_step_df <- updated_avg_step_df %>%
  mutate(avg_step = ifelse(sequence >= 1 & sequence <= split, aveg_step, avg_step)) %>%
  select(-aveg_step)  # Optionally remove the extra column after replacement


steps23 = updated_avg_step_df


# Plotting
ggplot(updated_avg_step_df, aes(x = sequence, y = avg_step, group = combined_group, color = combined_group)) +
  geom_line() +
  theme_minimal() +
  labs( x = "Sequence",
       y = "Average Step Size",
       color = "Combined Group")


steps_all = rbind(steps12,steps13,steps14,steps23,steps24,steps34)


avg_step_all <- steps_all %>%
  group_by(combined_group, sequence) %>%
  summarize(avg_step = mean(avg_step), .groups = 'drop')


# Plotting
p=ggplot(avg_step_all, aes(x = sequence, y = avg_step, group = combined_group, color = combined_group)) +
  geom_line() +
  theme_minimal() +
  labs( x = "Sequence",
        y = "Average Step Size",
        color = "Combined Group")+
  theme(axis.text = element_text(size = 26))+
  scale_x_continuous( limits = c(0, 990))


ggsave(p,height = 10,width = 22,filename = "divergence.svg")




























tail(steps12)










library(picante)


x = match.phylo.data(phy = tree,data = data)$data
groupvec=x$opt2

name = rownames(match.phylo.data(phy = tree,data = data)$data)



samp <- matrix(0, nrow = 2, ncol = length(name), dimnames = list(c("A", "B"), name))

# Populate the community matrix
for (i in 1:length(name)) {
  group <- groupvec[i]
  samp[group, name[i]] <- 1
}


dist.tree <- cophenetic.phylo(tree)

###B=superorganism
result = ses.mpd(samp = samp, dis = dist.tree)


n_A = choose(n = result$ntaxa[1],k = 2)
n_B = choose(n = result$ntaxa[2],k = 2)

param_A = result$mpd.obs[1]/n_A
param_B = result$mpd.obs[2]/n_B



#superorganisms
0.68/param_B
#simple
0.87/param_A






























######get angles of trajectories
# read data
data = read.csv("data.csv",row.names = "name")
dat = read.csv("pca_scores.csv",row.names = "X")
dat = dat[,1:4]
tree = read.newick("tree.nwk")

# match data with phylo
demo = match.phylo.data(phy = tree,data = dat)
dat = demo$data
phy = demo$phy
name.check(phy,dat)

# define function for column pairs
angle_calculator <- function(dat, index1, index2, phy, species = "M.rotundata"){
  
  # create vectors for axes
  vec_creator <- function(dat, index, phy, species = "M.rotundata"){
    demo = dat[,index]
    names(demo) = rownames(dat)
    x = as.vector(fastAnc(phy ,demo ,anc.states=setNames(c(dat[species,index]), c("78"))))
    return(x)
  }
  
  x = vec_creator(dat, index1, phy)
  y = vec_creator(dat, index2, phy)
  
  # combine vectors to df
  df <- data.frame(x=c(dat[,index1],x), y=c(dat[,index2],y))
  
  # extract the nodes each species go through and create path for each species for x and y
  path = nodepath(phy)
  a = lapply(path, function(i) df[i,1])
  b = lapply(path, function(i) df[i,2])
  
  # combine results
  a = data.frame(v1 = unlist(a), v2 = rep(seq(length(a)), lengths(a)))
  b = data.frame(v1 = unlist(b), v2 = rep(seq(length(b)), lengths(b)))
  demo = cbind(a, b[,1])
  colnames(demo) = c("x","species","y")
  
  # make values numeric
  demo$x = as.numeric(demo$x)
  demo$y = as.numeric(demo$y)
  
  # calculate angle
  groups = as.data.frame(data$angles)
  rownames(groups) = rownames(data)
  groups = match.phylo.data(phy = tree,data = groups)$data
  groups = as.vector(groups$`data[res$phy$tip.label, ]`)
  group_vec = groups[demo$species]
  
  # take absolute value
  demo = abs(demo)
  
  # get vectors
  df = mapply(c, demo$x, demo$y, SIMPLIFY = F)
  
  # Define a function to calculate angles between pairs of vectors
  calc_angle <- function(i, vec_list) {
    vec1 = vec_list[[i-1]]
    vec2 = vec_list[[i]]
    return(getAngle(vec1, vec2))
  }
  
  # Calculate angles for each pair of vectors
  angles <- sapply(2:length(df), calc_angle, vec_list = df)
  
  # Create data frame
  angle <- data.frame(cos = angles,
                      species = demo$species[-1],
                      state = group_vec[-1])
  
  # take complementary angle
  angle$cos = pi - angle$cos
  
  # remove invalid rows 
  remove = duplicated(angle$species);remove[1]=TRUE
  angle = angle[remove,]
  
  # add columns indexes
  angle$column_1 = index1
  angle$column_2 = index2
  
  return(angle)
}

# define column indexes
column_indexes = combn(1:ncol(dat), 2, simplify = FALSE)

# apply function for each pair of columns
result = bind_rows(lapply(column_indexes, function(index) {
  angle_calculator(dat, index[1], index[2], phy)
}))


# Replace "B" with "A" and "D" with "C" in 'state' column
result <- result %>%
  mutate(state = ifelse(state == "B", "A", 
                        ifelse(state == "D", "C", state)))

ticks=c(pi*0.5,pi*0.75,pi,pi*1.25)

angles = ggplot(data=result, aes(x=cos, group=state, fill=state)) +
  geom_density(adjust=3.5, alpha=.2) +
  theme_bw()+
  scale_x_continuous(breaks=ticks,limits = c(pi*0.5,pi*1.25))

ggsave("angles.svg", angles, width = 18,height = 10,dpi = 300)


####T.test COMPARE MEANS
t.test(cos ~ state, result)

#####F test to compare variances

var.test(data = result,cos~state)

var(result[result$state== 'A',"cos"])
var(result[result$state== 'C',"cos"])
