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

index = c(3,4)
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


























# demo12=demo
# demo34=demo

demo = cbind(demo12[,c(1,3)],demo34[,c(1,3)],demo12[,c(2,4)])
colnames(demo) = c("PC1","PC2","PC3","PC4","species","group")

demo <- demo %>%
  arrange(species) %>%
  group_by(species) %>%
  mutate(combined_group = ifelse(group %in% c("api", "bom", "mel"), "super", "simple"))

# Function to calculate the Euclidean distance between two points
calc_distance <- function(row1, row2) {
  sqrt(sum((row1 - row2) ^ 2))
}


# Function to calculate sequential distances within a species group
calc_sequential_distances <- function(species_data) {
  # Initialize an empty data frame to store results
  results <- tibble(FromRow = integer(),
                    ToRow = integer(),
                    Distance = numeric(),
                    Group = character(),
                    CombinedGroup = character())
  
  if(nrow(species_data) > 1) {
    for (i in 1:(nrow(species_data) - 1)) {
      distance <- calc_distance(species_data[i, c("PC1", "PC2", "PC3", "PC4")], species_data[i+1, c("PC1", "PC2", "PC3", "PC4")])
      new_row <- tibble(FromRow = i,
                        ToRow = i + 1,
                        Distance = distance,
                        Group = as.character(species_data$group[i]),
                        CombinedGroup = ifelse("combined_group" %in% names(species_data), as.character(species_data$combined_group[i]), NA_character_))
      results <- bind_rows(results, new_row)
    }
  }
  
  return(results)
}

# Calculate sequential distances for each species, ensuring group values are tracked
results <- demo %>%
  group_by(species) %>%
  do(calc_sequential_distances(.)) %>%
  ungroup()

# Display the results
print(results)





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
new <- cbind(results, all_lengths_vector)



# Convert species to character and rounding edge_length
names(new)[names(new) == "all_lengths_vector"] <- "edge_length"
new$species <- as.character(new$species)
# new$edge_length <- round(new$edge_length,digits = 1)

# Renaming and creating new columns
new <- transform(new, step = Distance / edge_length)


# Optimized replicate_steps_with_group function
replicate_steps_with_group <- function(species_id, df) {
  species_data <- subset(df, species == species_id)
  return(with(species_data, data.frame(
    species = species_id,
    group = rep(Group, edge_length*10),
    combined_group = rep(CombinedGroup, edge_length*10),
    step = rep(step, edge_length*10)
  )))
}

# Create the new DataFrame using lapply and do.call
new_df <- do.call(rbind, lapply(unique(new$species), function(species_id) replicate_steps_with_group(species_id, new)))

#adjust to node number
split=(round(100-branching.times(tree)["88"][1],1))*10

# Calculate sequence and average step
new_df <- new_df %>%
  group_by(species) %>%
  mutate(sequence = row_number()) %>%
  ungroup()


filtered_df <- new_df %>%
  group_by(sequence) %>%
  filter(!duplicated(step)) %>%
  ungroup()


all_species <- filtered_df %>%
  filter(sequence <= split) %>%
  group_by(sequence) %>%
  summarize(average_step = mean(step, na.rm = TRUE))

# Join the average steps with the original data frame and update the step values
result_df <- filtered_df %>%
  left_join(all_species, by = "sequence") %>%
  mutate(step = ifelse(sequence <= split, average_step, step)) %>%
  select(-average_step)  # Remove the average_step column after updating



avg_step_df <- result_df %>%
  group_by(combined_group, sequence) %>%
  summarize(avg_step = mean(step), .groups = 'drop')




# Plotting
p=ggplot(avg_step_df, aes(x = sequence, y = avg_step, group = combined_group, color = combined_group)) +
  geom_line() +
  theme_minimal() +
  labs( x = "Sequence",
        y = "Average Step Size",
        color = "Combined Group")


ggsave(p,height = 10,width = 22,filename = "diversification.svg")




######get angles of trajectories



# Function to calculate the dot product of two vectors
calc_dot_product <- function(v1, v2) {
  sum(v1 * v2)
}

# Function to calculate the magnitude of a vector
calc_magnitude <- function(v) {
  sqrt(sum(v^2))
}

# Function to calculate the angle between two vectors in radians
calc_angle <- function(v1, v2) {
  dot_product <- calc_dot_product(v1, v2)
  magnitude_v1 <- calc_magnitude(v1)
  magnitude_v2 <- calc_magnitude(v2)
  cos_theta <- dot_product / (magnitude_v1 * magnitude_v2)
  angle <- acos(cos_theta) # Angle in radians
  return(angle)
}

# Adjusted function to calculate sequential distances and angles within a species group
calc_sequential_distances_and_angles <- function(species_data) {
  results <- tibble(FromRow = integer(),
                    ToRow = integer(),
                    Distance = numeric(),
                    Angle = numeric(),
                    Group = character(),
                    CombinedGroup = character())
  
  if(nrow(species_data) > 1) {
    for (i in 1:(nrow(species_data) - 1)) {
      v1 <- species_data[i, c("PC1", "PC2", "PC3", "PC4")]
      v2 <- species_data[i+1, c("PC1", "PC2", "PC3", "PC4")]
      distance <- calc_distance(v1, v2)
      angle <- calc_angle(v1, v2) # Calculate angle in radians
      new_row <- tibble(FromRow = i,
                        ToRow = i + 1,
                        Distance = distance,
                        Angle = angle, # Store angle in radians
                        Group = as.character(species_data$group[i]),
                        CombinedGroup = ifelse("combined_group" %in% names(species_data), as.character(species_data$combined_group[i]), NA_character_))
      results <- bind_rows(results, new_row)
    }
  }
  
  return(results)
}

# Calculate sequential distances and angles for each species, ensuring group values are tracked
angles <- demo %>%
  group_by(species) %>%
  do(calc_sequential_distances_and_angles(.)) %>%
  ungroup()

# Display the results
print(angles)


new_angles <- cbind(angles, all_lengths_vector)

split=(100-branching.times(tree)["88"][1])

compute_cumulative_time <- function(df) {
  df %>%
    group_by(species) %>%
    mutate(time = cumsum(all_lengths_vector)) %>%
    ungroup()
}

# Apply the function to your DataFrame
new_angles <- compute_cumulative_time(new_angles)

new_angles <- new_angles %>%
  mutate(newGroup = ifelse(time < 30, 'all', CombinedGroup))


angles_filtered <- new_angles %>%
  group_by(newGroup) %>%
  filter(!duplicated(Angle)) %>%
  ungroup()

angles_filtered$Angle = pi - angles_filtered$Angle




ticks=c((pi-pi*1.5),(pi-pi*1.25),0,pi*0.25,pi,pi*0.5, pi*0.75,pi)

ticks=c((2*pi/6),(4*pi/6),pi,(pi*8)/6)

p=ggplot(data=angles_filtered, aes(x=Angle, group=newGroup, fill=newGroup)) +
  geom_density(adjust=3.5, alpha=.2) +
  theme_bw()+
  scale_x_continuous(breaks=ticks,limits = c((2*pi/6),(pi*8)/6))

ggsave("angles.svg", p, width = 18,height = 10,dpi = 300)


####T.test COMPARE MEANS
t.test(Angle ~ newGroup, angles_filtered)

group_all <- angles_filtered$Angle[angles_filtered$newGroup == 'simple']
group_other <- angles_filtered$Angle[angles_filtered$newGroup == 'super'] 

# Perform a t-test to compare the means
t_test_result <- t.test(group_all, group_other, alternative = "two.sided", var.equal = FALSE)
print(t_test_result)

sd(group_other)


#####F test to compare variances

var.test(group_all, group_other)

var(angles_filtered[angles_filtered$newGroup== 'super',"Angle"])
var(angles_filtered[angles_filtered$newGroup== 'simple',"Angle"])
var(angles_filtered[angles_filtered$newGroup== 'all',"Angle"])





