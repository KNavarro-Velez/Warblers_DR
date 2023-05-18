#                             Warblers DR Foraging Behavior
#     K. Navarro-Velez, A. Eppedio, E. Heiser, M. Gilbert, M. Reinoso-Perez, A. Dhondt 
#                                    May 10th 2023
#______________________________________________________________________________

# Set WD
setwd("~/github_files/Warblers_DR")

# Load packages
library(vegan)
library(tidyverse)
library(ggplot2)
library(pairwiseAdonis)
library(janitor)
library(cluster)
install.packages("cluster")

# set date ------
todays_date <- as.character(Sys.Date())

# Load the dataset ----

data <- read.csv("Target_species_and_variables.csv")

# organize dataset for the analysis
data2 <- data %>% 
  dplyr::select(ID.., Species, Foraging.Behavior) %>% 
  clean_names() %>% 
  arrange(id) %>% 
  mutate(count= 1) %>% 
  group_by(id, foraging_behavior, species) %>% 
  summarize(ct= sum(count)) %>%
  ungroup() %>% 
  pivot_wider(names_from = foraging_behavior, values_from = ct) %>% 
  rename("Fl"=`F`)

data3_max <- data2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(data2[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums) %>% 
  group_by(id, species) %>% 
  summarise(max_p= max(P), max_g= max(G), max_fl= max(Fl), max_h= max(H), max_s= max(S),max_c= max(C), max_k= max(K), mean_d= mean(D)) # instead of max we can also use the sum or the mean count.

data3_mean <- data2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(data2[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums) %>% 
  group_by(id, species) %>% 
  summarise(mean_p= mean(P), mean_g= mean(G), mean_fl= mean(Fl), mean_h= mean(H), mean_s= mean(S),mean_c= mean(C), mean_k= mean(K), mean_d= mean(D)) # instead of mean we can also use the sum or the mean count.


# calculate dissimilarity matrix----
data4_max <- data2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(data2[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums) %>% 
  group_by(species) %>% 
  summarise(max_p= max(P), max_g= max(G), max_fl= max(Fl), max_h= max(H), max_s= max(S),max_c= max(C), max_k= max(K), mean_d= mean(D)) 


#mean!!!!
data4_mean <- data2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(data2[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums) %>% 
  group_by(species) %>% 
  summarise(mean_p= mean(P), mean_g= mean(G), mean_fl= mean(Fl), mean_h= mean(H), mean_s= mean(S),mean_c= mean(C), mean_k= mean(K), mean_d= mean(D)) 



diss_matrix_max <- vegdist(data4_max[,-1], method = "bray")
print(diss_matrix_max)
max(diss_matrix_max)
min(diss_matrix_max)

diss_matrix_mean <- vegdist(data4_mean[,-1], method = "bray")
print(diss_matrix_mean)
max(diss_matrix_mean)
min(diss_matrix_mean)



species_10 <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA")
species_permanova <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA", "BBTO", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "OVEN", "PAWA", "PRAW", "YRWA", "PAWA", "PRAW", "YRWA", "PRAW", "YRWA")


# Perform PERMANOVA----

permanova_results <- adonis2(data3_max[, 3:10] ~ species,
        data = data2, )
# Species are significantly different in their behaviors.

# Print the PERMANOVA results
print(permanova_results)



# NMDS ----
nmds <- metaMDS(diss_matrix_max, k = 3)  # Set k = 3 for a 3-dimensional NMDS


# Extract NMDS coordinates----
nmds_coords <- data.frame(nmds$points)
nmds_coords$Species <- factor(species_10)

# Plot NMDS using ggplot2
ggplot(nmds_coords, aes(x = MDS1, y = MDS2, color = Species)) +
  geom_point(size = 3) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_minimal()+
  geom_encircle(aes(fill = Species), color = "black", expand = 0.1, alpha = 0.2)




#PLOT WITH MAX
# Perform k-means clustering on the dissimilarity matrix
k <- 4  # Set the number of clusters
kmeans_result <- kmeans(diss_matrix_max, centers = k)

# Add the cluster assignments to the NMDS coordinates data frame
nmds_coords$Cluster <- factor(kmeans_result$cluster)

# Create the NMDS plot with circles based on clusters
ggplot(nmds_coords, aes(x = MDS1, y = MDS2, color = Cluster)) +
  geom_point(size = 3) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_minimal() +
  geom_encircle(aes(fill = Cluster), color = "black", expand = 0.1, alpha = 0.2) +
  geom_text(aes(label = Species), nudge_x = 0.1)
  
  


#PLOT WITH MEAN
# Perform k-means clustering on the dissimilarity matrix
k <- 4  # Set the number of clusters
kmeans_result <- kmeans(diss_matrix_mean, centers = k)

# Add the cluster assignments to the NMDS coordinates data frame
nmds_coords$Cluster <- factor(kmeans_result$cluster)

# Create the NMDS plot with circles based on clusters
ggplot(nmds_coords, aes(x = MDS1, y = MDS2, color = Cluster)) +
  geom_point(size = 3) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_minimal() +
  geom_encircle(aes(fill = Cluster), color = "black", expand = 0.1, alpha = 0.2) +
  geom_text(aes(label = Species), nudge_x = 0.1)







#PLOT WITH MEAN
# Perform k-means clustering on the dissimilarity matrix
k <- 4  # Set the number of clusters
kmeans_result <- kmeans(diss_matrix_mean, centers = k)

# Add the cluster assignments to the NMDS coordinates data frame
nmds_coords$Cluster <- factor(kmeans_result$cluster)

# Create the NMDS plot with circles based on clusters
ggplot(nmds_coords, aes(x = MDS1, y = MDS2, color = Cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_encircle(aes(fill = Cluster), color = "black", expand = 0.1, alpha = 0.2) +
  geom_text(aes(label = Species), nudge_x = 0.1)

nmds_coords


#labs(x = "NMDS1", y = "NMDS3")


