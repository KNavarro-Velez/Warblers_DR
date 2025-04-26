#                         Warblers DR Foraging Behavior
# M. Gilbert,  A. Eppedio, E. Heiser, K. Navarro-Velez & A.A. Dhondt 
#                                July 14th 2023
#______________________________________________________________________________

#                                   SETUP

# Don't forget to set working directory
new.data = read.csv("Dom Republic 2023 Data Sheets - July (1).csv", header=T)

# Load packages
library(vegan)
library(tidyverse)
library(ggplot2)
library(pairwiseAdonis)
library(janitor)
library(cluster)
library(ggalt)
library(ggrepel)

# set date
todays_date <- as.character(Sys.Date())

# Load the dataset
data <- read.csv("Dom Republic 2023 Data Sheets - July.csv", header=T)


#______________________________________________________________________________

#                     ORANIZE FORAGING BEHAVIOR DATA


#Make data set to count number of different types of foraging behavior for each individual
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

#Edit previous data set to be the percentage of foraging tempts of the behavior type for that individual
data3_max <- data2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(data2[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums) %>% 
  group_by(id, species) %>% 
  summarise(max_p= max(P), max_g= max(G), max_fl= max(Fl), max_h= max(H), max_s= max(S),max_c= max(C), max_k= max(K), max_d= max(D)) # instead of max we can also use the sum or the mean count.

#Calculate maximum percentage of foraging attempts were this behavior for an 
#individual for species
data4_max <- data2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(data2[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums) %>% 
  group_by(species) %>% 
  summarise(max_p= max(P), max_g= max(G), max_fl= max(Fl), max_h= max(H), max_s= max(S),max_c= max(C), max_k= max(K), max_d= max(D)) 

#calculate what percentage of foraging attempts were this behavior for each species
data4_mean <- data2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(data2[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums) %>% 
  group_by(species) %>% 
  summarise(mean_p= mean(P), mean_g= mean(G), mean_fl= mean(Fl), mean_h= mean(H), mean_s= mean(S),mean_c= mean(C), mean_k= mean(K), mean_d= mean(D)) 


# _____________________________________________________________________________

#                     ORGANIZE FORAGING LOCATION DATA


#Make data set to count number of different types of foraging location for each individual
fata2 <- data %>% 
  dplyr::select(ID.., Species, Mixed.Location) %>% 
  clean_names() %>% 
  arrange(id) %>% 
  mutate(count= 1) %>% 
  group_by(id, mixed_location, species) %>% 
  summarize(ct= sum(count)) %>%
  ungroup() %>% 
  pivot_wider(names_from = mixed_location, values_from = ct) %>% 
  rename("FW"=`F`) %>%
  rename("DM"=`D`) %>%
  rename("GR"=`G`)
  
#Edit previous Data set to be the percentage of foraging attempts of the behavior type for that individual
fata3_max <- fata2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(fata2[,c(3,4,5,6,7,8,9,10,11,12,13,14)], na.rm = TRUE)) %>% 
  mutate(B= B/rowsums, L= L/rowsums, FW= FW/rowsums) %>% 
  mutate(OL= OL/rowsums, OF= OF/rowsums, UF= UF/rowsums) %>% 
  mutate(UL= UL/rowsums, DM= DM/rowsums) %>% 
  mutate(GR= GR/rowsums, T=T/rowsums, SW= SW/rowsums) %>% 
  mutate(W= W/rowsums) %>% 
  group_by(id, species) %>% 
  summarise(max_b= max(B), max_l= max(L), max_fw= max(FW), max_ol= max(OL), max_of= max(OF),max_uf= max(UF), max_ul= max(UL), max_dm= max(DM), max_gr= max(GR),max_t= max(T), max_sw= max(SW), max_w= max(W))

#Calculate maximum percentage of foraging attempts were this behavior for an 
#individual for species
fata4_max <- fata2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(fata2[,c(3,4,5,6,7,8,9,10,11,12,13,14)], na.rm = TRUE)) %>% 
  mutate(B= B/rowsums, L= L/rowsums, FW= FW/rowsums) %>% 
  mutate(OL= OL/rowsums, OF= OF/rowsums, UF= UF/rowsums) %>% 
  mutate(UL= UL/rowsums, DM= DM/rowsums) %>% 
  mutate(GR= GR/rowsums, T=T/rowsums, SW= SW/rowsums) %>% 
  mutate(W= W/rowsums) %>% 
  group_by(species) %>% 
  summarise(max_b= max(B), max_l= max(L), max_fw= max(FW), max_ol= max(OL), max_of= max(OF),max_uf= max(UF), max_ul= max(UL), max_dm= max(DM), max_gr= max(GR),max_t= max(T), max_sw= max(SW), max_w= max(W)) 

#calculate what percentage of foraging attempts were this behavior for each species
fata4_mean <- fata2 %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(fata2[,c(3,4,5,6,7,8,9,10,11,12,13,14)], na.rm = TRUE)) %>% 
  mutate(B= B/rowsums, L= L/rowsums, FW= FW/rowsums) %>% 
  mutate(OL= OL/rowsums, OF= OF/rowsums, UF= UF/rowsums) %>% 
  mutate(UL= UL/rowsums, DM= DM/rowsums) %>% 
  mutate(GR= GR/rowsums, T=T/rowsums, SW= SW/rowsums) %>% 
  mutate(W= W/rowsums) %>% 
  group_by(species) %>% 
  summarise(mean_b= mean(B), mean_l= mean(L), mean_fw= mean(FW), mean_ol= mean(OL), mean_of= mean(OF),mean_uf= mean(UF), mean_ul= mean(UL), mean_dm= mean(DM), mean_gr= mean(GR),mean_t= mean(T), mean_sw= mean(SW), mean_w= mean(W)) 


#______________________________________________________________________________

#             COMBINE FORAGING AND LOCATION DATA INTO SINGLE DATA FRAME


#Left_join, using the species names as keys
combined.data4_max = left_join(data4_max, fata4_max, by = "species")
combined.data4_mean = left_join(data4_mean, fata4_mean, by = "species")
data2.1 = left_join(data2, fata2, by = "id", relationship="many-to-many")
  combined.data2 <- data2.1[-c(11)]
  colnames(combined.data2)[2] = "species"
combined.data3_max.1 = left_join(data3_max, fata3_max, by = "id", relationship="many-to-many")
  combined.data3_max <- combined.data3_max.1[-c(11)]
  colnames(combined.data3_max)[2] = "species"


#______________________________________________________________________________

#                             CALCULATE DISSIMILARITY (foraging behavior)

#create dissimilarity matrix for maximums
diss_matrix_max <- vegdist(data4_max[,-1], method = "bray")
print(diss_matrix_max)
max(diss_matrix_max)
min(diss_matrix_max)

#create dissimilarity matrix for means
diss_matrix_mean <- vegdist(data4_mean[,-1], method = "bray")
print(diss_matrix_mean)
max(diss_matrix_mean)
min(diss_matrix_mean)


#______________________________________________________________________________

#                               RUN PERMANOVA (foraging behavior)
 

species_10 <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA")
species_permanova <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA", "BBTO", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "OVEN", "PAWA", "PRAW", "YRWA", "PAWA", "PRAW", "YRWA", "PRAW", "YRWA")
permanova_results <- adonis2(combined.data3_max[, 3:10] ~ species,
                             data = combined.data2, )
print(permanova_results)


#______________________________________________________________________________

#                                 RUN NMDS (FORAGING BEHAVIOR)

nmds <- metaMDS(diss_matrix_mean, k = 2, try = 50, trymax = 50)  # Set k = n for a n-dimensional NMDS
nmds_coords <- data.frame(nmds$points)
nmds_coords$Species <- factor(species_10)

k <- 4  # Set the number of clusters
kmeans_result <- kmeans(diss_matrix_mean, centers = k)

# Add the cluster assignments to the NMDS coordinates data frame
nmds_coords$Cluster <- factor(kmeans_result$cluster)

# Create the NMDS plot with circles based on clusters
ggplot(nmds_coords, aes(x = MDS1, y = MDS2, color = Cluster)) +
  geom_point(size = 4) +
  labs(x = "NMDS1", y = "NMDS2", title = "Dissimilarity of Foraging Attempt Type") +
  theme_minimal() +
  geom_text_repel(aes(label = Species), nudge_x = 0.07)
  #geom_segment(aes(x=0.25, y=0, xend=0.9, yend=0), arrow = arrow(length=unit(0.5, 'cm')), col="brown1") +
  #geom_segment(aes(x=0.25, y=0, xend=0, yend=0), arrow = arrow(length=unit(0.5, 'cm')), col="black")+
  #geom_segment(aes(x=0.25, y=0, xend=0.25, yend=-4e-04), arrow = arrow(length=unit(0.5, 'cm')), col="chartreuse3") +
  #geom_segment(aes(x=0.25, y=0, xend=0.25, yend=6e-04), arrow = arrow(length=unit(0.5, 'cm')), col="cornflowerblue")
#geom_encircle(aes(fill = Cluster), color = "black", expand = -0.05, alpha = 0.2)



#______________________________________________________________________________

#                             CALCULATE DISSIMILARITY (foraging location)

#create dissimilarity matrix for maximums
diss_matrix_maxl <-vegdist(fata4_max[,-1], method = "bray")
print(diss_matrix_maxl)
max(diss_matrix_maxl)
min(diss_matrix_maxl)

#create dissimilarity matrix for means
diss_matrix_meanl <- vegdist(fata4_mean[,-1], method = "bray")
print(diss_matrix_meanl)
max(diss_matrix_meanl)
min(diss_matrix_meanl)


#______________________________________________________________________________

#                               RUN PERMANOVA (foraging Location)


species_10l <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA")
species_permanoval <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA", "BBTO", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "OVEN", "PAWA", "PRAW", "YRWA", "PAWA", "PRAW", "YRWA", "PRAW", "YRWA")
permanova_resultsl <- adonis2(fata3_max[, 3:14] ~ species,
                             data = fata2, )
print(permanova_resultsl)


#______________________________________________________________________________

#                                 RUN NMDS (FORAGING Location)

nmdsl <- metaMDS(diss_matrix_meanl, k = 2, try = 50, trymax = 50)  # Set k = n for a n-dimensional NMDS
nmds_coordsl <- data.frame(nmdsl$points)
nmds_coordsl$Species <- factor(species_10l)

kl <- 4# Set the number of clusters
kmeans_resultl <- kmeans(diss_matrix_meanl, centers = kl)

# Add the cluster assignments to the NMDS coordinates data frame
nmds_coordsl$Cluster <- factor(kmeans_resultl$cluster)

# Create the NMDS plot with circles based on clusters
ggplot(nmds_coordsl, aes(x = MDS1, y = MDS2, color = Cluster)) +
  geom_point(size = 4) +
  labs(x = "NMDS1", y = "NMDS2", title = "Dissimilarity of Foraging Substrate Type") +
  theme_minimal() +
  geom_text_repel(aes(label = Species), nudge_x = 0.14)
  #geom_segment(aes(x=0.25, y=0, xend=0.9, yend=0), arrow = arrow(length=unit(0.5, 'cm')), col="brown1") +
  #geom_segment(aes(x=0.25, y=0, xend=0, yend=0), arrow = arrow(length=unit(0.5, 'cm')), col="black")+
  #geom_segment(aes(x=0.25, y=0, xend=0.25, yend=-4e-04), arrow = arrow(length=unit(0.5, 'cm')), col="chartreuse3") +
  #geom_segment(aes(x=0.25, y=0, xend=0.25, yend=6e-04), arrow = arrow(length=unit(0.5, 'cm')), col="cornflowerblue")
#geom_encircle(aes(fill = Cluster), color = "black", expand = -0.05, alpha = 0.2)



#______________________________________________________________________________

#                             CALCULATE DISSIMILARITY (combined)

#create dissimilarity matrix for maximums
diss_matrix_maxc <- vegdist(combined.data4_max[,-1], method = "bray")
print(diss_matrix_maxc)
max(diss_matrix_maxc)
min(diss_matrix_maxc)

#create dissimilarity matrix for means
diss_matrix_meanc <- vegdist(combined.data4_mean[,-1], method = "bray")
print(diss_matrix_meanc)
max(diss_matrix_meanc)
min(diss_matrix_meanc)


#______________________________________________________________________________

#                               RUN PERMANOVA (combined)


species_10c <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA")
species_permanovac <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA", "BBTO", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "OVEN", "PAWA", "PRAW", "YRWA", "PAWA", "PRAW", "YRWA", "PRAW", "YRWA")
permanova_resultsc <- adonis2(combined.data3_max[, 3:22] ~ species,
                              data = combined.data2, )
print(permanova_resultsc)


#______________________________________________________________________________

#                                 RUN NMDS (combined)

nmdsc <- metaMDS(diss_matrix_meanc, k = 2, try = 50, trymax = 50)  # Set k = n for a n-dimensional NMDS
nmdsc
nmds_coordsc <- data.frame(nmdsc$points)
nmds_coordsc$Species <- factor(species_10c)

kc <- 6# Set the number of clusters

kmeans_resultc <- kmeans(diss_matrix_meanc, centers = kc)

# Add the cluster assignments to the NMDS coordinates data frame
nmds_coordsc$Cluster <- factor(kmeans_resultc$cluster)
nmds_coordsc$Cluster = c("Leaf/Twig Gleaner", "Trunk Gleaner", "Flycatcher", "Berry Eater", "Leaf/Twig Gleaner", "Leaf/Twig Gleaner", "Ground Forager", "Insectivorous Generalist", "Leaf/Twig Gleaner", "Insectivorous Generalist")


# Create the NMDS plot with circles based on clusters
ggplot(nmds_coordsc, aes(x = MDS1, y = MDS2, color = Cluster)) +
  geom_point(size = 7) +
  labs(x = "NMDS1", y = "NMDS2", title = "Similarity in Foraging Attack and Substrate Types") +
  theme_minimal() +
  theme(text=element_text(size=17, color="black")) +
  theme(axis.text.y =element_text(size=15, color="black")) +
  theme(axis.text.x =element_text(size=13, color="black")) +
  geom_text_repel(aes(label = Species), size=7, nudge_x = 0.16) 

 # geom_text_repel(aes(label = Species), nudge_x = 0.14)
#geom_segment(aes(x=0.25, y=0, xend=0.9, yend=0), arrow = arrow(length=unit(0.5, 'cm')), col="brown1") +
#geom_segment(aes(x=0.25, y=0, xend=0, yend=0), arrow = arrow(length=unit(0.5, 'cm')), col="black")+
#geom_segment(aes(x=0.25, y=0, xend=0.25, yend=-4e-04), arrow = arrow(length=unit(0.5, 'cm')), col="chartreuse3") +
#geom_segment(aes(x=0.25, y=0, xend=0.25, yend=6e-04), arrow = arrow(length=unit(0.5, 'cm')), col="cornflowerblue")
#geom_encircle(aes(fill = Cluster), color = "black", expand = -0.05, alpha = 0.2)


#______________________________________________________________________________

#                    NMDS OF RAW HEIGHT and CANOPY HEIGHT


#read second data files with observations with no canopy height cut out
new.data <- read.csv("Dom Republic 2023 Data Sheets - July (1).csv", header=T)

height.stuff <- new.data %>% group_by(Species) %>%
  summarize(mean_height=0.0001*mean(Height.in.Canopy), .groups = "drop") %>%
  rename("species"=`Species`)

#Left_join, using the species names as keys
#full.data4_mean = left_join(combined.data4_mean, height.stuff, by = "species")

#______ add canopy height to it

canopy.height.data <- new.data %>% group_by(Species) %>%
  summarize(canopy_height=0.1*mean(Canopy.height), .groups = "drop") %>%
  rename("species"=`Species`)

#Left_join, using the species names as keys
height.stuff = left_join(height.stuff, canopy.height.data, by = "species")


#______ add raw bird height to it

bird.height.data <- new.data %>% group_by(Species) %>%
  summarize(bird_height=0.15*mean(Bird.Height), .groups = "drop") %>%
  rename("species"=`Species`)

#Left_join, using the species names as keys
height.stuff = left_join(height.stuff, bird.height.data, by = "species")


#______ create matrix

#create dissimilarity matrix for means
diss_matrix_meanf <- vegdist(height.stuff[,-1], method = "bray")
print(diss_matrix_meanf)
max(diss_matrix_meanf)
min(diss_matrix_meanf)

#______ calculate nmds and plot

nmdsf <- metaMDS(diss_matrix_meanf, k = 2, try = 50, trymax = 50)  # Set k = n for a n-dimensional NMDS
nmdsf
nmds_coordsf <- data.frame(nmdsf$points)
nmds_coordsf$Species <- factor(species_10c)

kf <- 5# Set the number of clusters
kmeans_resultf <- kmeans(diss_matrix_meanf, centers = kf)

# Add the cluster assignments to the NMDS coordinates data frame
nmds_coordsf$Cluster <- factor(kmeans_resultf$cluster)

# Create the NMDS plot with circles based on clusters
ggplot(nmds_coordsf, aes(x = MDS1, y = MDS2, color = Cluster)) +
  geom_point(size = 7) +
  labs(x = "NMDS1", y = "NMDS2", title = "Similarity in Foraging Height and Canopy Height") +
  theme_minimal() +
  theme(text=element_text(size=17, color="black")) +
  theme(axis.text.y =element_text(size=15, color="black")) +
  theme(axis.text.x =element_text(size=13, color="black")) +
  #theme(text =element_text(size=13)) +
  theme(axis.title.y = element_text(margin = margin(t = -50))) +
  geom_text_repel(aes(label = Species), size=7, nudge_x = 0.07)


#_______________________________________________________________________________

#                  PERMANOVA FOR FORAGING HEIGHT and CANOPY HEIGHT


species_10f <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA")
species_permanovaf <- c("AMRE", "BAWW", "BBTO", "BWVI", "CMWA", "NOPA",  "OVEN", "PAWA", "PRAW", "YRWA", "BBTO", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "BWVI", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "CMWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "NOPA", "OVEN", "PAWA", "PRAW", "YRWA", "OVEN", "PAWA", "PRAW", "YRWA", "PAWA", "PRAW", "YRWA", "PRAW", "YRWA")
hp.data <- new.data[ , c(3, 4, 10)]    
permanova_resultsf <- adonis2(hp.data[, 2:3] ~ Species,
                              data = hp.data) #I don't actually know if this is statistically valid
print(permanova_resultsf)


#______________________________________________________________________________

#             GRAPH OF FORAGING HEIGHT AND CANOPY HEIGHT WITH OVALS


#make species subset of the entire dataset
AMRE = subset(new.data, new.data$Species == "AMRE")
PRAW = subset(new.data, new.data$Species == "PRAW")
NOPA = subset(new.data, new.data$Species == "NOPA")
CMWA = subset(new.data, new.data$Species == "CMWA")
OVEN = subset(new.data, new.data$Species == "OVEN")
PAWA = subset(new.data, new.data$Species == "PAWA")
BWVI = subset(new.data, new.data$Species == "BWVI")
YRWA = subset(new.data, new.data$Species == "YRWA")
BBTO = subset(new.data, new.data$Species == "BBTO")
BAWW = subset(new.data, new.data$Species == "BAWW")


overlap <- rbind(AMRE, PRAW, NOPA, CMWA) #use only the data from these four species for the graph

ggplot(overlap, aes(y = Bird.Height, x = Canopy.height, color = Species)) +
         theme_minimal() +
         labs(x = "Canopy Height (m)", y = "Bird Height (m)", title = "") + 
         #geom_point(size = 1) + 
         xlim(0, 18) + ylim(0, 14) +
         theme(text = element_text(size = 20))  +
         geom_jitter(width = 0.5, height=0.5) +
         stat_ellipse (geom = "polygon", level=0.75, size=0.7, aes(fill=Species), alpha=0.2)
         #geom_encircle(mapping = NULL, data = NULL, stat = "identity",
                #position = "identity", na.rm = F, show.legend = NA,
                #inherit.aes = TRUE, spread=-5, alpha = 1, size=2, expand = 0.05) 
         #geom_encircle(aes(fill = Species), s_shape=0, color = "black", expand = 0.05, alpha = 0.2)
       

#______________________________________________________________________________

#                     Box Plot of HEIGHT WITHIN CANOPY


sorted_data <- data
sorted_data$Species <- factor(sorted_data$Species,     # Reorder factor levels
                         c("OVEN", "PAWA", "YRWA", "AMRE", "PRAW", "BBTO", "BAWW", "NOPA", "CMWA", "BWVI"))

cols <- c("OVEN" = "#00BF7D", "PAWA" = "#00BF7D", "YRWA" = "#E76BF3", "AMRE" = "#A3A500", "PRAW" = "#A3A500", "NOPA"="#00B0F6", "BAWW"="#00B0F6", "BBTO"="#00B0F6", "CMWA"="#F8766D", "BWVI"="#F8766D")

ggplot(sorted_data, aes(x=Species, y=X..Height.in.Canopy, fill=Species, alpha = 0.5)) +
  geom_boxplot(size=0.8) +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  theme(text=element_text(size=20, color="black")) +
  theme(axis.text.y =element_text(size=15, color="black")) +
  theme(axis.text.x =element_text(size=13, color="black")) +
  theme(axis.title.y = element_text(margin = margin(t = -50))) +
  labs(x = "\nSpecies", y = "Proportion of Foraging Height to Canopy Height\n") +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2.5, alpha=1)


#______________________________________________________________________________

#            Box Plot of HEIGHT WITHIN CANOPY

#actually, no, see other script (1_DR_Stacked_Bar_Plots.R)


AMRE = subset(data, data$Species == "AMRE")
PRAW = subset(data, data$Species == "PRAW")
NOPA = subset(data, data$Species == "NOPA")
CMWA = subset(data, data$Species == "CMWA")
OVEN = subset(data, data$Species == "OVEN")
PAWA = subset(data, data$Species == "PAWA")
BWVI = subset(data, data$Species == "BWVI")
YRWA = subset(data, data$Species == "YRWA")
BBTO = subset(data, data$Species == "BBTO")
BAWW = subset(data, data$Species == "BAWW")
BAWW = subset(data, data$Species == "BAWW")

table(YRWA$Foraging.Behavior)
table(PRAW$Foraging.Behavior)
table(NOPA$Foraging.Behavior)
table(AMRE$Foraging.Behavior)
table(CMWA$Foraging.Behavior)
table(OVEN$Foraging.Behavior)
table(PAWA$Foraging.Behavior)
table(BWVI$Foraging.Behavior)
table(BBTO$Foraging.Behavior)
table(BAWW$Foraging.Behavior)

mean(OVEN$Canopy.height)
mean(PAWA$Canopy.height)
mean(PRAW$Canopy.height)
mean(AMRE$Canopy.height)
mean(NOPA$Canopy.height)
mean(CMWA$Canopy.height)


flock = read.csv("1234.csv", header = T)
YRWA = subset(flock, flock$Species=="YRWA")
table(YRWA$Flock.)
YRWAheight = na.omit(YRWA$X..individuals)
YRWAheight
mean(YRWAheight)







