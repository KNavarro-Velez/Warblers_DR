#                             Warblers DR Foraging Behavior
#     K. Navarro-Velez, A. Eppedio, E. Heiser, M. Gilbert, M. Reinoso-Perez, A. Dhondt 
#                                    April 18th 2023
#______________________________________________________________________________


#Set wd
setwd("~/Downloads")

library(tidyverse) 
library(MASS)
library(vegan)
library(janitor) # install.packages("janitor")
library(pairwiseAdonis) #install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# NOTE: 
# TO INSTALL pairwiseAdonis you need to follow these steps:

#1. install.packages("devtools")
#2. library(devtools)
#3. install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#4. library(pairwiseAdonis)


#select variables data
data <- read_csv("Target_species_and_variables.csv")

#raw data
Dom_Rep_data <- read.csv("Dom Republic.csv")

# Data cleaning ----
data <- data %>%
  rename("index" = ...1, "for_be" = "Foraging Behavior")



# 1. Chi square----

# This is statistical test is very simple but it compares the frequency of the behaviors between species  from a random sample and determines if what we found is the same we would found with random data or it is significantly different s.


table1 <- table(data$Species, data$for_be)

chisq.test(table1)
fisher.test(table1, simulate.p.value = TRUE) # since expected counts are low I did a Fisher test as well and it is still significant



# 2. PERMANOVA test ----
# This test is an analysis of variance. It makes a geometric partitioning of the variation across the species using a distance matrix (similarity or dissimilarity matrix) in response to the behaviors

tmp_data <-data %>% 
  dplyr::select(Species, for_be, `ID #`) %>%
  mutate(ct = 1)

tmp_data <- tmp_data %>% 
  group_by(`ID #`, for_be, Species) %>%
  summarize(ct = sum(ct)) %>% 
  ungroup()

tmp_data <- tmp_data %>%
  pivot_wider(names_from = for_be, values_from = ct) %>% 
  rename("Fl"=`F`)


# standardize counts
tmp_data2 <- tmp_data %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(rowsums= rowSums(tmp_data[,c(3,4,5,6,7,8,9,10)], na.rm = TRUE)) %>% 
  mutate(P= P/rowsums, G= G/rowsums, H= H/rowsums) %>% 
  mutate(S= S/rowsums, C= C/rowsums, K= K/rowsums) %>% 
  mutate(D= D/rowsums, Fl= Fl/rowsums)

adonis2(tmp_data2[, 3:10] ~ Species,
        data = tmp_data, )

pairwise_adonis <- pairwise.adonis(tmp_data2[, 3:10], factors = as.factor(tmp_data$Species))




# 3. NMDS  with all data values----

data_behavior <- tmp_data2 %>%
  dplyr::select(Species:D) %>% 
  mutate(Species= as.factor(Species)) ####THIS IS MY LAST CHANGE I WAS TRYING TO USE ALL THE DATA FOR IT. I WANT TO TRANSPONSE SPECIES FOR BEHAVIORS AND RUN IT AGAIN

unique(data_behavior$Species)

nmds_be <- metaMDS(data_behavior[,-1], trace = FALSE, distance = "bray", maxit=99, trymax = 300, wascores = TRUE)# for final results increase maxit to 9999
plot(nmds_be, type = "t")
nmds_be$data # this is telling us the type of transformation made by the NMDS function

#As a rule of thumb literature has identified the following cut-off values for stress-level:

#Higher than 0.2 is poor (risks for false interpretation).
#0.1 - 0.2 is fair (some distances can be misleading for interpretation).
#0.05 - 0.1 is good (can be confident in inferences from plot).
#Less than 0.05 is excellent (this can be rare).

# check for good fit.
goodness(nmds_be) # this show us how well the ordination plots represent real data for each species

stressplot(nmds_be) # This shows how closely our ordination fits real world plot dissimilarities and how well we can interpret the ordination

colvec <- c("aliceblue","chocolate1", "darksalmon", "gray4", "lavenderblush2", "plum2", "royalblue3", "turquoise3")   # Identifies colors for group assignments
pchvec <- c(20, 20,20,20,20,20,20,20)   # Identifies character symbols for group assignments

nmds_be_groups <- unique(data_behavior[,1])

plot(nmds_be)
with(nmds_be_groups,
     points(nmds_be,
            display = "species", #speceis are the behaviors
            col = "red",
            pch = pchvec,
            bg = colvec))

#Create convex hulls that highlight point clusters based on grouping dataframe
ordihull(
  nmds_be,nmds_be_groups$Species,
  display = "sites", #species is behaviors
  draw = c("polygon"),
  col = NULL,
  border = c("aliceblue","chocolate1", "darksalmon", "gray4", "lavenderblush2", "plum2", "royalblue3", "turquoise3"),
  lty = c(1, 2, 1, 2,1,2,1,2),
  lwd = 2.5
)

# Calculating centroids 

# You can calculate centroids for your groups which can be viewed as the average position of observations in ordination space. 
# Calculating and plotting centroids of NMDS Result

scrs <-
  scores(nmds_be, display = "sites", "species")
cent <-
  aggregate(scrs ~ Species, data = nmds_be_groups, FUN = "mean")
names(cent) [-1] <- colnames(scrs)
points(cent [,-1],
       pch = c( 8 , 8 , 8, 8),
       col = c("aliceblue","chocolate1", "darksalmon", "gray4", "lavenderblush2", "plum2", "royalblue3", "turquoise3", "yellow4", "darkorange4"),
       bg = c("black"),
       lwd = 3.0,
       cex = 2.0 # Plots centroids as points on ordination
)
















#________________________________________________________________________________________#
#___________________________________________________________________________

# 4. NMDS with frequencies ----

data_behavior <- data %>%
  dplyr::select(Species, for_be, index)

unique(data_behavior$for_be)

Total_records_per_sp <- data_behavior %>%
  group_by(Species) %>%
  summarise(total_per_sp = n())

Total_records_per_sp_be <- data_behavior %>%  
  group_by(Species, for_be) %>% 
  summarise(behav_records= n()) %>% 
  pivot_wider(names_from = for_be, values_from = behav_records) %>% 
  replace(is.na(.), 0) %>% 
  mutate(Species= as.factor(Species))

nmds1 <- metaMDS(Total_records_per_sp_be[,-1], trace = FALSE, distance = "bray", maxit=99, trymax = 500, wascores = TRUE)
plot(nmds1, type = "t")
nmds1$data # this is telling us the type of transformation made by the NMDS function

#As a rule of thumb literature has identified the following cut-off values for stress-level:

#Higher than 0.2 is poor (risks for false interpretation).
#0.1 - 0.2 is fair (some distances can be misleading for interpretation).
#0.05 - 0.1 is good (can be confident in inferences from plot).
#Less than 0.05 is excellent (this can be rare).

# check for good fit.
goodness(nmds1) # this show us how well the ordination plots represent real data for each species

stressplot(nmds1) # This shows how closely our ordination fits real world plot dissimilarities and how well we can interpret the ordination



# next step: make it pretty.

plot(nmds1, "sites")   # this plots the species
orditorp(nmds1, "sites") # this pulls the names for species


colvec <- c("aliceblue","chocolate1", "darksalmon", "gray4", "lavenderblush2", "plum2", "royalblue3", "turquoise3", "yellow4", "darkorange4")   # Identifies colors for group assignments
pchvec <- c(20, 20,20,20,20,20,20,20,20,20)   # Identifies character symbols for group assignments

nmds1_groups <- Total_records_per_sp_be[,1]

plot(nmds1)
with(nmds1_groups,
     points(nmds1,
            display = "sites",
            col = "black",
            pch = pchvec[Species],
            bg = colvec[Species]))

#Create convex hulls that highlight point clusters based on grouping dataframe
ordihull(
  nmds1,
  nmds1_groups$Species,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("aliceblue","chocolate1", "darksalmon", "gray4", "lavenderblush2", "plum2", "royalblue3", "turquoise3", "yellow4", "darkorange4"),
  lty = c(1, 2, 1, 2),
  lwd = 2.5
)

# Calculating centroids 

# You can calculate centroids for your groups which can be viewed as the average position of observations in ordination space. 
# Calculating and plotting centroids of NMDS Result

scrs <-
  scores(nmds1, display = "sites", "species")
cent <-
  aggregate(scrs ~ Species, data = nmds1_groups, FUN = "mean")
names(cent) [-1] <- colnames(scrs)
points(cent [,-1],
       pch = c( 8 , 8 , 8, 8),
       col = c("aliceblue","chocolate1", "darksalmon", "gray4", "lavenderblush2", "plum2", "royalblue3", "turquoise3", "yellow4", "darkorange4"),
       bg = c("black"),
       lwd = 3.0,
       cex = 2.0 # Plots centroids as points on ordination
)





