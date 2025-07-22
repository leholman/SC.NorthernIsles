## scrappy analysis for fun! 

library(tidyverse)
library(reshape2)
library(ggridges)
library(vegan)
library(readr)

# Here we have made oiur input inot smaller 100MB files and now read them in
setwd("rawdata/MetagenomicFirstRound/")
con <- pipe("cat part.gz_* | gunzip -c", "rb")
data <- read_csv(con)
close(con)
data2 <-read.csv("results.csv")
setwd("../../")
#metadata
metadata <- read.csv("metadata.csv")



# Minimum amount of damage filter
DamMin = 0.0
D_Max = 1.0
#Lambda Likelihood Ratio
LR = 0
# Minimum reads for parsing taxa
MinRead = 100
# Minimum mean read length
MinLength = 35


### first let's make some observations with the metadata 

plot(metadata$Dct,metadata$duplication.rate,pch=16,ylab="% duplication")
abline(v=0,col="darkred",lwd=3)
text(0.2,60,labels="Underamplified", adj=0)
text(-0.2,60,labels="Overamplified", adj=1)
text(metadata$Dct,metadata$duplication.rate+1.5,labels = metadata$age,cex=0.7)

plot(metadata$Dct,metadata$low.complexity,pch=16,ylab="% low complexity")
abline(v=0,col="darkred",lwd=3)
text(0.2,0.32,labels="Underamplified", adj=0)
text(-0.2,0.32,labels="Overamplified", adj=1)

plot(metadata$Dct,metadata$dedupReads,pch=16,ylab="Post fastp reads")
abline(v=0,col="darkred",lwd=3)
text(0.2,2.2*10^8,labels="Underamplified", adj=0)
text(-0.2,2.2*10^8,labels="Overamplified", adj=1)
text(metadata$Dct,metadata$dedupReads+5*10^6,labels = metadata$age,cex=0.7)

# Let's make an animal subset of the data 
df1 <- data %>% filter(MAP_damage > DamMin, N_reads >= MinRead, mean_L > MinLength, MAP_significance  > LR,  grepl("Metazoa",tax_path), grepl("genus", tax_rank), grepl("", sample))
inputdata <- df1
inputdata$sample2 <- factor(inputdata$sample,levels= sort(unique(inputdata$sample),decreasing = TRUE))
inputdata$age <- metadata$age[match(inputdata$sample2,metadata$Sample)]
inputdata$site <- as.factor(metadata$CoreID[match(inputdata$sample2,metadata$Sample)])

raw.metazoa <- inputdata

# plants
df1 <- data %>% filter(MAP_damage > DamMin, N_reads >= MinRead, mean_L > MinLength, MAP_significance  > LR,  grepl("Viridiplantae",tax_path), grepl("genus", tax_rank), grepl("", sample))
inputdata <- df1
inputdata$sample2 <- factor(inputdata$sample,levels= sort(unique(inputdata$sample),decreasing = TRUE))
inputdata$age <- metadata$age[match(inputdata$sample2,metadata$Sample)]
inputdata$site <- as.factor(metadata$CoreID[match(inputdata$sample2,metadata$Sample)])

raw.plants <- inputdata

# bacteria 

df1 <- data %>% filter(MAP_damage > DamMin, N_reads >= MinRead, mean_L > MinLength, MAP_significance  > LR,  grepl("Bacteria",tax_path), grepl("genus", tax_rank), grepl("", sample))
inputdata <- df1
inputdata$sample2 <- factor(inputdata$sample,levels= sort(unique(inputdata$sample),decreasing = TRUE))
inputdata$age <- metadata$age[match(inputdata$sample2,metadata$Sample)]
inputdata$site <- as.factor(metadata$CoreID[match(inputdata$sample2,metadata$Sample)])

raw.bacteria <- inputdata

ggplot(raw.metazoa, aes(x = MAP_damage, y = sample2, group = sample2)) +
  geom_density_ridges2(scale = 3, alpha = 0.55) +  # Plot density ridges first
  geom_point(position = position_jitter(height = 0.1, width = 0), alpha = 0.6) + # Add points above ridges
  facet_wrap(~ site, scales = "free_y") + # Split into facets by coreID
  labs(y = "Samples (ordered by age)", x = "MAP Damage",title = "Metazoa MAP Damage")+
  geom_text(aes(x = max(MAP_damage), label = age), hjust = -0.1, size = 3,nudge_y = 0.4)

ggplot(raw.metazoa, aes(x = mean_L, y = sample2, group = sample2)) +
  geom_density_ridges2(scale = 3, alpha = 0.55) +  # Plot density ridges first
  geom_point(position = position_jitter(height = 0.1, width = 0), alpha = 0.6) + # Add points above ridges
  facet_wrap(~ site, scales = "free_y") + # Split into facets by coreID
  labs(y = "Samples (ordered by age)", x = "MAP Damage",title = "Metazoa Mean Length")+
  geom_text(aes(x = max(MAP_damage), label = age), hjust = -0.1, size = 3,nudge_y = 0.4)

ggplot(raw.bacteria, aes(x = MAP_damage, y = sample2, group = sample2)) +
  geom_density_ridges2(scale = 3, alpha = 0.55) +  # Plot density ridges first
  geom_point(position = position_jitter(height = 0.1, width = 0), alpha = 0.6) + # Add points above ridges
  facet_wrap(~ site, scales = "free_y") + # Split into facets by coreID
  labs(y = "Samples (ordered by age)", x = "MAP Damage",title = "Bacteria MAP Damage")+
  geom_text(aes(x = max(MAP_damage), label = age), hjust = -0.1, size = 3,nudge_y = 0.4)

ggplot(raw.bacteria, aes(x = mean_L, y = sample2, group = sample2)) +
  geom_density_ridges2(scale = 3, alpha = 0.55) +  # Plot density ridges first
  geom_point(position = position_jitter(height = 0.1, width = 0), alpha = 0.6) + # Add points above ridges
  facet_wrap(~ site, scales = "free_y") + # Split into facets by coreID
  labs(y = "Samples (ordered by age)", x = "MAP Damage",title = "Bacteria Mean Length")+
  geom_text(aes(x = max(MAP_damage), label = age), hjust = -0.1, size = 3,nudge_y = 0.4)

ggplot(raw.plants, aes(x = MAP_damage, y = sample2, group = sample2)) +
  geom_density_ridges2(scale = 3, alpha = 0.55) +  # Plot density ridges first
  geom_point(position = position_jitter(height = 0.1, width = 0), alpha = 0.6) + # Add points above ridges
  facet_wrap(~ site, scales = "free_y") + # Split into facets by coreID
  labs(y = "Samples (ordered by age)", x = "MAP Damage",title = "Plants MAP Damage")+
  geom_text(aes(x = max(MAP_damage), label = age), hjust = -0.1, size = 3,nudge_y = 0.4)

ggplot(raw.plants, aes(x = mean_L, y = sample2, group = sample2)) +
  geom_density_ridges2(scale = 3, alpha = 0.55) +  # Plot density ridges first
  geom_point(position = position_jitter(height = 0.1, width = 0), alpha = 0.6) + # Add points above ridges
  facet_wrap(~ site, scales = "free_y") + # Split into facets by coreID
  labs(y = "Samples (ordered by age)", x = "MAP Damage",title = "Plants Mean Length")+
  geom_text(aes(x = max(MAP_damage), label = age), hjust = -0.1, size = 3,nudge_y = 0.4)




