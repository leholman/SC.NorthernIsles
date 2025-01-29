## scrappy analysis for fun! 

library(tidyverse)
library(reshape2)
library(ggridges)
library(vegan)


metadata <- read.csv("metadata.csv")
data <- read_csv("results.csv")

# Minimum amount of damage filter
DamMin = 0.0
D_Max = 1.0
#Lambda Likelihood Ratio
LR = 0
# Minimum reads for parsing taxa
MinRead = 100
# Minimum mean read length
MinLength = 35


# Let's make an animal subset of the data 
df1 <- data %>% filter(MAP_damage > DamMin, N_reads >= MinRead, mean_L > MinLength, MAP_significance  > LR,  grepl("Metazoa",tax_path), grepl("", tax_rank), grepl("", sample))
inputdata <- df1
inputdata$sample2 <- factor(inputdata$sample,levels= sort(unique(inputdata$sample),decreasing = TRUE))
inputdata.genus <- df1[df1$tax_rank=="genus",]
genus.inputdata.wide <- dcast(inputdata.genus, tax_name ~ sample, value.var="N_reads")
genus.inputdata.wide[is.na(genus.inputdata.wide)] <- 0

rownames(genus.inputdata.wide) <- genus.inputdata.wide$tax_name
genus.inputdata.wide <- genus.inputdata.wide[,-1]

write.csv(genus.inputdata.wide,"~/Desktop/QucikLook/datawide.animal.csv")

  
  # Let's make a fish subset of the data 
df1 <- data %>% filter(MAP_damage > DamMin, N_reads >= MinRead, mean_L > MinLength, MAP_significance  > LR,  grepl("Actinopterygii",tax_path), grepl("", tax_rank), grepl("", sample))
inputdata < - df1
inputdata$sample2 <- factor(inputdata$sample2,levels= sort(unique(inputdata$sample2),decreasing = TRUE))
inputdata.genus <- df1[df1$tax_rank=="genus",]
genus.inputdata.wide <- dcast(inputdata.genus, tax_name ~ sample, value.var="N_reads")
genus.inputdata.wide[is.na(genus.inputdata.wide)] <- 0

genus.inputdata.wide[is.na(genus.inputdata.wide)] <- 0
rownames(genus.inputdata.wide) <- genus.inputdata.wide$tax_name
genus.inputdata.wide <- genus.inputdata.wide[,-1]


write.csv(genus.inputdata.wide,"~/Desktop/QucikLook/datawide.fish.csv")



df1 <- data %>% filter(MAP_damage > DamMin, N_reads >= MinRead, mean_L > MinLength, MAP_significance  > LR,  grepl("Actinopterygii",tax_path), grepl("", tax_rank), grepl("", sample))
inputdata < - df1
inputdata$sample2 <- factor(inputdata$sample2,levels= sort(unique(inputdata$sample2),decreasing = TRUE))
inputdata.genus <- df1[df1$tax_rank=="genus",]
genus.inputdata.wide <- dcast(inputdata.genus, tax_name ~ sample, value.var="N_reads")
genus.inputdata.wide[is.na(genus.inputdata.wide)] <- 0

genus.inputdata.wide[is.na(genus.inputdata.wide)] <- 0
rownames(genus.inputdata.wide) <- genus.inputdata.wide$tax_name
genus.inputdata.wide <- genus.inputdata.wide[,-1]


write.csv(genus.inputdata.wide,"~/Desktop/QucikLook/datawide.fish.csv")












