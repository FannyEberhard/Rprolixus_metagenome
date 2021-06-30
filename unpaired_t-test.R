#######################################################################################################################
########################################### unpaired t-test for Enterococcus spp. #####################################

library(data.table)
library(tidyverse)
library(dplyr)
library(compareGroups)
setwd("C:/path/to/kaiju/results")

# prepare data
file.list <- list.files("C:/path/to/kaiju/results")
csv.files <- grep("^reworked",file.list,value = T)

k <- c(1:13)
for (file in csv.files){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("data_set")){
    data_set <- read.csv(file, header=F, sep=",")
    data_set <- cbind(sample = file, data_set)
    data_set <- data_set[,k]
  }
  
  # if the merged dataset does exist, append to it
  if (exists("data_set")){
    temp_dataset <- read.csv(file, header=F, sep=",")
    temp_dataset <- cbind(sample = file, temp_dataset)
    temp_dataset <- temp_dataset[,k]
    data_set <- rbind(data_set, temp_dataset)
    rm(temp_dataset)
  }
}

#remove duplicated Mouseblood rows generated in the loop before
data_set2 <- data_set %>% distinct()

#select only rows with "Enterococcus"
data_sub <- data_set2[apply(data_set2,1,function(x) {any(c("Enterococcus") %in% x)}),]

#sum the reads of each sample
data_sub2 <- aggregate(data_sub$V1, by = list(sample = data_sub$sample), FUN = sum)

# add total number of reads taxonomically assigned by Kaiju
data_sub2$total_reads <- c(19736914,10085720,5271180,6214590,1712660,11135028,1300392,4754789,4754789,6195056,3188448,12992736,1300392,
                           4754789,1070773,5772637,2653998,8582582,1642134,1052187,1052187,2207211,2996442,7557866,4416528,5100985,1901755,
                           8285328,2207211,5277561,2101574)
#percentage of Enterococcus reads compared with all respectively assigned reads
data_sub3 <- data_sub2 %>% group_by(total_reads) %>% mutate(x2= x*100/total_reads)
barplot(data_sub3$x2, names.arg = data_sub3$sample,las = 2, col = c("yellow","blue","blue","red","red","orange","orange"))

#keep only sample and scaled reads columns; remove Mouseblood sample
p <- c(1,4)
data_comp <- data_sub3[,p]
data_comp <- data_comp[-1,]
# split cruzi and rangeli
cru <- c(1,2,3,4,7,8,9,10,13,14,15,16,19,20,21,22,25,26,27,28)
data_cruzi <- data_comp[cru,]
ran <- c(1,2,5,6,7,8,11,12,13,14,17,18,19,20,23,24,25,26,29,30)
data_rang <- data_comp[ran,]

# have a look at the density
plot(density(data_cruzi$x2))
plot(density(data_rang$x2))

# remove timepoints T0 and T1 as infection might not be established yet
g <- c(1:8)
data_cruzi2 <- data_cruzi[-g,]
data_rang2 <- data_rang[-g,]

## Perform shapiro-wilk test
shapiro.test(data_cruzi2$x2) # p-value > 0.05 -> normally distributed
shapiro.test(data_rang2$x2) # p-value > 0.05 -> normally distributed

## Plot distribution
qqnorm(data_cruzi2$x2)
qqline(data_cruzi2$x2, col = 2)
qqnorm(data_rang2$x2)
qqline(data_rang2$x2, col = 2)
# accept null hypothesis -> data is normally distributed

# Bartlett test of homogeneity of variances
library(graphics)
bartlett.test(data_cruzi2$x2, data_cruzi2$group) # p-value above 0.05 -> accept hypothesis that variances are homogene
bartlett.test(data_rang2$x2, data_rang2$group) # p-value above 0.05 -> accept hypothesis that variances are homogene

# unpaired t-Test
# need numeric vectors for both groups
e <- c(1,2,5,6,9,10)
data_control <- data_cruzi2[e,]
data_control <- data_control[,2]
data_control_v <- as.vector(t(data_control))
ic <- c(3,4,7,8,11,12)
data_cruzi2_inf <- data_cruzi2[ic,]
data_cruzi2_inf <- data_cruzi2_inf[,2]
data_cruzi2_inf_v <- as.vector(t(data_cruzi2_inf))
data_rang2_inf <- data_rang2[ic,]
data_rang2_inf <- data_rang2_inf[,2]
data_rang2_inf_v <- as.vector(t(data_rang2_inf))

# t-test
res_cruzi <- t.test(data_control_v, data_cruzi2_inf_v, var.equal = TRUE)
res_cruzi
res_rang <- t.test(data_control_v, data_rang2_inf_v, var.equal = TRUE)
res_rang
# -> p-value both below 0.05, which means the two groups (cont - inf) are significantly different

# obtain standard deviations of mean
sd(data_control_v)
sd(data_cruzi2_inf_v)
sd(data_rang2_inf_v)