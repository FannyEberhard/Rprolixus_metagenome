########################################################################################################################
######################################## Alpha diversity and species evenness ##########################################
library(vegan)
library(data.table)
library(tidyverse)
library(ggplot2)

setwd("C:/path/to/kaiju/results")

# prepare data
file.list <- list.files("C:/path/to/kaiju/results")
csv.files <- grep("Out",file.list,value = T)

for (file in csv.files){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("data_set")){
    data_set <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("data_set")){
    temp_dataset <- read.table(file, header=TRUE, sep=",")
    data_set <- rbind(data_set, temp_dataset)
    rm(temp_dataset)
  }
}

l <- c(1:25)
data_set <- data_set[-l,]
data_wider <- pivot_wider(data_set, names_from = order, values_from = reads, values_fill = 0, values_fn = sum)

# normalise rows
data_df <- as.data.frame(data_wider)
rownames(data_df) <- data_df[,1]
data_df <- data_df[,-1]
data_norm <- t(apply(data_df, 1, function(x)(x-min(x))/(max(x)-min(x))))

# alpha diversity
shannon <- diversity(data_norm[,-1], index="shannon") # 0 = complete uniformity, 5 = complete diversity
simpson <- diversity(data_norm[,-1], index = "simpson") # 0 = complete uniformity, 1 = complete diversity

## Species richness (S) and Pielou's evenness (J):
S <- specnumber(data_norm[,-1])
J <- shannon/log(S)
view(shannon)
plot(shannon)
plot(J)

######## Shannon Index
# make datatable shannon
shannon_table <-as.data.frame(shannon)
shannon_table$sample <- rownames(shannon_table)
shannon_table$no <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31)

# split shannon into T. cruzi and T. rangeli
c <- c(1,2,3,4,5,8,9,10,11,14,15,16,17,20,21,22,23,26,27,28,29)
shannon_cruzi <- shannon_table[c,-3]
shannon_cruzi$no <- c(1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,9,9,9,9)
r <- c(1,2,3,6,7,8,9,12,13,14,15,18,19,20,21,24,25,26,27,30,31)
shannon_rang <- shannon_table[r,-3]
shannon_rang$no <- c(1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,9,9,9,9)

# scatterplot shannon cruzi
z <- c(2,3,6,7,10,11,14,15,18,19)
w <- c(4,5,8,9,12,13,16,17,20,21)

gg <- ggplot(shannon_cruzi, aes (x = no, y = shannon)) +
  geom_point(aes(x = no, y = shannon), size = 3, shape = c(15,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17), 
             colour = c("goldenrod1", "midnightblue","midnightblue","firebrick","firebrick","midnightblue","midnightblue",
                        "firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick","midnightblue",
                        "midnightblue","firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick")) +
  geom_smooth(data=shannon_cruzi[z,],method='lm', formula= y~x, se = F, colour = "navy") +
  geom_smooth(data=shannon_cruzi[w,],method='lm', formula= y~x, se = F, colour = "indianred") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels=as.character(shannon_cruzi$no),breaks=shannon_cruzi$no)
plot(gg)
# scatterplot shannon rangeli
ggg <- ggplot(shannon_rang, aes (x = no, y = shannon)) +
  geom_point(aes(x = no, y = shannon), size = 3, shape = c(15,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17), 
             colour = c("goldenrod1", "midnightblue","midnightblue","firebrick","firebrick","midnightblue","midnightblue",
                        "firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick","midnightblue",
                        "midnightblue","firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick")) +
  geom_smooth(data=shannon_rang[z,],method='lm', formula= y~x, se = F, colour = "navy") +
  geom_smooth(data=shannon_rang[w,],method='lm', formula= y~x, se = F, colour = "indianred") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels=as.character(shannon_rang$no),breaks=shannon_rang$no)
plot(ggg)

######## Pielou's evenness
# make datatable pielou
pielou_table <-as.data.frame(J)
pielou_table$sample <- rownames(J)
pielou_table$no <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31)

# split pielou into T. cruzi and T. rangeli
c <- c(1,2,3,4,5,8,9,10,11,14,15,16,17,20,21,22,23,26,27,28,29)
pielou_cruzi <- pielou_table[c,]
pielou_cruzi$no <- c(1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,9,9,9,9)
r <- c(1,2,3,6,7,8,9,12,13,14,15,18,19,20,21,24,25,26,27,30,31)
pielou_rang <- pielou_table[r,]
pielou_rang$no <- c(1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,9,9,9,9)

# scatterplot pielou cruzi
ggp <- ggplot(pielou_cruzi, aes (x = no, y = J)) +
  geom_point(aes(x = no, y = J), size = 3, shape = c(15,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17), 
             colour = c("goldenrod1", "midnightblue","midnightblue","firebrick","firebrick","midnightblue","midnightblue",
                        "firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick","midnightblue",
                        "midnightblue","firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick")) +
  geom_smooth(data=pielou_cruzi[z,],method='lm', formula= y~x, se = F, colour = "navy") +
  geom_smooth(data=pielou_cruzi[w,],method='lm', formula= y~x, se = F, colour = "indianred") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels=as.character(pielou_cruzi$no),breaks=pielou_cruzi$no)
plot(ggp)
# scatterplot pielou rangeli
ggpp <- ggplot(pielou_rang, aes (x = no, y = J)) +
  geom_point(aes(x = no, y = J), size = 3, shape = c(15,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17,16,17), 
             colour = c("goldenrod1", "midnightblue","midnightblue","firebrick","firebrick","midnightblue","midnightblue",
                        "firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick","midnightblue",
                        "midnightblue","firebrick","firebrick","midnightblue","midnightblue","firebrick","firebrick")) +
  geom_smooth(data=pielou_rang[z,],method='lm', formula= y~x, se = F, colour = "navy") +
  geom_smooth(data=pielou_rang[w,],method='lm', formula= y~x, se = F, colour = "indianred") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels=as.character(pielou_rang$no),breaks=pielou_rang$no)
plot(ggpp)

# save as jpeg and pdf files
ggsave("Shannon_index_cruzi.jpeg", plot = gg, dpi = 300)
ggsave("Shannon_index_cruzi.pdf", plot = gg)
ggsave("Shannon_index_rangeli.jpeg", plot = ggg, dpi = 300)
ggsave("Shannon_index_rangeli.pdf", plot = ggg)
ggsave("Pielou_index_cruzi.jpeg", plot = ggp, dpi = 300)
ggsave("Pielou_index_cruzi.pdf", plot = ggp)
ggsave("Pielou_index_rangeli.jpeg", plot = ggpp, dpi = 300)
ggsave("Pielou_index_rangeli.pdf", plot = ggpp)