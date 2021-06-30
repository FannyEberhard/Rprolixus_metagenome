#######################################################################################################################
############################################# PCA #####################################################################

library(data.table)
library(tidyverse)
setwd("C:/path/to/kaiju/results")

file.list <- list.files("C:/path/to/kaiju/results")
csv.files <- grep("Out_all_order_reworked",file.list,value = T)

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

data_PCA <- prcomp(data_wider[,-1], center = TRUE, scale = TRUE)
summary(data_PCA)

plot(data_PCA$x[,1],data_PCA$x[,2], xlab="PC1 (60.91%)", ylab = "PC2 (10.52%)", main = "PC1 / PC2 - plot")

# include sample and group names
PCA_groups <- read.csv("C:/path/to/kaiju/results/PCA_groups.csv", header = T, sep = ",")
data_set_groups <- cbind(PCA_groups,data_wider)

# plot PC1 and PC2
library("factoextra")
fviz_pca_ind(data_PCA, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = data_set_groups$intestine, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = T,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Samples") +
  ggtitle("PC1 / PC2 - plot") +
  theme(plot.title = element_text(hjust = 0.5))

# 3D plot of PCA with first, second and third dimension
gr <- factor(data_set_groups[,3])
summary(gr)

install.packages("pca3d")
library(pca3d)
library(rgl)

rgl.open()

plot_3D <- pca3d(data_PCA,
                 group = gr,
                 axes.color = "black",
                 shape = "sphere",
                 radius = 0.9,
                 palette = c("#FFFFFF","#FF0000","#00FF00","#0000FF","#FF00FF","#000000"),# white (MB), red 0, green 1, 
                 #blue 2, purple3 , black 7
                 show.shadows = F,
                 show.scale = F,
                 fancy = F, 
                 bg = "white",
                 show.labels = T)

snapshotPCA3d(file="PCA_intestine.png")
rgl.snapshot(filename = "2d_plot.png")
snapshot3d(filename = "3plot.png", width = 4000, height = 4000)
rgl.postscript("plot.pdf",fmt="pdf")