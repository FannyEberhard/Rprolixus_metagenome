#######################################################################################################################
#################################### Visualisation of Kaiju results ###################################################
# Data prep
setwd("C:/path/to/kaiju/results")

# sort table so that taxonomic group 'order' is in one column (remove Terrbacteria group, PVC group and so on)
file.list <- list.files("C:/path/to/kaiju/results")
csv.files <- grep("kaiju",file.list,value = T)

for (file in csv.files){
  data <- read.csv(file , header = F)
  colnames(data) <- c("reads", "root", "life_form","kingdom", "phylum", "class", "order", "family","genus", "species", "subspecies")
  # keep only Bacteria and remove Viruses, Fungi...
  data <- data[str_detect(data$kingdom, "Bacteria"), ]
  
  # move rows with specific cells to the left
  data[data == "Terrabacteria group"] <- NA
  data[] <-  t(apply(data, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
  
  data[data == "PVC group"] <- NA
  data[] <-  t(apply(data, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
  
  data[data == "FCB group"] <- NA
  data[] <-  t(apply(data, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
  
  data[data == "Bacteroidetes/Chlorobi group"] <- NA
  data[] <-  t(apply(data, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
  
  data[data == "delta/epsilon subdivisions"] <- NA
  data[] <-  t(apply(data, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
  
  data[data == "unclassified Bacteria"] <- NA
  data[] <-  t(apply(data, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
  
  names(data) <- NULL
  
  # save table
  write.csv(data, paste0("reworked_",file), row.names = F, quote = F)
}

# reduce taxon table
file.list2 <- list.files("C:/path/to/kaiju/results")
reworked.files <- grep("reworked",file.list2,value = T)

for(i in reworked.files){
  # read csv file with unknown number of columns, firstly 100 columns are used, then all with NAs are discarded
  # otherwise last column is weirdly used as new row
  kaiju_results <- read.csv(paste("C:/path/to/kaiju/results/",i, sep = ""), 
                            col.names = paste("V",1:100), fill = T, sep = ",")
  kaiju_results <- kaiju_results[,which(!is.na(kaiju_results[1,]))]
  
  
  # keep only columns reads,kingdom and order
  kaiju_results_pp <- kaiju_results[,c(1,4,7)]
  # insert column names
  colnames(kaiju_results_pp)[1] <- "reads"
  colnames(kaiju_results_pp)[2] <- "kingdom"
  colnames(kaiju_results_pp)[3] <- "order"
  # keep only rows with Bacteria in kingdom
  kaiju_results_pp <- kaiju_results_pp[str_detect(kaiju_results_pp$kingdom, "Bacteria"), ]
  
  # replace NAs in column order with "unclassified"
  kaiju_results_pp$order[kaiju_results_pp$order==""]<-"unclassified"
  
  kaiju_results_ppp <- kaiju_results_pp
  
  # replace "Bacteria candidate order" with "unclassified"
  grep("Bacteria candidate order", kaiju_results_ppp)
  kaiju_results_pppp <- data.frame(lapply(kaiju_results_ppp, function(x) {gsub("Bacteria candidate order",
                                                                               "unclassified", x)}))
  
  # remove kingdom column
  kaiju_results_ppppp <- kaiju_results_pppp[,-2]
  
  # transform reads variable to numeric and sum reads value of all rows with the same string in column "order"
  kaiju_results_pppppp <- transform(kaiju_results_ppppp, reads = as.numeric(reads))
  kaiju_results_ppppppp <- aggregate(kaiju_results_pppppp[,-2], kaiju_results_pppppp["order"], sum)
  
  # give column "reads" colname back and insert column with sample name
  colnames(kaiju_results_ppppppp)[2] <- "reads"
  
  kaiju_results_ppppppp['sample'] = i
  sample_names <- kaiju_results_ppppppp[,ncol(kaiju_results_ppppppp)]
  kaiju_results_ppppppp <- kaiju_results_ppppppp[,-ncol(kaiju_results_ppppppp)]
  kaiju_results_ppppppp <- cbind(sample_names,kaiju_results_ppppppp)
  # and remove "kaiju_taxonpaths_" and ".csv" from sample name
  kaiju_results_pppppppp <- data.frame(lapply(kaiju_results_ppppppp, function(x) {gsub("reworked_kaiju.taxonpaths_",
                                                                                       "", x)}))
  kaiju_results_pppppppp <- data.frame(lapply(kaiju_results_pppppppp, function(x) {gsub(".csv",
                                                                                        "", x)}))
  
  # order by number of reads in column "reads"
  kaiju_results_pppppppp$reads <- as.numeric(as.character(kaiju_results_pppppppp$reads))
  kaiju_results_ppppppppp <- kaiju_results_pppppppp[order(kaiju_results_pppppppp$reads, decreasing = TRUE),]
  
  # sum of all reads and remove all bacterial orders which mapped less than 0.7% of all reads
  summe <- sum(kaiju_results_ppppppppp$reads)
  one_percent <- (summe*0.7)/100
  kaiju_results_finish <- kaiju_results_ppppppppp[kaiju_results_ppppppppp$reads >= one_percent,]
  
  write.csv(kaiju_results_finish, paste0("Out_all_order_",i), row.names = F, quote = F)
}

####################################### combine csv-files into four data.frames ########################################
file.list3 <- list.files("C:/path/to/kaiju/results")
reworked.files2 <- grep("Out_all_order_reworked",file.list3,value = T)

y <- c("ContAM","cruzAM")
x <- c("ContPM","cruzPM")
z <- c("ContAM","rangAM")
w <- c("ContPM","rangPM")

cruz_AM <- grep(paste(y, collapse = "|"), reworked.files2,value = T)
cruz_PM <- grep(paste(x, collapse = "|"), reworked.files2,value = T)
rang_AM <- grep(paste(z, collapse = "|"), reworked.files2,value = T)
rang_PM <- grep(paste(w, collapse = "|"), reworked.files2,value = T)

# T. cruzi AM
for (file in cruz_AM){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("T_cruz_AM")){
    T_cruz_AM <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("T_cruz_AM")){
    temp_dataset <- read.table(file, header=TRUE, sep=",")
    T_cruz_AM <- rbind(T_cruz_AM, temp_dataset)
    rm(temp_dataset)
  }
}

# remove first sample's rows as they are doubled because of the loop before
l <- c(1:23)
T_cruz_AM <- T_cruz_AM[-l,]

# T.cruzi PM
for (file in cruz_PM){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("T_cruz_PM")){
    T_cruz_PM <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("T_cruz_PM")){
    temp_dataset <- read.table(file, header=TRUE, sep=",")
    T_cruz_PM <- rbind(T_cruz_PM, temp_dataset)
    rm(temp_dataset)
  }
}

# remove first sample's rows as they are doubled because of the loop before
l <- c(1:13)
T_cruz_PM <- T_cruz_PM[-l,]

# T. rangeli AM
for (file in rang_AM){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("T_rang_AM")){
    T_rang_AM <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("T_rang_AM")){
    temp_dataset <- read.table(file, header=TRUE, sep=",")
    T_rang_AM <- rbind(T_rang_AM, temp_dataset)
    rm(temp_dataset)
  }
}

# remove first sample's rows as they are doubled because of the loop before
l <- c(1:23)
T_rang_AM <- T_rang_AM[-l,]

# T. rangeli PM
for (file in rang_PM){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("T_rang_PM")){
    T_rang_PM <- read.table(file, header=TRUE, sep=",")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("T_rang_PM")){
    temp_dataset <- read.table(file, header=TRUE, sep=",")
    T_rang_PM <- rbind(T_rang_PM, temp_dataset)
    rm(temp_dataset)
  }
}

# remove first sample's rows as they are doubled because of the loop before
l <- c(1:13)
T_rang_PM <- T_rang_PM[-l,]


# change column names
colnames(T_cruz_AM) <- c("Sample", "Order","Reads")
colnames(T_cruz_PM) <- c("Sample", "Order","Reads")
colnames(T_rang_AM) <- c("Sample", "Order","Reads")
colnames(T_rang_PM) <- c("Sample", "Order","Reads")

#################################### create stacked barplots for all four ##############################################
library(ggplot2)
q <- c("midnightblue","red4")
# T.cruzi AM
plot_T_cruzi_AM <- ggplot(T_cruz_AM, aes(fill=Order, y=Reads, x=Sample)) + 
  geom_bar(position="fill", stat="identity", width = 0.9)+ 
  scale_fill_manual(values = c("darkkhaki","lightgoldenrod","darkseagreen1","mediumturquoise","tomato2","maroon","darkolivegreen4",
                               "darkolivegreen3","firebrick","goldenrod1","darkseagreen3","darkseagreen4","darkred","goldenrod3","steelblue2",
                               "goldenrod4","slategray1","orchid4","steelblue4","orangered3","moccasin","steelblue1","coral4",
                               "tan4", "gray68", "midnightblue", "forestgreen")) + 
  labs(x = "Sample", y = "Bacterial reads mapped [%]") +
  theme(axis.text.x = element_text(angle = -90, colour = q)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
plot(plot_T_cruzi_AM)

# T.cruzi PM
plot_T_cruzi_PM <- ggplot(T_cruz_PM, aes(fill=Order, y=Reads, x=Sample)) + 
  geom_bar(position="fill", stat="identity", width = 0.9)+ 
  scale_fill_manual(values = c("darkkhaki","lightgoldenrod","darkseagreen1","mediumturquoise","tomato2","maroon","darkolivegreen4",
                               "darkolivegreen3","goldenrod1","darkseagreen3","darkseagreen4","darkred","goldenrod3","steelblue2",
                               "slategray1","steelblue4","orangered3","steelblue1","coral4",
                               "tan4", "gray68", "midnightblue")) + 
  labs(x = "Sample", y = "Bacterial reads mapped [%]") +
  theme(axis.text.x = element_text(angle = -90, colour = q)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
plot(plot_T_cruzi_PM)

# T.rangeli AM
plot_T_rang_AM <- ggplot(T_rang_AM, aes(fill=Order, y=Reads, x=Sample)) + 
  geom_bar(position="fill", stat="identity", width = 0.9)+ 
  scale_fill_manual(values = c("darkkhaki","lightgoldenrod","darkseagreen1","mediumturquoise","tomato2","maroon","darkolivegreen4",
                               "darkolivegreen3","goldenrod1","darkseagreen3","darkseagreen4","darkred","goldenrod3","steelblue2",
                               "slategray1","orchid4","steelblue4","orangered3","moccasin","steelblue1","coral4",
                               "tan4", "gray68", "midnightblue", "forestgreen")) + 
  labs(x = "Sample", y = "Bacterial reads mapped [%]") +
  theme(axis.text.x = element_text(angle = -90, colour = q)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
plot(plot_T_rang_AM)

# T.rangeli PM
plot_T_rang_PM <- ggplot(T_rang_PM, aes(fill=Order, y=Reads, x=Sample)) + 
  geom_bar(position="fill", stat="identity", width = 0.9)+ 
  scale_fill_manual(values = c("darkkhaki","lightgoldenrod","darkseagreen1","mediumturquoise","tomato2","maroon","darkolivegreen4",
                               "darkolivegreen3","goldenrod1","darkseagreen3","darkseagreen4","darkred","goldenrod3","steelblue2",
                               "slategray1","steelblue4","orangered3","steelblue1","coral4",
                               "tan4", "gray68", "midnightblue")) + 
  labs(x = "Sample", y = "Bacterial reads mapped [%]") +
  theme(axis.text.x = element_text(angle = -90, colour = q)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
plot(plot_T_rang_PM)

# save as jpeg and pdf files
ggsave("T.cruzi_AM_kaiju.jpeg", plot = plot_T_cruzi_AM, dpi = 300)
ggsave("T.cruzi_AM_kaiju.pdf", plot = plot_T_cruzi_AM)
ggsave("T.cruzi_PM_kaiju.jpeg", plot = plot_T_cruzi_PM, dpi = 300)
ggsave("T.cruzi_PM_kaiju.pdf", plot = plot_T_cruzi_PM)
ggsave("T.rang_AM_kaiju.jpeg", plot = plot_T_rang_AM, dpi = 300)
ggsave("T.rang_AM_kaiju.pdf", plot = plot_T_rang_AM)
ggsave("T.rang_PM_kaiju.jpeg", plot = plot_T_rang_PM, dpi = 300)
ggsave("T.rang_PM_kaiju.pdf", plot = plot_T_rang_PM)