#######################################################################################################################
##################################### heatmap of functional analysis ##################################################
install.packages("gplots")
library(gplots)
library(RColorBrewer)
install.packages("tidytable")
library(tidytable)

setwd("C:/path/to/KEGG/results")
# Documentation heatmap.2: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2

# data preparation
data <- read.csv("kegg-metabolism_modules.csv", header = T)
a <- c(2,3,9)
data <- data[,a]
data_wider <- pivot_wider.(data, names_from = "kegg_module", values_from = "module_completeness", values_fill = 0)

# keep only data for MAGs
b <- c(24,25,26,27,28,30)
data_wider <- data_wider[b,]

# assign labels in column 1 to "rnames"
rnames <- data_wider[,1]
# convert rnames into matrix and remove header, otherwise error in 'rownames(mat_data) <- rnames2'
rnames2 <- as.matrix(rnames)
rnames2 <- matrix(rnames2, dimnames = NULL)

# transform columns into a matrix
mat_data <- data.matrix(data_wider[,2:ncol(data_wider)])
# assign row names
rownames(mat_data) <- rnames2

# also create matrix which has sorted module_category
write.csv(mat_data, "mat_data.csv", sep = ",", dec = ".", quote = F)
# after sorting the columns (modules) by their module_category in Excel
mat_data2 <- read.csv("mat_data.csv", header = T)
row.names(mat_data2) <- mat_data2[,1]
mat_data2 <- mat_data2[,-1]
mat_data3 <- as.matrix(mat_data2)

install.packages("pheatmap")
library(pheatmap)

# insert column annotation (module_categories)
data2 <- read.csv("kegg-metabolism_modules.csv", header = T)
d <- c(3,6,7)
data2 <- data2[,d]
data3 <- unique(data2)
data3 <- data3[order(data3$module_category),]
row.names(data3) <- data3[,1]
data3 <- data3[,-1]
data3 <- data3[,1, drop = F]
colnames(data3) <- c("Module_category","Module_subcategory")
data3$Module_category<-gsub(" ", "_", data3$Module_category) # replace space in cells with underscore
data3$Module_subcategory<-gsub(" ", "_", data3$Module_subcategory)
data3$Module_subcategory<-gsub("-", "_", data3$Module_subcategory)

# insert row annotations (bacterial families)
Bin <- rownames(mat_data)
Family <- c("Enterobacterales", "Lactobacillales","Actinomycetales","Enterobacterales","Actinomycetales","Rickettsiales")
my_row_anno <- data.frame(Bin,Family)
rownames(my_row_anno) <- my_row_anno[,1]
my_row_anno <- my_row_anno[,2, drop = F]

# annotation colors
aka3 = list(Module_category = c(Amino_acid_metabolism = "darkorange", Biosynthesis_of_other_secondary_metabolites = "darkgreen",
                                Biosynthesis_of_terpenoids_and_polyketides = "goldenrod4",
                                Carbohydrate_metabolism = "khaki3", Energy_metabolism = "orangered", Gene_set = "gray25",
                                Glycan_metabolism = "darkolivegreen3", Lipid_metabolism = "goldenrod1", 
                                Metabolism_of_cofactors_and_vitamins = "indianred3", Module_set = "wheat1", 
                                Nucleotide_metabolism = "darkred", Xenobiotics_biodegradation = "yellow4"),
            Family = c(Actinomycetales = "orangered3", Enterobacterales = "red4", Lactobacillales = "royalblue4", 
                       Rickettsiales = "lightblue3"),
            Module_subcategory = c(Arginine_and_proline_metabolism = "lightgoldenrodyellow", Aromatic_amino_acid_metabolism = "orangered",
                                   Branched_chain_amino_acid_metabolism = "orange", Cysteine_and_methionine_metabolism = "sandybrown",
                                   Histidine_metabolism = "peachpuff2", Lysine_metabolism = "peru", Other_amino_acid_metabolism = "tan1",
                                   Polyamine_biosynthesis = "orange3", Serine_and_threonine_metabolism = "khaki2",
                                   Biosynthesis_of_beta_lactams = "darkolivegreen3", Biosynthesis_of_other_antibiotics = "gold4",
                                   Biosynthesis_of_other_bacterial_compounds = "darkgreen", Biosynthesis_of_phytochemical_compounds = "limegreen",
                                   Enediyne_biosynthesis = "orange4", Plant_terpenoid_biosynthesis = "navajowhite3", 
                                   Polyketide_sugar_unit_biosynthesis = "saddlebrown", Terpenoid_backbone_biosynthesis = "sienna3",
                                   Type_II_polyketide_biosynthesis = "tan", Central_carbohydrate_metabolism = "seashell", 
                                   Other_carbohydrate_metabolism = "peachpuff3", ATP_synthesis = "lightsalmon1", Carbon_fixation = "red2",
                                   Methane_metabolism = "orangered4", Nitrogen_metabolism = "hotpink3", Photosynthesis = "palevioletred", 
                                   Sulfur_metabolism = "firebrick2", Drug_resistance = "gray39", Pathogenicity = "gray10", Plant_pathogenicity = "gray85",
                                   Glycan_biosynthesis = "darkseagreen3", Glycosaminoglycan_metabolism = "goldenrod1", Lipopolysaccharide_metabolism = "darkolivegreen1",
                                   Fatty_acid_metabolism = "darkorange2", Lipid_metabolism = "gold", Sterol_biosynthesis = "khaki1", 
                                   Cofactor_and_vitamin_metabolism = "lightpink2", Metabolic_capacity = "ivory", Purine_metabolism = "tomato",
                                   Pyrimidine_metabolism = "tomato4", Aromatics_degradation = "yellowgreen"))

my_pheatmap <- pheatmap(mat_data3, 
                        annotation_col = data3,
                        annotation_row = my_row_anno,
                        annotation_colors = aka3,
                        annotation_names_col = T,
                        show_colnames = F,
                        cluster_cols = F, # if TRUE, columns are clustered by module completion
                        gaps_col = c(42,50,62,100,148,163,183,207,251,259,267),
                        border_color = "gray",
                        cellwidth = 3,
                        cellheight = 35,
                        treeheight_row = 40,
                        drop_levels = T)


# saving the heatmap
save_pheatmap_pdf <- function(x, filename, width=19, height=21) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(my_pheatmap, "heatmap.pdf")