# Library
#install.packages('tidyverse')
library('tidyverse')
library("RColorBrewer")
library(pvclust)
#library(ComplexHeatmap)

library('dplyr')
library(tidyr)
library('readr')
library('tibble')
library('ggplot2')
library('purrr')

library(stringr)
library(forcats)
library(lubridate)


# read in data
pca <- read_table2("Mangoum_PCA.eigenvec", col_names = FALSE)
eigenval <- scan("Mangoum_PCA.eigenval")


# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


# sort out the individual species and pops
# Phenotype
Phenotype <- rep(NA, length(pca$ind))
Phenotype[grep("Kisumu_Susceptible", pca$ind)] <- "Kis.Sus"
Phenotype[grep("Mangoum_Permethrin1x", pca$ind)] <- "Man.P1x"
Phenotype[grep("Mangoum_Permethrin5x", pca$ind)] <- "Man.P5x"
Phenotype[grep("Mangoum_Permethrin10x", pca$ind)] <- "Man.P10x"
Phenotype[grep("Mangoum_Unexposed", pca$ind)] <- "Man.Unx"
Phenotype[grep("MK_Permethrin1x", pca$ind)] <- "MK.P1x"
Phenotype[grep("MK_Unexposed", pca$ind)] <- "MK.Unx"
# location
##loc <- rep(NA, length(pca$ind))
##loc[grep("Mak", pca$ind)] <- "makobe"
##loc[grep("Pyt", pca$ind)] <- "python"
# combine - if you want to plot each in different colours
##Phenotype_loc <- paste0(Phenotype, "_", loc)

# remake data.frame
## pca <- as.tibble(data.frame(pca, Phenotype, loc, Phenotype_loc))

pca <- as.tibble(data.frame(pca, Phenotype))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:7, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = Phenotype)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = 
                               c(
                                 "Kis.Sus" = "#E41A1C",    # Red
                                 "Man.P10x" = "#FF7F00",   # Orange
                                 "Man.P1x" = "#4DAF4A",    # Green
                                 "Man.P5x" = "#984EA3",    # Purple
                                 "Man.Unx" = "#377EB8",    # Blue
                                 "MK.P1x" = "#FFFF33",     # Yellow
                                 "MK.Unx" = "#A65628"      # Brown
                               ))
b <- b + coord_equal() + theme_classic()


png(file="PCA.snp.png", width = 2000, height = 2000, res = 300)
b + xlab(paste0("PC1: (", signif(pve$pve[1], 3), "% variance)")) + ylab(paste0("PC2: (", signif(pve$pve[2], 3), "% variance)"))
dev.off()


sampleDists <- read_table2("Mangoum_distance.dist", col_names = FALSE)
id <- read_table2("Mangoum_distance.dist.id", col_names = FALSE)

# Redifine ID
Phenotype <- rep(NA, length(id$X2))
Phenotype[grep("Kisumu_Susceptible", id$X2)] <- "Kis.Sus"
Phenotype[grep("Mangoum_Permethrin1x", id$X2)] <- "Man.P1x"
Phenotype[grep("Mangoum_Permethrin5x", id$X2)] <- "Man.P5x"
Phenotype[grep("Mangoum_Permethrin10x", id$X2)] <- "Man.P10x"
Phenotype[grep("Mangoum_Unexposed", id$X2)] <- "Man.Unx"
Phenotype[grep("MK_Permethrin1x", id$X2)] <- "MK.P1x"
Phenotype[grep("MK_Unexposed", id$X2)] <- "MK.Unx"


# Add Phenotype column to id tibble
id$Phenotype <- Phenotype


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- id$Phenotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(file="sample_to_sample_distances.plot.png", width = 1600, height = 1200, pointsize = 12, res = 100)
heatmap(sampleDistMatrix, labRow=rownames(data), labCol = FALSE, scale="column", cexRow=1.5, col= colors)
dev.off()


# Calculate hierarchical clustering 

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- id$Phenotype
colnames(sampleDistMatrix) <- NULL


# Calculate hierarchical clustering  

# Convert tibble to matrix
dist_matrix <- as.matrix(sampleDists)

# Convert to 'dist' object
dist_obj <- as.dist(dist_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_obj)


# Plotting the dendrogram  
png(file="sample_Dendrogram.plot.png", width = 4000, height = 2500, res = 300)
plot(hc, labels = rownames(sampleDistMatrix), main = "Sample Dendrogram", xlab = "Samples", ylab = "Distance", hang = -1)
dev.off()


# Convert to matrix
dist_matrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- id$Phenotype
colnames(sampleDistMatrix) <- NULL

# Perform bootstrap hierarchical clustering
# Note: pvclust expects raw data, not a distance matrix, so we use the matrix directly
pv_result <- pvclust(dist_matrix, method.hclust = "average", method.dist = "euclidean", nboot = 1000)

# Save plot
png(file = "bootstrap_dendrogram.png", width = 4000, height = 2500, res = 400)
plot(pv_result, labels = rownames(sampleDistMatrix), print.pv = "bp", print.edge = FALSE)
pvrect(pv_result, alpha = 0.95) # Highlight clusters with AU p-value > 95%
dev.off()

