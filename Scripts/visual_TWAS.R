# PLOTTING FST Values in ggplot2
# Loading library

library(tidyverse)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(grid)
library(circlize)


## FST analysis

Fst <- read.delim("transcript_Mangoum_RNAseqfst.snp.txt", header = T, sep = "\t", dec = ".")

Fst <- Fst %>%
  filter(chrom %in% c("2L", "2R", "3L", "3R", "X"))

Fst[is.na(Fst)] <- 0
Fst[Fst < 0] <- 0

Fst$start <- as.integer(Fst$start)
Fst$end <- as.integer(Fst$end)

names(Fst)
str(Fst)


## Diversity analysis

div <- read.delim("Mangoum_diversity.snp_cleaned.txt", header = T, sep = "\t", dec = ".")

names(div)
div$start <- as.integer(div$start)
div$end <- as.integer(div$end)

Pi <- div[,c(1,2,3,4,7,10, 13, 16, 19, 22)]
Taj <- div[,c(1,2,3,6,9,12, 15, 18, 21, 24)]

Taj[is.na(Taj)] <- 0
Pi[is.na(Pi)] <- 0

## tajimas_d analysis

names(Taj)
str(Taj)
head(Taj)


library(dplyr)

Annotation <- read.delim("Annotation.txt", header = TRUE, sep = "\t", dec = ".") %>%
  mutate(
    Start = as.numeric(Start),
    End   = as.numeric(End),
    start = Start + 1,
    end   = End,
    chrom = Chromosome
  ) %>%
  select(chrom, start, end, GeneID, Gene_Name)


Annotation$start <- as.integer(Annotation$start)
Annotation$end <- as.integer(Annotation$end)

str(Annotation)
head(Annotation, 2)


Div_params <- Annotation %>%
  left_join(Fst, by = c("chrom", "start", "end"))%>%
  left_join(Taj, by = c("chrom", "start", "end")) %>%
  left_join(Pi, by = c("chrom", "start", "end")) 
  

names(Div_params)
head(Div_params,2)

write.xlsx(Div_params, "Div_params.xlsx")


GOI <- read_excel("GOI.xlsx")

Div_params_GOI <- Div_params %>%
  filter(Gene_Name %in% GOI$Gene_Name) %>%
  select(GeneID, Gene_Name, chrom, start, end, Fst_Man_CS = "Kisumu_Susceptible.MK_Unexposed", 
         Fst_Man_RS = "Kisumu_Susceptible.Mangoum_Permethrin1x", 
         Fst_Man_RC = "Mangoum_Unexposed.MK_Permethrin1x", 
         Fst_Man_10x.1x = "Mangoum_Permethrin10x.Mangoum_Permethrin1x",
         Fst_MK_RC = "MK_Permethrin1x.MK_Unexposed",
         Taj_Susceptible = "Kisumu_Susceptible.tajimas_d",
         Taj_Hybrid = "MK_Unexposed.tajimas_d",
         Taj_Field = "Mangoum_Unexposed.tajimas_d",
         Pi_Susceptible = "Kisumu_Susceptible.theta_pi",
         Pi_Hybrid = "MK_Unexposed.theta_pi",
         Pi_Field = "Mangoum_Unexposed.theta_pi")


head(Div_params_GOI)


write.xlsx(Div_params_GOI, "Div_params_GOI.xlsx")


plot_Fst_Man_RS <- ggplot(Div_params_GOI, aes(x = Gene_Name, y = Fst_Man_RS)) +
  geom_col(fill = "steelblue") +   # bar plot
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate labels
  ) +
  labs(
    x = "Gene Name",
    y = "FST",
    title = "Field Res. vs Sus. FST values per Gene "
  )+
  ylim(0, 1) +
  theme(
    plot.title = element_text(face = "italic", colour = "black"),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  )



png("plot_Fst_Man_RS.png", width = 6000, height = 3000, res = 400)
plot_Fst_Man_RS
dev.off()


plot_Fst_Man_CS <- ggplot(Div_params_GOI, aes(x = Gene_Name, y = Fst_Man_CS)) +
  geom_col(fill = "steelblue") +   # bar plot
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate labels
  ) +
  labs(
    x = "Gene Name",
    y = "FST",
    title = "Field Cntr. vs Sus. FST values per Gene "
  )+
  ylim(0, 1) +
  theme(
    plot.title = element_text(face = "italic", colour = "black"),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  )



png("plot_Fst_Man_CS.png", width = 6000, height = 3000, res = 400)
plot_Fst_Man_CS
dev.off()



plot_Fst_Man_RC <- ggplot(Div_params_GOI, aes(x = Gene_Name, y = Fst_Man_RC)) +
  geom_col(fill = "steelblue") +   # bar plot
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate labels
  ) +
  labs(
    x = "Gene Name",
    y = "FST",
    title = "Field Res. vs Cntr. FST values per Gene "
  )+
  ylim(0, 1) +
  theme(
    plot.title = element_text(face = "italic", colour = "black"),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  )



png("plot_Fst_Man_RC.png", width = 6000, height = 3000, res = 400)
plot_Fst_Man_RC
dev.off()


plot_Fst_Man_10x.1x <- ggplot(Div_params_GOI, aes(x = Gene_Name, y = Fst_Man_10x.1x)) +
  geom_col(fill = "steelblue") +   # bar plot
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate labels
  ) +
  labs(
    x = "Gene Name",
    y = "FST",
    title = "Field Perm10x. vs Perm1x. FST values per Gene "
  )+
  ylim(0, 1) +
  theme(
    plot.title = element_text(face = "italic", colour = "black"),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  )



png("plot_Fst_Man_10x.1x.png", width = 6000, height = 3000, res = 400)
plot_Fst_Man_10x.1x
dev.off()



plot_Fst_MK_RC <- ggplot(Div_params_GOI, aes(x = Gene_Name, y = Fst_MK_RC)) +
  geom_col(fill = "steelblue") +   # bar plot
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate labels
  ) +
  labs(
    x = "Gene Name",
    y = "FST",
    title = "Hybrid Res. vs Cntr. FST values per Gene "
  )+
  ylim(0, 1) +
  theme(
    plot.title = element_text(face = "italic", colour = "black"),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  )



png("plot_Fst_MK_RC.png", width = 6000, height = 3000, res = 400)
plot_Fst_MK_RC
dev.off()




library(gridExtra)
library(grid) 

# Function to add plot labels  
add_plot_label <- function(plot, label) {  
  plot +  
    ggplot2::labs(tag = label) +  
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 10, face = "bold"))  
}  

# Apply the function to your plots  
plot_Fst_Man_RS_labeled <- add_plot_label(plot_Fst_Man_RS, "A")  
plot_Fst_Man_CS_labeled <- add_plot_label(plot_Fst_Man_CS, "B")  
plot_Fst_Man_RC_labeled <- add_plot_label(plot_Fst_Man_RC, "C")
plot_Fst_Man_10x.1x_labeled <- add_plot_label(plot_Fst_Man_10x.1x, "D")
plot_Fst_MK_RC_labeled <- add_plot_label(plot_Fst_MK_RC, "E")

library(ggpubr)


# Arrange the plots using grid.arrange  

Fst_plot <- ggarrange(plot_Fst_Man_RS_labeled, 
                      plot_Fst_Man_CS_labeled,
                      plot_Fst_Man_RC_labeled,
                      plot_Fst_Man_10x.1x_labeled,
                      plot_Fst_MK_RC_labeled,
                      ncol=1, nrow=5, legend = FALSE)

Fst_plot


png("Fst_plot.png", width = 7000, height = 8000, res = 500)
Fst_plot
dev.off()













Taj_params_GOI <- Div_params %>%
  filter(Gene_Name %in% GOI$Gene_Name) %>%
  select(GeneID, Gene_Name, chrom, start, end, Sus. = "Kisumu_Susceptible.tajimas_d",
         Hybrid = "MK_Unexposed.tajimas_d",
         Field = "Mangoum_Unexposed.tajimas_d")

library(reshape2)

Taj_params_GOI.melt <- reshape2::melt(
  data = Taj_params_GOI,
  measure.vars = c("Sus.", "Hybrid", "Field"),
  value.name = "Taj",
  variable.name = "Comparison"
) 

head(Taj_params_GOI.melt, 2)

Taj_plot <- ggplot(Taj_params_GOI.melt, aes(x = Gene_Name, y = Taj)) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = c(-1, -2, -3, -4), linetype = "dashed", color = "gray50") +
  facet_wrap(~ Comparison, nrow = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene Name", y = "Tajima's D", title = "Tajima's D per Gene by Population") +
  theme(
    plot.title = element_text(face = "italic", colour = "black"),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  ) 

png("Taj_plot.png", width = 6000, height = 3000, res = 400)
Taj_plot
dev.off()




Pi_params_GOI <- Div_params %>%
  filter(Gene_Name %in% GOI$Gene_Name) %>%
  select(GeneID, Gene_Name, chrom, start, end, Sus. = "Kisumu_Susceptible.theta_pi",
         Hybrid = "MK_Unexposed.theta_pi",
         Field = "Mangoum_Unexposed.theta_pi")

library(reshape2)

Pi_params_GOI.melt <- reshape2::melt(
  data = Pi_params_GOI,
  measure.vars = c("Sus.", "Hybrid", "Field"),
  value.name = "Pi",
  variable.name = "Comparison"
) 

head(Pi_params_GOI.melt, 2)

Pi_plot <- ggplot(Pi_params_GOI.melt, aes(x = Gene_Name, y = Pi)) +
  geom_point(color = "steelblue", size = 2) +
#  geom_hline(yintercept = c(-1, -2, -3, -4), linetype = "dashed", color = "gray50") +
  facet_wrap(~ Comparison, nrow = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene Name", y = "Nucleotide diversity", title = "Diversity per Gene by Population") +
  theme(
    plot.title = element_text(face = "italic", colour = "black"),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  ) 

png("Pi_plot.png", width = 6000, height = 3000, res = 400)
Pi_plot
dev.off()



library(gridExtra)
library(grid) 

# Function to add plot labels  
add_plot_label <- function(plot, label) {  
  plot +  
    ggplot2::labs(tag = label) +  
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 10, face = "bold"))  
}  

# Apply the function to your plots  
Taj_plot_labeled <- add_plot_label(Taj_plot, "A")  
Pi_plot_labeled <- add_plot_label(Pi_plot, "B")  


library(ggpubr)


# Arrange the plots using grid.arrange  

Div_plot <- ggarrange(Taj_plot_labeled, 
                      Pi_plot_labeled,
                      ncol=1, nrow=2, legend = FALSE)

Div_plot


png("Div_plot.png", width = 7000, height = 8000, res = 600)
Div_plot
dev.off()





# Plotting heatmap


R


library("tidyverse")
library(qqman)
library(viridis)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(grid)
library(circlize)


data <- read.table("transcript_Mangoum_RNAseq.ann.tsv", header = TRUE) 


# Remove "ANN...." prefix, "GEN." prefix, and ".FREQ" suffix  
colnames(data) <- gsub("^ANN\\.\\.\\.\\.|^GEN\\.|\\.\\.FREQ$", "", colnames(data))  

# Display the modified column names  
print(colnames(data)) 

# Step 1: Keep the first 10 columns as character  
data[1:10] <- lapply(data[1:10], as.character)  


# Step 2: Convert the rest of the columns to numeric  
numeric_columns <- colnames(data)[11:ncol(data)]  
data[numeric_columns] <- lapply(data[numeric_columns], function(x) {  
  x <- gsub("%", "", x)  # Remove percentage signs  
  as.numeric(ifelse(x == ".", NA, x))  # Convert to numeric, setting '.' to NA  
})  


# Step 3: Calculate Means for Each Group Based on the Basename  
mean_replicates_all_groups <- function(data) {  
  # Extract unique bases for grouping  
  group_basenames <- unique(sub("[0-9]+$", "", grep("[0-9]+$", colnames(data), value = TRUE)))  
  
  for (group in group_basenames) {  
    # Find columns that correspond to this group  
    cols_to_mean <- grep(paste0("^", group, "[0-9]+$"), colnames(data), value = TRUE)  
    
    if (length(cols_to_mean) > 0) {  
      # Calculate the mean across these columns and store it in a new column  
      data[[paste0(group, "_Mean")]] <- rowMeans(data[cols_to_mean], na.rm = TRUE)  # Compute mean  
    }  
  }  
  return(data)  
}  

# Step 4: Calculate Means for All Replicate Groups  
results <- mean_replicates_all_groups(data)  

results <- results %>% 
  distinct(CHROM,POS,HGVS_C, .keep_all = TRUE)


rownames(results) <- paste(results$CHROM,":",results$POS,":",results$HGVS_C,":",results$HGVS_P,":",results$GENEID) 
colnames (results)

results$SNP <- rownames(results)



# Read the text file into R  
file_path <- "groups.txt"  # Update the path to your text file  
group_data <- readLines(file_path)  # Read all lines from the file  

# Initialize an empty list to store the groups  
groups <- list()  

# Loop through each line to populate the list  
for (line in group_data) {  
  # Split the line into group name and members  
  parts <- strsplit(line, ":")[[1]]  
  
  # Check for the correct number of parts  
  if (length(parts) == 2) {  
    group_name <- parts[1]  # Group name  
    members <- unlist(strsplit(parts[2], ","))  # Members list  
    groups[[group_name]] <- members  # Add to the list  
  }  
}  


colnames(results) <- gsub("\\.merged$", "", colnames(results))
results$CHROM <- gsub("^AgamP4_", "", results$CHROM)
results$HGVS_C <- gsub("^c.", "", results$HGVS_C)
results$HGVS_P <- gsub("^p.", "", results$HGVS_P)
head(results, 2)


Annotation <- read_excel("Annotation.xlsx")

Annotation <- Annotation[!duplicated(Annotation$GeneID), ]

Annotation <- Annotation %>%
  select(GeneID, Gene_Name, Chromosome, Start, End, Description, Annotation)

Annotation <- Annotation %>%
  rename(CHROM = Chromosome, START = Start, END = End)

library(data.table)


# Convert to data.table by reference
setDT(results)
setDT(Annotation)


# Convert START and END to numeric in both tables
results[, `:=`(START = as.numeric(POS), END = as.numeric(POS))]
Annotation[, `:=`(START = as.numeric(START), END = as.numeric(END))]

# Set keys again (optional but safe)
setkey(Annotation, CHROM, START, END)
setkey(results, CHROM, START, END)

# Now perform the overlap join
joined <- foverlaps(results, Annotation, type = "within", nomatch = 0L)

combined_results <- joined

# Collapse GeneID per original row
#annotated <- joined[, .(
#  GeneID = paste(unique(GeneID), collapse = ";"),
#  Gene_Name = paste(unique(Gene_Name), collapse = ";"),
#  Description = paste(unique(Description), collapse = ";"),
#  Annotation = paste(unique(Annotation), collapse = ";")
#), by = .(CHROM, POS)]


# Merge back with original results
#combined_results <- merge(results, joined, by = c("CHROM", "POS"), all.x = TRUE)

#combined_results <- results%>%
#  left_join(Annotation, by = c('GENEID' = 'GeneID'))



combined_results$cSNP <- paste(combined_results$Gene_Name,"|",
                               combined_results$CHROM,"|",
                               combined_results$HGVS_C,"|", combined_results$HGVS_P)

combined_results <- combined_results[complete.cases(combined_results$EFFECT), ]

summary(combined_results)
head(combined_results,2)



# Read the text file into R  
#file_path <- "CHR_of_interest.txt"  # Update the path to your text file  
#CHR_of_interest <- readLines(file_path)  # Read all lines from the file 

CHR_of_interest <- c("2L", "2R", "3L", "3R", "X")

combined_results <- combined_results %>%
  filter(CHROM %in% CHR_of_interest)

head(combined_results,2)

# Mutate the CHROM column to numeric values based on their factor levels  
combined_results <- combined_results %>%  
  mutate(CHROM_NUM = as.numeric(factor(CHROM)))  

summary(combined_results)
head(combined_results,2)

file_path <- "GOI.txt"  # Update the path to your text file  
GOI <- readLines(file_path)  # Read all lines from the file 



Fisher <- read.delim("Mangoum_RNAseq_fisher_results_cleaned_with_new_header.txt", header = T, sep = "\t", dec = ".")

combined_results$POS <- as.numeric(combined_results$POS)

Fisher$POS <- as.numeric(Fisher$POS)


combined_results <- combined_results%>%
  left_join(Fisher, by = c("CHROM", "POS"))


# Columns to convert
cols_to_convert <- c("Kisumu_Susceptible.Mangoum_Permethrin10x",
                     "Kisumu_Susceptible.Mangoum_Permethrin1x",
                     "Kisumu_Susceptible.Mangoum_Permethrin5x",
                     "Kisumu_Susceptible.Mangoum_Unexposed",
                     "Kisumu_Susceptible.MK_Permethrin1x",
                     "Kisumu_Susceptible.MK_Unexposed",
                     "Mangoum_Permethrin10x.Mangoum_Permethrin1x",
                     "Mangoum_Permethrin10x.Mangoum_Permethrin5x",
                     "Mangoum_Permethrin10x.Mangoum_Unexposed",
                     "Mangoum_Permethrin10x.MK_Permethrin1x",
                     "Mangoum_Permethrin10x.MK_Unexposed",
                     "Mangoum_Permethrin1x.Mangoum_Permethrin5x",
                     "Mangoum_Permethrin1x.Mangoum_Unexposed",
                     "Mangoum_Permethrin1x.MK_Permethrin1x",
                     "Mangoum_Permethrin1x.MK_Unexposed",
                     "Mangoum_Permethrin5x.Mangoum_Unexposed",
                     "Mangoum_Permethrin5x.MK_Permethrin1x",
                     "Mangoum_Permethrin5x.MK_Unexposed",
                     "Mangoum_Unexposed.MK_Permethrin1x",
                     "Mangoum_Unexposed.MK_Unexposed",
                     "MK_Permethrin1x.MK_Unexposed"
)

# Convert to numeric, replacing "na" with NA
combined_results[, (cols_to_convert) := lapply(.SD, function(x) as.numeric(ifelse(x == "na", NA, x))), .SDcols = cols_to_convert]

head(combined_results,2)


#write.xlsx(combined_results, "combined_results.xlsx")

Misenses_variant <- combined_results %>%
  filter(EFFECT == "missense_variant")

head(Misenses_variant, 2)



combined_results_GOI <- combined_results%>%
  filter(Gene_Name %in% GOI)

head(combined_results_GOI, 2)

#test <- na.omit(combined_results_GOI)

#head(test)

# Get current column names
#col_names <- colnames(combined_results_GOI)

# Replace "CC_ALL" with "Fisher_"
#col_names <- gsub("CC_ALL", "Fisher_", col_names)


# Assign the modified names back to the data frame
#colnames(combined_results_GOI) <- col_names

## Save all snps

# 1. Remove columns with suffixes like 1, 2, ..., 9
cols_to_remove <- grep("[0-9]+$", colnames(combined_results_GOI))
cleaned_data <- combined_results_GOI[, -cols_to_remove]

# 2. Split data by 'Gene_Name'
groups <- split(cleaned_data, cleaned_data$Gene_Name)

# 3. Create a new Excel workbook and add each group as a sheet
wb <- createWorkbook()

for (gene in names(groups)) {
  addWorksheet(wb, sheetName = gene)
  writeData(wb, sheet = gene, x = groups[[gene]])
}

# 4. Save the workbook
saveWorkbook(wb, "Variant_Calling_GOI_selected.xlsx", overwrite = TRUE)



## Save all only missense_variant snps

# 1. Remove columns with suffixes like 1, 2, ..., 9
cols_to_remove <- grep("[0-9]+$", colnames(combined_results_GOI))
cleaned_data <- combined_results_GOI[, -cols_to_remove]

cleaned_data <- combined_results_GOI%>%
  filter(EFFECT == "missense_variant")

# 2. Split data by 'Gene_Name'
groups <- split(cleaned_data, cleaned_data$Gene_Name)

# 3. Create a new Excel workbook and add each group as a sheet
wb <- createWorkbook()

for (gene in names(groups)) {
  addWorksheet(wb, sheetName = gene)
  writeData(wb, sheet = gene, x = groups[[gene]])
}

# 4. Save the workbook
saveWorkbook(wb, "Variant_Calling_GOI_missense.xlsx", overwrite = TRUE)


## Save all only significant missense_variant snps

# 1. Remove columns with suffixes like 1, 2, ..., 9
cols_to_remove <- grep("[0-9]+$", colnames(combined_results_GOI))
cleaned_data <- combined_results_GOI[, -cols_to_remove]

cleaned_data <- combined_results_GOI%>%
  filter(EFFECT == "missense_variant" & 
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              if_any(starts_with("MK_Permethrin1x.MK_Unexposed"), ~ !is.na(.) & . > 1.30103)) |
           (Mangoum_Unexposed > Kisumu_Susceptible & Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("Mangoum_Permethrin1x.Mangoum_Unexposed"), ~ !is.na(.) & . > 1.30103)))

# 2. Split data by 'Gene_Name'
groups <- split(cleaned_data, cleaned_data$Gene_Name)

# 3. Create a new Excel workbook and add each group as a sheet
wb <- createWorkbook()

for (gene in names(groups)) {
  addWorksheet(wb, sheetName = gene)
  writeData(wb, sheet = gene, x = groups[[gene]])
}

# 4. Save the workbook
saveWorkbook(wb, "Variant_Calling_GOI_signif_missense.xlsx", overwrite = TRUE)



combined_results_GOI <- na.omit(combined_results_GOI)


## Heatmap of supporting reads :

Candidate_misenses_FREQ <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant")

Candidate_misenses_FREQ <- na.omit(Candidate_misenses_FREQ)
head(Candidate_misenses_FREQ, 1)

write.xlsx(Candidate_misenses_FREQ, "Candidate_misenses_FREQ.xlsx")



#Candidate_misenses_FREQ <- combined_results_GOI %>%
#  filter(EFFECT == "missense_variant")%>%
#  select(Chr, cSNP, all_of(ends_with("__Mean")))%>%
#  rename_with(
#    ~ gsub("__Mean$", "", .),
#    .cols = ends_with("__Mean")
#  )



Candidate_misenses_FREQ <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant")%>%
  select(Chr = CHROM, cSNP, Kisumu_Susceptible, Mangoum_Unexposed, 
         Mangoum_Permethrin1x, Mangoum_Permethrin5x,
         Mangoum_Permethrin10x, MK_Unexposed, 
         MK_Permethrin1x)


Candidate_misenses_FREQ <- na.omit(Candidate_misenses_FREQ)
head(Candidate_misenses_FREQ, 1)

write.xlsx(Candidate_misenses_FREQ, "Candidate_misenses_FREQ_heat.xlsx")



# Add a new column 'pos_number' by extracting the numeric position from 'HGVS_C'
combined_results_GOI <- combined_results_GOI %>%
  mutate(pos_number = as.numeric(str_extract(HGVS_C, "-?\\d+$")))

Candidate_upstream_FREQ <- combined_results_GOI %>%
  filter(EFFECT == "upstream_gene_variant" & (is.na(pos_number) | pos_number >= -1000))%>%
  select(Chr = CHROM, cSNP, Kisumu_Susceptible, Mangoum_Unexposed, 
         Mangoum_Permethrin1x, Mangoum_Permethrin5x,
         Mangoum_Permethrin10x, MK_Unexposed, 
         MK_Permethrin1x)%>%
  mutate(cSNP = gsub("\\| NA", "", cSNP))



Candidate_upstream_FREQ <- na.omit(Candidate_upstream_FREQ)
head(Candidate_upstream_FREQ, 4)

write.xlsx(Candidate_upstream_FREQ, "Candidate_upstream_FREQ_heat.xlsx")



###Melting datasets for misenses

Candidate_misenses_FREQ <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant" & 
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              if_any(starts_with("MK_Permethrin1x.MK_Unexposed"), ~ !is.na(.) & . > 1.30103)) |
              (Mangoum_Unexposed > Kisumu_Susceptible & Mangoum_Permethrin1x > Mangoum_Unexposed & 
                 if_any(starts_with("Kisumu_Susceptible.Mangoum_Permethrin1x"), ~ !is.na(.) & . > 1.30103)))%>%
  select(Chr = CHROM, cSNP, Kis.Sus = "Kisumu_Susceptible", Man.Unx = "Mangoum_Unexposed", 
         Man.Perm1x = "Mangoum_Permethrin1x", Man.Perm5x ="Mangoum_Permethrin5x",
         Man.Perm10x = "Mangoum_Permethrin10x", MK.Unx = "MK_Unexposed", 
         MK.Perm1x = "MK_Permethrin1x", Annotation, Description, 
         Man_CS = Kisumu_Susceptible.Mangoum_Unexposed, 
         Man_RC = Mangoum_Permethrin1x.Mangoum_Unexposed, 
         Man_RS = Kisumu_Susceptible.Mangoum_Permethrin1x, 
         MK_RC = MK_Permethrin1x.MK_Unexposed, 
         Man_10X.1X = Mangoum_Permethrin10x.Mangoum_Permethrin1x)

Candidate_misenses_FREQ <- na.omit(Candidate_misenses_FREQ)
head(Candidate_misenses_FREQ, 1)

write.xlsx(Candidate_misenses_FREQ, "Candidate_misenses_FREQ_heat_signif.xlsx")


names(Candidate_misenses_FREQ)
str(Candidate_misenses_FREQ)
summary(Candidate_misenses_FREQ)



Candidate_misenses <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant") %>%
  mutate(
    Man_CS = 10^-(Kisumu_Susceptible.Mangoum_Unexposed),
    Man_RC = 10^-(Mangoum_Permethrin1x.Mangoum_Unexposed),
    MK_CS = 10^-(Kisumu_Susceptible.MK_Unexposed),
    MK_RC = 10^-(MK_Permethrin1x.MK_Unexposed), 
    Man_10X.1X = 10^-(Mangoum_Permethrin10x.Mangoum_Permethrin1x)
  ) %>%
  select(Chr = "CHROM", cSNP, Kis.Sus = "Kisumu_Susceptible", Man.Unx = "Mangoum_Unexposed", 
         Man.Perm1x = "Mangoum_Permethrin1x", Man.Perm5x ="Mangoum_Permethrin5x",
         Man.Perm10x = "Mangoum_Permethrin10x", MK.Unx = "MK_Unexposed", 
         MK.Perm1x = "MK_Permethrin1x", Man_CS, Man_RC, MK_CS, MK_RC, Man_10X.1X, EFFECT, Annotation, Description)


write.xlsx(Candidate_misenses, "Candidate_misenses.xlsx")


# Melt the dataset 

library(reshape2) 

Candidate_misenses_FREQ_sub<-melt(data = Candidate_misenses_FREQ,  id=1:2, value.name = "Frequency", variable.name = "Population")


head(Candidate_misenses_FREQ_sub)



# Filter the dataset for chromosome  
Candidate_misenses_FREQ_GOI <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant") %>%
  mutate(
    Man_CS = 10^-(Kisumu_Susceptible.Mangoum_Unexposed),
    Man_RC = 10^-(Mangoum_Permethrin1x.Mangoum_Unexposed),
    Man_RS = 10^-(Kisumu_Susceptible.Mangoum_Permethrin1x),
    MK_RC = 10^-(MK_Permethrin1x.MK_Unexposed), 
    Man_10X.1X = 10^-(Mangoum_Permethrin10x.Mangoum_Permethrin1x)
  ) %>%
  filter(EFFECT == "missense_variant" & 
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("MK_RC"), ~ !is.na(.) & . < 0.05)) |
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("Man_RC"), ~ !is.na(.) & . < 0.05)))%>%
  select(Chr = "CHROM", cSNP, Kis.Sus = "Kisumu_Susceptible", Man.Unx = "Mangoum_Unexposed", 
         Man.Perm1x = "Mangoum_Permethrin1x", Man.Perm5x ="Mangoum_Permethrin5x",
         Man.Perm10x = "Mangoum_Permethrin10x", MK.Unx = "MK_Unexposed", 
         MK.Perm1x = "MK_Permethrin1x", Man_RS, Man_CS, Man_RC, MK_RC, Man_10X.1X, EFFECT, Annotation, Description)

write.xlsx(Candidate_misenses_FREQ_GOI, "Candidate_misenses_FREQ_GOI.xlsx")




### GOI

# Filter the dataset for chromosome  
Candidate_misenses_FREQ_GOI <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant") %>%
  mutate(
    Man_CS = 10^-(Kisumu_Susceptible.Mangoum_Unexposed),
    Man_RC = 10^-(Mangoum_Permethrin1x.Mangoum_Unexposed),
    Man_RS = 10^-(Kisumu_Susceptible.Mangoum_Permethrin1x),
    MK_RC = 10^-(MK_Permethrin1x.MK_Unexposed), 
    Man_10X.1X = 10^-(Mangoum_Permethrin10x.Mangoum_Permethrin1x)
  ) %>%
  filter(EFFECT == "missense_variant" & 
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("MK_RC"), ~ !is.na(.) & . < 0.05)) |
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("Man_RC"), ~ !is.na(.) & . < 0.05)))%>%
  select(Chr = "CHROM", cSNP, Kis.Sus = "Kisumu_Susceptible", Man.Unx = "Mangoum_Unexposed", 
         Man.Perm1x = "Mangoum_Permethrin1x", Man.Perm5x ="Mangoum_Permethrin5x",
         Man.Perm10x = "Mangoum_Permethrin10x", MK.Unx = "MK_Unexposed", 
         MK.Perm1x = "MK_Permethrin1x", Man_RS, Man_CS, Man_RC, MK_RC, Man_10X.1X, EFFECT, Annotation, Description)



# Melt frequency columns
freq_melt_GOI <- melt(
  Candidate_misenses_FREQ_GOI,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Kis.Sus", "Man.Unx", "Man.Perm1x", "Man.Perm5x",
                   "Man.Perm10x", "MK.Unx", "MK.Perm1x"),
  variable.name = "Population",
  value.name = "Frequency"
)

# Melt p-value columns
pval_melt_GOI <- melt(
  Candidate_misenses_FREQ_GOI,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Man_RS", "Man_CS", "Man_RC", "MK_RC", "Man_10X.1X"),
  variable.name = "Comparison",
  value.name = "Pvalue"
)

# Merge frequency and p-value data by Chr and cSNP
merged_data_GOI <- left_join(freq_melt_GOI, pval_melt_GOI, by = c("Chr", "cSNP"))


# Reshape Frequency and Pvalue into long format with correct group names
long_plot_data_GOI <- merged_data_GOI %>%
  select(cSNP, Population, Frequency, Comparison, Pvalue) %>%
  pivot_longer(
    cols = c(Frequency, Pvalue),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      Metric == "Frequency" ~ Population,
      Metric == "Pvalue" ~ Comparison
    )
  ) %>%
  select(cSNP, Metric, Group, Value)


plot_GOI <- ggplot(long_plot_data_GOI, aes(x = Group, y = cSNP, fill = Value)) +
  geom_tile(color = "black") +
  facet_wrap(~Metric, scales = "free_x") +
  scale_fill_viridis(
    option = "F",
    alpha = 1,
    direction = -1,
    begin = 0.45,
    end = 1,
    limits = c(0, 100)  # Set the limits of the fill scale to 0 and 100
  ) +
  geom_text(aes(label = round(Value, 2)), size = 5, fontface = "bold", colour = "black") +
  labs(
    title = "Gene of interest: Frequency and P-values",
    x = "",
    y = "Variant",
    fill = "Frequency (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold.italic"),
    legend.position = "top",
    axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")
  )


png(file="Candidate_misenses_FREQ_sub_GOI.png", width = 4500, height = 6000, pointsize = 30, res = 300)
plot_GOI 
dev.off()  



### VGSC

# Filter the dataset for chromosome  
Candidate_misenses_FREQ_VGSC <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant") %>%
  filter(CHROM == "2L", POS > 1967083, POS < 2714472) %>%
  mutate(
    Man_CS = 10^-(Kisumu_Susceptible.Mangoum_Unexposed),
    Man_RC = 10^-(Mangoum_Permethrin1x.Mangoum_Unexposed),
    Man_RS = 10^-(Kisumu_Susceptible.Mangoum_Permethrin1x),
    MK_RC = 10^-(MK_Permethrin1x.MK_Unexposed), 
    Man_10X.1X = 10^-(Mangoum_Permethrin10x.Mangoum_Permethrin1x)
  ) %>%
  filter(EFFECT == "missense_variant" & 
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("MK_RC"), ~ !is.na(.) & . < 0.05)) |
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("Man_RC"), ~ !is.na(.) & . < 0.05)))%>%
  select(Chr = "CHROM", cSNP, Kis.Sus = "Kisumu_Susceptible", Man.Unx = "Mangoum_Unexposed", 
         Man.Perm1x = "Mangoum_Permethrin1x", Man.Perm5x ="Mangoum_Permethrin5x",
         Man.Perm10x = "Mangoum_Permethrin10x", MK.Unx = "MK_Unexposed", 
         MK.Perm1x = "MK_Permethrin1x", Man_RS, Man_CS, Man_RC, MK_RC, Man_10X.1X, EFFECT, Annotation, Description)



# Melt frequency columns
freq_melt_VGSC <- melt(
  Candidate_misenses_FREQ_VGSC,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Kis.Sus", "Man.Unx", "Man.Perm1x", "Man.Perm5x",
                   "Man.Perm10x", "MK.Unx", "MK.Perm1x"),
  variable.name = "Population",
  value.name = "Frequency"
)

# Melt p-value columns
pval_melt_VGSC <- melt(
  Candidate_misenses_FREQ_VGSC,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Man_RS", "Man_CS", "Man_RC", "MK_RC", "Man_10X.1X"),
  variable.name = "Comparison",
  value.name = "Pvalue"
)

# Merge frequency and p-value data by Chr and cSNP
merged_data_VGSC <- left_join(freq_melt_VGSC, pval_melt_VGSC, by = c("Chr", "cSNP"))


# Reshape Frequency and Pvalue into long format with correct group names
long_plot_data_VGSC <- merged_data_VGSC %>%
  select(cSNP, Population, Frequency, Comparison, Pvalue) %>%
  pivot_longer(
    cols = c(Frequency, Pvalue),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      Metric == "Frequency" ~ Population,
      Metric == "Pvalue" ~ Comparison
    )
  ) %>%
  select(cSNP, Metric, Group, Value)


plot_VGSC <- ggplot(long_plot_data_VGSC, aes(x = Group, y = cSNP, fill = Value)) +
  geom_tile(color = "black") +
  facet_wrap(~Metric, scales = "free_x") +
  scale_fill_viridis(
    option = "F",
    alpha = 1,
    direction = -1,
    begin = 0.45,
    end = 1,
    limits = c(0, 100)  # Set the limits of the fill scale to 0 and 100
  ) +
  geom_text(aes(label = round(Value, 2)), size = 5, fontface = "bold", colour = "black") +
  labs(
    title = "A. VGSC: Frequency and P-values",
    x = "Group",
    y = "Variant",
    fill = "Frequency (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold.italic"),
    legend.position = "top",
    axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")
  )


png(file="Candidate_misenses_FREQ_sub_VGSC.png", width = 4500, height = 6000, pointsize = 30, res = 300)
plot_VGSC 
dev.off()  


### Others

# Filter the dataset for chromosome  
Candidate_misenses_FREQ_Others <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant") %>%
  filter(!(CHROM == "2L" & POS > 1967083 & POS < 2714472)) %>%
  mutate(
    Man_CS = 10^-(Kisumu_Susceptible.Mangoum_Unexposed),
    Man_RC = 10^-(Mangoum_Permethrin1x.Mangoum_Unexposed),
    Man_RS = 10^-(Kisumu_Susceptible.Mangoum_Permethrin1x),
    MK_RC = 10^-(MK_Permethrin1x.MK_Unexposed), 
    Man_10X.1X = 10^-(Mangoum_Permethrin10x.Mangoum_Permethrin1x)
  ) %>%
  filter(EFFECT == "missense_variant" & 
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("MK_RC"), ~ !is.na(.) & . < 0.05)) |
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("Man_RC"), ~ !is.na(.) & . < 0.05)))%>%
  select(Chr = "CHROM", cSNP, Kis.Sus = "Kisumu_Susceptible", Man.Unx = "Mangoum_Unexposed", 
         Man.Perm1x = "Mangoum_Permethrin1x", Man.Perm5x ="Mangoum_Permethrin5x",
         Man.Perm10x = "Mangoum_Permethrin10x", MK.Unx = "MK_Unexposed", 
         MK.Perm1x = "MK_Permethrin1x", Man_RS, Man_CS, Man_RC, MK_RC, Man_10X.1X, EFFECT, Annotation, Description)



# Melt frequency columns
freq_melt_Others <- melt(
  Candidate_misenses_FREQ_Others,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Kis.Sus", "Man.Unx", "Man.Perm1x", "Man.Perm5x",
                   "Man.Perm10x", "MK.Unx", "MK.Perm1x"),
  variable.name = "Population",
  value.name = "Frequency"
)

# Melt p-value columns
pval_melt_Others <- melt(
  Candidate_misenses_FREQ_Others,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Man_RS", "Man_CS", "Man_RC", "MK_RC", "Man_10X.1X"),
  variable.name = "Comparison",
  value.name = "Pvalue"
)

# Merge frequency and p-value data by Chr and cSNP
merged_data_Others <- left_join(freq_melt_Others, pval_melt_Others, by = c("Chr", "cSNP"))


# Reshape Frequency and Pvalue into long format with correct group names
long_plot_data_Others <- merged_data_Others %>%
  select(cSNP, Population, Frequency, Comparison, Pvalue) %>%
  pivot_longer(
    cols = c(Frequency, Pvalue),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      Metric == "Frequency" ~ Population,
      Metric == "Pvalue" ~ Comparison
    )
  ) %>%
  select(cSNP, Metric, Group, Value)


plot_Others <- ggplot(long_plot_data_Others, aes(x = Group, y = cSNP, fill = Value)) +
  geom_tile(color = "black") +
  facet_wrap(~Metric, scales = "free_x") +
  scale_fill_viridis(
    option = "F",
    alpha = 1,
    direction = -1,
    begin = 0.45,
    end = 1,
    limits = c(0, 100)  # Set the limits of the fill scale to 0 and 100
  ) +
  geom_text(aes(label = round(Value, 2)), size = 5, fontface = "bold", colour = "black") +
  labs(
    title = "B. Others candidates: Frequency and P-values",
    x = "Group",
    y = "Variant",
    fill = "Frequency (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold.italic"),
    legend.position = "top",
    axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")
  )


png(file="Candidate_misenses_FREQ_sub_Others.png", width = 4500, height = 6000, pointsize = 30, res = 300)
plot_Others 
dev.off()  


### CytochromeP450

# Filter the dataset for chromosome  
Candidate_misenses_FREQ_CytochromeP450 <- combined_results_GOI %>%
  filter(EFFECT == "missense_variant") %>%
  filter(!(CHROM == "2L" & POS > 1967083 & POS < 2714472) & Annotation == "Cytochrome P450") %>%
  mutate(
    Man_CS = 10^-(Kisumu_Susceptible.Mangoum_Unexposed),
    Man_RC = 10^-(Mangoum_Permethrin1x.Mangoum_Unexposed),
    Man_RS = 10^-(Kisumu_Susceptible.Mangoum_Permethrin1x),
    MK_RC = 10^-(MK_Permethrin1x.MK_Unexposed), 
    Man_10X.1X = 10^-(Mangoum_Permethrin10x.Mangoum_Permethrin1x)
  ) %>%
  filter(EFFECT == "missense_variant" & 
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("MK_RC"), ~ !is.na(.) & . < 0.05)) |
           (Mangoum_Unexposed > Kisumu_Susceptible & MK_Permethrin1x > MK_Unexposed & 
              Mangoum_Permethrin1x > Mangoum_Unexposed & 
              if_any(starts_with("Man_RC"), ~ !is.na(.) & . < 0.05)))%>%
  select(Chr = "CHROM", cSNP, Kis.Sus = "Kisumu_Susceptible", Man.Unx = "Mangoum_Unexposed", 
         Man.Perm1x = "Mangoum_Permethrin1x", Man.Perm5x ="Mangoum_Permethrin5x",
         Man.Perm10x = "Mangoum_Permethrin10x", MK.Unx = "MK_Unexposed", 
         MK.Perm1x = "MK_Permethrin1x", Man_RS, Man_CS, Man_RC, MK_RC, Man_10X.1X, EFFECT, Annotation, Description)



# Melt frequency columns
freq_melt_CytochromeP450 <- melt(
  Candidate_misenses_FREQ_CytochromeP450,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Kis.Sus", "Man.Unx", "Man.Perm1x", "Man.Perm5x",
                   "Man.Perm10x", "MK.Unx", "MK.Perm1x"),
  variable.name = "Population",
  value.name = "Frequency"
)

# Melt p-value columns
pval_melt_CytochromeP450 <- melt(
  Candidate_misenses_FREQ_CytochromeP450,
  id.vars = c("Chr", "cSNP"),
  measure.vars = c("Man_RS", "Man_CS", "Man_RC", "MK_RC", "Man_10X.1X"),
  variable.name = "Comparison",
  value.name = "Pvalue"
)

# Merge frequency and p-value data by Chr and cSNP
merged_data_CytochromeP450 <- left_join(freq_melt_CytochromeP450, pval_melt_CytochromeP450, by = c("Chr", "cSNP"))


# Reshape Frequency and Pvalue into long format with correct group names
long_plot_data_CytochromeP450 <- merged_data_CytochromeP450 %>%
  select(cSNP, Population, Frequency, Comparison, Pvalue) %>%
  pivot_longer(
    cols = c(Frequency, Pvalue),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      Metric == "Frequency" ~ Population,
      Metric == "Pvalue" ~ Comparison
    )
  ) %>%
  select(cSNP, Metric, Group, Value)


plot_CytochromeP450 <- ggplot(long_plot_data_CytochromeP450, aes(x = Group, y = cSNP, fill = Value)) +
  geom_tile(color = "black") +
  facet_wrap(~Metric, scales = "free_x") +
  scale_fill_viridis(
    option = "F",
    alpha = 1,
    direction = -1,
    begin = 0.45,
    end = 1,
    limits = c(0, 100)  # Set the limits of the fill scale to 0 and 100
  ) +
  geom_text(aes(label = round(Value, 2)), size = 5, fontface = "bold", colour = "black") +
  labs(
    title = "B. Cytochrome P450: Frequency and P-values",
    x = "Group",
    y = "Variant",
    fill = "Frequency (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold.italic"),
    legend.position = "top",
    axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")
  )


png(file="Candidate_misenses_FREQ_sub_CytochromeP450.png", width = 4500, height = 6000, pointsize = 30, res = 300)
plot_CytochromeP450 
dev.off()  


