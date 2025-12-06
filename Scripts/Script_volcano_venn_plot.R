# Loading library

library(tidyverse)
library(reshape2)
library(openxlsx)
library(readxl)
library(RColorBrewer)
library(ggrepel) 

# Data importation

setwd("D:/OneDrive - Centre for Research in Infectous Diseases/Analysis/Ngousso/")


data <- read_excel("Merged_DESeq2_Results.xlsx")  


# Read the comparisons from a text file without showing warnings  
comparisons_file <- "comparisons.txt"  # Specify your file path here  
specific_comparisons <- suppressWarnings(readLines(comparisons_file))


# Read the combinaisons for venndiagram  
combinaisons_file <- "combinaisons.txt"  # Specify your file path here  
specific_combinaisons <- suppressWarnings(readLines(combinaisons_file))

normalized_counts_averages <- read_excel("normalized_counts_averages.xlsx")

# Read and deduplicate annotation file  
Annotation <- read_excel("Annotation.xlsx")  
Annotation <- Annotation[!duplicated(Annotation$GeneID), ] 

str(data)
print(specific_combinaisons)
str(specific_comparisons)


## filter all Down regulated genes 

# Initialize an empty list to store created dataframes (optional, if needed)  
comparison_list <- list()  

# Create a new Excel workbook  
wb <- createWorkbook()  

# Step 1: Filter datasets for each comparison, join with annotation and counts, and save to sheets  
for (comp in specific_comparisons) {  
  # Select relevant columns  
  selected_cols <- c("GeneID", grep(paste0(comp, "$"), names(data), value = TRUE))  
  selected_data <- data %>%  
    select(all_of(selected_cols))  
  
  # Determine regulation column  
  if (grepl("vs.Ng.Sus$", comp)) {  
    reg_col <- paste0("regulation_2_", comp)  
  } else {  
    reg_col <- paste0("regulation_1.5_", comp)  
  }  
  
  # Filter for "Down" regulation  
  if (reg_col %in% names(data)) {  
    filtered_data <- selected_data %>%  
      filter(.data[[reg_col]] == "Down")  
  } else {  
    message("Column not found for: ", comp)  
    next  
  }  
  
  # Join with normalized counts and annotation  
  joined_data <- filtered_data %>%  
    left_join(normalized_counts_averages, by = c("GeneID" = "Gene")) %>%  
    left_join(Annotation, by = "GeneID")  
  
  # Store the dataframe with a variable named after the comparison  
  assign(comp, joined_data, envir = .GlobalEnv)  
  
  # Write to Excel worksheet  
  sheet_name <- substr(comp, 1, 31)  
  addWorksheet(wb, sheetName = sheet_name)  
  writeData(wb, sheet = sheet_name, joined_data)  
}  

# Step 2: Merge datasets based on combinations, clean, and save to sheets  
for (comb in specific_combinaisons) {  
  # Split combination string into dataset names  
  datasets_names <- strsplit(comb, "_")[[1]]  
  
  # Retrieve existing datasets  
  datasets_list <- lapply(datasets_names, function(name) {  
    if (exists(name, envir = .GlobalEnv)) {  
      get(name, envir = .GlobalEnv)  
    } else {  
      warning("Dataset ", name, " does not exist.")  
      NULL  
    }  
  })  
  
  # Remove NULL entries  
  datasets_list <- Filter(Negate(is.null), datasets_list)  
  
  # Merge datasets  
  if (length(datasets_list) > 1) {  
    merged_data <- Reduce(function(x, y) full_join(x, y, by = "GeneID"), datasets_list)  
  } else if (length(datasets_list) == 1) {  
    merged_data <- datasets_list[[1]]  
  } else {  
    warning("No datasets to merge for: ", comb)  
    next  
  }  
  
  # Append annotation columns (if not included)  
  # (Optional: perform join if needed; here, assuming already included or not necessary)  
  
  # Clean duplicate columns: keep only one version, preferring `.x`  
  remove_duplicate_columns <- function(df) {
    col_names <- colnames(df)
    base_names <- gsub("(\\.x|\\.y)$", "", col_names)
    keep_col <- rep(FALSE, length(col_names))
    
    for (base in unique(base_names)) {
      matches <- grep(paste0("^", base, "(\\.x|\\.y)?$"), col_names)
      col_x_match <- grep(paste0("^", base, "\\.x$"), col_names[matches])
      
      if (length(col_x_match) > 0) {
        # Keep the .x version
        keep_idx <- matches[col_x_match[1]]
      } else {
        # Keep the first match (which might be the version without suffix)
        keep_idx <- matches[1]
      }
      keep_col[keep_idx] <- TRUE
    }
    df_clean <- df[, keep_col, drop = FALSE]
    return(df_clean)
  } 
  
  # Apply column cleanup  
  cleaned_data <- remove_duplicate_columns(merged_data)  
  
  # Rename columns to remove only the suffix .x or .y
  names(cleaned_data) <- gsub("(\\.x|\\.y)$", "", names(cleaned_data))
  
  # Save to Excel sheet
  sheet_name <- substr(comb, 1, 31)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, cleaned_data)
}

# Save the workbook after all sheets are added
saveWorkbook(wb, "Analysis_Results_Down.xlsx", overwrite = TRUE)



## filter all Up regulated genes 

# Initialize an empty list to store created dataframes (optional, if needed)  
comparison_list <- list()  

# Create a new Excel workbook  
wb <- createWorkbook()  

# Step 1: Filter datasets for each comparison, join with annotation and counts, and save to sheets  
for (comp in specific_comparisons) {  
  # Select relevant columns  
  selected_cols <- c("GeneID", grep(paste0(comp, "$"), names(data), value = TRUE))  
  selected_data <- data %>%  
    select(all_of(selected_cols))  
  
  # Determine regulation column  
  if (grepl("vs.Ng.Sus$", comp)) {  
    reg_col <- paste0("regulation_2_", comp)  
  } else {  
    reg_col <- paste0("regulation_1.5_", comp)  
  }  
  
  # Filter for "Down" regulation  
  if (reg_col %in% names(data)) {  
    filtered_data <- selected_data %>%  
      filter(.data[[reg_col]] == "Up")  
  } else {  
    message("Column not found for: ", comp)  
    next  
  }  
  
  # Join with normalized counts and annotation  
  joined_data <- filtered_data %>%  
    left_join(normalized_counts_averages, by = c("GeneID" = "Gene")) %>%  
    left_join(Annotation, by = "GeneID")  
  
  # Store the dataframe with a variable named after the comparison  
  assign(comp, joined_data, envir = .GlobalEnv)  
  
  # Write to Excel worksheet  
  sheet_name <- substr(comp, 1, 31)  
  addWorksheet(wb, sheetName = sheet_name)  
  writeData(wb, sheet = sheet_name, joined_data)  
}  

# Step 2: Merge datasets based on combinations, clean, and save to sheets  
for (comb in specific_combinaisons) {  
  # Split combination string into dataset names  
  datasets_names <- strsplit(comb, "_")[[1]]  
  
  # Retrieve existing datasets  
  datasets_list <- lapply(datasets_names, function(name) {  
    if (exists(name, envir = .GlobalEnv)) {  
      get(name, envir = .GlobalEnv)  
    } else {  
      warning("Dataset ", name, " does not exist.")  
      NULL  
    }  
  })  
  
  # Remove NULL entries  
  datasets_list <- Filter(Negate(is.null), datasets_list)  
  
  # Merge datasets  
  if (length(datasets_list) > 1) {  
    merged_data <- Reduce(function(x, y) full_join(x, y, by = "GeneID"), datasets_list)  
  } else if (length(datasets_list) == 1) {  
    merged_data <- datasets_list[[1]]  
  } else {  
    warning("No datasets to merge for: ", comb)  
    next  
  }  
  
  # Append annotation columns (if not included)  
  # (Optional: perform join if needed; here, assuming already included or not necessary)  
  
  # Clean duplicate columns: keep only one version, preferring `.x`  
  remove_duplicate_columns <- function(df) {
    col_names <- colnames(df)
    base_names <- gsub("(\\.x|\\.y)$", "", col_names)
    keep_col <- rep(FALSE, length(col_names))
    
    for (base in unique(base_names)) {
      matches <- grep(paste0("^", base, "(\\.x|\\.y)?$"), col_names)
      col_x_match <- grep(paste0("^", base, "\\.x$"), col_names[matches])
      
      if (length(col_x_match) > 0) {
        # Keep the .x version
        keep_idx <- matches[col_x_match[1]]
      } else {
        # Keep the first match (which might be the version without suffix)
        keep_idx <- matches[1]
      }
      keep_col[keep_idx] <- TRUE
    }
    df_clean <- df[, keep_col, drop = FALSE]
    return(df_clean)
  } 
  
  # Apply column cleanup  
  cleaned_data <- remove_duplicate_columns(merged_data)  
  
  # Rename columns to remove only the suffix .x or .y
  names(cleaned_data) <- gsub("(\\.x|\\.y)$", "", names(cleaned_data))
  
  # Save to Excel sheet
  sheet_name <- substr(comb, 1, 31)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, cleaned_data)
}

# Save the workbook after all sheets are added
saveWorkbook(wb, "Analysis_Results_Up.xlsx", overwrite = TRUE)



# Volcano plot Field




# Volcano plot for Ng.Perm1x.vs.Ng.Sus

res_Ng.Perm1x.vs.Ng.Sus <- data %>%
  select(GeneID, Gene_Name, "log2FoldChange" = log2FoldChange_Ng.Perm1x.vs.Ng.Sus, "FoldChange" = FC_Ng.Perm1x.vs.Ng.Sus,
         "padj" = padj_Ng.Perm1x.vs.Ng.Sus, "gene_type" = regulation_2_Ng.Perm1x.vs.Ng.Sus, 
         Description, Product_Description, Function, Annotation)


# Obtain gene_type counts ------------------------------------------------------           
#res_Ng.Perm1x.vs.Ng.Sus %>%
#  count(gene_type)

# Modify legend labels by re-ordering gene_type levels -------------------------
res_Ng.Perm1x.vs.Ng.Sus <- res_Ng.Perm1x.vs.Ng.Sus %>%
  mutate(gene_type = fct_relevel(gene_type, "Up", "Down")) 


# Import Gene of interest table -----------------------------------
GOI_annotation <- read_excel("GOI.xlsx") 

# Get the unique values in the 'Code' column  
unique_codes <- unique(GOI_annotation$Code)  

# Determine the number of unique codes  
num_codes <- length(unique_codes)  

# Generate a color palette using RColorBrewer.  Adjust the palette name if needed.  

#my_colors <- c('cyan4', 'darkred', 'darkblue', 
#               'purple', 'darkgreen', 'magenta', 'orange', 
#               'yellow', 'green', 'black')#brewer.pal(n = num_codes, name = "Set3")

my_colors <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00"#, # orange
  #"black", "gold1",
  #"skyblue2", "#FB9A99", # lt pink
  #"palegreen2",
  #"#CAB2D6", # lt purple
  #"#FDBF6F", # lt orange
  #"khaki2",
  #"maroon"#, "orchid1", "steelblue4", "deeppink1", "blue1",
  #"darkturquoise", "green1", "yellow4", "yellow3",
  #"darkorange4", "brown"
)
# Create a named vector where names are the unique codes and values are the colors  
names(my_colors) <- unique_codes  

GOI_annotation$color <- my_colors[GOI_annotation$Code]

GOI_Ng.Perm1x.vs.Ng.Sus <- res_Ng.Perm1x.vs.Ng.Sus %>%
  filter(gene_type == "Up" | gene_type == "Down") %>%
  left_join(GOI_annotation, by = c('Annotation' = 'Annotation'))

GOI_Ng.Perm1x.vs.Ng.Sus <- GOI_Ng.Perm1x.vs.Ng.Sus %>% 
  drop_na(Code)

GOI_Ng.Perm1x.vs.Ng.Sus$Code <- factor(GOI_Ng.Perm1x.vs.Ng.Sus$Code) 

GOI_Ng.Perm1x.vs.Ng.Sus <- GOI_Ng.Perm1x.vs.Ng.Sus[order(GOI_Ng.Perm1x.vs.Ng.Sus$FoldChange, decreasing=TRUE),]


write.xlsx(GOI_Ng.Perm1x.vs.Ng.Sus, "GOI_Ng.Perm1x.vs.Ng.Sus.xlsx")

# Filter selected candidates  

# Filter top 10 Up-regulated genes
top_up_Ng.Perm1x.vs.Ng.Sus <- GOI_Ng.Perm1x.vs.Ng.Sus %>%
  filter(gene_type == "Up") %>% 
  slice_head(n = 10)

# Filter top 10 Down-regulated genes
top_down_Ng.Perm1x.vs.Ng.Sus <- GOI_Ng.Perm1x.vs.Ng.Sus %>%
  filter(gene_type == "Down") %>%
  slice_tail(n = 10)

# Combine the results if needed
selected_candidates <- bind_rows(top_up_Ng.Perm1x.vs.Ng.Sus, top_down_Ng.Perm1x.vs.Ng.Sus)


volcano_GOI_Ng.Perm1x.vs.Ng.Sus <- ggplot(data = res_Ng.Perm1x.vs.Ng.Sus,
                                          aes(x = log2FoldChange,
                                              y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.5)  +
  geom_point(data = GOI_Ng.Perm1x.vs.Ng.Sus,
             aes(fill = Code),
             shape = 21,
             size = 4,
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-log2(2), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = GOI_Ng.Perm1x.vs.Ng.Sus$color) + 
  scale_y_continuous(breaks = c(seq(0, 200, 50)),     
                     limits = c(0, 200)) +
  scale_x_continuous(breaks = c(seq(-12, 12, 6)),     
                     limits = c(-12, 12)) +
  labs(title = "Volcano Plot Res vs Sus \nFDR < 0.05 and absolute FC >= 2",
       x = "log2(fold change)",
       y = "-log10(FDR)",
       fill= "") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_label_repel(data          = subset(selected_candidates, log2FoldChange > 0),
                   aes(label = Gene_Name),
                   nudge_x       = 12 - subset(selected_candidates, log2FoldChange > 0)$log2FoldChange,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y",
                   hjust         = 0) +
  geom_label_repel(data          = subset(selected_candidates, log2FoldChange < 0),
                   aes(label = Gene_Name),
                   nudge_x       = -12 - subset(selected_candidates, log2FoldChange < 0)$log2FoldChange,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y",
                   hjust         = 1)

png(file="volcano_GOI_Ng.Perm1x.vs.Ng.Sus.png", width = 3000, height = 4000, res = 400)  
volcano_GOI_Ng.Perm1x.vs.Ng.Sus
dev.off()  





# Volcano plot for Ng.Unx.vs.Ng.Sus

res_Ng.Unx.vs.Ng.Sus <- data %>%
  select(GeneID, Gene_Name, "log2FoldChange" = log2FoldChange_Ng.Unx.vs.Ng.Sus, "FoldChange" = FC_Ng.Unx.vs.Ng.Sus,
         "padj" = padj_Ng.Unx.vs.Ng.Sus, "gene_type" = regulation_2_Ng.Unx.vs.Ng.Sus, 
         Description, Product_Description, Function, Annotation)


# Obtain gene_type counts ------------------------------------------------------           
#res_Ng.Unx.vs.Ng.Sus %>%
#  count(gene_type)

# Modify legend labels by re-ordering gene_type levels -------------------------
res_Ng.Unx.vs.Ng.Sus <- res_Ng.Unx.vs.Ng.Sus %>%
  mutate(gene_type = fct_relevel(gene_type, "Up", "Down")) 


# Import Gene of interest table -----------------------------------
GOI_annotation <- read_excel("GOI.xlsx") 

# Get the unique values in the 'Code' column  
unique_codes <- unique(GOI_annotation$Code)  

# Determine the number of unique codes  
num_codes <- length(unique_codes)  

# Generate a color palette using RColorBrewer.  Adjust the palette name if needed.  

#my_colors <- c('cyan4', 'darkred', 'darkblue', 
#               'purple', 'darkgreen', 'magenta', 'orange', 
#               'yellow', 'green', 'black')#brewer.pal(n = num_codes, name = "Set3")

my_colors <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00"#, # orange
  #"black", "gold1",
  #"skyblue2", "#FB9A99", # lt pink
  #"palegreen2",
  #"#CAB2D6", # lt purple
  #"#FDBF6F", # lt orange
  #"khaki2",
  #"maroon"#, "orchid1", "steelblue4", "deeppink1", "blue1",
  #"darkturquoise", "green1", "yellow4", "yellow3",
  #"darkorange4", "brown"
)
# Create a named vector where names are the unique codes and values are the colors  
names(my_colors) <- unique_codes  

GOI_annotation$color <- my_colors[GOI_annotation$Code]

GOI_Ng.Unx.vs.Ng.Sus <- res_Ng.Unx.vs.Ng.Sus %>%
  filter(gene_type == "Up" | gene_type == "Down") %>%
  left_join(GOI_annotation, by = c('Annotation' = 'Annotation'))

GOI_Ng.Unx.vs.Ng.Sus <- GOI_Ng.Unx.vs.Ng.Sus %>% 
  drop_na(Code)

GOI_Ng.Unx.vs.Ng.Sus$Code <- factor(GOI_Ng.Unx.vs.Ng.Sus$Code) 

GOI_Ng.Unx.vs.Ng.Sus <- GOI_Ng.Unx.vs.Ng.Sus[order(GOI_Ng.Unx.vs.Ng.Sus$FoldChange, decreasing=TRUE),]


write.xlsx(GOI_Ng.Unx.vs.Ng.Sus, "GOI_Ng.Unx.vs.Ng.Sus.xlsx")

# Filter selected candidates  

# Filter top 10 Up-regulated genes
top_up_Ng.Unx.vs.Ng.Sus <- GOI_Ng.Unx.vs.Ng.Sus %>%
  filter(gene_type == "Up") %>% 
  slice_head(n = 10)

# Filter top 10 Down-regulated genes
top_down_Ng.Unx.vs.Ng.Sus <- GOI_Ng.Unx.vs.Ng.Sus %>%
  filter(gene_type == "Down") %>%
  slice_tail(n = 10)

# Combine the results if needed
selected_candidates <- bind_rows(top_up_Ng.Unx.vs.Ng.Sus, top_down_Ng.Unx.vs.Ng.Sus)


volcano_GOI_Ng.Unx.vs.Ng.Sus <- ggplot(data = res_Ng.Unx.vs.Ng.Sus,
                                       aes(x = log2FoldChange,
                                           y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.5)  +
  geom_point(data = GOI_Ng.Unx.vs.Ng.Sus,
             aes(fill = Code),
             shape = 21,
             size = 4,
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-log2(2), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = GOI_Ng.Unx.vs.Ng.Sus$color) + 
  scale_y_continuous(breaks = c(seq(0, 200, 50)),     
                     limits = c(0, 200)) +
  scale_x_continuous(breaks = c(seq(-12, 12, 6)),     
                     limits = c(-12, 12)) +
  labs(title = "Volcano Plot Cntr vs Sus \nFDR < 0.05 and absolute FC >= 2",
       x = "log2(fold change)",
       y = "-log10(FDR)",
       fill= "") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_label_repel(data          = subset(selected_candidates, log2FoldChange > 0),
                   aes(label = Gene_Name),
                   nudge_x       = 12 - subset(selected_candidates, log2FoldChange > 0)$log2FoldChange,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y",
                   hjust         = 0) +
  geom_label_repel(data          = subset(selected_candidates, log2FoldChange < 0),
                   aes(label = Gene_Name),
                   nudge_x       = -12 - subset(selected_candidates, log2FoldChange < 0)$log2FoldChange,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y",
                   hjust         = 1)

png(file="volcano_GOI_Ng.Unx.vs.Ng.Sus.png", width = 3000, height = 4000, res = 400)  
volcano_GOI_Ng.Unx.vs.Ng.Sus
dev.off()  







# Volcano plot for Ng.Perm1x.vs.Ng.Unx

res_Ng.Perm1x.vs.Ng.Unx <- data %>%
  select(GeneID, Gene_Name, "log2FoldChange" = log2FoldChange_Ng.Perm1x.vs.Ng.Unx, "FoldChange" = FC_Ng.Perm1x.vs.Ng.Unx,
         "padj" = padj_Ng.Perm1x.vs.Ng.Unx, "gene_type" = regulation_1.5_Ng.Perm1x.vs.Ng.Unx, 
         Description, Product_Description, Function, Annotation)


# Obtain gene_type counts ------------------------------------------------------           
#res_Ng.Perm1x.vs.Ng.Unx %>%
#  count(gene_type)

# Modify legend labels by re-ordering gene_type levels -------------------------
res_Ng.Perm1x.vs.Ng.Unx <- res_Ng.Perm1x.vs.Ng.Unx %>%
  mutate(gene_type = fct_relevel(gene_type, "Up", "Down")) 


# Import Gene of interest table -----------------------------------
GOI_annotation <- read_excel("GOI.xlsx") 

# Get the unique values in the 'Code' column  
unique_codes <- unique(GOI_annotation$Code)  

# Determine the number of unique codes  
num_codes <- length(unique_codes)  

# Generate a color palette using RColorBrewer.  Adjust the palette name if needed.  

#my_colors <- c('cyan4', 'darkred', 'darkblue', 
#               'purple', 'darkgreen', 'magenta', 'orange', 
#               'yellow', 'green', 'black')#brewer.pal(n = num_codes, name = "Set3")


my_colors <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00"#, # orange
  #"black", "gold1",
  #"skyblue2", "#FB9A99", # lt pink
  #"palegreen2",
  #"#CAB2D6", # lt purple
  #"#FDBF6F", # lt orange
  #"khaki2",
  #"maroon"#, "orchid1", "steelblue4", "deeppink1", "blue1",
  #"darkturquoise", "green1", "yellow4", "yellow3",
  #"darkorange4", "brown"
)
# Create a named vector where names are the unique codes and values are the colors  
names(my_colors) <- unique_codes  

GOI_annotation$color <- my_colors[GOI_annotation$Code]

GOI_Ng.Perm1x.vs.Ng.Unx <- res_Ng.Perm1x.vs.Ng.Unx %>%
  filter(gene_type == "Up" | gene_type == "Down") %>%
  left_join(GOI_annotation, by = c('Annotation' = 'Annotation'))

GOI_Ng.Perm1x.vs.Ng.Unx <- GOI_Ng.Perm1x.vs.Ng.Unx %>% 
  drop_na(Code)

GOI_Ng.Perm1x.vs.Ng.Unx$Code <- factor(GOI_Ng.Perm1x.vs.Ng.Unx$Code) 

GOI_Ng.Perm1x.vs.Ng.Unx <- GOI_Ng.Perm1x.vs.Ng.Unx[order(GOI_Ng.Perm1x.vs.Ng.Unx$FoldChange, decreasing=TRUE),]


write.xlsx(GOI_Ng.Perm1x.vs.Ng.Unx, "GOI_Ng.Perm1x.vs.Ng.Unx.xlsx")

# Filter selected candidates  

# Filter top 10 Up-regulated genes
top_up_Ng.Perm1x.vs.Ng.Unx <- GOI_Ng.Perm1x.vs.Ng.Unx %>%
  filter(gene_type == "Up") %>% 
  slice_head(n = 10)

# Filter top 10 Down-regulated genes
top_down_Ng.Perm1x.vs.Ng.Unx <- GOI_Ng.Perm1x.vs.Ng.Unx %>%
  filter(gene_type == "Down") %>%
  slice_tail(n = 10)

# Combine the results if needed
selected_candidates <- bind_rows(top_up_Ng.Perm1x.vs.Ng.Unx, top_down_Ng.Perm1x.vs.Ng.Unx)


volcano_GOI_Ng.Perm1x.vs.Ng.Unx <- ggplot(data = res_Ng.Perm1x.vs.Ng.Unx,
                                          aes(x = log2FoldChange,
                                              y = -log10(padj))) + 
  geom_point(colour = "grey", alpha = 0.5)  +
  geom_point(data = GOI_Ng.Perm1x.vs.Ng.Unx,
             aes(fill = Code),
             shape = 21,
             size = 4,
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
             linetype = "dashed") +
  scale_fill_manual(values = GOI_Ng.Perm1x.vs.Ng.Unx$color) + 
  scale_y_continuous(breaks = c(seq(0, 30, 5)),     
                     limits = c(0, 30)) +
  scale_x_continuous(breaks = c(seq(-12, 12, 6)),     
                     limits = c(-12, 12)) +
  labs(title = "Volcano Plot Res vs Cntr \nFDR < 0.05 and absolute FC >= 1.5",
       x = "log2(fold change)",
       y = "-log10(FDR)",
       fill= "") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_label_repel(data          = subset(selected_candidates, log2FoldChange > 0),
                   aes(label = Gene_Name),
                   nudge_x       = 12 - subset(selected_candidates, log2FoldChange > 0)$log2FoldChange,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y",
                   hjust         = 0) +
  geom_label_repel(data          = subset(selected_candidates, log2FoldChange < 0),
                   aes(label = Gene_Name),
                   nudge_x       = -12 - subset(selected_candidates, log2FoldChange < 0)$log2FoldChange,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y",
                   hjust         = 1)

png(file="volcano_GOI_Ng.Perm1x.vs.Ng.Unx.png", width = 3000, height = 4000, res = 400)  
volcano_GOI_Ng.Perm1x.vs.Ng.Unx
dev.off()  




library(ggpubr)


Volcano <- ggarrange(volcano_GOI_Ng.Perm1x.vs.Ng.Sus, 
                     volcano_GOI_Ng.Unx.vs.Ng.Sus, 
                     volcano_GOI_Ng.Perm1x.vs.Ng.Unx,  
                     nrow = 1, common.legend = TRUE, legend="bottom")

# Arrange the plots using grid.arrange  
png(file="DEG_Ngousso.png", width = 6000, height = 4000, res = 400)  
Volcano
dev.off()  





## Similar entities analysis

# Load necessary libraries
library(readxl)
library(dplyr)


# Read target genes from Excel file
target_genes <- GOI_Ng.Unx.vs.Ng.Sus %>%
  filter(gene_type == "Up") %>% 
  slice_head(n = 10)
# Assuming the column with gene IDs is named "GeneID" or similar
target_gene_list <- target_genes$GeneID

# Read expression data (CSV format)
expr_data <- read_excel("normalized_counts_averages.xlsx")

# Check if your expr_data is numeric
# If not, convert variables (columns) to numeric
expr_data <- as.data.frame(expr_data)

row.names(expr_data) <- expr_data$Gene
expr_data_numeric <- expr_data[ , -1]
expr_data_numeric <- data.frame(lapply(expr_data_numeric, as.numeric))
rownames(expr_data_numeric) <- expr_data$Gene



# Suppose target_gene_list contains gene IDs matching the rownames
results_list <- list()

for (gene in target_gene_list) {
  if (gene %in% rownames(expr_data_numeric)) {
    target_expr <- as.numeric(expr_data_numeric[gene, ])  # ensure numeric vector
    
    # Calculate correlations
    cor_scores <- apply(expr_data_numeric, 1, function(x) {
      x <- as.numeric(x)
      cor(x, target_expr, method = "pearson", use = "pairwise.complete.obs")
    })
    
    # Remove NA
    cor_scores <- na.omit(cor_scores)
    
    # Get top 20
    top_hits <- sort(cor_scores, decreasing = TRUE)[1:20]
    
    # Save results
    results_list[[gene]] <- data.frame(
      Target_Gene = gene,
      Similar_Gene = names(top_hits),
      Correlation = as.numeric(top_hits)
    )
  } else {
    warning(paste("Gene", gene, "not found in expression data"))
  }
}

# Combine all results into one data frame
all_results <- do.call(rbind, results_list) %>%
  left_join(normalized_counts_averages, by = c('Similar_Gene' = 'Gene'))%>%  
  left_join(Annotation, by = c('Similar_Gene' = 'GeneID'))

Similar_entities <- all_results %>%
  select(Target_Gene, "Gene ID"=Gene_Name, "Similarity"=Correlation, Ng.Perm1x, Ng.Unx, Ng.Sus, Annotation)

write.xlsx(Similar_entities, "Similar_entities_analysis.xlsx")

library(openxlsx)
write.xlsx(Similar_entities, "C:/Users/ARNAUD/OneDrive - Centre for Research in Infectous Diseases/Analysis/Ngousso/Similar_entities_analysis.xlsx")




# Venn diagram

library('VennDiagram')
library(limma)
library(gridExtra)

DEG_All <- data %>%
  select(GeneID, Gene_Name, 
         "RS" = regulation_2_Ng.Perm1x.vs.Ng.Sus , 
         "RC" = regulation_1.5_Ng.Perm1x.vs.Ng.Unx , 
         "CS" = regulation_2_Ng.Unx.vs.Ng.Sus)

tail(DEG_All, 10)

DEG_All <- DEG_All %>%  
  mutate(  
    RS = ifelse(RS == "Up", 1, ifelse(RS == "Down", -1, 0)),  
    RC = ifelse(RC == "Up", 1, ifelse(RC == "Down", -1, 0)),  
    CS = ifelse(CS == "Up", 1, ifelse(CS == "Down", -1, 0)))

DEG_All[is.na(DEG_All)] <- 0  


# All comparison
DEG_All_matrix<- data.matrix(DEG_All[,3:5])
head(DEG_All_matrix)

Venn_DEG_All <- grid::grobTree(
  grid.draw(vennDiagram(DEG_All_matrix, 
                        names = c("R-S", "R-C", "C-S"),
                        include=c("up", "down"),
                        cex=1.2,
                        counts.col=c("red", "blue"),
                        circle.col = c("#440154ff", '#21908dff', '#fde725ff')))
)



# Arrange the plots using grid.arrange  
png(file="Venn_DEG_All.png", width = 2000, height = 2000, res = 300)  
grid::grobTree(
  grid.draw(vennDiagram(DEG_All_matrix, 
                        names = c("R-S", "R-C", "C-S"),
                        include=c("up", "down"),
                        cex=1.2,
                        counts.col=c("red", "blue"),
                        circle.col = c("#440154ff", '#21908dff', '#fde725ff')))
)
dev.off()  

