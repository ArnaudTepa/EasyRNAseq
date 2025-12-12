## Define cores and work directory

options(mc.cores = 6L)
dir <- setwd("C:/Users/ARNAUD/OneDrive - Centre for Research in Infectous Diseases/Analysis/Ngousso/")


list.files(dir)



library(tximport)  
library(tximportData)  
library(rhdf5)  
library(DESeq2)  
library(apeglm)  
library(limma)  
library(pheatmap)  
library(edgeR)  
library(EnhancedVolcano)  
library(VennDiagram)  
library(readxl)  
library(tidyverse)  
library(ggplot2)  
library(reshape2)  
library(knitr)  
library(affy)  
library(openxlsx)  
library(ComplexHeatmap)  
library(grid)  
library(circlize)  
library(RColorBrewer)  
library(DGEobj.utils)  
library(chromoMap)  
library(wesanderson)  
library(viridisLite)  
library(ClusterR)   
library(cluster)   
library(stats)  
library(gplots)  
library(dendextend)  
library(factoextra)  
library(NbClust)  
library(venn)
library(ggdendro)

# Read and deduplicate annotation file  
Annotation <- read_excel(file.path(dir, "Annotation.xlsx"))  
Annotation <- Annotation[!duplicated(Annotation$GeneID), ] 

# Read the comparisons from a text file without showing warnings  
comparisons_file <- "comparisons.txt"  # Specify your file path here  
specific_comparisons <- suppressWarnings(readLines(comparisons_file))


# Read the combinaisons for venndiagram  
combinaisons_file <- "combinaisons.txt"  # Specify your file path here  
specific_combinaisons <- suppressWarnings(readLines(combinaisons_file))


# Load data  
cts <- read.table(file.path(dir, "counts.txt"), header = TRUE)  
coldata <- read.table(file.path(dir, "design.txt"), header = TRUE)  
coldata$condition <- factor(coldata$condition)  
coldata$Strain <- factor(coldata$Strain)  

# Check integrity of data  
cts <- cts[, rownames(coldata)]  
all(rownames(coldata) == colnames(cts))

# Descriptive statistics  
# Print dimensions of cts  
print(dim(cts))  

# Print summary of all columns dynamically  
print(head(summary(cts[, 1:ncol(cts)])))   

# Calculate statistics per sample  
stats.per.sample <- data.frame(t(do.call(cbind, lapply(cts, summary))))  
print(head(stats.per.sample))  

# Automatically generate color codes based on unique conditions  
unique_conditions <- unique(coldata$condition)  

# Generate a color palette for the number of unique conditions  
condition_colors <- brewer.pal(n = length(unique_conditions), name = "Set1")  

# Create a named vector to assign colors to conditions  
named_colors <- setNames(condition_colors, unique_conditions)  

# Assign colors to the coldata dataframe  
coldata$color <- named_colors[coldata$condition]  

# Display the updated coldata  
print(coldata)  

  
epsilon <- 1 # pseudo-count to avoid problems with log(0)  

# Define function for plotting  
plot_data <- function(data, color_col, file_name, plot_type, main_title) {  
  png(file = file_name, width = 2000, height = 800, res = 200)  
  if (plot_type == "boxplot") {  
    boxplot(log2(data + epsilon), col = color_col, pch = ".",   
            horizontal = FALSE, cex.axis = 0.5, las = 2,  
            xlab = "Samples", ylab = "log2(Counts +1)", main = main_title)  
  } else if (plot_type == "density") {  
    plotDensity(log2(data + epsilon), lty = 1, col = color_col, lwd = 2,  
                ylab="Density", xlab="log2(Counts +1)", main = main_title)  
    grid()  
    legend("topright", legend = names(named_colors), col = named_colors, lwd = 2)  
  }  
  dev.off()  
}  

# Boxplots and Density plots for raw counts  
plot_data(cts, coldata$color, "Boxplot.raw.count.png", "boxplot", "Raw counts boxplot")  
plot_data(cts, coldata$color, "Densityplot.raw.count.png", "density", "Raw counts density plot")  



# DESeq2 processing  
dds_raw <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)  
dds_raw <- estimateSizeFactors(dds_raw)  
normalized_counts <- counts(dds_raw, normalized=TRUE)  

# Boxplots and Density plots for normalized counts  
plot_data(normalized_counts, coldata$color, "Boxplot.Normalized.count.png", "boxplot", "Normalized counts boxplot")  
plot_data(normalized_counts, coldata$color, "Densityplot.Normalized.count.png", "density", "Normalized counts density plot")  

## Estimate average read counts
# Transpose the normalized counts to have samples in rows  
normalized_counts_t <- as.data.frame(t(normalized_counts))  

# Combine the conditions from coldata with the transposed normalized counts  
# Remove both the "Strain" and "color" columns  
combined_data <- cbind(condition = coldata[,-c(1, 3)], normalized_counts_t) # Removing columns 1 (Strain) and 3 (color)  


# Calculate the average for each condition  
averages <- combined_data %>%  
  group_by(condition) %>%  
  summarise(across(everything(), mean, na.rm = TRUE))  

# Transpose the averages data frame  
normalized_counts_averages <- as.data.frame(t(averages))  

# Set the first row as the column names for better readability  
# This step takes the first row and sets it as the columns  
colnames(normalized_counts_averages) <- as.character(normalized_counts_averages[1, ])  
# Remove the first row now that it's set as column names  
normalized_counts_averages <- normalized_counts_averages[-1, ]  


# Add the row names (gene names) as the first column  
normalized_counts_averages <- cbind(Gene = rownames(normalized_counts_averages), normalized_counts_averages)  


# Save the transposed averages to an Excel file  
write.xlsx(normalized_counts_averages, "normalized_counts_averages.xlsx")  

# Optionally print a message confirming the save  
print("The averages have been saved to 'normalized_counts_averages.xlsx'")  


 

# PCA analysis 
dds <- DESeqDataSetFromMatrix(countData = round(normalized_counts), colData = coldata, design = ~ condition)  
vsd <- vst(dds, blind = FALSE)  
png(file = "PCA.plot.png", width = 4000, height = 4000, pointsize = 20, res = 400)  
pca_plot <- plotPCA(vsd, ntop = nrow(cts), intgroup = "condition") +    
  theme_classic() +   
  theme(  
    legend.position = "right",  
    legend.text = element_text(size = 15),  
    axis.text = element_text(size = 15),  
    legend.title = element_text(size = 18, face = "bold"),  # This line sets the legend title style  
    title = element_text(size = 18, face = "bold")  
  ) +    
  labs(color = "Phenotype") +  # This sets the legend title to "Phenotype"  
  scale_color_manual(values = named_colors)  # Define the colors based on the color column in coldata  
print(pca_plot)
dev.off()  

# Sample distance heatmap

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$Strain, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(file="sample_to_sample_distances.plot.png", width = 1600, height = 1200, pointsize = 12, res = 144)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# Calculate hierarchical clustering 
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- NULL


# Calculate hierarchical clustering  
hc <- hclust(sampleDists)  


# Plotting the dendrogram  
png(file="sample_Dendrogram.plot.png", width = 4000, height = 2500, pointsize = 20, res = 300)
plot(hc, labels = rownames(sampleDistMatrix), main = "Sample Dendrogram", xlab = "Samples", ylab = "Distance", hang = -1)
dev.off()


# DESeq2 analysis for DEG   


# Create DESeqDataSet  
dds <- DESeqDataSetFromMatrix(countData = round(normalized_counts), colData = coldata, design = ~ condition)  

# Run the DESeq analysis initially for the dataset  
dds <- DESeq(dds)  


# Initialize a list to store data frames for each comparison result  
result_list <- list()  
gene_lists <- list()  # A new list to store the up and down genes for Venn diagrams  

# Print current levels to debug  
print("Current levels of coldata$condition:")  
print(levels(coldata$condition))  

# Loop through the specified comparisons directly from the text file  
for (contrast_name in specific_comparisons) {  
  comp <- strsplit(contrast_name, ".vs.")[[1]]  
  print(paste("Checking comparison:", contrast_name))  # Debug print  
  
  # Validate the comparison against the results names  
  if (length(comp) != 2 || !(comp[1] %in% levels(coldata$condition)) || !(comp[2] %in% levels(coldata$condition))) {  
    warning(paste("Invalid comparison specified:", contrast_name))  
    next  # Skip to the next comparison  
  }  
  
  # Conduct the differential analysis  
  res <- results(dds, contrast = c("condition", comp[1], comp[2]))  
  
  # Convert results to a data frame  
  res_df <- as.data.frame(res)  
  res_df <- cbind(GeneID = rownames(res_df), res_df)  
  
  # Calculate FC: 2 raised to the power of the absolute value of log2FoldChange  
  res_df$FC <- ifelse(res_df$log2FoldChange < 0,  
                      -2^abs(res_df$log2FoldChange),  
                      2^res_df$log2FoldChange)  
  
  # Classify based on cutoffs for regulation_2  
  res_df$regulation_2 <- "Not Significant"  
  res_df$regulation_2[res_df$FC > 2 & res_df$padj < 0.05] <- "Up"  
  res_df$regulation_2[res_df$FC < -2 & res_df$padj < 0.05] <- "Down"  
  
  # Classify based on cutoffs for regulation_1.5  
  res_df$regulation_1.5 <- "Not Significant"  
  res_df$regulation_1.5[res_df$FC > 1.5 & res_df$padj < 0.05] <- "Up"  
  res_df$regulation_1.5[res_df$FC < -1.5 & res_df$padj < 0.05] <- "Down"  
  
  # Determine which regulation cutoff to use for Venn diagrams  
  if (grepl("Kis.Sus", contrast_name)) {  
    regulation_type <- "regulation_2"  
  } else {  
    regulation_type <- "regulation_1.5"  
  }  
  
  # Update column names to include the comparison name  
  colnames(res_df) <- gsub("baseMean", paste0("baseMean_", contrast_name), colnames(res_df))  
  colnames(res_df) <- gsub("log2FoldChange", paste0("log2FoldChange_", contrast_name), colnames(res_df))  
  colnames(res_df) <- gsub("lfcSE", paste0("lfcSE_", contrast_name), colnames(res_df))  
  colnames(res_df) <- gsub("stat", paste0("stat_", contrast_name), colnames(res_df))  
  colnames(res_df) <- gsub("pvalue", paste0("pvalue_", contrast_name), colnames(res_df))  
  colnames(res_df) <- gsub("padj", paste0("padj_", contrast_name), colnames(res_df))  
  
  # Rename combined additional columns to include comparison name  
  colnames(res_df)[which(colnames(res_df) == "FC")] <- paste0("FC_", contrast_name)  
  colnames(res_df)[which(colnames(res_df) == "regulation_2")] <- paste0("regulation_2_", contrast_name)  
  colnames(res_df)[which(colnames(res_df) == "regulation_1.5")] <- paste0("regulation_1.5_", contrast_name)  
  
  # Add result data frame to the list  
  result_list[[contrast_name]] <- res_df  # Storing results by contrast name  
  
  # Save each DESeq2 result to the environment with the name `deseq2_<contrast_name>`  
  assign(paste0(contrast_name), res_df, envir = .GlobalEnv)  
  
  # Create an output filename based on the comparison  
  output_filename <- paste0("DESeq2_results_", contrast_name, ".xlsx")  
  
  # Save the results to the xlsx file  
  write.xlsx(res_df, file = output_filename, rowNames = FALSE)  
  
  # Print the comparison and the output filename  
  print(paste("Comparison:", contrast_name, "Output saved to:", output_filename))  
}  


# Initialize an empty data frame to hold combined results  
combined_results <- NULL  

# Loop through each result data frame in the result list  
for (name in names(result_list)) {  
  df <- result_list[[name]]  
  
  # Perform the join, leaving the first one as the base  
  if (is.null(combined_results)) {  
    combined_results <- df  # Initialize combined_results with the first data frame  
  } else {  
    combined_results <- left_join(combined_results, df, by = "GeneID", suffix = c("", paste0("_", name)))  
  }  
}  


# Join with Annotation and Normalized Counts  
combined_results <- left_join(combined_results, normalized_counts_averages, by = c("GeneID" = "Gene"))  
combined_results <- left_join(combined_results, Annotation, by = "GeneID")  

# Filter the data for glycosyltransferase annotations  
#combined_results <- combined_results[combined_results$Annotation == "glucosyltransferase", ]  


# Define the output filename for merged results  
merged_output_file <- "Merged_DESeq2_Results.xlsx"  

# Write the combined results to a single Excel file  
write.xlsx(combined_results, file = merged_output_file, row.names = FALSE)  

# Print confirmation message for merged results  
print(paste("Merged results saved to:", merged_output_file)) 


# Print the available contrast names in result_list for debugging  
print("Available contrasts in result_list:")  
print(names(result_list))  

left_join_combination <- function(df1, df2, df3) {  
  df1 %>%  
    left_join(df2, by = "GeneID") %>%  
    left_join(df3, by = "GeneID") %>% 
    left_join(normalized_counts_averages, by = c("GeneID" = "Gene")) %>%
    left_join(Annotation, by = "GeneID")
}  

# Create a list to store results  
results_list <- list()  

# Loop through each combination  
for (combination in specific_combinaisons) {  
  # Split the combination string into individual data frame names  
  df_names <- strsplit(combination, "_")[[1]]  
  
  # Ensure all required data frames are available in the environment  
  df_list <- lapply(df_names, get)  
  
  # Perform the left joins  
  result <- do.call(left_join_combination, df_list)  
  
  # Store the result in the results list, using the combination name as the key  
  results_list[[combination]] <- result  
}  

# Check the results  
results_list  

# Save each result to an Excel file  
for (combination in names(results_list)) {  
  # Create a filename  
  filename <- paste0(gsub("_", "-", combination), ".xlsx")  
  
  # Save the corresponding data frame to the Excel file  
  write.xlsx(results_list[[combination]], file = filename, row.names = FALSE)  
} 


left_join_combination_filtered <- function(df1, df2, df3) {  
  # Perform left joins  
  combined_df <- df1 %>%  
    left_join(df2, by = "GeneID") %>%  
    left_join(df3, by = "GeneID") %>%  
    left_join(normalized_counts_averages, by = c("GeneID" = "Gene")) %>%  
    left_join(Annotation, by = "GeneID") %>%
    filter(str_detect(GO_Names, paste(Selected_function, collapse = "|")))  # Modifying your filter condition
  
 
  # Initialize a column for upregulated contrasts  
  combined_df <- combined_df %>% mutate(upregulated_contrasts = 0)  
  
  # Loop through the contrast names in result_list to create the regulation counts  
  for (name in names(result_list)) {  
    # Check if the current name includes "Susceptible"  
    if (grepl("Kis.Sus", name)) {  
      regulation_type <- paste0("regulation_2_", name)  
    } else {  
      regulation_type <- paste0("regulation_1.5_", name)  
    }  
    
    # Ensure the regulation_type column exists before using it  
    if (regulation_type %in% names(combined_df)) {  
      combined_df <- combined_df %>%  
        mutate(upregulated_contrasts = upregulated_contrasts + ifelse(get(regulation_type) == "Up", 1, 0))  
    }  
  }  
  
  # Filter to keep rows with up-regulated genes in at least 2 contrasts  
  combined_df <- combined_df %>% filter(upregulated_contrasts >= 2)  
  
  return(combined_df)  
}
 
# Create a list to store results  
results_list <- list()  

# Loop through each combination  
for (combination in specific_combinaisons) {  
  # Split the combination string into individual data frame names  
  df_names <- strsplit(combination, "_")[[1]]  
  
  # Ensure all required data frames are available in the environment  
  df_list <- lapply(df_names, get)  
  
  # Perform the left joins  
  result <- do.call(left_join_combination_filtered, df_list)  
  
  # Store the result in the results list, using the combination name as the key  
  results_list[[combination]] <- result  
}  

# Check the results  
results_list  


# Function to order data frames by the first FC (Fold Change) column it finds  
order_by_fc <- function(df) {  
  # Find the relevant FC columns  
  fc_column_names <- names(df)[grepl("FC_", names(df))]  
  
  # Sort the data frame by the first found FC column, if it exists  
  if (length(fc_column_names) > 0) {  
    # Use the first FC column for sorting  
    fc_column_name <- fc_column_names[1]  
    
    df <- df %>%  
      arrange(desc(.data[[fc_column_name]]))  # Sorting in descending order  
  }  
  
  return(df)  
}  

# Apply the ordering function to each data frame in results_list  
ordered_results_list <- lapply(results_list, order_by_fc)  

# Check the top entry of one of the sorted data frames  
head(ordered_results_list[[1]], 1)  # Adjust the index according to your needs

# Save each result to an Excel file  
for (combination in names(ordered_results_list)) {  
  # Create a filename  
  filename <- paste0(gsub("_", "-", combination), "_filtered", ".xlsx")  
  
  # Save the corresponding data frame to the Excel file  
  write.xlsx(ordered_results_list[[combination]], file = filename, row.names = FALSE)  
} 





# Create an empty list to store filtered data frames for each contrast  
filtered_df_list <- list()  

# Loop through the combined results to select the appropriate regulation column  
for (name in names(result_list)) {  
  # Check if the current name includes "Susceptible"  
  if (grepl("Kis.Sus", name)) {  
    regulation_type <- paste0("regulation_2_", name)  
  } else {  
    regulation_type <- paste0("regulation_1.5_", name)  
  }  
  
  # Check if the regulation column exists in the combined_results  
  if (regulation_type %in% colnames(combined_results)) {  
    # Create the filtered dataframe with GeneID and the relevant regulation column  
    filtered_df <- combined_results %>%  
      select(GeneID, all_of(regulation_type))  
    
    # Add to the list of filtered data frames  
    filtered_df_list[[name]] <- filtered_df  
  } else {  
    print(paste("Regulation column not found for:", regulation_type))  
  }  
}  

# Initialize final_filtered_results with the first data frame if it exists  
if (length(filtered_df_list) > 0) {  
  final_filtered_results <- filtered_df_list[[1]]  # Start with the first filtered dataframe  
  
  # Loop through the remaining filtered data frames and merge them by GeneID  
  for (i in 2:length(filtered_df_list)) {  
    final_filtered_results <- left_join(final_filtered_results, filtered_df_list[[i]], by = "GeneID", suffix = c("", paste0("_", names(filtered_df_list)[i])))  
  }  
  
  # Optionally, view the combined results  
  print(head(final_filtered_results))  
  
  # Define an output filename for the merged results  
  merged_output_file <- "Merged_Filtered_DESeq2_Results.xlsx"  
  
  # Write the merged results to an Excel file  
  write.xlsx(combined_results, file = merged_output_file, row.names = FALSE)  
  
  # Print confirmation message for merged results  
  print(paste("Merged filtered results saved to:", merged_output_file))  
} else {  
  print("No filtered results to merge.")  
}  



# Gene Ontology enrichment


# Read the combinaisons from a text file without showing warnings  
combinaisons_file <- "comparisons2.txt"  # Specify your file path here  
specific_combinaisons <- suppressWarnings(readLines(combinaisons_file))


candidate_genes_id <- read_excel(file.path(dir, "candidate_genes.xlsx"))
candidate_genes_id <- as.list(candidate_genes_id)


# Load annotation file generated with blast2GO and remove duplicated rows  
Annotation <- read_excel(file.path(dir, "Annotation.xlsx")) %>%  
  distinct(GeneID, .keep_all = TRUE)  # Remove duplicated GeneIDs  

# Function to extract GO terms and IDs based on prefix  
extract_GO_terms <- function(go_ids, go_names, prefix) {  
  id_pattern <- paste0(prefix, ":GO:")  
  terms_pattern <- paste0(prefix, ":")  
  
  go_id_result <- sapply(go_ids, function(x) {  
    go_terms <- unlist(strsplit(x, ";"))  
    go_terms <- go_terms[grep(id_pattern, go_terms)]  
    gsub(paste0(prefix, ":"), "", go_terms)  
  }) %>%   
    sapply(paste, collapse = ";")  # Combine go terms if multiple  
  
  go_terms_result <- sapply(go_names, function(x) {  
    go_terms <- unlist(strsplit(x, ";"))  
    go_terms <- go_terms[grep(paste0(prefix, ":"), go_terms)]  
    gsub(paste0(prefix, ":"), "", go_terms)  
  }) %>%  
    sapply(paste, collapse = ";")  # Combine go terms if multiple  
  
  return(list(go_ids = go_id_result, go_terms = go_terms_result))  
}  

# Extract GO terms for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC)  
bp_results <- extract_GO_terms(Annotation$GO_IDs, Annotation$GO_Names, "P")  
mf_results <- extract_GO_terms(Annotation$GO_IDs, Annotation$GO_Names, "F")  
cc_results <- extract_GO_terms(Annotation$GO_IDs, Annotation$GO_Names, "C")  

# Add extracted GO terms to the Annotation dataframe  
Annotation <- Annotation %>%  
  mutate(  
    BP_GO_ID = bp_results$go_ids,  
    BP_GO_terms = bp_results$go_terms,  
    MF_GO_ID = mf_results$go_ids,  
    MF_GO_terms = mf_results$go_terms,  
    CC_GO_ID = cc_results$go_ids,  
    CC_GO_terms = cc_results$go_terms  
  )  

# Clean up and filter only relevant columns  
GO_file <- Annotation %>%  
  select(GeneID, BP_GO_ID, BP_GO_terms, MF_GO_ID, MF_GO_terms, CC_GO_ID, CC_GO_terms, Enzyme_Codes, Enzyme_Names)  

# Display the first few rows of the final data frame  
head(GO_file) 


# Create a function to separate GO IDs and their corresponding terms while maintaining correspondence  
separate_go_terms <- function(data, id_col, term_col) {  
  data %>%  
    mutate(temp_IDs = strsplit(!!sym(id_col), ";"),  
           temp_terms = strsplit(!!sym(term_col), ";")) %>%  
    # Create lists of data frames, aligning IDs and terms  
    mutate(data_pairs = map2(temp_IDs, temp_terms, ~ {  
      # Create a data frame from the pairs  
      # Length of the IDs and terms lists  
      len <- max(length(.x), length(.y))  
      # Fill with NA where lists are shorter  
      df <- data.frame(  
        ID = c(.x, rep(NA, len - length(.x))),  
        Term = c(.y, rep(NA, len - length(.y))),  
        stringsAsFactors = FALSE  
      )  
      return(df)  
    })) %>%  
    # Unnest the data pairs into long format  
    unnest(data_pairs) %>%  
    # Select and rename to return only the original column names  
    select(GeneID, ID, Term) %>%  
    rename_with(~ id_col, .cols = ID) %>%
    rename_with(~ term_col, .cols = Term)
}  


# Create subsets for further analysis  
GO_file_MF <- separate_go_terms(GO_file, "MF_GO_ID", "MF_GO_terms")  
head(GO_file_MF)  # View the results  

GO_file_BP <- separate_go_terms(GO_file, "BP_GO_ID", "BP_GO_terms")  
head(GO_file_BP)  # View the results  

GO_file_CC <- separate_go_terms(GO_file, "CC_GO_ID", "CC_GO_terms")  
head(GO_file_CC)  # View the results  

GO_file_Enzyme <- separate_go_terms(GO_file, "Enzyme_Codes", "Enzyme_Names")  
head(GO_file_Enzyme)  # View the results  


# Function to truncate labels  
truncate_labels <- function(labels, max_length = 50) {  
  sapply(labels, function(x) {  
    if (is.na(x)) {  # Check for NA values  
      return(NA)     # Return NA if the value is NA  
    } else if (nchar(x) > max_length) {  
      paste0(substr(x, 1, max_length - 3), "...")  # Truncate and add ellipses  
    } else {  
      x  
    }  
  })  
}  

# Initialize a list to store all Result DataFrames if desired
results_list <- list()

# Prepare your gene lists: Define your genes of interest and the background gene list.  
background_genes <- unique(GO_file$GeneID)  # All genes from GO data  

# Function to perform GO enrichment analysis  
perform_go_enrichment <- function(candidate_genes, go_data, go_id_col, go_terms_col, output_prefix, list_name, background_genes) {  
  
  # Create contingency tables for GO IDs  
  go_counts_id <- go_data %>%  
    group_by(!!sym(go_id_col)) %>%  
    summarize(total_genes_in_GO = n(), .groups = 'drop')  
  
  cand_go_counts_id <- go_data %>%  
    filter(GeneID %in% candidate_genes) %>%  
    group_by(!!sym(go_id_col)) %>%  
    summarize(Cands_in_GO = n(), .groups = 'drop')  
  
  # Merge the counts for GO IDs  
  enrichment_data_id <- full_join(go_counts_id, cand_go_counts_id, by = go_id_col) %>%  
    mutate(  
      Cands_in_GO = replace_na(Cands_in_GO, 0),  
      total_genes_in_GO = replace_na(total_genes_in_GO, 0),  
      total_genes_in_background = length(background_genes),  
      GeneRatio = Cands_in_GO / total_genes_in_GO,  
      fold_enrichment = (Cands_in_GO / length(candidate_genes)) /  
        (total_genes_in_GO / total_genes_in_background)  
    )  
  
  # Calculate p-values using Fisher's exact test for GO IDs  
  enrichment_data_id <- enrichment_data_id %>%  
    rowwise() %>%  
    mutate(  
      p_value = fisher.test(matrix(c(Cands_in_GO,  
                                     length(candidate_genes) - Cands_in_GO,  
                                     total_genes_in_GO - Cands_in_GO,  
                                     total_genes_in_background - length(candidate_genes)),  
                                   nrow = 2))$p.value  
    )   
  
  # Adjust p-values  
  enrichment_data_id <- enrichment_data_id %>%  
    mutate(adj_p_value = p.adjust(p_value, method = "bonferroni"))  
  
  all_significant_results_id <- enrichment_data_id %>%  
    filter(adj_p_value < 0.05, fold_enrichment > 2) %>%  
    arrange(desc(Cands_in_GO))
  
  # Save all results for GO IDs as environment object
  assign(paste0(output_prefix, "_", list_name, "_all_significant_results_id"), all_significant_results_id, envir = .GlobalEnv)
  # Or, alternatively, save to the results list:
  results_list[[paste0(output_prefix, "_", list_name, "_all_significant_results_id")]] <- all_significant_results_id
  write.xlsx(all_significant_results_id, paste0(output_prefix, "_", list_name, "_all_significant_results_id.xlsx")) 
  
  # Filter results for GO IDs and keep only the top 30 enriched
  significant_results_id <- enrichment_data_id %>%  
    filter(adj_p_value < 0.05, fold_enrichment > 2) %>%  
    arrange(desc(Cands_in_GO)) %>%{
      rbind(head(., 30))
    }  # Select top 30 enriched GO IDs
  
  # Save results for top 30 GO IDs as environment object
  assign(paste0(output_prefix, "_", list_name, "_top_30_significant_results_id"), significant_results_id, envir = .GlobalEnv)
  results_list[[paste0(output_prefix, "_", list_name, "_top_30_significant_results_id")]] <- significant_results_id
  write.xlsx(significant_results_id, paste0(output_prefix, "_", list_name, "_top_30_significant_results_id.xlsx"))
  
  # Apply truncation to GO IDs for plotting
  significant_results_id[[go_id_col]] <- truncate_labels(significant_results_id[[go_id_col]])

  # Plot for GO IDs  
  go_id_plot <- ggplot(significant_results_id, aes(x = fold_enrichment, y = reorder(!!sym(go_id_col), GeneRatio), size = Cands_in_GO, color = adj_p_value)) +  
    geom_point(alpha = 0.8) +  
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +  
    labs(x = "Fold enrichment", y = "GO ID", size = "N genes", color = "FDR",    
         title = str_wrap(paste(output_prefix, "ID for", list_name), width = 50)) +  
    theme_minimal() +  
    theme(legend.position = "right")  
  
  # Save the plot as PNG  
  ggsave(paste0(output_prefix, "_", list_name, "_GO_ID_plot.png"), plot = go_id_plot, width = 1500, height = 3000, unit = "px", dpi = 400)  
  
  
  # Create contingency tables for GO Terms
  go_counts_terms <- go_data %>%  
    group_by(!!sym(go_terms_col)) %>%  
    summarize(total_genes_in_GO = n(), .groups = 'drop')  
  
  cand_go_counts_terms <- go_data %>%  
    filter(GeneID %in% candidate_genes) %>%  
    group_by(!!sym(go_terms_col)) %>%  
    summarize(Cands_in_GO = n(), .groups = 'drop')  
  
  # Merge the counts for GO Terms
  enrichment_data_terms <- full_join(go_counts_terms, cand_go_counts_terms, by = go_terms_col) %>%  
    mutate(  
      Cands_in_GO = replace_na(Cands_in_GO, 0),  
      total_genes_in_GO = replace_na(total_genes_in_GO, 0),  
      total_genes_in_background = length(background_genes), 
      GeneRatio = Cands_in_GO / total_genes_in_GO,  
      fold_enrichment = (Cands_in_GO / length(candidate_genes)) /  
        (total_genes_in_GO / total_genes_in_background)  
    )  
  
  # Calculate p-values using Fisher's exact test for GO Terms
  enrichment_data_terms <- enrichment_data_terms %>%  
    rowwise() %>%  
    mutate(  
      p_value = fisher.test(matrix(c(Cands_in_GO,  
                                     length(candidate_genes) - Cands_in_GO,  
                                     total_genes_in_GO - Cands_in_GO,  
                                     total_genes_in_background - length(candidate_genes)),  
                                   nrow = 2))$p.value  
    )   
  
  # Adjust p-values
  enrichment_data_terms <- enrichment_data_terms %>%  
    mutate(adj_p_value = p.adjust(p_value, method = "bonferroni"))  
  
  all_significant_results_terms <- enrichment_data_terms %>%  
    filter(adj_p_value < 0.05, fold_enrichment > 2) %>%  
    arrange(desc(Cands_in_GO))
  
  # Save all results for GO terms as environment object
  assign(paste0(output_prefix, "_", list_name, "_all_significant_results_terms"), all_significant_results_terms, envir = .GlobalEnv)
  # Or, alternatively, save to the results list:
  results_list[[paste0(output_prefix, "_", list_name, "_all_significant_results_terms")]] <- all_significant_results_terms
  write.xlsx(all_significant_results_terms, paste0(output_prefix, "_", list_name, "_all_significant_results_terms.xlsx")) 
  
  # Filter results for GO Terms and keep only the top 30 enriched
  significant_results_terms <- enrichment_data_terms %>%   
    filter(adj_p_value < 0.05, fold_enrichment > 2) %>%  
    arrange(desc(Cands_in_GO)) %>%{
      rbind( head(., 30))
    }  # Select top 20 enriched GO Terms
  
  # Save results for top 20 GO Terms as environment object
  assign(paste0(output_prefix, "_", list_name, "_top_30_significant_results_terms"), significant_results_terms, envir = .GlobalEnv)
  results_list[[paste0(output_prefix, "_", list_name, "_top_30_significant_results_terms")]] <- significant_results_terms
  write.xlsx(significant_results_terms, paste0(output_prefix, "_", list_name, "_top_30_significant_results_terms.xlsx"))
  
  # Apply truncation to GO Terms for plotting
  significant_results_terms[[go_terms_col]] <- truncate_labels(significant_results_terms[[go_terms_col]])

  # Plot for GO terms  
  go_terms_plot <- ggplot(significant_results_terms, aes(x = fold_enrichment, y = reorder(!!sym(go_terms_col), GeneRatio), size = Cands_in_GO, color = adj_p_value)) +  
    geom_point(alpha = 0.8) +  
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +  
    labs(x = "Fold enrichment", y = "GO terms", size = "N genes", color = "FDR",    
         title = str_wrap(paste(output_prefix, "terms for", list_name), width = 50)) +  
    theme_minimal() +  
    theme(legend.position = "right")  
  
  # Save the plot as PNG  
  ggsave(paste0(output_prefix, "_", list_name, "_GO_terms_plot.png"), plot = go_terms_plot, width = 2500, height = 3000, unit = "px", dpi = 400)  

  # Combine Terms and IDs in same file
  
  
  all_significant_results <- significant_results_id %>%
    left_join(
      all_significant_results_terms,
      by = c(
        "total_genes_in_GO",
        "Cands_in_GO",
        "total_genes_in_background",
        "GeneRatio",
        "fold_enrichment",
        "p_value",
        "adj_p_value"
      )
    ) %>%
    rename(
      "Genes in GO" = total_genes_in_GO,
      "Genes in background" = total_genes_in_background,
      "p value" = p_value,
      "N genes" = Cands_in_GO,
      "Gene Ratio" = GeneRatio,
      "Fold enrichment" = fold_enrichment,
      "FDR" = adj_p_value
    )
  
  write.xlsx(
    all_significant_results,
    paste0(output_prefix, "_", list_name, "_all_significant_results.xlsx")
  )
  
  
}

# Initialize an empty list to store candidate genes for each combination  
candidate_genes_list <- list()  

# Iterate through each name in candidate_genes_id  
for (combination_name in names(candidate_genes_id)) {  
  
  # Extract the candidate genes for the current combination  
  candidate_genes <- candidate_genes_id[[combination_name]]  
  
  # Store the result in the candidate_genes_list with the name as the current combination  
  candidate_genes_list[[combination_name]] <- candidate_genes  
  
  # Use the combination name as the list name for GO enrichment analysis  
  list_name <- combination_name  # Directly use the combination name  
  
  # Use the background genes in your function calls  
  perform_go_enrichment(candidate_genes, GO_file_MF, "MF_GO_ID", "MF_GO_terms", "MF", list_name, background_genes)  
  perform_go_enrichment(candidate_genes, GO_file_BP, "BP_GO_ID", "BP_GO_terms", "BP", list_name, background_genes)  
  perform_go_enrichment(candidate_genes, GO_file_CC, "CC_GO_ID", "CC_GO_terms", "CC", list_name, background_genes)  
  perform_go_enrichment(candidate_genes, GO_file_Enzyme, "Enzyme_Codes", "Enzyme_Names", "Enzyme", list_name, background_genes)  
}



