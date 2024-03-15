#!/usr/bin/env Rscript

## Script name: DEAnalysis.R
##
## Purpose of script: Differential Expression Analysis
##
## Author: Klaudia Adamowicz
## Email: klaudia.adamowicz@uni-hamburg.de
##
## Date Created: 2024-02-29
##
## Copyright (c) Dr. Tanja Laske, 2023

# Load required libraries -----
suppressPackageStartupMessages({
  required_packages <- c("optparse","data.table","dplyr","stringr","ggrepel", "rjson", "tidyr")
  for(package in required_packages){
    #if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
    library(package, character.only = TRUE, quietly = TRUE)
  }
})

## Parse arguments -----

# set up arguments
parser <- OptionParser()
parser <- add_option(parser, c("-i","--in_dir"), help="Directory with files from DE analysis", default="")
parser <- add_option(parser, c("-o","--out_dir"), help="Directory for output files from Summarizer", default="")

# Adding new options for thresholds with defaults
parser <- add_option(parser, c("--logFC"), help = "Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--logFC_up"), help = "Upper log2 fold change threshold (dividing into upregulated)", type = "numeric", default = 1)
parser <- add_option(parser, c("--logFC_down"), help = "Lower log2 fold change threshold (dividing into downregulated)", type = "numeric", default = -1)
parser <- add_option(parser, c("--p_adj"), help = "Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--alpha"), help = "Threshold for adjusted p-values or p-values", type = "numeric", default = 0.05)

# get command line options, if help option encountered print help and exit
args <- parse_args(parser)

in_dir <- args$in_dir
out_dir <- args$out_dir

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) #stops warnings if folder already exists

## Read all DE files ----

#out_dir <- args$out_dir
de_data_dir <- file.path(in_dir, "de_data")

# Get the list of .txt files in the directory
file_list <- list.files(de_data_dir, pattern = "\\.txt$", full.names = TRUE)

# Function to read file and create a list with dataframe and file name
read_file_to_list <- function(file_path) {
  # Read the file into a dataframe
  de_results <- read.table(file_path, sep = " ", header = TRUE)
  
  # Extract the file name without extension and prefix in one line
  comparison <- sub("^DEdata.", "", sub("\\.txt$", "", basename(file_path)))
  
  # Return a list with dataframe and file name
  list(de_results = de_results, file_name = basename(file_path), comparison = comparison)
}

# Apply the function to all files and create the list
de_list <- lapply(file_list, read_file_to_list)

names(de_list) <- sapply(de_list, function(x) x$comparison)

## Update DE results change by thresholds ----

# Function to update the Change column
update_change <- function(df) {
  if (args$logFC) {
    if (args$p_adj) {
      # Check both logFC thresholds and adjusted p-value
      df$Change <- ifelse(df$adj.P.Val < args$alpha & df$logFC > args$logFC_up, "Up regulated",
                          ifelse(df$adj.P.Val < args$alpha & df$logFC < args$logFC_down, "Down regulated", "No Change"))
    } else {
      # Check logFC thresholds and raw p-value
      df$Change <- ifelse(df$P.Value < args$alpha & df$logFC > args$logFC_up, "Up regulated",
                          ifelse(df$P.Value < args$alpha & df$logFC < args$logFC_down, "Down regulated", "No Change"))
    }
  } else {
    # If logFC is FALSE, just check alpha (p-value threshold)
    p_col <- ifelse(args$p_adj, "adj.P.Val", "P.Value")
    df$Change <- ifelse(df[[p_col]] < args$alpha, "Significant", "No Change")
  }
  return(df)
}

# Apply the function to all dataframes in the list
de_list <- lapply(de_list, function(x) {
  x$de_results <- update_change(x$de_results)
  return(x)
})

## Calculate statistics ----

### Top comparisons ----

# Initialize an empty list to store counts
comparison_counts <- list()

# Loop through each comparison in de_list
for (comparison in names(de_list)) {
  df <- de_list[[comparison]]$de_results
  
  if (args$logFC) {
    # Count the occurrences for "Up regulated" and "Down regulated"
    up_count <- sum(df$Change == "Up regulated", na.rm = TRUE)
    down_count <- sum(df$Change == "Down regulated", na.rm = TRUE)
    significant_count <- up_count + down_count  # Up and Down are both considered significant
  } else {
    # When logFC is FALSE, count only "Significant"
    up_count <- 0
    down_count <- 0
    significant_count <- sum(df$Change == "Significant", na.rm = TRUE)
  }
  
  # Store the counts in a list
  comparison_counts[[comparison]] <- c(up = up_count, down = down_count, total = significant_count)
}

# Convert the list to a dataframe
top_comparisons <- do.call(rbind, comparison_counts)
top_comparisons <- data.frame(comparison = rownames(top_comparisons), top_comparisons, row.names = NULL)

# Sort the dataframe by the 'total' column
top_comparisons <- top_comparisons[order(-top_comparisons$total), ]

write.table(top_comparisons, file.path(out_dir,"top-comparisons.tsv"), sep="\t", row.names=FALSE)

### Top IDs ----

# Function to process each dataframe and count gene occurrences
process_genes <- function(df, logFC) {
  # Function to get unique genes and count them once per category
  count_unique_genes <- function(genes, label) {
    unique_genes <- unique(unlist(strsplit(genes, ";")))
    if (length(unique_genes) == 0) {
      return(data.frame(gene=character(), label=factor(levels = c("up", "down", "total")), count=integer()))
    }
    data.frame(gene = unique_genes, label = label, count = 1)
  }
  
  # Prepare gene lists for each category
  up_genes <- if(logFC) count_unique_genes(df$gene[which(df$Change == "Up regulated")], "up") else 0
  down_genes <- if(logFC) count_unique_genes(df$gene[which(df$Change == "Down regulated")], "down") else 0
  significant_genes <- count_unique_genes(df$gene[which(df$Change %in% c("Up regulated", "Down regulated", "Significant"))], "total")
  
  # Combine all gene data
  all_genes <- rbind(up_genes, down_genes, significant_genes)
  
  # Return processed data
  if (!is.null(all_genes) && nrow(all_genes) > 0) {
    all_genes %>% 
      group_by(gene, label) %>% 
      summarise(count = sum(count)) %>% 
      spread(label, count, fill = 0)
  } else {
    return(data.frame(gene=character(), up=integer(), down=integer(), total=integer()))
  }
}

# Apply function to each comparison and bind rows
gene_counts <- bind_rows(lapply(de_list, function(x) process_genes(x$de_results, args$logFC)))

# Final aggregation across all comparisons
final_gene_counts <- gene_counts %>% 
  group_by(gene) %>% 
  summarise(up = sum(up, na.rm = TRUE), 
            down = sum(down, na.rm = TRUE), 
            total = sum(total, na.rm = TRUE)) %>% 
  as.data.frame()

# Sort by the total column
final_gene_counts <- final_gene_counts[order(-final_gene_counts$total), ]

write.table(final_gene_counts, file.path(out_dir,"top-ids.tsv"), sep="\t", row.names=FALSE)

