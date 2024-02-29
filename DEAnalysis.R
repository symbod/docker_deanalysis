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
  required_packages <- c("optparse","data.table","dplyr","stringr","ggrepel", "rjson")
  for(package in required_packages){
    if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
    library(package, character.only = TRUE, quietly = TRUE)
  }
  if(!require("limma",character.only = TRUE, quietly = TRUE)) BiocManager::install("limma")
})

# Methods --------------------------

# Method performing limma
# input:
# ----- set_condition: vector of experimental design specifying the condition(s)
# ----- set_counts: counts (rows = proteins, columns = samples)
# ----- set_comparisons: vector of comparisons
perform_de <- function(set_condition, set_counts, set_comparisons){
  # create design matrix
  groupsM <- as.factor(set_condition)
  designM <- model.matrix(~0+groupsM) 
  colnames(designM) <-levels(groupsM) 
  fit <- lmFit(set_counts, designM)
  
  # create contrasts
  contr <- makeContrasts(contrasts = set_comparisons, levels = colnames(coef(fit)))
  
  fit2 <- contrasts.fit(fit, contr)
  ebfit <- eBayes(fit2, trend = TRUE)
  return(ebfit)
}


# Method to extract results of fit object
# input:
# ----- set_fit: fit object of the perform_de method
# ----- set_comparisons: vector of comparisons
# ----- out_dir: output directory
# ----- lfc_up: log fold change threshold
# ----- lfc_down: log fold change threshold
# ----- alpha: p-value or p.adjust threshold (depending on padj)
# ----- strict: TRUE if < and >, FALLSE if <= and >=
# ----- padj: TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
# ----- logFC_thr: TRUE if threshold for logFC should be set, FALSE if not
compare_de_expr <- function(set_fit, set_comparisons, out_dir, lfc_up = args$logFC_up, lfc_down = args$logFC_down, alpha = args$alpha, 
                            strict = FALSE, padj = args$p_adj, logFC_thr=args$logFC, write_output = TRUE){
  de_results <- list()
  for (i in 1:length(set_comparisons)){
    # get results for each comparison
    top.table <- topTable(set_fit, sort.by = "P", number=Inf, coef=c(i))
    gene_reg <- setDT(top.table, keep.rownames = "gene") # save row names as column
    
    # different threshold settings
    if (logFC_thr){
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$adj.P.Val < alpha, "Up regulated", ifelse(gene_reg$logFC < lfc_down & gene_reg$adj.P.Val < alpha, "Down regulated", "No change"))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$P.Value < alpha, "Up regulated", ifelse(gene_reg$logFC < lfc_down & gene_reg$P.Value < alpha, "Down regulated", "No change"))
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$adj.P.Val <= alpha, "Up regulated", ifelse(gene_reg$logFC <= lfc_down & gene_reg$adj.P.Val <= alpha, "Down regulated", "No change"))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$P.Value <= alpha, "Up regulated", ifelse(gene_reg$logFC <= lfc_down & gene_reg$P.Value <= alpha, "Down regulated", "No change"))
        }
      }
    } else {
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val < alpha, "Up regulated", "No change")
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value < alpha, "Up regulated", "No change")
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val <= alpha, "Up regulated", "No change")
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value <= alpha, "Up regulated", "No change")
        }
      }
    }
    gene_reg$Change <- as.factor(gene_reg$Change)
    if (write_output){
      write.table(gene_reg, file = file.path(out_dir,"de_data",paste0("DEdata.", set_comparisons[[i]],".txt")),sep=" ",row.names = FALSE) 
    }
    de_results[[set_comparisons[[i]]]] <- gene_reg
  }
  return(de_results)
}


# Method to do a differential expression analysis
# input:
# ----- counts: counts (rows = proteins, columns = samples)
# ----- md: experimental design data.frame
# ----- condition_name: colname of vector of experimental design representing the condition to analyse
# ----- comparisons: vector of comparisons
# ----- out_dir: output directory
# ----- lfc_up: log fold change threshold
# ----- lfc_down: log fold change threshold
# ----- alpha: p-value or p.adjust threshold (depending on padj)
# ----- strict: TRUE if < and >, FALLSE if <= and >=
# ----- padj: TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
# ----- logFC_thr: TRUE if threshold for logFC should be set, FALSE if not
# ----- plot: TRUE if plots should be rendered, FALSE if not
# ----- write_output: TRUE if files with results should be written, else FALSE
run_DE <- function(counts, md, condition_name, comparisons, out_dir, lfc_up = args$logFC_up, lfc_down = args$logFC_down, alpha = args$alpha, 
                   strict = FALSE, padj = args$p_adj, logFC_thr=args$logFC, plot = TRUE, write_output = TRUE){
  if (write_output){
    # create output directories
    dir.create(out_dir, showWarnings = FALSE) #stops warnings if folder already exists
    dir.create(file.path(out_dir,"de_data"), showWarnings = FALSE) #stops warnings if folder already exists
    
  }
  
  # get condition vector
  condition <- md[, condition_name]
  
  # perform DE analysis pairwise
  pairwise_de <- perform_de(set_condition = condition, set_counts = counts, set_comparisons = comparisons)
  # extract results
  de_results <- compare_de_expr(set_fit = pairwise_de, set_comparisons = comparisons, out_dir = out_dir, lfc_up = lfc_up, 
                                lfc_down = lfc_down, alpha = alpha, strict = strict, padj = padj, logFC_thr = logFC_thr, write_output = write_output)
  
  # save overview values
  overview <- save_de_results(de_results = de_results, out_dir = out_dir, cond = condition_name, write_output = write_output)
  
  return(list(de_results = de_results, cond_sum = overview$sum, cond_updown = overview$updown))
}


# Method to save overview data.frame of the results of compare_de_expr
# input:
# ----- de_results: results of the method compare_de_expr
# ----- out_dir: output directory
save_de_results <- function(de_results, out_dir, cond, write_output = TRUE){
  comparisons <- names(de_results)
  elements <- unique(unlist(strsplit(comparisons, "-")))
  sum <- data.frame(matrix(ncol=length(elements),nrow=length(elements), dimnames=list(elements, elements)))
  updown <- data.frame(matrix(ncol=length(elements), nrow=length(elements), dimnames=list(elements, elements)))
  for (comp in comparisons){
    de_result <- de_results[[comp]]
    sample_one <- strsplit(comp, split="-")[[1]][1]
    sample_two <- strsplit(comp, split="-")[[1]][2]
    
    # save meta data into json 
    sample_group <- if (cond == "TimeCond") list("Timepoint", "Condition") else if (cond == "TimeCondLoc") list("Timepoint", "Condition", "Location") else list("Condition")
    meta = list(metadata=list(
      samples_group_A = meta_data[meta_data[[cond]] %in% c(sample_one)]$Column_name,
      samples_group_B = meta_data[meta_data[[cond]] %in% c(sample_two)]$Column_name,
      group_A = str_split(sample_one,"_")[[1]],
      group_B = str_split(sample_two,"_")[[1]],
      sample_group = sample_group
    ), DE_file = file.path(out_dir,"de_data",paste0("DEdata.",sample_one,'-',sample_two,".txt"))
    )
    write(toJSON(meta), file.path(out_dir,"de_data",paste0("DEdata.",sample_one,'-',sample_two,".txt.json")))
    
    # save number of DEs
    sum[sample_one, sample_two] = nrow(de_result[de_result$Change != "No change"])
    updown[sample_one, sample_two] = nrow(de_result[de_result$Change == "Up regulated"])
    updown[sample_two, sample_one] = nrow(de_result[de_result$Change == "Down regulated"])
  }
  
  if (write_output){
    write.csv(sum, file = file.path(out_dir, paste0("DE_overview_", cond, "-sum.csv")),row.names = TRUE) 
    write.csv(updown, file = file.path(out_dir, paste0("DE_overview_", cond, "-updown.csv")),row.names = TRUE)
  }
  return(list(sum=sum, updown=updown))
}


# function to make unsupervised clustering of the samples
# and creating plots for visualization. 
# Two types of clusterings are supported: 
# -- common: clustering based solely on genes shared over all samples (PCA)
# -- pairwise: clustering based on shared genes between each sample pair 
#            individually (MDS)
#
# input:
# ----- conditions: factor containing sample identificators
# ----- dge_obj: DGEObject containing count data
# ----- filename: name for the plot file
# ----- set_type: type of plot "common" for PCA and "pairwise" for MDS
# ----- out_dir: output directory for the plot
clustering_plot <- function(conditions, data, filename, set_type, out_dir){
  name_list <- list(common='pca', pairwise='mds')
  # create clustering
  #cond <- factor(plyr::mapvalues(names(data), from = meta_data$Sample_name, to = meta_data[[by]], warn_missing=FALSE))
  coord <- plotMDS(data, pch=19, gene.selection=set_type)
  dev.off()
  write.table(data.frame(meta_data$Column_name, coord$x, coord$y), file = file.path(out_dir, paste0(name_list[[set_type]],".csv")), sep=",", col.names=TRUE, row.names = FALSE)
}



# Differential Expression Analysis -----------------

## Parse arguments -----
# set up arguments
parser <- OptionParser()
parser <- add_option(parser, c("-m","--meta_file"), help="Meta data description file")
parser <- add_option(parser, c("-c","--count_file"), help="Preprocessed count file")
parser <- add_option(parser, c("-o","--out_dir"), help="Directory for output files", default="")

# Adding new options for thresholds with defaults
parser <- add_option(parser, c("--logFC"), help = "Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--logFC_up"), help = "Upper log2 fold change threshold (dividing into upregulated)", type = "numeric", default = 1)
parser <- add_option(parser, c("--logFC_down"), help = "Lower log2 fold change threshold (dividing into downregulated)", type = "numeric", default = -1)
parser <- add_option(parser, c("--p_adj"), help = "Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--alpha"), help = "Threshold for adjusted p-values or p-values", type = "numeric", default = 0.05)


# get command line options, if help option encountered print help and exit
args <- parse_args(parser)
# check if mandatory options are set
check_options <- function(tags){
  for (tag in tags){
    if (is.null(args[tag])){
      print_help(parser)
      stop("Missing mandatory option.", call.=FALSE)
    }
  }
}
check_options(c('meta_file','count_file'))

# save arguments
meta_file_path <- args$meta_file 
count_file_path <- args$count_file 
out_dir <- args$out_dir

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) #stops warnings if folder already exists
dir.create(file.path(out_dir,"de_data"), showWarnings = FALSE, recursive = TRUE) # for de data

## Load data ----
meta_data <- fread(meta_file_path)
count_data <- fread(count_file_path)

## Prepare data ----

### Correct data ----
# remove ref
meta_data <- meta_data[meta_data$ID != "ref",]

# convert timepoint column
meta_data[, Timepoint := as.numeric(Timepoint)]

### Rename Columns ---

# rename from file name to sample name
names(count_data) <- plyr::mapvalues(names(count_data), from = meta_data$Column_name, to = meta_data$Sample_name, warn_missing=FALSE) 
meta_data <- subset(meta_data, meta_data$Sample_name %in% names(count_data))

# remove columns that are not in meta_data
columns_to_keep <- c("Protein.IDs", "Gene.Names", meta_data$Sample_name)
existing_columns <- columns_to_keep[columns_to_keep %in% names(count_data)]
count_data <- count_data[, existing_columns, with = FALSE]

### Encode Columns ----

mappings <- list(
  Condition = list(setNames(paste0("C", seq_along(unique(meta_data$Condition))), unique(meta_data$Condition))),
  Location = if("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1) 
    list(setNames(paste0("L", seq_along(unique(meta_data$Location))), unique(meta_data$Location))),
  Timepoint = list(setNames(paste0("T", seq_along(sort(unique(meta_data$Timepoint)))), sort(unique(meta_data$Timepoint))))
)

write(toJSON(mappings), file.path(out_dir, "Metadata_encoding.txt.json")) 

reverse_mappings <- list(
  Condition = if (!is.null(mappings$Condition)) setNames(names(mappings$Condition[[1]]), unlist(mappings$Condition[[1]])) else NULL,
  Location = if (!is.null(mappings$Location)) setNames(names(mappings$Location[[1]]), unlist(mappings$Location[[1]])) else NULL,
  Timepoint = if (!is.null(mappings$Timepoint)) setNames(names(mappings$Timepoint[[1]]), unlist(mappings$Timepoint[[1]])) else NULL
)

write(toJSON(reverse_mappings), file.path(out_dir, "Metadata_reverse_encoding.txt.json")) 

# rename condition
meta_data[, Condition := mappings$Condition[[1]][Condition]]
# rename timepoint
meta_data[, Timepoint := as.character(Timepoint)]
meta_data[, Timepoint := mappings$Timepoint[[1]][Timepoint]]
# rename location
if(!is.null(mappings$Location)) {
  meta_data[, Location := mappings$Location[[1]][Location]]
}

### Create matrix and samples ----

# create count matrix 
counts <- as.data.frame(count_data[,c(meta_data$Sample_name), with=F], )
rownames(counts) <- count_data$Protein.IDs
counts <- counts[rowSums(is.na(counts)) != ncol(counts), ] # remove where full row is NA


## Save PCA and MDS ----
clustering_plot(data = counts, filename = paste0("mds_plot_",by_type), set_type = "pairwise", out_dir = out_dir)
clustering_plot(data = counts, filename = paste0("pca_plot_",by_type), set_type = "common", out_dir = out_dir)

## Pairwise comparison of each condition x time point ----
# save sample groups
target_column <- ifelse("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1, 
                        "TimeCondLoc", 
                        "TimeCond")
# Create the target column based on the presence of 'Location'
meta_data[, (target_column) := if(target_column == "TimeCondLoc") {
  paste(meta_data$Timepoint, meta_data$Condition, meta_data$Location, sep = "_")} else {
  paste(meta_data$Timepoint, meta_data$Condition, sep = "_")}]

### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data[[target_column]]))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    # only if time point or condition are the same!
    split_a = str_split(elements[index_a],"_")[[1]]
    split_b = str_split(elements[index_b],"_")[[1]]
    # Count the number of differences
    differences <- sum(split_a != split_b)
    # Only add to comparisons if exactly one of timepoint, condition, or location is different
    if (differences == 1){
      comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
    }
  }
}
### Perform DE analysis ----
results <- run_DE(counts=counts, md=as.data.frame(meta_data), condition_name=target_column, comparisons=comparisons, out_dir=out_dir, write_output = TRUE)

## Pairwise comparison of each condition ----
### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data$Condition))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
  }
}
### Perform DE analysis ----
results <- run_DE(counts=counts, md=as.data.frame(meta_data), condition_name="Condition", comparisons=comparisons, out_dir=out_dir, write_output = TRUE)

