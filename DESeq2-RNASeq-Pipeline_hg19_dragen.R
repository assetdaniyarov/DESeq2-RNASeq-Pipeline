# DESeq2 Analysis Script

library(DESeq2)
library(tximport)
library(rtracklayer)
library(dplyr)

# file paths
gtf_file <- "gencode.v19.annotation.gtf"
input_folder <- "./input/"  # Input CSV
output_folder <- "./output/"  # Output CSV
annotated_folder <- "./output/annotated/"  # Annotated CSV
mt_genes_folder <- "./output/mt_genes/"  # MT genes

# output directories
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(annotated_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(mt_genes_folder, showWarnings = FALSE, recursive = TRUE)

# GTF
gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

transcripts <- gtf_df %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id, gene_name)

# files
process_data <- function(dataset_id) {
  cat("Processing dataset:", dataset_id, "\n")
  
  files <- list.files(path = input_folder, pattern = paste0(dataset_id, ".quant.sf$"), full.names = TRUE)
  sample_names <- gsub(".quant.sf", "", basename(files))
  names(files) <- sample_names
  
txi <- tximport(files, type = "salmon", txOut = TRUE)
  
  # metadata
group_labels <- sub("_\\d+$", "", sample_names)
  metadata <- data.frame(sample = sample_names, group = factor(group_labels))
  
  cat("Number of samples:", nrow(metadata), "\n")
  cat("Unique groups:", length(unique(metadata$group)), "\n")
  
  # DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~group)
dds <- dds[rowSums(counts(dds)) >= 1, ]
  
  # DESeq2
dds <- DESeq(dds, fitType = "local")
  
  # PCA
  vsd <- vst(dds, blind = FALSE, fitType = "local")
  pca_plot <- plotPCA(vsd, intgroup = "group")
  print(pca_plot)
  
  # comparisons
  comparisons <- list(
    c("GM2_80S_" %s+% dataset_id, "0_80S_" %s+% dataset_id),
    c("GM2_pol_" %s+% dataset_id, "0_pol_" %s+% dataset_id),
    c("GM_total_" %s+% dataset_id, "0_total_" %s+% dataset_id)
  )
  
  analyze_comparison <- function(dds, group2, group1, name) {
    res <- results(dds, contrast = c("group", group2, group1), alpha = 0.1)
    res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
    res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.1, ]
    
    res_df <- as.data.frame(res)
    res_df <- cbind(transcript_id = rownames(res_df), res_df)
    output_file <- file.path(output_folder, paste0(name, "_results.csv"))
    write.csv(res_df, output_file, row.names = FALSE)
    cat("Results saved:", output_file, "\n")
  }
  
  for (comp in comparisons) {
    name <- paste0(comp[2], "_vs_", comp[1])
    analyze_comparison(dds, comp[2], comp[1], name)
  }
}

process_data("549")
process_data("231")

annotate_results <- function() {
  csv_files <- list.files(output_folder, pattern = "\.csv$", full.names = TRUE)
  
  for (csv_file in csv_files) {
    cat("Annotating:", csv_file, "\n")
    input_data <- read.csv(csv_file, header = TRUE)
    
    if (!"transcript_id" %in% colnames(input_data)) {
      cat("Skipping file, missing 'transcript_id':", csv_file, "\n")
      next
    }
    
    annotated_data <- input_data %>%
      inner_join(transcripts, by = "transcript_id") %>%
      select(transcript_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
    
    output_file <- file.path(annotated_folder, paste0(basename(tools::file_path_sans_ext(csv_file)), "_annotated.csv"))
    write.csv(annotated_data, output_file, row.names = FALSE)
    cat("Annotated file saved:", output_file, "\n")
  }
}

annotate_results()

# MT genes
extract_mt_genes <- function() {
  input_files <- list.files(annotated_folder, pattern = "\.csv$", full.names = TRUE)
  
  for (file in input_files) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    filtered_data <- subset(data, grepl("^MT-", gene_name))
    if (nrow(filtered_data) > 0) {
      output_file <- file.path(mt_genes_folder, paste0(tools::file_path_sans_ext(basename(file)), "_mt.csv"))
      write.csv(filtered_data, output_file, row.names = FALSE)
      cat("MT gene file saved:", output_file, "\n")
    } else {
      cat("No MT genes found in:", basename(file), "\n")
    }
  }
}

extract_mt_genes()

cat("All processing complete!\n")
