library(DESeq2)
library(tximport)
library(dplyr)
library(rtracklayer)

# 549

files <- c(
  "0_80S_549_1.quant.sf", "0_80S_549_2.quant.sf", "0_80S_549_3.quant.sf",
  "0_pol_549_1.quant.sf", "0_pol_549_2.quant.sf", "0_pol_549_3.quant.sf",
  "0_total_549_1.quant.sf", "0_total_549_2.quant.sf", "0_total_549_3.quant.sf",
  "AA_80S_549_1.quant.sf", "AA_80S_549_2.quant.sf", "AA_80S_549_3.quant.sf",
  "AA_pol_549_1.quant.sf", "AA_pol_549_2.quant.sf", "AA_pol_549_3.quant.sf",
  "AA_total_549_1.quant.sf", "AA_total_549_2.quant.sf", "AA_total_549_3.quant.sf",
  "Gl_80S_549_1.quant.sf", "Gl_80S_549_2.quant.sf", "Gl_80S_549_3.quant.sf",
  "Gl_pol_549_1.quant.sf", "Gl_pol_549_2.quant.sf", "Gl_pol_549_3.quant.sf",
  "Gl_total_549_1.quant.sf", "Gl_total_549_2.quant.sf", "Gl_total_549_3.quant.sf",
  "GM1_80S_549_1.quant.sf", "GM1_80S_549_2.quant.sf", "GM1_80S_549_3.quant.sf",
  "GM1_pol_549_1.quant.sf", "GM1_pol_549_2.quant.sf", "GM1_pol_549_3.quant.sf",
  "GM2_80S_549_1.quant.sf", "GM2_80S_549_2.quant.sf", "GM2_80S_549_3.quant.sf",
  "GM2_pol_549_1.quant.sf", "GM2_pol_549_2.quant.sf", "GM2_pol_549_3.quant.sf",
  "GM_total_549_1.quant.sf", "GM_total_549_2.quant.sf", "GM_total_549_3.quant.sf",
  "Tor_80S_549_1.quant.sf", "Tor_80S_549_2.quant.sf", "Tor_80S_549_3.quant.sf",
  "Tor_pol_549_1.quant.sf", "Tor_pol_549_2.quant.sf", "Tor_pol_549_3.quant.sf",
  "Tor_total_549_1.quant.sf", "Tor_total_549_2.quant.sf", "Tor_total_549_3.quant.sf"
)

sample_names <- gsub(".quant.sf", "", files)
names(files) <- sample_names

txi <- tximport(files, type = "salmon", txOut = TRUE)

group_labels <- c(
  rep("0_80S_549", 3), rep("0_pol_549", 3), rep("0_total_549", 3),
  rep("AA_80S_549", 3), rep("AA_pol_549", 3), rep("AA_total_549", 3),
  rep("Gl_80S_549", 3), rep("Gl_pol_549", 3), rep("Gl_total_549", 3),
  rep("GM1_80S_549", 3), rep("GM1_pol_549", 3),
  rep("GM2_80S_549", 3), rep("GM2_pol_549", 3), rep("GM_total_549", 3),
  rep("Tor_80S_549", 3), rep("Tor_pol_549", 3), rep("Tor_total_549", 3)
)

metadata <- data.frame(
  sample = sample_names,
  group = factor(group_labels)
)

cat("Number of samples:", nrow(metadata), "\n")
cat("Unique groups:", length(unique(metadata$group)), "\n")

dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~group)
dds <- dds[rowSums(counts(dds)) >= 1, ]

comparisons <- list(
  c("GM1_80S_549", "0_80S_549"), c("GM2_80S_549", "0_80S_549"),
  c("GM1_80S_549", "AA_80S_549"), c("GM2_80S_549", "AA_80S_549"),
  c("GM1_80S_549", "Gl_80S_549"), c("GM2_80S_549", "Gl_80S_549"),
  c("GM1_80S_549", "Tor_80S_549"), c("GM2_80S_549", "Tor_80S_549"),
  c("GM1_pol_549", "0_pol_549"), c("GM2_pol_549", "0_pol_549"),
  c("GM1_pol_549", "AA_pol_549"), c("GM2_pol_549", "AA_pol_549"),
  c("GM1_pol_549", "Gl_pol_549"), c("GM2_pol_549", "Gl_pol_549"),
  c("GM1_pol_549", "Tor_pol_549"), c("GM2_pol_549", "Tor_pol_549"),
  c("GM_total_549", "0_total_549"), c("GM_total_549", "AA_total_549"),
  c("GM_total_549", "Gl_total_549"), c("GM_total_549", "Tor_total_549")
)

analyze_comparison <- function(dds, group2, group1, name) {
  res <- results(dds, contrast = c("group", group2, group1), alpha = 0.1)
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.1, ]
  
  cat(paste0(name, ":\n\n"))
  summary(res)
  
  res_df <- as.data.frame(res)
  res_df <- cbind(transcript_id = rownames(res_df), res_df)
  write.csv(res_df, paste0(name, "_results_hg19_549.csv"), row.names = FALSE)
}

dds <- DESeq(dds, fitType = "local")

# PCA
vsd <- vst(dds, blind = FALSE, fitType = "local")
pca_plot <- plotPCA(vsd, intgroup = "group")
print(pca_plot)

for (i in seq_along(comparisons)) {
  group1 <- comparisons[[i]][1]
  group2 <- comparisons[[i]][2]
  name <- paste0(group2, "_vs_", group1)
  analyze_comparison(dds, group2, group1, name)
}

# 231

files <- c(
  "0_80S_231_1.quant.sf", "0_80S_231_2.quant.sf", "0_80S_231_3.quant.sf",
  "0_pol_231_1.quant.sf", "0_pol_231_2.quant.sf", "0_pol_231_3.quant.sf",
  "0_total_231_1.quant.sf", "0_total_231_2.quant.sf", "0_total_231_3.quant.sf",
  "AA_80S_231_1.quant.sf", "AA_80S_231_2.quant.sf", "AA_80S_231_3.quant.sf",
  "AA_pol_231_1.quant.sf", "AA_pol_231_2.quant.sf", "AA_pol_231_3.quant.sf",
  "AA_total_231_1.quant.sf", "AA_total_231_2.quant.sf", "AA_total_231_3.quant.sf",
  "Gl_80S_231_1.quant.sf", "Gl_80S_231_2.quant.sf", "Gl_80S_231_3.quant.sf",
  "Gl_pol_231_1.quant.sf", "Gl_pol_231_2.quant.sf", "Gl_pol_231_3.quant.sf",
  "Gl_total_231_1.quant.sf", "Gl_total_231_2.quant.sf", "Gl_total_231_3.quant.sf",
  "GM2_80S_231_1.quant.sf", "GM2_80S_231_2.quant.sf", "GM2_80S_231_3.quant.sf",
  "GM2_pol_231_1.quant.sf", "GM2_pol_231_2.quant.sf", "GM2_pol_231_3.quant.sf",
  "GM_total_231_1.quant.sf", "GM_total_231_2.quant.sf", "GM_total_231_3.quant.sf",
  "Tor_80S_231_1.quant.sf", "Tor_80S_231_2.quant.sf", "Tor_80S_231_3.quant.sf",
  "Tor_pol_231_1.quant.sf", "Tor_pol_231_2.quant.sf", "Tor_pol_231_3.quant.sf",
  "Tor_total_231_1.quant.sf", "Tor_total_231_2.quant.sf", "Tor_total_231_3.quant.sf"
)

sample_names <- gsub(".quant.sf", "", files)
names(files) <- sample_names

txi <- tximport(files, type = "salmon", txOut = TRUE)

group_labels <- c(
  rep("0_80S_231", 3), rep("0_pol_231", 3), rep("0_total_231", 3),
  rep("AA_80S_231", 3), rep("AA_pol_231", 3), rep("AA_total_231", 3),
  rep("Gl_80S_231", 3), rep("Gl_pol_231", 3), rep("Gl_total_231", 3),
  rep("GM2_80S_231", 3), rep("GM2_pol_231", 3), rep("GM_total_231", 3),
  rep("Tor_80S_231", 3), rep("Tor_pol_231", 3), rep("Tor_total_231", 3)
)

metadata <- data.frame(
  sample = sample_names,
  group = factor(group_labels)
)

cat("Number of samples:", nrow(metadata), "\n")
cat("Unique groups:", length(unique(metadata$group)), "\n")

dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~group)
dds <- dds[rowSums(counts(dds)) >= 1, ]

comparisons <- list(
  c("GM2_80S_231", "0_80S_231"), c("GM2_80S_231", "AA_80S_231"), 
  c("GM2_80S_231", "Gl_80S_231"), c("GM2_80S_231", "Tor_80S_231"),
  c("GM2_pol_231", "0_pol_231"), c("GM2_pol_231", "AA_pol_231"),
  c("GM2_pol_231", "Gl_pol_231"), c("GM2_pol_231", "Tor_pol_231"),
  c("GM_total_231", "0_total_231"), c("GM_total_231", "AA_total_231"),
  c("GM_total_231", "Gl_total_231"), c("GM_total_231", "Tor_total_231")
)

analyze_comparison <- function(dds, group2, group1, name) {
  res <- results(dds, contrast = c("group", group2, group1), alpha = 0.1)
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.1, ]
  
  cat(paste0(name, ":\n\n"))
  summary(res)
  
  res_df <- as.data.frame(res)
  res_df <- cbind(transcript_id = rownames(res_df), res_df)
  write.csv(res_df, paste0(name, "_results_hg19_231.csv"), row.names = FALSE)
}

dds <- DESeq(dds, fitType = "local")

# PCA
vsd <- vst(dds, blind = FALSE, fitType = "local")
pca_plot <- plotPCA(vsd, intgroup = "group")
print(pca_plot)

for (i in seq_along(comparisons)) {
  group1 <- comparisons[[i]][1]
  group2 <- comparisons[[i]][2]
  name <- paste0(group2, "_vs_", group1)
  analyze_comparison(dds, group2, group1, name)
}

# ===========================
# annotation
# ===========================

gtf_file <- "gencode.v19.annotation.gtf"
input_folder <- "./hg19/"  # input csv
output_folder <- "./hg19/csv/"  # output csv

gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

transcripts <- gtf_df %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id, gene_name)

csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

for (csv_file in csv_files) {
  cat("Processing file:", csv_file, "\n")
  
  input_data <- read.csv(csv_file, header = TRUE, sep = ",", quote = "\"", fill = TRUE, comment.char = "")
  
  if (!"transcript_id" %in% colnames(input_data)) {
    cat("Skipping file, 'transcript_id' column is missing:", csv_file, "\n")
    next
  }
  
  annotated_data <- input_data %>%
    inner_join(transcripts, by = "transcript_id") %>%
    select(transcript_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
  
  output_file <- file.path(output_folder, paste0(basename(tools::file_path_sans_ext(csv_file)), "_annotated.csv"))
  
  write.csv(annotated_data, output_file, row.names = FALSE, quote = FALSE)
  cat("Annotated file saved as:", output_file, "\n")
}


# ===========================
# MT genes
# ===========================

input_folder <- "./hg19/csv_annotated/"  # input csv
output_folder <- "./hg19/csv_annotated_mt_genes/"  # output csv

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

input_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

for (file in input_files) {
  data <- read.csv(file, stringsAsFactors = FALSE)
  filtered_data <- subset(data, grepl("^MT-", gene_name))
  if (nrow(filtered_data) > 0) {
    output_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(file)), "_mt.csv"))
    write.csv(filtered_data, output_file, row.names = FALSE)
    cat("File saved:", basename(output_file), "\n")
  } else {
    cat("File missing (no MT genes):", basename(file), "\n")
  }
}

cat("All files processed!\n")
