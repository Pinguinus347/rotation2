if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

BiocManager::install("DESeq2")

BiocManager::install("statmod")

# Attempting DESeq2 pipeline

library("DESeq2")

# List all files (adjust pattern if needed)
files <- list.files(path = "/home/callum/rotation2/transcription_data/raw_contact", pattern = "*.txt", full.names = TRUE)
head(files)
# Read each file and extract the counts
# Assume gene IDs are in the first column, counts in the second
count_list <- lapply(files, function(file) {
  read.table(file, header = FALSE, sep = "", stringsAsFactors = FALSE)
})
head(count_list)
# Merge all data frames by the first column (gene names)
merged_counts <- Reduce(function(x, y) merge(x, y, by = 1), count_list)
head(merged_counts)
# Set gene names as rownames
rownames(merged_counts) <- merged_counts[, 1]
merged_counts <- merged_counts[, -1]  # remove first column

# Optionally rename columns based on filenames
colnames(merged_counts) <- gsub("\\.txt$", "", basename(files))
head(merged_counts)

cts <- as.matrix(merged_counts)

head(count_matrix, n = 50)

# Getting condition information
pasAnno <- read.csv("/home/callum/rotation2/transcription_data/raw_contact/conditions.csv")
# Renaming rows for sample name, and removing this column
rownames(pasAnno) <- pasAnno[,1]
coldata <- pasAnno[,-1]
head(coldata)

all(rownames(coldata) == colnames(cts))


# Initialising DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

dds
# Load edgeR
library(edgeR)

# List all files (adjust pattern if needed)
files <- list.files(path = "/home/callum/rotation2/transcription_data/raw_contact", pattern = "*.txt", full.names = TRUE)
head(files)
# Read each file and extract the counts
# Assume gene IDs are in the first column, counts in the second
count_list <- lapply(files, function(file) {
  read.table(file, header = FALSE, sep = "", stringsAsFactors = FALSE)
})
head(count_list)
# Merge all data frames by the first column (gene names)
merged_counts <- Reduce(function(x, y) merge(x, y, by = 1), count_list)
head(merged_counts)
# Set gene names as rownames
rownames(merged_counts) <- merged_counts[, 1]
merged_counts <- merged_counts[, -1]  # remove first column

# Optionally rename columns based on filenames
colnames(merged_counts) <- gsub("\\.txt$", "", basename(files))
head(merged_counts)

count_matrix <- as.matrix(merged_counts)

head(count_matrix)

group <- factor(c("filter","filter", "filter", "liquid", "liquid", "liquid"))  # modify as needed


# 3. Create DGEList object
dge <- DGEList(counts = count_matrix, group = group)
head(dge)

# Removing genes which are not expressed, I think min default is 10
keep <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]
# 3. Normalize with TMM
dge <- calcNormFactors(dge, method = "TMM")
dge$samples

plotMDS(dge)
# Batch effects are not mentioned in methods or metadata
library(statmod)
# Estimate dispersion, without design optimisation
# Making robust true doesnt seem to change the common dispersion or BCV
dge <- estimateDisp(dge,robust=FALSE)
dge$common.dispersion
plotBCV(dge)

fit <- glmQLFit(dge)
plotQLDisp(fit)

qlf <- glmQLFTest(fit)
topTags(qlf)

top <- rownames(topTags(qlf))
cpm(dge)[top,]
summary(decideTests(qlf))

plotMD(qlf)
abline(h = 0, col = "red")

# Get results
results <- topTags(qlf, n = Inf)$table
head(results)

write.csv(results, file = "edgeR_results.csv", row.names = TRUE)


# 4. Estimate dispersions (i think same as doing this in one step)
dge <- estimateCommonDisp(dge)     # common dispersion
dge <- estimateTagwiseDisp(dge)    # tagwise (per-gene) dispersion



# 5. Differential expression testing
et <- exactTest(dge)
topTags(et)
# 6. Get top results
results <- topTags(et, n = Inf)$table
head(results)

write.csv(results, file = "edgeR_results.csv", row.names = TRUE)

# 4. Normalize using TMM (default in edgeR)
dge <- calcNormFactors(dge, method = "TMM")

# 5. Estimate dispersion (empirical Bayes)
dge <- estimateDisp(dge)

# 6. Fit a negative binomial GLM model
fit <- glmFit(dge)

# 7. Conduct likelihood ratio test
lrt <- glmLRT(fit)
head(lrt)
# 8. Extract top differentially expressed genes
results <- topTags(lrt, n = Inf)$table
head(results)
# 9. Filter for FDR < 0.05
significant_genes <- results[results$FDR < 0.05, ]
results <- topTags(significant_genes, n = 50)$table
head(significant_genes, n = 50)



# 10. Further filter by p-value < 0.001 (optional, stricter)
sig_genes_strict <- significant_genes[significant_genes$PValue < 0.001, ]

library(ggplot2)

# Add -log10 FDR and significance flags
results$negLogFDR <- -log10(results$FDR)
results$significant <- with(results, FDR < 0.05)  # adjust threshold as needed

ggplot(results, aes(x = -logFC, y = negLogFDR)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10 FDR"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

library(ggrepel)

top_genes <- results[results$significant, ][1:50, ]  # top 10 significant
print(top_genes)
ggplot(results, aes(x = -logFC, y = negLogFDR)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes))) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal()
