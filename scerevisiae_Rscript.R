# load the essential libraries
library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")
library("vioplot")
library("RColorBrewer")
library("biomaRt")

# Setup the Directory 
setwd("D:/R/Pan_NGS_data/Rna_seq/genome")

database <- makeTxDbFromGFF("Saccharomyces_cerevisiae.R64-1-1.110.gtf", 
                      format = "gtf", organism = "Saccharomyces cerevisiae", 
                      dataSource = "https://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/")


# Get the exons per gene, and compute bp lengths of all genes
exons <- exonsBy(database, by = "gene")
gene.lengths <- lapply(exons, function(x){ sum(width(reduce(x))) })

setwd("D:/R/Pan_NGS_data/Rna_seq/Results")

# Define the samples and their read IDs (SRA run accessions)
samples_name <- c("AMB + LF rep 1", "AMB + LF rep 2", "AMB + LF rep 3","control 1", "control 2", "control 3")
names(samples_name) <- c( "SRR3396388", "SRR3396389", "SRR3396390", "SRR3396391", "SRR3396392", "SRR3396392")

# Construct the file path - 
files <- c()

for (s in samples) {
  fp <- paste0(s, ".bam")
  if (file.exists(fp)) files <- c(files, fp)
}

# Load in the BAM files
bams_files <- BamFileList(files, yieldSize = 100000, asMates=TRUE)
# Overlap BAM reads and genes
overlap <- summarizeOverlaps(exons, bams_files, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
overlap


# Perform PCA on log2-transformed RPKM data
library(ggplot2)
library(ggfortify)

pca_data <- prcomp(t(geneRPKM_log2))
autoplot(pca_data, data = overview, colour = "AMB + LF")





# Extract the raw-reads per gene
geneRawReads <- assay(overlap)
geneRawReads[1:10,]
colnames(geneRawReads) <- gsub(".bam", "", colnames(geneRawReads))
write.table(geneRawReads, "gene_raw_reads.txt", sep = "\t", quote = FALSE)

# Calculate the RPKM values per gene
# RPKM = (10^9 * C)/(N * L)
# C = Number of reads mapped to a gene
# N = Total mapped reads in the sample
# L = gene length in base-pairs for a gene

# Get the total number of reads per samples
totalMappedReads <- apply(geneRawReads, 2, sum)

# Loop through all genes, compute RPKM
geneIdx <- 1
geneRPKM <- t(apply(geneRawReads, 1, function(C){
  geneLength <- as.numeric(gene.lengths[geneIdx])
  RPKM <- (10 ^ 9 * C) / (totalMappedReads * geneLength)
  geneIdx <<- geneIdx + 1
  return(round(RPKM, d = 1))
}))

geneRPKM[1:10,]

op <- par(mar = c(8,4,2,2))

#Violin distribution plot
# Define the sample names and colors
sample_names <- c("Control", "AMB + LF")
sample_colors <- c("red", "yellow")

# Create a vioplot for RPKM with the updated sample names and colors
vioplot(geneRPKM, col = sample_colors, ylab = "Reads", las = 2)
legend("topright", inset = c(-0.07, -0.12), legend = sample_names, fill = sample_colors, xpd = TRUE)


write.table(geneRPKM, "RPKM.txt", sep = "\t", quote = FALSE)

# Quantile normalization of RPKM values
geneRPKM_norm <- round(normalize.quantiles(as.matrix(geneRPKM)), d = 1)
colnames(geneRPKM_norm) <- colnames(geneRPKM)
rownames(geneRPKM_norm) <- rownames(geneRPKM)

par(mar = c(8,4,2,2))
# Violin distribution plot
vioplot(geneRPKM_norm, col = sample_colors, ylab = "Normalized", las = 2)
legend("topright", inset = c(-0.07, -0.14), legend = sample_names, fill = sample_colors, xpd = TRUE)

write.table(geneRPKM_norm, "gene_RPKM_normalized.txt", sep = "\t", quote = FALSE)

# Log2 transform normalized RPKM (microarray-like)
geneRPKM_log2 <- round(log2(geneRPKM_norm), d = 1)
geneRPKM_log2[geneRPKM_log2 < 0] <- 0

par(mar = c(8,4,2,2))
# Violin distribution plot
vioplot(geneRPKM_log2, col = sample_colors, ylab = "log2(Normalized)", las = 2)
legend("topright", inset = c(-0.07, -0.14), legend = sample_names, fill = sample_colors, xpd = TRUE)

write.table(geneRPKM_log2, "gene_RPKM_normalized_log2.txt", sep = "\t", quote = FALSE)

plot(geneRPKM_log2[,1], geneRPKM_norm[,1])
plot(geneRPKM_log2[,1], geneRPKM[,1])
plot(geneRPKM_log2[,1], readcount[,1])
plot(geneRPKM_log2[,1], geneRPKM_norm[,1])

hist(geneRPKM_log2[,1])

# P-values and Log2 fold change
# Compute p-values and fold changes for each individual sample
pvals <- apply(geneRPKM_log2, 1, function(x){
  tryCatch(t.test(x[1:3], x[4:6])$p.value, 
           error = function(x){return(NA);})
})
plot(pvals)
hist(pvals)

fc <- apply(geneRPKM_log2, 1, function(x){
  tryCatch(log2(mean(x[1:3]) / mean(x[4:6])), 
           error = function(x){return(NA);})
})

hist(fc)

## Assign colors based on P-values
colz <- rep("black", length(pvals))
colz[which(pvals < 5e-2)] <- "red"
colz[which(pvals < 1e-2)] <- "gold"
colz[which(pvals < 1e-3)] <- "blue"

# Volcano plot (x = fc, y = -log10(P-values))
plot(fc, -log10(pvals), col=colz, pch=18, main ="Volcano plot", xlab="Fold Change")
legend("topleft", pch=18, c("<0.05", "<0.01", "<0.001"), col = c("red", "gold", "blue"))

# Down & Up regulated genes
down <- geneRPKM_log2[which(pvals < 5e-2 & fc < -0.3),]
dclust <- down[hclust(dist(down))$order,]
up <- geneRPKM_log2[which(pvals < 5e-2 & fc > 0.3),]
uclust <- up[hclust(dist(up))$order,]

# Gene IDs of up/Down regulated genes
geneIDs <- c(rownames(dclust), rownames(uclust))

write.table(geneIDs, "gene.txt", sep = "\t", quote = FALSE)

## Custom heatmap using the spectral colors
heatmap_data <- rbind(dclust, uclust)
heatmap_colnames <- colnames(geneRPKM_log2)[1:6]

heatmap(x = heatmap_data, Rowv = NA, Colv = NA, 
        col = brewer.pal(11, "Blues"), main = "Differentially Expressed Genes Heatmap",
        xlab = "Samples", ylab = "Genes", margins = c(8, 6))

axis(2, at = 1:length(geneIDs), labels = geneIDs, las = 2, cex.axis = 0.7)
axis(1, at = 1:3, labels = heatmap_colnames, las = 2, col.axis = "black")
axis(1, at = 3:6, labels = heatmap_colnames, las = 2, col.axis = "black")

# Compute mean expression and standard deviations
means <- t(apply(geneRPKM_log2, 1, function(x){
  tryCatch(round(c(mean(x[1:3]),mean(x[4:6])),1), error = function(x){return(NA);})
}))

sds <- t(apply(geneRPKM_log2, 1, function(x){
  tryCatch(round(c(sd(x[1:3]),sd(x[4:6])),1), error = function(x){return(NA);})
}))


colnames(means) <- c("AMB + LF", "Control")
colnames(sds) <- c("AMB + LF", "Control")

# Create an overview table
overview <- cbind("Control" = means[, "Control"], 
                  "Control(SD)" = sds[, "Control"], 
                  "AMB + LF" = means[, "AMB + LF"], 
                  "AMB + LF(SD)" = sds[, "AMB + LF"], 
                  FC = round(fc,1),
                  P = round(pvals,6))
overview[1:10,]

# Load the biomaRt library

library(biomaRt)

bio.mart <- useMart("ensembl", "scerevisiae_gene_ensembl")



mattr <- c("ensembl_gene_id", "external_gene_name", 
           "chromosome_name", "start_position", "end_position",
           "description")

res.bm <- getBM(attributes = mattr, 
                filters = c("ensembl_gene_id"), 
                values = geneIDs, mart = bio.mart)
rownames(res.bm) <- res.bm[, "ensembl_gene_id"]

# Merge biomaRt results with the overview
p1 <- res.bm[geneIDs, c("external_gene_name", "chromosome_name", "start_position", "end_position")]
overview <- cbind(p1, overview[geneIDs,], res.bm[geneIDs, "description"])
colnames(overview)[1:4] <- c("GeneName", "Chr", "Start", "End")
colnames(overview)[11] <- c("Description")

overview[1:10,1:10]

# Write out the table
write.table(overview, "overview_annotation.txt", sep = "\t", quote = FALSE)