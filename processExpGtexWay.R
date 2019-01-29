#!/usr/bin/env Rscript
# Process expression file in way recommended by GTEx for eQTL

require(getopt)
require(preprocessCore)

opt <- getopt(matrix(c("inTable", "i", 1, "character",
		       "minRpkm", "r", 1, "numeric",
		       "maxSmp", "m", 1, "numeric",
		       "outFile", "o", 1, "character"),
		     ncol=4,
		     byrow=TRUE))
print(opt)

# Load data
p.rpkm <- read.table(opt$inTable, sep="\t", header=T, row.names=1, check.names=FALSE)

# Processing, following GTEx guidelines:

# GTEx recommends filtering on >= 10 individuals having >0.1 RPKM
ExpSamples <- apply(p.rpkm, 1, function(tt) { sum(tt >= opt$minRpkm) })
p.rpkm <- p.rpkm[which(ExpSamples >= opt$maxSmp), ]

cat("Genes with expression of <=", opt$minRpkm, "in >=", opt$maxSmp, "samples filtered out\n")
cat("Number of gene IDs passing filtering:", nrow(p.rpkm), "\n")

# Log2
p.rpkm <- log2(p.rpkm)
geneIds <- rownames(p.rpkm)
sampleIds <- names(p.rpkm)

# Values with RPKM of 0 are set to -Inf in log2 step
# To avoid problems when scaling, they are set to minimum observed value above 0
minObs <- min(p.rpkm[p.rpkm > -Inf])
p.rpkm[p.rpkm == -Inf] <- minObs

# Normalize across all samples
# p.rpkm <- (p.rpkm - mean(unlist(p.rpkm))) / sd(unlist(p.rpkm))
q.rpkm <- normalize.quantiles(as.matrix(p.rpkm))  

# For outlier correction, perform Inverse Normal Transformation: rank values across samples then map to a standard normal
# Inverse Normal Transformation for each gene
q.rpkm <- t(apply(q.rpkm, 1, function(tt) qnorm( (rank(tt)-0.5)/length(tt))))

q.rpkm <- as.data.frame(q.rpkm)
names(q.rpkm) <- sampleIds
ready.rpkm <- cbind(geneIds = geneIds, q.rpkm)

# Print output with genes as rows and samples as columns
write.table(ready.rpkm, file=opt$outFile, quote=F, sep="\t", row.names=FALSE)
