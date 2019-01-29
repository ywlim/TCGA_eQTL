#!/usr/bin/env Rscript

# Run PEER (https://www.sanger.ac.uk/resources/software/peer/) on expression data
# The input to PEER are the post-normalization expression values

require(getopt)
require(peer)

opt <- getopt(matrix(c("inExp", "i", 1, "character",
		       "outCov", "o", 1, "character"),
		       ncol=4, byrow=TRUE))
print(opt)

numFactors <- 15

expr <- read.table(opt$inExp, row.names=1, header=TRUE)
expr <- t(expr)

model <- PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setNk(model, numFactors)
PEER_getNk(model)
PEER_update(model)
factors <- PEER_getX(model)

factors <- t(factors)
colnames(factors) <- rownames(expr)
factors <- cbind(paste("PEER", 1:numFactors,sep="_"), factors)

write.table(factors, opt$outCov, row.names=F, quote=F, sep="\t")
