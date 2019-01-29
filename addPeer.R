#!/usr/bin/env Rscript

# Add PEER factors to EVs

require(getopt)

opt <- getopt(matrix(c("covarIn", "i", 1, "character",
		       "covarOut", "o", 1, "character",
		       "peerIn", "a", 1, "character"),
		     ncol=4,
		     byrow=T))

print(opt)

covars <- read.table(opt$covarIn, header=T, row.names=1, sep="\t")
colnames(covars) <- sub("\\.", "-", colnames(covars))

peers <- read.table(opt$peerIn, header=T, row.names=1, sep="\t")
colnames(peers) <- sub("\\.", "-", colnames(peers))

commonSamples <- intersect(colnames(covars), colnames(peers))
covars <- rbind(covars[, commonSamples], peers[, commonSamples])
			  
rNames <- rownames(covars)
write.table(cbind(rNames, covars), file=opt$covarOut,
	                quote=F, sep="\t", row.names=FALSE, col.names=c("id", colnames(covars)))
