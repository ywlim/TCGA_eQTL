#!/usr/bin/env Rscript

# Add covariate other that EVs

require(getopt)

opt <- getopt(matrix(c("covarIn", "i", 1, "character",
		       "covarOut", "o", 1, "character",
		       "addCovarFile", "a", 1, "character",
		       "sampIdHeader", "s", 1, "character",
		       "addCovar", "l", 1, "character"),
		     ncol=4,
		     byrow=T))

print(opt)

# opt$covarIn <- "/gnet/is6/HumGenet/analysis/HumGenet_dbGap/GTEx/eQTLs/tissues/stomach/stomach.evOnly.cov"
# opt$covarOut <- "/gnet/is6/HumGenet/analysis/HumGenet_dbGap/GTEx/eQTLs/tissues/stomach/stomach.cov"
# opt$addCovarFile <- "/gnet/is6/HumGenet/analysis/HumGenet_dbGap/GTEx/eQTLs/pheno/gtexMetadata.tsv"
# opt$sampIdHeader <- "SUBJID"
# opt$addCovar <- "sex"

covarIn <- opt$covarIn
covarOut <- opt$covarOut
addCovarFile <- opt$addCovarFile
addCovar <- strsplit(opt$addCovar, ",")[[1]] #take comma separate input
message("addCovar1 is: ", addCovar[1])
message("addCovar2 is: ", addCovar[2])

covars <- read.table(covarIn, header=T, row.names=1, sep="\t",check.names = FALSE)
colnames(covars) <- sub("\\.", "-", colnames(covars))  # change TCGA-string.string to TCGA-string-string
dic <- read.table(addCovarFile, header=T, sep="\t")

toAddCols <- as.character(c(opt$sampIdHeader, addCovar))
toAdd <- dic[which(dic[, opt$sampIdHeader] %in% colnames(covars)), toAddCols]
toAdd <- unique(toAdd)
rownames(toAdd) <- toAdd[, opt$sampIdHeader]

# set gender factor to numeric value
toAdd$gender <- plyr::mapvalues(toAdd$gender, from=levels(toAdd$gender), to=c(1,2))

covars <- rbind(covars, t(toAdd[colnames(covars), addCovar, drop=FALSE]))
# Get rid of samples where no covars could be found

			  
rNames <- rownames(covars)
write.table(cbind(rNames, covars), file=opt$covarOut,
	                quote=F, sep="\t", row.names=FALSE, col.names=c("id", colnames(covars)))
