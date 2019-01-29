#!/usr/bin/env Rscript
# harmonizes samples in meg, mee and cnv

# Scipt makes sure same sample IDs used for everything
# Also gets rid of samples for which either expression or SNP array data isn't available
# meg is SNP array data, mee is expression data, cnv is copy number data

require(getopt)

opt <- getopt(matrix(c("mee", "e", 1, "character",
	  		   "meg", "g", 1, "character",
			   "cnv", "v", 1, "character",
		       "dicTable", "d", 1, "character",
		       "outMee", "o", 1, "character",
		       "outMeg", "p", 1, "character",
			   "outCnv", "q", 1, "character",
		       "gdsFile", "h", 1, "character",
		       "outCov", "c", 1, "character"),
		       ncol=4,
		       byrow=TRUE))
print(opt)

mee <- read.table(opt$mee, header=T, row.names=1)

# read in cnv file
cnv <- read.table(opt$cnv, header=T, row.names=1)

# Read in .meg whilst ignoring rows with duplicate row names (e.g. ".")
meg <- read.table(opt$meg, header=T)
uniqRsIds <- names(which(table(meg[, 1]) == 1))
meg <- meg[which(meg[, 1] %in% uniqRsIds), ]
rownames(meg) <- meg[, 1]
meg <- meg[, 2:dim(meg)[2]]

cat("Read in", opt$meg, "with", dim(meg)[1], "unique SNPs and", dim(meg)[2], "samples\n")

# Read in table on how different sample IDs correspond with each other
# Required format: expression_ID \t genotype_ID
dicTable <- read.table(opt$dicTable, fill=T, header=F, stringsAsFactors=FALSE)
dicTable <- subset(dicTable, V2 != "")
dic <- as.character(dicTable$V2)
names(dic) <- as.character(dicTable$V1)

# R substitutes "-" with ".". Reverse this
# Names starting with numbers are prefixed by "X" by R. Reverse this
names(meg) <- gsub("^X", "", names(meg))
names(meg) <- gsub("\\.", "-", names(meg))

names(mee) <- gsub("^X", "", names(mee))
names(mee) <- gsub("\\.", "-", names(mee))

names(cnv) <- gsub("^X", "", names(cnv))
names(cnv) <- gsub("\\.", "-", names(cnv))

# Subset data so that only samples in dic are kept
if (length(intersect(names(mee), names(dic))) == 0) { cat("Problem: No overlap between dic and mee names\n") }
if (length(intersect(names(meg), dic)) == 0) { cat("Problem: No overlap between dic and meg names\n") }
if (length(intersect(names(cnv), dic)) == 0) { cat("Problem: No overlap between dic and cnv names\n") }

mee <- mee[, intersect(names(dic), names(mee))]
meg <- meg[, intersect(names(meg), dic)]
cnv <- cnv[, intersect(names(cnv), dic)]

# Translate mee sample IDs
colnames(mee) <- dic[colnames(mee)]
cat("Sample IDs in RNAseq data changed so that they are the same than in genotyping data.\n")

# Get overlap
bothSets <- intersect(colnames(mee), colnames(meg))
threeSets <- intersect(colnames(cnv), bothSets)

# Keep only columns in intersect
mee <- mee[, threeSets]
meg <- meg[which(!is.na(rownames(meg))), threeSets]
cnv <- cnv[, threeSets]
cat("Removed samples that don't have all genotyping and RNAseq, and CNV data.", length(threeSets), "samples remaining\n")

# Write the files
rNames <- rownames(mee)
write.table(cbind(rNames, mee), file=opt$outMee,
	    quote=F, sep="\t", row.names=FALSE, col.names=c("id", threeSets))
cat("Wrote cleaned RNAseq data to", opt$outMee, "\n")

rNames <- rownames(meg)
write.table(cbind(rNames, meg), file=opt$outMeg,
	    quote=F, sep="\t", row.names=FALSE, col.names=c("id", threeSets))
cat("Wrote cleaned genotyping data to", opt$outMeg, "\n")

rNames <- rownames(cnv)
write.table(cbind(rNames, cnv), file=opt$outCnv,
            quote=F, sep="\t", row.names=FALSE, col.names=c("id", threeSets))
cat("Wrote cleaned cnv data to", opt$outCnv, "\n")

# Convert to GDS and do PCA
require(gdsfmt)
require(SNPRelate)

# Get 150k random SNPs
rand150k <- as.matrix(meg[sample(dim(meg)[1], 150000), ])

snpgdsCreateGeno(opt$gdsFile, 
		 genmat = rand150k,
		 sample.id = colnames(rand150k),
		 snp.id = rownames(rand150k),
		 snpfirstdim = TRUE,
		 compress.annotation = "")
cat("Created GDS file with 150,000 random sites for PCA in", opt$gdsFile, "\n")

# Read in GDS
gds <- openfn.gds(opt$gdsFile)

# Do PCA and print EVs
pca <- snpgdsPCA(gds)

evs <- t(pca$eigenvect)
evs <- evs[1:3, ] # Only the first three EVs likely to be relevant
rNames <- paste("EV", 1:dim(evs)[1], sep=".")
cNames <- pca$sample.id
write.table(cbind(rNames, evs), file=opt$outCov, quote=F, sep="\t", row.names=FALSE, col.names=c("id", cNames))
cat("Wrote eigenvectors to file", opt$outCov, "\n")
cat("Complete.\n")
