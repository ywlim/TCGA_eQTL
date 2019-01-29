#!/usr/bin/env Rscript
# additional preparation step to process cnv files
# after the "make normalize" step
# so that the genes are the same in the mee and cnv files

require(getopt)

opt <- getopt(matrix(c("mee", "e", 1, "character",
                        "cnv", "v", 1, "character",
			"outCnv", "q", 1, "character"),
			ncol=4,
			byrow=TRUE))
print(opt)

mee <- read.table(opt$mee, header=T, row.names=1)
cnv <- read.table(opt$cnv, header=T, row.names=1, check.names=FALSE)
names(cnv) <- gsub("\\.", "-", names(cnv))

# arrange genes in cnv file so that it's the same as the order in mee file
library(dplyr)
key <- as.data.frame(rownames(mee))
names(key) <- "id"
cnv <- mutate(cnv, id = rownames(cnv))
cnv <- left_join(key, cnv, by = "id") %>% select(starts_with("TCGA"))
rownames(cnv) <- key$id

# Write the file

rNames <- as.data.frame(rownames(cnv))
names(rNames) <- "id"
write.table(cbind(rNames, cnv), file=opt$outCnv, , quote=F, sep="\t", row.names=FALSE)
