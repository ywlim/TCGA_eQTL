#!/usr/bin/env Rscript
require(getopt)
opt <- getopt(matrix(c(
        'geno','g',1,"character",
	    'expr','e',1,"character",
	    'covar','c',1,"character",
	    'snploc','s',1,"character",
	    'geneloc','l',1,"character",
	    'cisP','u',1,"numeric",
	    'traP','v',1,"numeric",
	    'dist','d',1,"numeric",
	    'cnv', 'n',1,"character",
	    'methyl','m',1,"character",
	    'outprefix','o',1,"character")
	,ncol=4,byrow=T));
print(opt);

# opt <- list();
# opt$outprefix = 'COAD'
# opt$cisP = 1e-5
# opt$traP = 0
# opt$dist = 1e6


# U: pvOutputThreshold.cis  (cis-pvalue) 
# V: pvOutputThreshold
# 
# Cis-mode:
# U = 1e-5, V = 0
# pvOutputThreshold.cis = 1e-5, 
# pvOutputThreshold = 0
# 
# Trans-mode
# U = 0, V = 1e-10
# pvOutputThreshold.cis = 0
# pvOutputThreshold = 1e-10

library(MatrixEQTL)

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR <- CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR <- CROSS

SNP_file_name <- opt$geno;

expression_file_name <- opt$expr;

covariates_file_name <- opt$covar;

cnv_file_name <- opt$cnv;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Output file name
output_file_name_cis = paste(opt$outprefix,".cis.txt",sep="");
output_file_name_tra = paste(opt$outprefix,".tra.txt",sep="");

# Only associations significant at this level will be saved
if(exists("cisP",opt)){
    pvOutputThreshold_cis = opt$cisP;
}else{
    pvOutputThreshold_cis = 1e-2;
}
if(exists("traP",opt)){
    pvOutputThreshold_tra = opt$traP;
}else{
    pvOutputThreshold_tra = 2e-2;
}

# Distance for local gene-SNP pairs
if(exists("dist",opt)){
    cisDist <- opt$dist;
}else{
    cisDist = 1e6;
}

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 200000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Load copy number variations

cnv = SlicedData$new();
cnv$fileDelimiter = "\t";      # the TAB character
cnv$fileOmitCharacters = "NA"; # denote missing values;
cnv$fileSkipRows = 1;          # one row of column labels
cnv$fileSkipColumns = 1;       # one column of row labels
cnv$fileSliceSize = 2000;      # read file in slices of 2,000 rows
if(length(covariates_file_name)>0) {
    cnv$LoadFile(cnv_file_name);
}

## Run the analysis
genepos = read.table(opt$geneloc, header = TRUE, stringsAsFactors = FALSE);

library(data.table)
snpspos = fread(opt$snploc, header = TRUE, stringsAsFactors = FALSE);
snpspos = as.data.frame(snpspos)

stopifnot (cnv$fileSliceSize == gene$fileSliceSize)
# cnv and gene must be harmonized row wise
stopifnot(rownames(as.matrix(cnv))==rownames(as.matrix(gene)))

# change path to Matrix_eQTL_mod_partition_SS.R here 
source("Matrix_eQTL_mod_partition_SS.R")

me = Matrix_eQTL_main_mod(
    snps = snps,
    gene = gene,
    cnv = cnv,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

pdf(paste(opt$outprefix,".pval.pdf",sep=""));
plot(me)
dev.off();


