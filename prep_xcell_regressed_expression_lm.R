# take RNA-seq, transform to log2, regress on 34 xcell gene signatures one by one
# then reverse log2, and write new expression values 
# so that PEERs can be estimated from them

library(dplyr)
library(tidyr)
library(readr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
indi <- args[1]
suffix <- args[2]

#exp_file <- paste0("/gne/research/data/dnaseq/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/pheno/hg38/", indi, ".hg38.eset.tumor.nrpkm.tsv")
# because running the entire rna-seq is too slow, i broke rna-seq into mini files that i can run in parallel
exp_file <- paste0("/gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/genesets/xcell_genesets_regressed_peers/Mini_exp/", indi, "_split_", suffix)
exp <- read_delim(exp_file, delim = "\t")
xcell_file <- paste0("/gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/genesets/xcell_genesets_lm/", indi, "/", indi, "_geneset_scores.txt")
xcell <- read_delim(xcell_file, delim = "\t")

exp_t <- data.frame(t(exp[,-1]))
# add g_ in front of geneId for easy selection later
names(exp_t) <- paste0("g_", exp$geneId)
exp_t <- rownames_to_column(exp_t, var = "sample")

xcell_t <- data.frame(t(xcell[,-1]))
names(xcell_t) <- xcell$X1
xcell_t <- rownames_to_column(xcell_t, var = "sample")

df <- inner_join(exp_t, xcell_t, by = "sample")

# log2 transform to linearize data (add a small pseudocount)
df_log <- cbind.data.frame(df$sample, log2(df[,-1] + 0.00001))
names(df_log) <- names(df)

all_genesets <- unlist(xcell$X1)
all_genes <- names(exp_t[,-1])

sample <- df_log$sample

residual2 <- lapply(all_genes, function(x) {
  gene <- x
  print(gene)
  
  mini <- select(df_log, gene, all_genesets)
  colnames(mini)[1] <- "gene"
  
  # for loop instead of lapply
  # loop all genesets, regressing one out at a time
  for (j in 1:length(all_genesets)) {
    geneset <- all_genesets[j]
    mini_slim <- select(mini, 1, one_of(geneset))
    names(mini_slim) <- c("gene", "geneset")
    model <- lm(gene ~ geneset, data = mini_slim)
    # replace original expression with residual
    mini$gene <- resid(model)
  }
  
  # return the final expression for the gene, after 34 rounds of regression
  return(mini$gene)
})

new <- as.data.frame(do.call(cbind, residual2))
names(new) <- all_genes

# sanity check
# plot(df_log$g_1, df_log$Smooth_muscle)
# plot(new$g_1, df_log$Smooth_muscle)

# remove log2
new_raw <- 2^new
new_raw <- cbind.data.frame(sample, new_raw)

new_raw_t <- as.data.frame(t(new_raw[,-1]))
names(new_raw_t) <- new_raw$sample
new_raw_t <- rownames_to_column(new_raw_t, var = "geneId")

# remove g_ from geneids
new_raw_t$geneId <- sub("g_", "", new_raw_t$geneId)

outfile = paste0("/gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/genesets/xcell_genesets_lm/Mini_exp_regressed/", indi, "_out_", suffix)
write.table(new_raw_t, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)




