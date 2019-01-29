# this script runs mediation test to see if a candidate gene 
# is mediating the association between a SNP and a gene signature
# it takes arguments so that i can run jobs in parallel
# the args should be file name in RDS and a number

library(tidyr)
library(ggplot2)
library(purrr)
library(mediation)
library(readr)
library(dplyr)

select <- dplyr::select

args = commandArgs(trailingOnly=TRUE)

igis <- read_delim("/gne/research/workspace/limy19/TCGAgermline/trunk/data/IGIS_hg38_id_symbol.tsv", delim = "\t")
igis <- igis %>% rename(id = Entrez)

# previously prepared RDS file
xcell <- readRDS(args[1])

med_clean <- function(mmm) {
  # tidy and extract results from mediation analysis
  mmm %>%
    capture.output() %>%
    discard(`==`, "") -> lines
  
  med_df <- lines[which(grepl("^ ", lines)):(which(grepl("^---", lines))-1)] %>%
    sub("^       ", "med", .) %>%
    gsub("95% ", "95_", .) %>%
    gsub("CI ", "ci_", .) %>%
    sub("p-", "p_", .) %>%
    sub("Total Effect", "Total_effect", ., fixed=TRUE)  %>%
    sub("Prop. Mediated", "Prop_mediated", ., fixed=TRUE)  %>%
    #sub("p_value", "p_value sig", ., fixed=TRUE)  %>%
    #gsub("*", "", ., fixed=TRUE) %>%
    textConnection() %>%
    read.table(header=FALSE, fill = TRUE, stringsAsFactors = FALSE) 
  # critical fill = true to avoid error when there is uneven column number
  
  # total effect of genotype on xcell score
  tot_effect_b <- med_df$V2[4]
  tot_effect_ci_lo <- med_df$V3[4]
  tot_effect_ci_hi <- med_df$V4[4]
  tot_effect_p <- med_df$V5[4]
  
  # direct effect of genotype on xcell score (after including potential mediator gene) - should become insignificant if gene is mediator
  dir_effect_b <- med_df$V2[3]
  dir_effect_ci_lo <- med_df$V3[3]
  dir_effect_ci_hi <- med_df$V4[3]
  dir_effect_p <- med_df$V5[3]
  
  # mediation effect (are the two models significantly different?)
  med_effect_b <- med_df$V2[2]
  med_effect_ci_lo <- med_df$V3[2]
  med_effect_ci_hi <- med_df$V4[2]
  med_effect_p <- med_df$V5[2]
  
  vec <- data.frame(tot_effect_b, tot_effect_ci_lo, tot_effect_ci_hi, tot_effect_p,
                    dir_effect_b, dir_effect_ci_lo, dir_effect_ci_hi, dir_effect_p,
                    med_effect_b, med_effect_ci_lo, med_effect_ci_hi, med_effect_p)
  
  return(vec)
}

all_indi <- unique(xcell$indication)

all2 <- do.call("rbind", sapply(1:length(all_indi), function(y) {
  indi <- all_indi[[y]]
  mini <- dplyr::filter(xcell, indication == indi)
  
  all_snps <- mini$SNP
  all_genes <- mini$gene
  all_genesets <- mini$geneset
  
  # get expression, geno and xcell files
  # expression
  exp_file <- paste0("/gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/tissues/Results_1001_2016/", indi, "/", indi, ".ready.mee")
  exp0 <- read_delim(exp_file, delim = "\t") %>%
    rename(id = geneIds)
  exp <- left_join(exp0, igis, by = "id") %>% dplyr::select(Symbol, 2:ncol(exp0))
  
  # xcell score
  xcell_file <- paste0("/gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/genesets/xcell_genesets/", indi, "/", indi, ".ready.mee")
  xcell <- read_delim(xcell_file, delim = "\t")
  
  # genotype
  geno_file <- paste0("/gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/genesets/xcell_genesets/Mini_geno/", indi, "_mini_geno.txt")
  geno <- read_delim(geno_file, delim = "\t")
  
  # make sure the gene in xcell_he is also in my RNAseq file
  all_genes2 <- all_genes[all_genes %in% exp$Symbol]
  all_snps2 <- all_snps[all_genes %in% exp$Symbol]
  all_genesets2 <- all_genesets[all_genes %in% exp$Symbol]

  # loop every snp-gene-geneset pair
  all <- do.call("rbind", sapply(1:length(all_snps2), function(x) {
    gene <- all_genes2[[x]]
    geneset <- all_genesets2[[x]]
    snp <- all_snps2[[x]]
    print(c(indi, gene, geneset, snp))
    
    mini_exp <- dplyr::filter(exp, Symbol == gene)
    mini_xcell <- dplyr::filter(xcell, geneIds == geneset)
    mini_geno <- dplyr::filter(geno, id == snp)
    
    mini_geno <- as.data.frame(t(mini_geno[,-1])) %>% tibble::rownames_to_column()
    names(mini_geno) <- c("sample", "geno")
    mini_exp <- as.data.frame(t(mini_exp[,-1])) %>% tibble::rownames_to_column()
    names(mini_exp) <- c("sample", "gene_exp")
    mini_xcell <- as.data.frame(t(mini_xcell[,-1])) %>% tibble::rownames_to_column()
    names(mini_xcell) <- c("sample", "xcell_score")
    
    temp <- inner_join(mini_exp, mini_geno, by = "sample")
    df <- inner_join(temp, mini_xcell, by = "sample")
    df <- df[!is.na(df$geno),]
    
    # make sure geno is dose incremental 0 1 2 and not categorical
    df$geno <- as.numeric(as.character(df$geno))
    
    # mediation test
    model1 <- lm(gene_exp ~ geno, df)
    model2 <- lm(xcell_score ~ geno + gene_exp, df)
    results <- mediate(model1, model2, treat="geno", mediator="gene_exp",
                       boot=TRUE, sims=500)
    med <- summary(results)
    
    # get tidy result
    # trycatch is to avoid error - will return NAs
    tidy_med <- tryCatch(med_clean(med), error = function(e) {
      foo <- data.frame(t(rep(NA, 12)))
      names(foo) <- c("tot_effect_b", "tot_effect_ci_lo", "tot_effect_ci_hi", "tot_effect_p",
      "dir_effect_b", "dir_effect_ci_lo", "dir_effect_ci_hi", "dir_effect_p",
      "med_effect_b", "med_effect_ci_lo", "med_effect_ci_hi", "med_effect_p")
      return(foo)
    })
    
    p <- data.frame(snp, gene, geneset, indi, tidy_med)
                
    return(p)
    
  }, simplify = FALSE))
}, simplify = FALSE))

# this will generate many output files to be concatenated later
write.table(all2, file = paste0("/gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/genesets/xcell_genesets/Mediation/mediation_out_", args[2], ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")






