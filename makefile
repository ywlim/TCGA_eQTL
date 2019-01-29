# makefile for processing RNA-seq and genotype data for cis-eQTL discovery 
# start in empty projectDir
# run the commands using make tagname (eg. make process)

##########
## SETUP
##########

# set project name here
project = THCA

# set directory containing snploc, gene expression files etc.
phenoDir = /gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/pheno/

# set path to RNA-seq, genotype, copy number files
# gene x samples (first column: geneId is EntrezID, not gene symbol)
exp = $(phenoDir)/$(project).hg38.eset.tumor.nrpkm.tsv
cnv = $(phenoDir)/$(project).ready.cnv
geno = /gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/geno/TCGA_fsrm_maf01

# sample level annotation (contains age, gender etc)
annot = $(phenoDir)/$(project).sample_level_annotations.tsv

# set path to scripts
scriptDir = /gne/research/data/dnaseq/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/scripts

# set up Rscript options for running MatrixEQTL
runMatrixEQTL = /gne/research/data/HumGenet/analysis/HumGenet_dbGap/TCGA/eQTLs/scripts/run_MatrixEQTL_mod_partition_SS.R 

##########
## RUN
##########

# prepare genotype tsv file
# previously removed SNPs with missing genotyping rate > 10% and maf < 10%
# keep only inviduals in the current indication
# keep only normal blood derived genotype (-10 in the sample name)
prepGeno:
	cut -f 1 sampleDic.$(project).tsv > TCGA_$(project).data 
	grep -f TCGA_$(project).data $(geno).fam | grep '\-10\s' | awk '{print $$1,$$2}' > TCGA_$(project).dbgap.samples
	plink --noweb --bfile $(geno) --keep TCGA_$(project).dbgap.samples --out $(project)_fsrm_maf01 --recodeA
	mkdir Geno_temp
	cd Geno_temp; \
	$(scriptDir)/snp_raw_to_tsv_tcga.pl ../$(project)_fsrm_maf01.raw; \
	cat out1 | sed 's/_/\t/g' | cut -f 1 > temp; \
	mv temp out1; \
	paste out* > ../$(project)_fsrm_maf01.tsv
	rm -rf Geno_temp

# This R script does several things to get input files ready. Specifically:
# It grabs only the samples which are listed in sampleDic
# It harmonizes sample IDs between genotyping, RNAseq and copy number data
# It removes samples for which all types of data aren't available
# It computes covariates for PCA
process:
	$(scriptDir)/prepFiles.R -e $(exp) -g $(project)_fsrm_maf01.tsv -v $(cnv) -o $(project).preProcessed.mee -p $(project).ready.meg -q $(project).ready.cnv -h $(project).150k.gds -c $(project).evOnly.cov > process.log

# Process expression data (normalization)
normalize:
	$(scriptDir)/processExpGtexWay.R -i $(project).preProcessed.mee -r .1 -m 10 -o $(project).ready.mee > normalize.log

# harmonize cnv file with normalized ready.mee file
processCnv:
	$(scriptDir_tcga)/prepCnv.R -e $(project).ready.mee -v $(project).ready.cnv -q $(project).ready.cnv

# Add PEER factors to covariates file
addPeer:
	$(scriptDir)/runPeer.R -i $(project).ready.mee -o $(project).peer.cov
	$(scriptDir)/addPeer.R -i $(project).evOnly.cov -a $(project).peer.cov -o $(project).withPeer.cov 

# Add gender and age to covariates file
addCovar:
	$(scriptDir_tcga)/addCovars.R -i $(project).withPeer.cov -o $(project).cov -a $(annot) -s Patient_ID -l gender,age

# Submit the cis analysis job 
# change -u to decide what p-value cutoff to print (eg. without cut-off is -u 1)
runCis:
	bsub -e $(project).cis.%J -o $(project).cis.%J -R "rusage[mem=32]" "$(runMatrixEQTL) -g $(project).ready.meg -s $(phenoDir)/snploc_GRCh38.txt -e $(project).ready.mee -c $(project).cov -n $(project).ready.cnv -l $(phenoDir)/geneloc.b38.tsv -d 1000000 -u 1 -v 0 -o $(project).cis.hg38"
