# Modified from Matrix_eQTL_engine.R

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

# http://cran.r-project.org/web/packages/policies.html
OutputSaver_FDR <- setRefClass("OutputSaver_FDR",
	fields = list(
		sdata = ".listBuilder",
		gdata = ".listBuilder",
		cdata = ".listBuilder",
		bdata = ".listBuilder",
		ss_total = ".listBuilder",
		ss_cov_residual = ".listBuilder",
		ss_global_residual = ".listBuilder",
		fraction_var_genotype = ".listBuilder",
		fid = "list",
		testfun1 = "list",
		pvfun1 = "list"
	),
	methods = list(
		initialize = function () {
			sdata <<- MatrixEQTL:::.listBuilder$new();  # <<- is for assignment to variable in global namespace
			gdata <<- MatrixEQTL:::.listBuilder$new();
			cdata <<- MatrixEQTL:::.listBuilder$new();
			bdata <<- MatrixEQTL:::.listBuilder$new();
			ss_total <<- MatrixEQTL:::.listBuilder$new();  			 # added by Haiyin
			ss_cov_residual <<- MatrixEQTL:::.listBuilder$new();		 # added by Haiyin: sum of squares of residual after regressing out all covariates 
			ss_global_residual <<- MatrixEQTL:::.listBuilder$new();		# sum of squares of residual after regressing out global covariate
			fraction_var_genotype <<- MatrixEQTL:::.listBuilder$new(); # added by Haiyin: fraction of variance explained by genotype

			fid <<- list(0);
			testfun1 <<- list(0);
			pvfun1 <<- list(0);
			return(.self);
		},
		start = function(filename, statistic_name, unused1, unused2, testfun, pvfun) {
			testfun1 <<- list(testfun);
			pvfun1 <<- list(pvfun);
			if(length(filename) > 0) {
				if(class(filename) == "character") {
					fid <<- list(file(description = filename, open = "wt", blocking = FALSE, raw = FALSE), TRUE);
				} else {
					fid <<- list(filename, FALSE)
				}
				writeLines( paste("SNP", 
								  "gene",
								  statistic_name,  # beta and t-stat
								  "p-value", 
								  "FDR", 
								  "SS_total", 
								  "SS_global_residual", 
								  "SS_allcov_residual", 						  
								  "Fraction_genotype", sep = "\t"), fid[[1]]);
			} else {
				fid <<- list();
			}
		},
		update = function(spos, gpos, sta, beta = NULL, SS_total, SS_global_residual, SS_cov_residual, fraction_geno) {
			if(length(sta)>0) {
				sdata$add(spos);
				gdata$add(gpos);
				cdata$add(sta );
				if(!is.null(beta ))
					bdata$add(beta );
				#####  added by Haiyin #####
				ss_total$add(SS_total[gpos]);
				ss_cov_residual$add(SS_cov_residual[gpos])
				ss_global_residual$add(SS_global_residual[gpos])
				fraction_var_genotype$add(fraction_geno)
			}
			return(.self);
		},
		getResults = function( gene, snps, FDR_total_count) {
			pvalues = NULL;
 			if(cdata$n > 0) {
 				tests = testfun1[[1]](cdata$unlist());  # t-test is done based on mat[[1]]
 				cdata <<- MatrixEQTL:::.listBuilder$new(); #why clear out memory after retreiving results ???
 				
 				pvalues = pvfun1[[1]](tests);
 				ord = order(pvalues);
 				
 				tests = tests[ord];
 				pvalues = pvalues[ord];
 				
 				FDR = pvalues * FDR_total_count / (1:length(pvalues));
 				FDR[length(FDR)] = min(FDR[length(FDR)], 1);
 				FDR = rev(cummin(rev(FDR)));
 				
 				snps_names  = rownames(snps)[sdata$unlist()[ord]];
 				sdata <<- MatrixEQTL:::.listBuilder$new();
				gene_names  = rownames(gene)[gdata$unlist()[ord]];
 				gdata <<- MatrixEQTL:::.listBuilder$new();
 				
 				##### added by Haiyin #####
 				#browser()
 				SS_total = ss_total$unlist()[ord];
 				SS_residual = ss_cov_residual$unlist()[ord];
 				SS_global_residual = ss_global_residual$unlist()[ord];
 				Fraction_geno = fraction_var_genotype$unlist()[ord];

 				beta = NULL;
 				if(bdata$n > 0)
 					beta = bdata$unlist()[ord];
				
 				if(length(fid)>0)	{	
					step = 1000; ########### 100000
					for( part in 1:ceiling(length(FDR)/step) ) {
	 					fr = (part-1)*step + 1;
	 					to = min(part*step, length(FDR));
						dump = data.frame(snps_names[fr:to],
										  gene_names[fr:to],
										  if(is.null(beta)) tests[fr:to] else list(beta[fr:to],tests[fr:to]),
										  pvalues[fr:to],
										  FDR[fr:to], 
										  ##### added by Haiyin #####
										  SS_total[fr:to],
										  SS_global_residual[fr:to],
										  SS_residual[fr:to],									  
										  Fraction_geno[fr:to],
										  row.names = NULL, 
										  check.rows = FALSE, 
										  check.names = FALSE, 
										  stringsAsFactors = FALSE);
						write.table(dump, file = fid[[1]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
					}
 				}
			} else {
				cat("No significant associations were found.\n", file = if(length(fid)>0){fid[[1]]}else{""});
			}
			if(length(fid)>0)	{	
				if(fid[[2]]) {
					close(fid[[1]]);
				}
	 			fid <<- list();
			}
 			
 			if(!is.null(pvalues)) {
 				eqtls = list( snps = snps_names,
							 				gene = gene_names,
											statistic = tests,
											pvalue = pvalues,
											FDR = FDR);
 				if(!is.null(beta))
 					eqtls$beta = beta;
 			} else {
 				eqtls = list( snps = character(),
							 				gene = character(),
				 							beta = numeric(),
											statistic = numeric(),
											pvalue = numeric(),
											FDR = numeric());
 			}
			return(list(eqtls = data.frame(eqtls)));
		}
	)
)


Matrix_eQTL_main_mod = function(	
						snps, 
						gene, 
					 	cnv,
						cvrt = SlicedData$new(), 
						output_file_name = "", 
						pvOutputThreshold = 1e-5,
						useModel = modelLINEAR, 
						errorCovariance = numeric(), 
						verbose = TRUE, 
						output_file_name.cis = "", 
						pvOutputThreshold.cis = 0,
						snpspos = NULL, 
						genepos = NULL,
						cisDist = 1e6,
 						pvalue.hist = FALSE,
						min.pv.by.genesnp = FALSE,
						noFDRsaveMemory = FALSE) {
	################################# Basic variable checks #################################
 	{				
		# status("Performing basic checks of the input variables");
		stopifnot( "SlicedData" %in% class(gene) );
		stopifnot( any(c("SlicedData","SlicedData.fmt") %in% class(snps)) );
		stopifnot( "SlicedData" %in% class(cvrt) );
		
		# Check dimensions
		if( min(snps$nRows(),snps$nCols()) == 0 )
			stop("Empty genotype dataset");
		if( min(gene$nRows(),gene$nCols()) == 0 )
			stop("Empty expression dataset");
		if( snps$nCols() != gene$nCols() )
			stop("Different number of samples in the genotype and gene expression files");
		if( cvrt$nRows()>0 ) {
			if( snps$nCols() != cvrt$nCols() )
				stop("Wrong number of samples in the matrix of covariates");
		}

		stopifnot( class(pvOutputThreshold) == "numeric" );
		stopifnot( length(pvOutputThreshold) == 1 );
		stopifnot( pvOutputThreshold >= 0 );
		stopifnot( pvOutputThreshold <= 1 );

		stopifnot(  class(noFDRsaveMemory) == "logical" );
		stopifnot( length(noFDRsaveMemory) == 1 );

		if( pvOutputThreshold > 0 ) {
			stopifnot( !((length(output_file_name) == 0) && noFDRsaveMemory) )
			stopifnot( length(output_file_name) <= 1 );
			if( length(output_file_name) == 1 ) {
				stopifnot( class(output_file_name) %in% c("character","connection") );
			}
		}
		
		stopifnot( class(pvOutputThreshold.cis) == "numeric" );
		stopifnot( length(pvOutputThreshold.cis) == 1 );
		stopifnot( pvOutputThreshold.cis >= 0 );
		stopifnot( pvOutputThreshold.cis <= 1 );
		stopifnot( !((pvOutputThreshold > 0) & (pvOutputThreshold.cis > 0) & (pvOutputThreshold > pvOutputThreshold.cis)) );
		stopifnot( (pvOutputThreshold > 0) | (pvOutputThreshold.cis > 0) );

		stopifnot( class(useModel) == class(modelLINEAR) );
		stopifnot( length(useModel) == 1 );
		stopifnot( useModel %in% c(modelLINEAR, modelANOVA, modelLINEAR_CROSS) );
		if( useModel %in%  c(modelLINEAR, modelLINEAR_CROSS) ) {
			if( snps$nCols() <= cvrt$nRows() + 1 + 1) {
				stop("The number of covariates exceeds the number of samples.\nLinear regression can not be fit.")
			}
		}
		if( useModel == modelLINEAR_CROSS ) {
			if( cvrt$nRows() == 0 ) {
				stop( "Model \"modelLINEAR_CROSS\" requires at least one covariate" );
			}
		}
		if( useModel == modelANOVA ) {
			n.anova.groups = getOption("MatrixEQTL.ANOVA.categories", 3);
			stopifnot( n.anova.groups == floor(n.anova.groups) );
			stopifnot( n.anova.groups >= 3 );
# 			stopifnot( n.anova.groups < snps$nCols() - cvrt$nRows() - 2 );
			if( snps$nCols() <= cvrt$nRows() + n.anova.groups) {
				stop("The number of covariates exceeds the number of samples.\nLinear regression (ANOVA) can not be fit.")
			}
		}
		
		stopifnot(  class(verbose) == "logical" );
		stopifnot( length(verbose) == 1 );

		stopifnot(  class(min.pv.by.genesnp) == "logical" );
		stopifnot( length(min.pv.by.genesnp) == 1 );
	
		if( pvOutputThreshold.cis > 0 ) {
			stopifnot( !((length(output_file_name.cis) == 0) && noFDRsaveMemory) )
			stopifnot( length(output_file_name.cis) <= 1 );
			if( length(output_file_name.cis) == 1 ) {
				stopifnot( class(output_file_name.cis) %in% c("character","connection") );
			}

# 			stopifnot( class(output_file_name.cis) == "character" );
# 			stopifnot( length(output_file_name.cis) == 1 );
			stopifnot( class(snpspos) == "data.frame" );
			stopifnot( ncol(snpspos) == 3 );
			stopifnot( nrow(snpspos) > 0 );
			stopifnot( class(snpspos[1,3]) %in% c("integer", "numeric") )
			stopifnot( !any(is.na(snpspos[,3])) )
			stopifnot( class(genepos) == "data.frame" );
			stopifnot( ncol(genepos) == 4 );
			stopifnot( nrow(genepos) > 0 );
			stopifnot( class(genepos[1,3]) %in% c("integer", "numeric") )
			stopifnot( class(genepos[1,4]) %in% c("integer", "numeric") )
			stopifnot( !any(is.na(genepos[,3])) )
			stopifnot( !any(is.na(genepos[,4])) )
			stopifnot( nzchar(output_file_name.cis) )
		}
		
		if( pvOutputThreshold > 0 ) {
			stopifnot( nzchar(output_file_name) )
		}
		
		stopifnot( class(errorCovariance) %in% c("numeric", "matrix") );
		errorCovariance = as.matrix(errorCovariance);
		if(length(errorCovariance)>0) {
			if( nrow(errorCovariance) != ncol(errorCovariance) ) {
				stop("The covariance matrix is not square");
			}	
			if( nrow(errorCovariance) != snps$nCols() ) {
				stop("The covariance matrix size does not match the number of samples");
			}
			if( !all(errorCovariance == t(errorCovariance)) ) {
				stop("The covariance matrix is not symmetric");
			}
		}
	}
	################################# Initial setup #########################################
	{
		gene.std = MatrixEQTL:::.listBuilder$new();
		snps.std = MatrixEQTL:::.listBuilder$new();
		
		dont.clone.gene = getOption("MatrixEQTL.dont.preserve.gene.object", FALSE)
		if(is.null(dont.clone.gene))
			dont.clone.gene = FALSE;
		
		if( !dont.clone.gene )
			gene = gene$Clone();
		# snps = snps$Clone(); # snps is read only
		cvrt = cvrt$Clone();

		params = list(
			output_file_name = output_file_name, 
			pvOutputThreshold = pvOutputThreshold,
			useModel = useModel, 
			errorCovariance = errorCovariance , 
			verbose = verbose, 
			output_file_name.cis = output_file_name.cis, 
			pvOutputThreshold.cis = pvOutputThreshold.cis,
			cisDist = cisDist ,
	 		pvalue.hist = pvalue.hist,
			min.pv.by.genesnp = min.pv.by.genesnp);

		if( verbose ) {
			lastTime = 0;
			status <- function(text) {
				# gc();
				newTime = proc.time()[3];
				if(lastTime != 0) {
					cat("Task finished in ", newTime-lastTime, " seconds\n");
				}
				cat(text,"\n");
				lastTime <<- newTime;
				unused = flush.console();
			}
		} else {
			status = function(text){}
		}
		start.time = proc.time()[3];
	}
	################################# Error covariance matrix processing ####################
	{
		if( length(errorCovariance) > 0 ) {
			status("Processing the error covariance matrix");
			eig = eigen(errorCovariance, symmetric = TRUE)
			d = eig$values;
			v = eig$vectors;
			#  errorCovariance == v %*% diag(d) %*% t(v)
			#  errorCovariance^0.5 == v*sqrt(d)*v" (no single quotes anymore)
			#  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v"
			if( any(d<=0) ) {
				stop("The covariance matrix is not positive definite");
			}
			correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v);
			rm( eig, v, d, errorCovariance )
		} else {
			rm( errorCovariance );
			correctionMatrix = numeric();
		}
	}
	################################# Matching gene and SNPs locations ######################
	if( pvOutputThreshold.cis > 0 ) {
		status("Matching data files and location files")
		
		# names in the input data	
		gene_names = rownames(gene);
		snps_names = rownames(snps);
		
		# gene range, set: left<right
		if(any(genepos[,3] > genepos[,4])) {
			temp3 = genepos[,3];
			temp4 = genepos[,4];
			genepos[,3] = pmin(temp3,temp4);
			genepos[,4] = pmax(temp3,temp4);
			rm(temp3, temp4);
		}
		
		# match with the location data
		genematch = match( gene_names, genepos[ ,1],  nomatch = 0L);
		usedgene = matrix(FALSE, nrow(genepos), 1); # genes in "genepos" that are matching  "gene_names"
		usedgene[ genematch ] = TRUE;
		if( !any(genematch) ) {
			stop("Gene names do not match those in the gene location file.");
		}
		cat( sum(genematch>0), "of", length(gene_names), " genes matched\n");
		
		
		snpsmatch = match( snps_names, snpspos[ ,1],  nomatch = 0L);
		usedsnps = matrix(FALSE, nrow(snpspos),1);
		usedsnps[ snpsmatch ] = TRUE;
		if( !any(snpsmatch) ) {
			stop("SNP names do not match those in the SNP location file.");
		}
		cat( sum(snpsmatch>0), "of", length(snps_names), " SNPs matched\n");
		
		# list used chr names
		chrNames = unique(c( as.character(unique(snpspos[usedsnps,2])), 
												 as.character(unique(genepos[usedgene,2])) ))
		chrNames = chrNames[ sort.list( suppressWarnings(as.integer(chrNames)), 
																		method = "radix", na.last = TRUE ) ];
		# match chr names
		genechr = match(genepos[,2],chrNames);
		snpschr = match(snpspos[,2],chrNames);
		
		# max length of a chromosome
		chrMax = max( snpspos[usedsnps, 3], genepos[usedgene, 4], na.rm = TRUE) + cisDist;
		
		# Single number location for all rows in "genepos" and "snpspos"
 		genepos2 = as.matrix(genepos[ ,3:4, drop = FALSE] + (genechr-1)*chrMax);
 		snpspos2 = as.matrix(snpspos[ ,3  , drop = FALSE] + (snpschr-1)*chrMax);
		
		# the final location arrays;
		snps_pos = matrix(0,length(snps_names),1);
		snps_pos[snpsmatch>0, ] = snpspos2[snpsmatch, , drop = FALSE];
		snps_pos[rowSums(is.na(snps_pos))>0, ] = 0;
		snps_pos[snps_pos==0] = (length(chrNames)+1) * (chrMax+cisDist);
		rm(snps_names, snpsmatch, usedsnps, snpschr, snpspos2)
		
		gene_pos = matrix(0,length(gene_names),2);
		gene_pos[genematch>0, ] = genepos2[genematch, , drop = FALSE];
		gene_pos[rowSums(is.na(gene_pos))>0, ] = 0;
		gene_pos[gene_pos==0] = (length(chrNames)+2) * (chrMax+cisDist);
		rm(gene_names, genematch, usedgene, genechr, genepos2)
		rm(chrNames, chrMax);

		if( is.unsorted(snps_pos) ) {
			status("Reordering SNPs\n");
			ordr = sort.list(snps_pos);
			snps$RowReorder(ordr);
			snps_pos = snps_pos[ordr, , drop = FALSE];
			rm(ordr);
		}
		if( is.unsorted(rowSums(gene_pos)) ) {
			status("Reordering genes\n");
			ordr = sort.list(rowSums(gene_pos));
			gene$RowReorder(ordr);
			gene_pos = gene_pos[ordr, , drop = FALSE];
			rm(ordr);
		}
		
		# Slice it back.
		geneloc = vector("list", gene$nSlices())
		gene_offset = 0;
		for(gc in 1:gene$nSlices()) {
			nr = gene$GetNRowsInSlice(gc);
			geneloc[[gc]] = gene_pos[gene_offset + (1:nr), , drop = FALSE];
			gene_offset = gene_offset + nr;	
		}
		rm(gc, gene_offset, gene_pos);
		
		snpsloc = vector("list", snps$nSlices())
		snps_offset = 0;
		for(sc in 1:snps$nSlices()) {
			nr = snps$GetNRowsInSlice(sc);
			snpsloc[[sc]] = snps_pos[snps_offset + (1:nr), , drop = FALSE];
			snps_offset = snps_offset + nr;	
		}
		rm(nr, sc, snps_offset, snps_pos);
	}
	################################# Covariates processing #################################
	{	
		status("Processing covariates");
		if( useModel == modelLINEAR_CROSS ) {
			last.covariate = as.vector(tail( cvrt$getSlice(cvrt$nSlices()), n = 1));
		}		
		if( cvrt$nRows()>0 ) {
			cvrt$SetNanRowMean();
			cvrt$CombineInOneSlice();
			cvrt = rbind(matrix(1,1,snps$nCols()),cvrt$getSlice(1));
		} else {
			cvrt = matrix(1,1,snps$nCols());
		}
		# Correct for the error covariance structure
		if( length(correctionMatrix)>0 ) {
			cvrt = cvrt %*% correctionMatrix;
		}
		# Orthonormalize covariates
		# status("Orthonormalizing covariates");
		q = qr(t(cvrt));
		if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * snps$nCols() ) {
			stop("Colinear or zero covariates detected");
		}
		cvrt = t( qr.Q(q) );

		q.ethnicity <- t( qr.Q( qr(t(cvrt[1:4,])) ) ) 		# QR decomposition of intercept and 3 PCs
		q.ethnicity.PEER <- t( qr.Q( qr(t(cvrt[1:19,])) ) ) # QR decomp of intercept + 3PCs + 15 PEER
		q.ethnicity.PEER.gender <- t( qr.Q( qr(t(cvrt[1:20,])) ) ) # QR decomp of interc + 3PCs + 15 PEER + gender
		#saveRDS(q, file = 'Global_QR.RDS')
		rm( q );
	}
	{	
        status("Define functions to augment covariate matrix")
        # Gram Schmidt orthognalization: X*Y_T*Y
        # tcrossprod(x,y) is x %*% t(y)

        orthoCVRT <- function(newcvrt, Q.base){
                new.cvrt.vect <- newcvrt - tcrossprod(newcvrt, Q.base) %*% Q.base # orthogonalized residual dimension
                new.cvrt.vect <- new.cvrt.vect/sqrt(sum(new.cvrt.vect^2)) # normalize
        }

        # Given a covariate vector for one gene, get orthonormalized vector w.r.t. t(Q) %*% Q
        orthoQ <- function(cvrt.vect, QQ) {
            ortho.cvrt.vect <- cvrt.vect - cvrt.vect %*% QQ           # orthogonalize
            ortho.cvrt.vect <- ortho.cvrt.vect/sqrt(sum(ortho.cvrt.vect ^ 2)) # make normal length
        }


	}


	################################# Gene expression processing ############################
	{
		status("Processing gene expression data (imputation, residualization, etc.)");
		# Impute gene expression: whenever there's missing value, replace with row mean of gene expression
		gene$SetNanRowMean();
		# cnv$SetNanRowMean(); # should this be done?
		# Correct for the error covariance structure
		if( length(correctionMatrix)>0 ) {
			gene$RowMatrixMultiply(correctionMatrix);
		}
		# Orthogonolize expression w.r.t. covariates
		# status("Orthogonolizing expression w.r.t. covariates");
		gene_offsets = double(gene$nSlices()+1);

        Q = cvrt  # QR orthonormalized column vectors of the sample-wide covariates
        QQ = t(Q) %*% Q  # sample number x sample number

        QQ.ethnicity = t(q.ethnicity) %*% q.ethnicity
        QQ.ethnicity.PEER = t(q.ethnicity.PEER) %*% q.ethnicity.PEER
        QQ.ethnicity.PEER.gender = t(q.ethnicity.PEER.gender) %*% q.ethnicity.PEER.gender


        status("Regressing out gene-level CNV, generating projection matrix for residualization, etc.)");
        SS_total = c()
        SS_cov_residual = c()
        SS_global_residual = c()
        SS_ethnicity_residual = c()
        SS_ethnicity_PEER_residual = c()
        SS_ethnicity_PEER_gender_residual = c()

		for( sl in 1:gene$nSlices() ) {
			slice = gene$getSlice(sl);
			#saveRDS(slice, file = paste('expression.slice', sl,'RDS', sep = '.'))

			gene_offsets[sl+1] = gene_offsets[sl] + nrow(slice);
			rowsq1 = rowSums(slice^2); #Expression data has been inverse normal transformed so it's centered at 0
			SS_total = c(SS_total, rowsq1)

			cnv_slice = cnv$getSlice(sl)
            stopifnot(cnv$fileSliceSize == gene$fileSliceSize)
            
            # we can orthogonalize CNV for all genes with the global covariates in one go, then normalize as follows
            # CNV_residual = cnv_slice - cnv_slice %*% qq  #  still need to normalize each row
            # but to extend to more than one gene-level covariates, say CNV + methylation, need to do this gene by gene

            projection <- lapply(1:nrow(slice), function(geneidx) {

            	project_ethnicity 			  <- slice[geneidx,] %*% QQ.ethnicity
                project_ethnicity.PEER 		  <- slice[geneidx,] %*% QQ.ethnicity.PEER
                project_ethnicity.PEER.gender <- slice[geneidx,] %*% QQ.ethnicity.PEER.gender
                project_global 				  <- slice[geneidx,] %*% QQ

                # for genes with missing CNV data (NA for any samples)
                # don't control for CNV
                if (any(is.na(cnv_slice[geneidx, ]))) {
                    project_all <- project_global   
                } else {
                	# to determine contribution only from CNV, first project only on global covariates
                	# project_global <- slice[geneidx,] %*% QQ  
                    # get orthonormalized vector of CNV for the gene in question (geneidx in the slice)
                    cnv_gene <- orthoQ(cnv_slice[geneidx, ], QQ)  
                    projection_matrix <- QQ + t(cnv_gene) %*% cnv_gene
                    project_all <- slice[geneidx,] %*% projection_matrix
                }
                return(list(project_all, 
                			project_global, 
                			project_ethnicity, 
                			project_ethnicity.PEER, 
                			project_ethnicity.PEER.gender))
            })

            # fit over all global covariates + CNV
            slice_cvrt_fit <- do.call(rbind, lapply(projection,'[[',1))
            rowsq2 = rowSums((slice-slice_cvrt_fit)^2);

            slice_global_fit <- do.call(rbind, lapply(projection,'[[',2))
            rowsq3 = rowSums((slice-slice_global_fit)^2)

			slice_ethnicity_fit <- do.call(rbind, lapply(projection,'[[',3))
			rowsq4 = rowSums((slice-slice_ethnicity_fit)^2)

			slice_ethnicity.PEER_fit <- do.call(rbind, lapply(projection,'[[',4))
			rowsq5 = rowSums((slice-slice_ethnicity.PEER_fit)^2)  

			slice_ethnicity.PEER.gender_fit <- do.call(rbind, lapply(projection,'[[',5))
			rowsq6 = rowSums((slice-slice_ethnicity.PEER.gender_fit)^2)  

			# residual after taking out all covariates
			SS_cov_residual = c(SS_cov_residual, rowsq2)  
			# residual after regressing out only global covariates
			SS_global_residual = c(SS_global_residual, rowsq3)
			# residual after fitting only intercept, gender and age
			SS_ethnicity_residual = c(SS_ethnicity_residual, rowsq4)
			# residual after regressing out only intercept, gender age, ethnicity PCs
			SS_ethnicity_PEER_residual = c(SS_ethnicity_PEER_residual, rowsq5)
			# residual after regressing out only intercept, gender age, ethnicity PCs
			SS_ethnicity_PEER_gender_residual = c(SS_ethnicity_PEER_gender_residual, rowsq6)

            # original matrix operation
			#slice = slice - tcrossprod(slice,cvrt) %*% cvrt;
			# residual after taking out projection onto all covariates (global and gene level)
            slice = slice - slice_cvrt_fit

			# kill rows colinear with the covariates
			delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
			slice[delete.rows,] = 0;
			rowsq2[delete.rows] = 1;
			div = sqrt(rowsq2); #sqrt( rowSums(slice^2) );
# 			div[ div == 0 ] = 1;
			gene.std$set(sl, div);
			gene$setSlice(sl, slice / div); # sum of squares of of each gene in new slice = 1
		}


		# save sum of squares regression statistics
		ss_stats <- data.frame( geneID = rownames(as.matrix(gene)),
								# sum of squares
								SS_total = SS_total, 					 # total residual, number is ~ number of samples b/c INT process								
								SS_ethnicity_residual = SS_ethnicity_residual, 
								SS_ethnicity_PEER_residual = SS_ethnicity_PEER_residual,
								SS_ethnicity_PEER_gender_residual = SS_ethnicity_PEER_gender_residual,  
								SS_global_residual = SS_global_residual, # residual after regressing out all global/sample level covariates
								SS_cov_residual =  SS_cov_residual,		 # residual after regressing out all residuals: global + CNV
								# fractions of residuals
								fraction_var_intercept_ethnicity = 1 - SS_ethnicity_residual/SS_total,	# fraction of total variance due to intercept + gender + age
								fraction_var_PEER = (SS_ethnicity_residual-SS_ethnicity_PEER_residual)/SS_total,  # fraction of total variance due to PEER
								fraction_var_gender = (SS_ethnicity_PEER_residual-SS_ethnicity_PEER_gender_residual)/SS_total,
								fraction_var_age = (SS_ethnicity_PEER_gender_residual- SS_global_residual)/SS_total,
								fraction_var_global_covariates = 1- SS_global_residual/SS_total,
								fraction_var_cnv = (SS_global_residual - SS_cov_residual) /SS_total,
								fraction_var_all_covariates = 1 - SS_cov_residual/SS_total)
		write.table(ss_stats,
					file = paste0(output_file_name.cis, '.partial.sum.of.squares.txt'),
					quote = FALSE, sep = "\t", row.names = FALSE)

		# rr <- sweep(rbind(SS_global_residual, SS_cov_residual), 1, SS_total, '/')
		# max(rr[,1]-rr[,2])  # maximum fraction of variance due to CNV
		rm(rowsq1, rowsq2, delete.rows, div);
		rm( sl, slice );
		#gene$RowRemoveZeroEps();
	}
	################################# snps_process, testfun, pvfun, threshfun, afun  ########
	{
		# snps_process - preprocess SNPs slice
		#
		# afun --- abs for signed stats, identity for non-negative
		# threshfun --- internal stat threshold for given p-value
		# testfun --- t or F statistic from the internal one
		# pvfun --- p-value from the t or F statistic
		
		nSamples = snps$nCols();
		nGenes = gene$nRows();
		nSnps  = snps$nRows();
		nCov = nrow(cvrt) + 1;
		# nVarTested = length(snps_list); # set in case(useModel)
		# dfNull = nSamples - nCov;
		# d.f. of the full model
		betafun = NULL;
		
		if( useModel == modelLINEAR ) {
			snps_process = function(x) {
				return( list(MatrixEQTL:::.SetNanRowMean(x)) );
			};
			nVarTested = 1;
			dfFull = nSamples - nCov - nVarTested;
			statistic.fun = function(mat_list) {
				return( mat_list[[1]] );
			}
			afun = function(x) {return(abs(x))};
			threshfun = function(pv) {
				thr = qt(pv/2, dfFull, lower.tail = FALSE);
				thr = thr^2;
				thr = sqrt(  thr / (dfFull + thr) );
				thr[pv >= 1] = 0;
				thr[pv <= 0] = 1;
				return( thr );
			}
			testfun = function(x) { return( x * sqrt( dfFull / (1 - MatrixEQTL:::.my.pmin(x^2,1))));	}
			pvfun = function(x) { return( MatrixEQTL:::.pv.nz(pt(-abs(x),dfFull)*2)); }
			thresh.cis = threshfun(pvOutputThreshold.cis);
			thresh = threshfun(pvOutputThreshold);
			betafun = function(stat, ss, gg, select) {
				return(stat * gene.std$get(gg)[select[,1]] / snps.std$get(ss)[select[,2]]);
			}
		} else 
		if( useModel == modelANOVA ) {
			snps_process = function(x).SNP_process_split_for_ANOVA(x,n.anova.groups);
			nVarTested = n.anova.groups - 1;
			dfFull = nSamples - nCov - nVarTested;
# 			statistic.fun = function(mat_list) {
# 				return( mat_list[[1]]^2 + mat_list[[2]]^2 );
# 			}
			statistic.fun = function(mat_list) {
				x = mat_list[[1]]^2;
				for( j in 2:length(mat_list) )
					x = x + mat_list[[j]]^2;
				return( x );
			}
			afun = identity;
			threshfun = function(pv) {
				thr = qf(pv, nVarTested, dfFull, lower.tail = FALSE);
				thr = thr / (dfFull/nVarTested + thr);
				thr[pv >= 1] = 0;
				thr[pv <= 0] = 1;
				return( thr );
			}
			testfun = function(x) { return( x / (1 - MatrixEQTL:::.my.pmin(x,1)) * (dfFull/nVarTested) ); }
			pvfun = function(x) { return( MatrixEQTL:::.pv.nz(pf(x, nVarTested, dfFull, lower.tail = FALSE)) ); }
			thresh.cis = threshfun(pvOutputThreshold.cis);
			thresh = threshfun(pvOutputThreshold);
		} else 
		if( useModel == modelLINEAR_CROSS ) {
			last.covariate = as.vector( last.covariate );
			snps_process = .SNP_process_split_for_LINEAR_CROSS = function(x) {
				out = vector("list", 2);
				out[[1]] = MatrixEQTL:::.SetNanRowMean(x);
				out[[2]] = t( t(out[[1]]) * last.covariate );
				return( out );
			};
			nVarTested = 1;
			dfFull = nSamples - nCov - nVarTested - 1;
			statistic.fun = function(mat_list) {
				return( mat_list[[2]] / sqrt(1 - mat_list[[1]]^2) );
			}
			afun = function(x) {return(abs(x))};
			threshfun = function(pv) {
				thr = qt(pv/2, dfFull, lower.tail = FALSE);
				thr = thr^2;
				thr = sqrt(  thr / (dfFull + thr) );
				thr[pv >= 1] = 0;
				thr[pv <= 0] = 1;
				return( thr );
			}
			testfun = function(x) { return( x * sqrt( dfFull / (1 - MatrixEQTL:::.my.pmin(x^2,1))));	}
			pvfun = function(x) { return( MatrixEQTL:::.pv.nz(pt(-abs(x),dfFull)*2 )); }		
			thresh.cis = threshfun(pvOutputThreshold.cis);
			thresh = threshfun(pvOutputThreshold);				
			betafun = function(stat, ss, gg, select) {
				return(stat * gene.std$get(gg)[select[,1]] / snps.std$get(ss)[select[,2]]);
			}
		}
		params$dfFull = dfFull;
	}
	################################# Saver class(es) creation ##############################
	{
		status("Creating output file(s)");
		if(noFDRsaveMemory) {
			if( pvOutputThreshold > 0 ) {
				saver.tra = MatrixEQTL:::.OutputSaver_direct$new();
			}
			if( pvOutputThreshold.cis > 0 ) {
				saver.cis = MatrixEQTL:::.OutputSaver_direct$new();
			}
		} else {
			if( pvOutputThreshold > 0 ) {
				saver.tra = OutputSaver_FDR$new();  # revised function by Haiyin
			}
			if( pvOutputThreshold.cis > 0 ) {
				saver.cis = OutputSaver_FDR$new();  # revised function by Haiyin
			}
		}
		if( pvOutputThreshold > 0 )
			if( pvOutputThreshold * gene$nRows() * snps$nRows() > 1000000 )
				if(!noFDRsaveMemory)
					cat("Warning: pvOutputThreshold may be too large.\nExpected number of findings > ", 
							pvOutputThreshold * gene$nRows() * snps$nRows(),"\n");
		if( (useModel == modelLINEAR) || (useModel == modelLINEAR_CROSS) ) {
			statistic_name = "t-stat";
		} else if( useModel == modelANOVA ) {
			statistic_name = "F-test";
		}
		if(!is.null(betafun))  # if there's beta values, print both beta and t-statistics
			statistic_name = paste("beta\t",statistic_name, sep="");
		if( pvOutputThreshold > 0 )
			saver.tra$start(output_file_name,     statistic_name, snps, gene, testfun, pvfun);
		if( pvOutputThreshold.cis > 0 )
			saver.cis$start(output_file_name.cis, statistic_name, snps, gene, testfun, pvfun);
		rm( statistic_name );
	}
	################################# Some useful functions #################################
	{
		orthonormalize.snps = function(cursnps, ss) {
			for(p in 1:length(cursnps)) {
				if(length(correctionMatrix)>0) {
					cursnps[[p]] = cursnps[[p]] %*% correctionMatrix;
				}
				rowsq1 = rowSums(cursnps[[p]]^2);
				cursnps[[p]] = cursnps[[p]] - tcrossprod(cursnps[[p]],cvrt) %*% cvrt;
				for(w in MatrixEQTL:::.seq(1L,p-1L))
					cursnps[[p]] = cursnps[[p]] - rowSums(cursnps[[p]]*cursnps[[w]]) * cursnps[[w]];
				rowsq2 = rowSums(cursnps[[p]]^2);
				delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
				cursnps[[p]][delete.rows,] = 0;
				div = sqrt( rowsq2 );
				div[ delete.rows ] = 1;
# 				show(c(rowsq2,rowsq1, div));
				cursnps[[p]] = cursnps[[p]]/div;
			}
			snps.std$set(ss, div);
			return(cursnps);
		}
# 		if( pvOutputThreshold.cis > 0 ) {
# 			is.cis.pair = function(gg,ss) {
# 				return(!( ( snpsloc[[ss]][1, 1] - tail( geneloc[[gg]][ , 2], n = 1L) > cisDist) |
# 					    ( geneloc[[gg]][1, 1] - tail( snpsloc[[ss]]      , n = 1L) > cisDist) ) );
# 			}
# 		}
 		if( pvOutputThreshold.cis > 0 ) {
# 			sn.l = sapply(snpsloc, function(x)x[1] );
# 			sn.r = sapply(snpsloc, function(x)tail(x,1) );
# 			ge.l = sapply(geneloc, function(x)x[1,1] );
# 			ge.r = sapply(geneloc, function(x)x[nrow(x) , 2] );
			sn.l = sapply(snpsloc, "[", 1 );
			sn.r = sapply(snpsloc, tail, 1 );
			ge.l = sapply(geneloc, "[", 1, 1 );
			ge.r = sapply( lapply(geneloc, tail.matrix, 1 ), "[", 2);
			gg.1 = findInterval( sn.l , ge.r + cisDist +1) + 1; # which gene slice the start of a snp slice map to, length 12
# 			cat(gg.1,"\n")
			gg.2 = findInterval( sn.r , ge.l - cisDist );		# which gene slice the end of a snp slice map to
# 			cat(gg.2,"\n")
			rm(sn.l, sn.r, ge.l, ge.r);
 		}

	}	
	################################# Prepare counters and histogram bins ###################
	{
		pvbins = NULL; # bin edges for p-values
		statbins = 0;  # bin edges for the test statistic (|t| or F)
		do.hist = FALSE;
		if( length(pvalue.hist) == 1 ) {
			if(pvalue.hist == "qqplot") {
				pvbins = c(0, 10^rev(seq(0, log10(.Machine$double.xmin)-1, -0.05)));
			} else
			if( is.numeric(pvalue.hist) ) {
				pvbins = seq(from = 0, to = 1, length.out = pvalue.hist+1);
			} else
			if( pvalue.hist == TRUE ) {
				pvbins = seq(from = 0, to = 1, length.out = 100+1);
			}
		} else
		if( is.numeric(pvalue.hist) && (length(pvalue.hist) > 1) ) {
			pvbins = pvalue.hist;
		}
		if( is.null(pvbins) && (pvalue.hist != FALSE) ) {
			stop("Wrong value of pvalue.hist. Must be FALSE, TRUE, \"qqplot\", or numerical");
		}
		do.hist = !is.null(pvbins);
		if( do.hist ) {
			pvbins = sort(pvbins);
			statbins = threshfun(pvbins);
			if( pvOutputThreshold > 0) {
				hist.all = MatrixEQTL:::.histogrammer$new(pvbins, statbins);
			}
			if( pvOutputThreshold.cis > 0) {
				hist.cis = MatrixEQTL:::.histogrammer$new(pvbins, statbins);
			}
		}
		rm( pvbins, statbins);
		if(min.pv.by.genesnp) {
			if( pvOutputThreshold > 0) {
				minpv.tra = MatrixEQTL:::.minpvalue$new(snps,gene);
			}
			if( pvOutputThreshold.cis > 0) {
				minpv.cis = MatrixEQTL:::.minpvalue$new(snps,gene);
			}
		}
	}
	################################# Main loop #############################################
	{
		beta = NULL;
		n.tests.all = 0;
		n.tests.cis = 0;
		n.eqtls.tra = 0;
		n.eqtls.cis = 0;
		
		status("Performing eQTL analysis");
		# ss = 1; gg = 1;
		# ss = snps$nSlices(); gg = gene$nSlices();
		
		snps_offset = 0;
		for(ss in 1:snps$nSlices()) {
# 		for(ss in 1:min(2,snps$nSlices())) { #for debug
			cursnps = NULL;
			nrcs = snps$GetNRowsInSlice(ss);
			
			# loop only through the useful stuff
			for(gg in if(pvOutputThreshold>0){1:gene$nSlices()}else{MatrixEQTL:::.seq(gg.1[ss],gg.2[ss])} ) {
				gene_offset = gene_offsets[gg];
				curgene = gene$getSlice(gg);  # 2000 x 183
				nrcg = nrow(curgene);
				if(nrcg == 0) next;
				
				rp = "";
				
				statistic = NULL;
				select.cis.raw = NULL;
				## do cis analysis
# 				if( (pvOutputThreshold.cis > 0) && ( is.cis.pair(gg, ss) ) ) {
				if( (pvOutputThreshold.cis > 0) && (gg >= gg.1[ss]) && (gg <= gg.2[ss]) ) {
					
					if( is.null( statistic ) ) {
						if( is.null( cursnps ) ) {
							# orthnormalize snps w.r.t. global covariates (cvrt)
							cursnps = orthonormalize.snps( snps_process( snps$getSlice(ss) ), ss ); # cursnps[[1]] is 200,000 x 183
						}					
						mat = vector("list", length(cursnps));
						for(d in 1:length(cursnps)) {
							mat[[d]] = tcrossprod(curgene, cursnps[[d]]);
						}
						statistic = statistic.fun( mat ); 	# identity function for modelLINEAR: mat[[1]], 2000 x 200,000
						astatistic = afun(statistic);		# absolute value of mat[[1]], 200 x 200,000
# 						rm(mat);
					}
					#browser()
# 					sn.l = findInterval(geneloc[[gg]][ ,1] - cisDist-1  +1   , snpsloc[[ss]]);
# 					sn.r = findInterval(geneloc[[gg]][ ,2] + cisDist    -1   , snpsloc[[ss]]);
					
					# index of snp in snpsloc[[ss]] 1Mb left of gene start
					# snpsloc[[ss]][sn.l]-geneloc[[gg]][,1] = 1Mb
					# geneloc[[gg]] is 2000 x 2
					# snpsloc[[ss]] is 200,000 x 1
					# when slice/ss 1, last snp is at 169,348,665, but last gene 2000 has right boundary of 225,978,168, way beyond las snp in current slice
					# find interval simply returns 200,000, the last SNP, which is >> 1Mb away, so cis-eQTLs of this genes (current gg slice) it will be handled at the next slice
					# how is chromosome managed? location of last snp on snploc[[22]] is 5,305,699,251
					sn.l = findInterval(geneloc[[gg]][ ,1] - cisDist-1, snpsloc[[ss]]); # vector of len 2000, index to left most snp within 1Mb of gene start
					# index of snp in snpsloc[[ss]] 1Mb right of gene end, max is 200,000
					sn.r = findInterval(geneloc[[gg]][ ,2] + cisDist, snpsloc[[ss]]);   # vecotr of len 2000, index to right most snp within 1Mb of gene end
					
					# vectorized index of those sn.l to sn.r, based on dimention of statistic (2000 x200,000)

					xx = unlist(lapply(which(sn.r>sn.l),FUN=function(x){(sn.l[x]:(sn.r[x]-1))*nrow(statistic)+x}))  
					# max index is 400M = 2Kx200K, for a slice, there could be 0 or thousands of snps passing threshold, 4138
					select.cis.raw = xx[ astatistic[xx] >= thresh.cis ];  # vector index of snps that pass thresh.cis
					select.cis = arrayInd(select.cis.raw, dim(statistic)) # array index of those SNPs where astatistic surpasses thresh.cis
					
					n.tests.cis = n.tests.cis + length(xx);
					n.eqtls.cis = n.eqtls.cis + length(select.cis.raw);
					
					if( do.hist )	
						hist.cis$update(astatistic[xx]);
					
					if( min.pv.by.genesnp ) {
	# 					minpv.cis$updatecis(ss, gg, arrayInd(xx, dim(statistic)), astatistic[xx])
						temp = double(length(astatistic));
						dim(temp) = dim(astatistic);
						temp[xx] = astatistic[xx];
						minpv.cis$update(ss, gg, temp);
					}
					
					if(!is.null(betafun))
						beta = betafun(mat[[length(mat)]][select.cis.raw], ss, gg, select.cis);
							
					# fitting R = beta * G + epsilon
					# beta or mat[[1]] will be betwen -1 and 1 given 
					# ||R|| = 1 (curgene normalized by sqrt(SS_residual))
					# ||G|| = 1 (rowSums(cursnps[[1]]^2) are all 1, done by orthonormalize.snps)

					# added by Haiyin
					# this value is between 0 and 1					
					SS_g = (mat[[1]][select.cis]^2) # fraction of variance explained by each genotype	

					# results to save
					saver.cis$update(snps_offset + select.cis[ , 2], # vec of len 4138: column index in astatistic
									 gene_offset + select.cis[ , 1], # vec of len 4138: row index to astatistic
									 statistic[select.cis.raw], # beta value in mat[[1]], those significant
									 beta,
									 ##### added by Haiyin #####
									 SS_total, 			 # length 21114, to be indexed by gene_offset + select.cis[ , 1] in $update
									 SS_global_residual, # length 21114, to be indexed as above
									 SS_cov_residual, 	 # length 21114, to be indexed as above
									 SS_g  # already the same length as beta
									 );

	# 				statistic.select.cis  = statistic[ select.cis ];
	# 				test = testfun( statistic.select.cis );
	# 				pv = pvfun(test);
	# 				Saver.cis$WriteBlock( cbind(snps_offset + select.cis[ , 2], gene_offset + select.cis[ , 1], test, pv) );
	# 				counter.cis$Update(gg, ss, select.cis, pv, n.tests = length(xx), if(do.hist) afun(statistic[xx]) )
					# Repor number of cis eQTL found in current gg-ss slice combination
					rp = paste(rp, ", ", formatC(n.eqtls.cis, big.mark=",", format = "f", digits = 0), " cis-eQTLs", sep = "");
				}
				## do trans/all analysis
				if(pvOutputThreshold>0) {
					if( is.null( statistic ) ) {
						if( is.null( cursnps ) ) {
							cursnps = orthonormalize.snps( snps_process( snps$getSlice(ss) ), ss ); # cursnps[[1]] is 200,000x183
						}
						mat = vector("list", length(cursnps));
						for(d in 1:length(cursnps)) {
							mat[[d]] = tcrossprod(curgene, cursnps[[d]]); # 1x200,000
						}
						statistic = statistic.fun( mat );
						astatistic = afun(statistic);    
# 						rm(mat);
					}
	
					if( do.hist )	
						hist.all$update(astatistic);
	
					if(!is.null(select.cis.raw)) 
						astatistic[xx] = -1;
	# 					select.tra.raw = select.tra.raw[!(select.tra.raw %in% select.cis.raw)];
					select.tra.raw = which( astatistic >= thresh);
					select.tra = arrayInd(select.tra.raw, dim(statistic))
					
					n.eqtls.tra = n.eqtls.tra + length(select.tra.raw);
					n.tests.all = n.tests.all + length(statistic);

					if(!is.null(betafun))
						beta = betafun(mat[[length(mat)]][select.tra.raw], ss, gg, select.tra);
										
					SS_g = (mat[[1]][select.tra]^2) 
					saver.tra$update( snps_offset + select.tra[ , 2],
									  gene_offset + select.tra[ , 1],
									  statistic[select.tra.raw],
									  beta,
									  ##### added by Haiyin #####
									  SS_total, 			 # length 21114, to be indexed by gene_offset + select.cis[ , 1] in $update
									  SS_global_residual, # length 21114, to be indexed as above
									  SS_cov_residual, 	 # length 21114, to be indexed as above
									  SS_g  # already the same length as beta
									);
									 
					
					if( min.pv.by.genesnp ) 
						minpv.tra$update(ss, gg, astatistic)
					
	# 				statistic.select.tra = statistic[ select.tra ];
	# 				test = testfun( statistic.select.tra );
	# 				pv = pvfun( test );
	# 				Saver$WriteBlock( cbind( snps_offset + select.tra[ , 2], gene_offset + select.tra[ , 1], test, pv) );
	# 				counter$Update(gg, ss, select.tra, pv, n.tests = nrcs*nrcg, if(do.hist) afun(statistic) )
					rp = paste(rp, ", ", formatC(n.eqtls.tra, big.mark=",", format = "f", digits = 0), if(pvOutputThreshold.cis > 0)" trans-"else" ","eQTLs", sep = "")
				}
	
				#gene_offset = gene_offset + nrcg;
				if( !is.null(statistic) ) {
					per = 100*(gg/gene$nSlices() + ss-1) / snps$nSlices();
					cat( formatC(floor(per*100)/100, format = "f", width = 5, digits = 2), "% done" , rp, "\n", sep = "");
	 				flush.console();
				}
			} # gg in 1:gene$nSlices()
			snps_offset = snps_offset + nrcs;
		} # ss in 1:snps$nSlices()
	}
	################################# Results collection ####################################
	{
		rez = list(time.in.sec = proc.time()[3] - start.time);
		rez$param = params;
		
		if(pvOutputThreshold.cis > 0) {
			rez.cis = list(ntests = n.tests.cis, neqtls = n.eqtls.cis);
			rez.cis = c(rez.cis, saver.cis$getResults( gene, snps, n.tests.cis) );
			if(do.hist)
				rez.cis = c(rez.cis, hist.cis$getResults() );
			if(min.pv.by.genesnp)
				rez.cis = c(rez.cis, minpv.cis$getResults(snps, gene, pvfun = function(x){pvfun(testfun(x))}) );
		}
		
		if(pvOutputThreshold>0) {
			rez.all = list(ntests = n.tests.all, neqtls = n.eqtls.tra + n.eqtls.cis);
			if(pvOutputThreshold.cis > 0) {
				rez.tra = list(ntests = n.tests.all - n.tests.cis, neqtls = n.eqtls.tra);
				rez.tra = c(rez.tra, saver.tra$getResults( gene, snps, n.tests.all - n.tests.cis) );
			} else {
				rez.all = c(rez.all, saver.tra$getResults( gene, snps, n.tests.all              ) );
			}
			if(do.hist) {
				rez.all = c(rez.all, hist.all$getResults() );
				if(pvOutputThreshold.cis > 0) {
					rez.tra$hist.bins = rez.all$hist.bins;
					rez.tra$hist.counts = rez.all$hist.counts - rez.cis$hist.counts;
				}
			}
			if(min.pv.by.genesnp) {
				if(pvOutputThreshold.cis > 0) {
					rez.tra = c(rez.tra, minpv.tra$getResults(snps, gene, pvfun = function(x){pvfun(testfun(x))}) );
				} else {
					rez.all = c(rez.all, minpv.tra$getResults(snps, gene, pvfun = function(x){pvfun(testfun(x))}) );
				}
			}
		}
	
		if(exists("rez.all")>0)
			rez$all = rez.all;
		if(exists("rez.tra")>0)
			rez$trans = rez.tra;
		if(exists("rez.cis")>0)
			rez$cis = rez.cis;	
		
		class(rez) = c(class(rez),"MatrixEQTL");
		status("");
	}
# 	cat("s std ",snps.std$get(1),"\n");
# 	cat("g std ",gene.std$get(1),"\n");
	################################# Results collection ####################################
	return(rez);
}

.histme = function(m, name1, name2, ...) {
	cnts = m$hist.counts;
	bins = m$hist.bins;
	ntst = m$ntests;
	centers = 0.5 * (bins[-1L] + bins[-length(bins)]);
	density = 0.5 / (bins[-1L] - centers) * cnts / ntst;
	ntext = paste("Histogram for ", name1, formatC(ntst, big.mark=",", format = "f", digits = 0), name2, " p-values ",sep="");
	r = structure(list(breaks = bins, counts = cnts, density = density,
	      mids = centers, equidist = FALSE), class = "histogram");
	plot(r, main = ntext, ylab = "Density", xlab = "P-values", ...)
	abline( h = 1, col = "blue");
	return(invisible());
}

.qqme = function(m, lcol, cex, pch, ...) {
	cnts = m$hist.counts;
	bins = m$hist.bins;
	ntst = m$ntests;
	
	cusu = cumsum(cnts) / ntst;
	ypos = bins[-1][is.finite(cusu)];
	xpos = cusu[is.finite(cusu)];
	lines(-log10(xpos), -log10(ypos), col = lcol, ...);
# 	lines(xpos, ypos, col = lcol, ...);
	if(length(m$eqtls$pvalue)==0)
		return();
	ypvs = -log10(m$eqtls$pvalue);
	xpvs = -log10(1:length(ypvs) / ntst);
	if(length(ypvs) > 1000) {
		# need to filter a bit, make the plotting faster
		levels = as.integer( xpvs/xpvs[1] * 1e3);
		keep = c(TRUE, diff(levels)!=0);
		levels = as.integer( ypvs/ypvs[1] * 1e3);
		keep = keep | c(TRUE, diff(levels)!=0);
		ypvs = ypvs[keep];
		xpvs = xpvs[keep];
		rm(keep)
	}
	points(xpvs, ypvs, col = lcol, pch = pch, cex = cex, ...);
}

plot.MatrixEQTL = function(x, cex = 0.5, pch = 19, xlim = NULL, ylim = NULL, main = NULL, ...) {
# 	cat(class(main),'\n')
	if( x$param$pvalue.hist == FALSE ) {
		warning("Cannot plot p-value distribution: the information was not recorded.\nUse pvalue.hist!=FALSE.");
		return(invisible());
	}
	if( x$param$pvalue.hist == "qqplot" ) {
		xmin = 1/max(x$cis$ntests, x$all$ntests);
		ymax = NULL;
		if(!is.null(ylim)) {
			ymax = ylim[2];
		} else {
			ymax = -log10(min( 
					x$cis$eqtls$pvalue[1],   x$cis$hist.bins[  c(FALSE,x$cis$hist.counts>0)][1],
					x$all$eqtls$pvalue[1],   x$all$hist.bins[  c(FALSE,x$all$hist.counts>0)][1],
					x$trans$eqtls$pvalue[1], x$trans$hist.bins[c(FALSE,x$trans$hist.counts>0)][1],
					na.rm = TRUE))+0.1;
		}
		if(ymax == 0) {
			ymax = -log10(.Machine$double.xmin)
		}
		if(!is.null(ymax))
			ylim = c(0,ymax);
		
		if(is.null(xlim))
			xlim =  c(0, -log10(xmin/1.5));
		
		plot(numeric(),numeric(), xlab = "-Log10(p-value), theoretical",
			ylab = "-Log10(p-value), observed",
			xlim = c(0, -log10(xmin/1.5)),
			ylim = ylim,
			xaxs="i", yaxs="i", ...);
		lines(c(0,1e3), c(0,1e3), col = "gray");
		if((x$param$pvOutputThreshold > 0) && (x$param$pvOutputThreshold.cis > 0)) {
			.qqme( x$cis, "red", cex, pch, ...);
			.qqme( x$trans, "blue", cex, pch, ...);
			if(is.null(main)) {
				main = paste("QQ-plot for",
					formatC(x$cis$ntests, big.mark=",", format = "f", digits = 0),
					"local and",
					formatC(x$trans$ntests, big.mark=",", format = "f", digits = 0),
					"distant p-values");
			}
			lset = c(1,2,4);
		} else
		if(x$param$pvOutputThreshold.cis > 0) {
			.qqme(x$cis, "red", cex, pch, ...);
			if(is.null(main)) {
				main = paste("QQ-plot for",
					formatC(x$cis$ntests, big.mark=",", format = "f", digits = 0),
					"local p-values");
			}
			lset = c(1,4);
		} else {
			.qqme(x$all, "blue", cex, pch, ...);
			if(is.null(main)) {
				main = paste("QQ-plot for all",
					formatC(x$all$ntests, big.mark=",", format = "f", digits = 0),
					"p-values");
			}
			lset = c(3,4);
		}
		title(main);

		legend("topleft",
			c("Local p-values","Distant p-values","All p-values","diagonal")[lset],
			col =      c("red","blue","blue","gray")[lset],
			text.col = c("red","blue","blue","gray")[lset],
			pch = 20, lwd = 1, pt.cex = c(1,1,1,0)[lset])
	} else {
		if((x$param$pvOutputThreshold > 0) && (x$param$pvOutputThreshold.cis > 0)) {
			par(mfrow=c(2,1));
			.histme(x$cis, "", " local", ...);
			tran = list(hist.counts = x$all$hist.counts - x$cis$hist.counts,
					hist.bins = x$all$hist.bins,
					ntests =  x$all$ntests - x$cis$ntests);
			.histme(x$trans,""," distant", ...);
			par(mfrow=c(1,1));
		} else
		if(x$param$pvOutputThreshold.cis > 0) {
			.histme(x$cis, "", " local", ...);
		} else {
			.histme(x$all, "all ", ""  , ...);
		}
	}
	return(invisible());
}

