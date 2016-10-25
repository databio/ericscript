### ensgene converter

toens <- function(ericscriptfolder, genename) {
	
	load(file.path(ericscriptfolder, "lib", "data", "EnsemblGene.GeneInfo.RData"))
	ensgene <- as.character(EnsemblGene.GeneInfo$EnsemblGene)[which(as.character(EnsemblGene.GeneInfo$GeneName) == genename)]
	if (length(ensgene) == 0) {ensgene <- NA}
	return(ensgene)
	
}


convertToComplement<-function(x) {
	
	bases=c("A","C","G","T")
	xx<-unlist(strsplit(toupper(x), NULL))
	paste(unlist(lapply(xx, function(bbb) {
						if(bbb=="A") compString <- "T"
						if(bbb=="C") compString <- "G"
						if(bbb=="G") compString <- "C"
						if(bbb=="T") compString <- "A"
						if(!bbb %in% bases) compString <- "N"
						return(compString)
						})),collapse="")
	
}

### import results from algorithm''s output

Import_ericscript <- function(outputpath) {
	
	filename <- grep(".results.total.tsv", list.files(outputpath), value = T)
	
	if (length(filename) == 1) {
		xx <- read.delim(file.path(outputpath, filename), sep = "\t", header = T)
		gene5 <- as.character(xx$EnsemblGene1)
		gene3 <- as.character(xx$EnsemblGene2)
		nreads <- as.numeric(as.character(xx$spanningreads))
		score <- as.numeric(as.character(xx$EricScore))
		seq <- as.character(xx$JunctionSequence)
		
		algout <- list()
		algout$gene5 <- gene5
		algout$gene3 <- gene3
		algout$nreads <- nreads
		algout$score <- score
		algout$seq <- seq
		
		return(algout)
	} else if (length(filename) == 0) {
		
		return(0)    
		
	} else if (length(filename) > 1) {
		
		return(length(filename))
		
	}
    
	
}



Import_defuse <- function(outputpath) {
	
	filename <- grep(".classify.tsv", list.files(outputpath), value = T)
    
	if (length(filename) == 1) {
		
		xx <- read.delim(file.path(outputpath, filename), sep = "\t", header = T)
		gene5 <- as.character(xx$gene1)
		gene3 <- as.character(xx$gene2)
		nreads <- as.numeric(as.character(xx$splitr_count))
		score <- as.numeric(as.character(xx$probability))
		seq <- rep("", dim(xx)[1])
		for (seqd in 1: dim(xx)[1]) {
			tmp <-  unlist(strsplit(as.character(xx$splitr_sequence[seqd]), "|", fixed = T))
			seq[seqd] <- paste(substr(tmp[1], (nchar(tmp[1])-29), nchar(tmp[1])), substr(tmp[2], 1, 30), sep = "")
		}
		
		
		algout <- list()
		algout$gene5 <- gene5
		algout$gene3 <- gene3
		algout$nreads <- nreads
		algout$score <- score
		algout$seq <- seq
		
		return(algout)
	} else if (length(filename) == 0) {
		
		return(0)    
		
	} else if (length(filename) > 1) {
		
		return(length(filename))
		
	}
	
}


Import_chimerascan <- function(outputpath) {
	
	filename <- grep("chimeras.bedpe", list.files(outputpath), value = T)
	
	if (length(filename) == 1) {
		
		xx <- read.delim(file.path(outputpath, filename), sep = "\t", header = T)
		gene1tmp <- as.character(xx$genes5p)
		gene2tmp <- as.character(xx$genes3p)
		nreadstmp <- as.numeric(as.character(xx$total_frags))
		scoretmp <- as.numeric(as.character(xx$score))
		
		gene1 <- c()
		gene2 <- c()
		nreads <- c()
		score <- c()
		
		for (i in 1: length(gene1tmp)) {
			
			if((length(grep(",", gene1tmp[i])) > 0) & (length(grep(",", gene2tmp[i])) == 0)) {
				
				gene1 <- c(gene1, unlist(strsplit(gene1tmp[i], ",")))
				myrep <- length(unlist(strsplit(gene1tmp[i], ",")))
				gene2 <- c(gene2, rep(gene2tmp[i], myrep))
				nreads <- c(nreads, rep(nreadstmp[i], myrep))
				score <- c(score, rep(scoretmp[i], myrep))  
				
			} else if ((length(grep(",", gene1tmp[i])) == 0) & (length(grep(",", gene2tmp[i])) > 0)) {
				
				gene2 <- c(gene2, unlist(strsplit(gene2tmp[i], ",")))
				myrep <- length(unlist(strsplit(gene2tmp[i], ",")))
				gene1 <- c(gene1, rep(gene1tmp[i], myrep))
				nreads <- c(nreads, rep(nreadstmp[i], myrep))
				score <- c(score, rep(scoretmp[i], myrep))  
				
			} else if ((length(grep(",", gene1tmp[i])) > 0) & (length(grep(",", gene2tmp[i])) > 0)) {
				
				gene1tmp1 <- unlist(strsplit(gene1tmp[i], ","))
				gene2tmp1 <- unlist(strsplit(gene2tmp[i], ","))
				myrep1 <- length(unlist(strsplit(gene1tmp[i], ",")))
				myrep2 <- length(unlist(strsplit(gene2tmp[i], ",")))     
				
				for (j in 1: myrep1) {
					
					gene1 <- c(gene1, rep(gene1tmp1[j], myrep2))
					gene2 <- c(gene2, gene2tmp1)
					nreads <- c(nreads, rep(nreadstmp[i], myrep2))
					score <- c(score, rep(scoretmp[i], myrep2))
					
				}
				
			} else {
				
				gene1 <- c(gene1, gene1tmp[i])
				gene2 <- c(gene2, gene2tmp[i])
				nreads <- c(nreads, nreadstmp[i])
				score <- c(score, scoretmp[i])
				
			}
			
		}
		
		seq <- rep(NA, length(gene1))
		gene5 <- rep(NA, length(gene1))
		gene3 <- rep(NA, length(gene1))
		
		for (i in 1: length(gene1)) {
			
			gene5[i] <- toens(ericscriptfolder, gene1[i])
			gene3[i] <- toens(ericscriptfolder, gene2[i])
			
		}  
		
		algout <- list()
		algout$gene5 <- gene5
		algout$gene3 <- gene3
		algout$nreads <- nreads
		algout$score <- score
		algout$seq <- seq
		
		return(algout)
	} else if (length(filename) == 0) {
		
		return(0)    
		
	} else if (length(filename) > 1) {
		
		return(length(filename))
		
	}
	
}



Import_shortfuse <- function(outputpath) {
	
	filename <- grep("fusion_counts.bedpe", list.files(outputpath), value = T)
	
	if (length(filename) == 1) {
		
		xx <- read.delim(file.path(outputpath, filename), sep = "\t", header = F)
		gene1tmp <- as.character(xx[, 11])
		gene2tmp <- as.character(xx[, 12])
		nreadstmp <- as.numeric(as.character(xx[, 8]))
		scoretmp <- as.numeric(as.character(xx[, 8]))
		
		gene12 <- paste(gene1tmp, gene2tmp)
		ixNOdupgene <- which(duplicated(gene12) == F)
		
		gene5 <- rep(NA, length(ixNOdupgene))
		gene3 <- rep(NA, length(ixNOdupgene))
		
		for (i in 1: length(ixNOdupgene)) {
			
			gene5[i] <- toens(ericscriptfolder, gene1tmp[ixNOdupgene[i]])
			gene3[i] <- toens(ericscriptfolder, gene2tmp[ixNOdupgene[i]])
			
		}  
		
		nreads <- nreadstmp[ixNOdupgene]
		score <- scoretmp[ixNOdupgene]
		seq <- rep(NA, length(nreads))
		
		algout <- list()
		algout$gene5 <- gene5
		algout$gene3 <- gene3
		algout$nreads <- nreads
		algout$score <- score
		algout$seq <- seq
		
		return(algout)
	} else if (length(filename) == 0) {
		
		return(0)    
		
	} else if (length(filename) > 1) {
		
		return(length(filename))
		
	}
	
}




Import_fusionmap <- function(outputpath) {
	
	filename <- grep("FusionReport.txt", list.files(outputpath), value = T)
	
	if (length(filename) == 1) {
		
		xx <- read.delim(file.path(outputpath, filename), sep = "\t", header = T)
		
		transcriptversetmp <- as.character(xx$FusionGene)
		genelisttmp <- unlist(strsplit(transcriptversetmp, "->"))
		gene1tmp <- genelisttmp[seq(1, length(genelisttmp), by = 2)]
		gene2tmp <- genelisttmp[seq(2, length(genelisttmp), by = 2)]
		gene1tmp1 <- as.character(xx$KnownGene1)
		gene2tmp1 <- as.character(xx$KnownGene2)
		seqtmp <- as.character(xx$FusionJunctionSequence)
		ix.reverse <- which(gene1tmp %in% gene1tmp1 == F)
		for (ii in 1: length(ix.reverse)) {  
			seqtmp[ix.reverse[ii]] <- convertToComplement(reverse(seqtmp[ix.reverse[ii]]))
		}
		nreadstmp <- as.numeric(as.character(xx[, 2]))
		scoretmp <- as.numeric(as.character(xx[, 2]))
		
		gene1 <- c()
		gene2 <- c()
		nreads <- c()
		score <- c()
		seq <- c()
		
		for (i in 1: length(gene1tmp)) {
			
			if((length(grep(",", gene1tmp[i])) > 0) & (length(grep(",", gene2tmp[i])) == 0)) {
				
				gene1 <- c(gene1, unlist(strsplit(gene1tmp[i], ",")))
				myrep <- length(unlist(strsplit(gene1tmp[i], ",")))
				gene2 <- c(gene2, rep(gene2tmp[i], myrep))
				nreads <- c(nreads, rep(nreadstmp[i], myrep))
				score <- c(score, rep(scoretmp[i], myrep))  
				seq <- c(seq, rep(seqtmp[i], myrep))
				
			} else if ((length(grep(",", gene1tmp[i])) == 0) & (length(grep(",", gene2tmp[i])) > 0)) {
				
				gene2 <- c(gene2, unlist(strsplit(gene2tmp[i], ",")))
				myrep <- length(unlist(strsplit(gene2tmp[i], ",")))
				gene1 <- c(gene1, rep(gene1tmp[i], myrep))
				nreads <- c(nreads, rep(nreadstmp[i], myrep))
				score <- c(score, rep(scoretmp[i], myrep))  
				seq <- c(seq, rep(seqtmp[i], myrep))
				
				
			} else if ((length(grep(",", gene1tmp[i])) > 0) & (length(grep(",", gene2tmp[i])) > 0)) {
				
				gene1tmp1 <- unlist(strsplit(gene1tmp[i], ","))
				gene2tmp1 <- unlist(strsplit(gene2tmp[i], ","))
				myrep1 <- length(unlist(strsplit(gene1tmp[i], ",")))
				myrep2 <- length(unlist(strsplit(gene2tmp[i], ",")))     
				
				for (j in 1: myrep1) {
					
					gene1 <- c(gene1, rep(gene1tmp1[j], myrep2))
					gene2 <- c(gene2, gene2tmp1)
					nreads <- c(nreads, rep(nreadstmp[i], myrep2))
					score <- c(score, rep(scoretmp[i], myrep2))     
					seq <- c(seq, rep(seqtmp[i], myrep2))
					
				}
				
			} else {
				
				gene1 <- c(gene1, gene1tmp[i])
				gene2 <- c(gene2, gene2tmp[i])
				nreads <- c(nreads, nreadstmp[i])
				score <- c(score, scoretmp[i])
				seq <- c(seq, seqtmp[i])
				
			}
			
		}
		
		gene5 <- rep(NA, length(gene1))
		gene3 <- rep(NA, length(gene1))
		
		for (i in 1: length(gene1)) {
			
			gene5[i] <- toens(ericscriptfolder, gene1[i])
			gene3[i] <- toens(ericscriptfolder, gene2[i])
			
		}  
		
		algout <- list()
		algout$gene5 <- gene5
		algout$gene3 <- gene3
		algout$nreads <- nreads
		algout$score <- score
		algout$seq <- seq
		
		return(algout)
	} else if (length(filename) == 0) {
		
		return(0)    
		
	} else if (length(filename) > 1){
		
		return(length(filename))
		
	}
	
}




Import_unknown <- function(outputpath) {
	
	filename <- grep("ericsim", list.files(outputpath), value = T)
	
	if (length(filename) == 1) {
		
		xx <- read.delim(file.path(outputpath, filename), sep = "\t", header = T)
		
		gene5 <- as.character(xx$gene5)
		gene3 <- as.character(xx$gene3)
		nreads <- as.numeric(as.character(xx$nread))
		score <- as.numeric(as.character(xx$score))
		seq <- as.character(x$seq)
		
		algout <- list()
		algout$gene5 <- gene5
		algout$gene3 <- gene3
		algout$nreads <- nreads
		algout$score <- score
		algout$seq <- seq
		
		return(algout)
		
	} else if (length(filename) == 0) {
		
		return(0)    
		
	} else if (length(filename) > 1) {
		
		return(length(filename))
		
	}
	
	
	
}












