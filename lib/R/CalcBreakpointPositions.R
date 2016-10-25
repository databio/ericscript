## re-calculate breakpoints positions for samples analysed with ericscript < 0.4.0
## and re-estimation of ericscore and blacklist

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
ericscriptfolder <- as.character(split.vars[1])
outputfolder <- split.vars [2]
dbfolder <- split.vars [3]
refid <- as.character(split.vars[4])
genomeref <- as.character(split.vars[5])

flag.ada <- require(ada, quietly = T)
if (flag.ada == F) {
  require(kernlab, quietly = T)  
}
load(file.path(ericscriptfolder, "lib","data", "_resources", "BlackList.RData"))
load(file.path(ericscriptfolder, "lib","data", "_resources", "DataModel.RData"))
load(file.path(dbfolder, "data", refid, "EnsemblGene.Structures.RData"))

myls <- list.files(outputfolder, pattern = "Summary.RData")
myls.tsv.total <- list.files(outputfolder, pattern = ".results.total.tsv")
myls.tsv.filt <- list.files(outputfolder, pattern = ".results.filtered.tsv")

if (length(myls) == 1 & length(myls.tsv.total) == 1 & length(myls.tsv.filt)) {
  samplename <- gsub(".results.total.tsv", "", myls.tsv.total)  
  load(file.path(outputfolder, myls))
  cat(paste("[EricScript] Re-estimating EricScore for sample ", samplename, "... ", sep = ""))
  
  ensgenename1 <- as.character(SummaryMat$EnsemblGene1)
  ensgenename2 <- as.character(SummaryMat$EnsemblGene2)
  genename1 <- as.character(SummaryMat$GeneName1)
  genename2 <- as.character(SummaryMat$GeneName2)
  gjs.score <- as.numeric(as.character(SummaryMat$GJS))
  edge.score <- as.numeric(as.character(SummaryMat$ES))
  nreads.score <- as.numeric(as.character(SummaryMat$US))
  cov.score <- as.numeric(as.character(SummaryMat$GeneExpr_Fused))
  myscores <- cbind(gjs.score, edge.score, nreads.score, cov.score)
  colnames(myscores) <- c("probs.gjs", "probs.es", "probs.us", "cov")
  myscores <- data.frame(myscores)
  if (flag.ada) {
    myada <- ada(control~., data = DataScores,loss="exponential", nu = 0.1)
    ericscore <- as.numeric(predict(myada, myscores, type = "probs")[,2])
  } else {
    sig <- sigest(control~., data = DataScores, frac = 1, na.action = na.omit, scaled = TRUE)[2]
    model <- ksvm(control~., data = DataScores, type = "C-svc", kernel = "rbfdot", kpar = list(sigma = sig), C = 1, prob.model = TRUE)
    ericscore <- predict(model, myscores, type = "probabilities")[,2]  
  }
  
  myblacklist <- rep("", length(genename1))
  ix.bl <- which((genename1 %in% gene.bl1 & genename2 %in% gene.bl2) | (genename1 %in% gene.bl2 & genename2 %in% gene.bl1))
  if (length(ix.bl) > 0) {
    for (bli in 1: length(ix.bl)) {
      ix.bli <- which((gene.bl1 == genename1[ix.bl[bli]] & gene.bl2 == genename2[ix.bl[bli]]) | (gene.bl2 == genename1[ix.bl[bli]] & gene.bl1 == genename2[ix.bl[bli]]))
      myblacklist[ix.bl[bli]] <- paste("Frequency:", sum(freq.bl[ix.bli]))
    }
  }
  
  cat("done. \n")
  cat(paste("[EricScript] Re-calculating breakpoint positions for sample ", samplename, "... ", sep = ""))
  
  myseq <- as.character(SummaryMat$JunctionSequence) 
  left_junction <- substr(myseq, 1, 50)
  right_junction <- substr(myseq, 51, 100)
  
  chr1 <- rep("", length(ensgenename1))
  chr2 <- rep("", length(ensgenename1))
  genestart1 <- rep("", length(ensgenename1))
  genestart2 <- rep("", length(ensgenename1))
  geneend1 <- rep("", length(ensgenename1))
  geneend2 <- rep("", length(ensgenename1))
  strand1 <- rep("", length(ensgenename1))
  strand2 <- rep("", length(ensgenename1))
  
  generef <- as.character(EnsemblGene.Structures$EnsemblGene)
  chrref <- as.character(EnsemblGene.Structures$Chromosome)
  genestartref <- as.character(EnsemblGene.Structures$geneStart)
  geneendref <- as.character(EnsemblGene.Structures$geneEnd)
  strandref <- as.character(EnsemblGene.Structures$Strand)
  for (i in 1: length(ensgenename1)) {
    
    ix.ref <- which(generef == ensgenename1[i])
    chr1[i] <- chrref[ix.ref]
    genestart1[i] <- genestartref[ix.ref]
    geneend1[i] <- geneendref[ix.ref]
    strand1[i] <- strandref[ix.ref]
    
    ix.ref <- which(generef == ensgenename2[i])
    chr2[i] <- chrref[ix.ref]
    genestart2[i] <- genestartref[ix.ref]
    geneend2[i] <- geneendref[ix.ref]
    strand2[i] <- strandref[ix.ref]
  }
  
  
  # NEW Find GenomicPosition (50nt)
  for (i in 1: length(ensgenename1)) {
    if (i == 1) {
      cat(paste("@", i, "_", 1, "\n", left_junction[i], "\n+\n", gsub(", ", "", toString(rep("I", nchar(left_junction[i])))), "\n", "@", i, "_", 2, "\n", right_junction[i],"\n+\n", gsub(", ", "", toString(rep("I", nchar(right_junction[i])))), sep = ""), sep = "\n", file = file.path(outputfolder, "out", "findgenomicpos.fq"), append = F)    
    } else {
      cat(paste("@", i, "_", 1, "\n", left_junction[i], "\n+\n", gsub(", ", "", toString(rep("I", nchar(left_junction[i])))), "\n", "@", i, "_", 2, "\n", right_junction[i],"\n+\n", gsub(", ", "", toString(rep("I", nchar(right_junction[i])))), sep = ""), sep = "\n", file = file.path(outputfolder, "out", "findgenomicpos.fq"), append = T)  
    }
  }
  system(paste("bwa aln", "-R 50", genomeref, file.path(outputfolder, "out", "findgenomicpos.fq"), ">", file.path(outputfolder, "out", "findgenomicpos.fq.sai"), "2>>", file.path(outputfolder, "out", ".ericscript.log")))
  system(paste("bwa samse", "-n 50", genomeref, file.path(outputfolder, "out", "findgenomicpos.fq.sai"), file.path(outputfolder, "out", "findgenomicpos.fq"), ">", file.path(outputfolder, "out", "findgenomicpos.fq.tmp"), "2>>", file.path(outputfolder, "out", ".ericscript.log")))
  system(paste("cat", file.path(outputfolder, "out", "findgenomicpos.fq.tmp"), "|", file.path(ericscriptfolder, "lib", "perl", "xa2multi.pl"), "-", "|","grep -v -e \'^\\@\' -",">", file.path(outputfolder, "out", "findgenomicpos.fq.sam")))
  xx.pos <- read.delim(file.path(outputfolder, "out", "findgenomicpos.fq.sam"), sep = "\t", header= F)
  genpos_1 <- rep(0,length(ensgenename1))
  genpos_2 <- rep(0,length(ensgenename1))
  id.pos <- as.character(xx.pos[[1]])
  flag.pos <- as.character(xx.pos[[2]])
  chr.pos <- as.character(xx.pos[[3]])
  if (length(grep("chr", chr.pos)) > 0) {
    chr.pos <- gsub("chr", "", chr.pos)
  }
  pos.pos <- as.numeric(as.character(xx.pos[[4]]))
  mapq.pos <- as.numeric(as.character(xx.pos[[5]]))
  
  for (i in 1: length(ensgenename1)) {
    
    ## for 5' gene
    ix.mypos <- which(id.pos == paste(i, "_1", sep = ""))
    chr.pos.ix <- chr.pos[ix.mypos]
    flag.pos.ix <- flag.pos[ix.mypos]
    pos.pos.ix <- pos.pos[ix.mypos]
    mapq.pos.ix <- mapq.pos[ix.mypos] 
    ix.okpos <- which(chr.pos.ix == chr1[i] & pos.pos.ix >= as.numeric(genestart1[i]) & pos.pos.ix <= as.numeric(geneend1[i]))
    if (length(ix.okpos) > 1) {
      ix.okpos <- ix.okpos[which.max(mapq.pos.ix[ix.okpos])]
    }
    if (length(ix.okpos) > 0) {
      if (flag.pos.ix[ix.okpos] == 16) {
        genpos_1[i] <- pos.pos.ix[ix.okpos] 
      } else {
        genpos_1[i] <- pos.pos.ix[ix.okpos] + 49
      }
    }
    ## for 3' gene
    ix.mypos <- which(id.pos == paste(i, "_2", sep = ""))
    chr.pos.ix <- chr.pos[ix.mypos]
    flag.pos.ix <- flag.pos[ix.mypos]
    pos.pos.ix <- pos.pos[ix.mypos]
    mapq.pos.ix <- mapq.pos[ix.mypos] 
    ix.okpos <- which(chr.pos.ix == chr2[i] & pos.pos.ix >= as.numeric(genestart2[i]) & pos.pos.ix <= as.numeric(geneend2[i]))
    if (length(ix.okpos) > 1) {
      ix.okpos <- ix.okpos[which.max(mapq.pos.ix[ix.okpos])]
    }
    if (length(ix.okpos) > 0) {
      if (flag.pos.ix[ix.okpos] == 16) {
        genpos_2[i] <- pos.pos.ix[ix.okpos] + 49
      } else {
        genpos_2[i] <- pos.pos.ix[ix.okpos] 
      }
    }
  }
  
  # NEW Find GenomicPosition (25nt_1)
  ix.na.pos_1 <- which(genpos_1 == 0)
  ix.na.pos_2 <- which(genpos_2 == 0)
  if (length(ix.na.pos_1 ) > 0 | length(ix.na.pos_2 ) > 0) {
  left_junction.trim <- substr(left_junction, 26, 50)
  right_junction.trim <- substr(right_junction, 1, 25)
  for (i in 1: length(ensgenename1)) {
    if (i == 1) {
      cat(paste("@", i, "_", 1, "\n", left_junction.trim[i], "\n+\n", gsub(", ", "", toString(rep("I", nchar(left_junction.trim[i])))), "\n", "@", i, "_", 2, "\n", right_junction.trim[i],"\n+\n", gsub(", ", "", toString(rep("I", nchar(right_junction.trim[i])))), sep = ""), sep = "\n", file = file.path(outputfolder, "out", "findgenomicpos.fq"), append = F)    
    } else {
      cat(paste("@", i, "_", 1, "\n", left_junction.trim[i], "\n+\n", gsub(", ", "", toString(rep("I", nchar(left_junction.trim[i])))), "\n", "@", i, "_", 2, "\n", right_junction.trim[i],"\n+\n", gsub(", ", "", toString(rep("I", nchar(right_junction.trim[i])))), sep = ""), sep = "\n", file = file.path(outputfolder, "out", "findgenomicpos.fq"), append = T)  
    }
  }
  system(paste("bwa aln", "-R 50", genomeref, file.path(outputfolder, "out", "findgenomicpos.fq"), ">", file.path(outputfolder, "out", "findgenomicpos.fq.sai"), "2>>", file.path(outputfolder, "out", ".ericscript.log")))
  system(paste("bwa samse", "-n 50", genomeref, file.path(outputfolder, "out", "findgenomicpos.fq.sai"), file.path(outputfolder, "out", "findgenomicpos.fq"), ">", file.path(outputfolder, "out", "findgenomicpos.fq.tmp"), "2>>", file.path(outputfolder, "out", ".ericscript.log")))
  system(paste("cat", file.path(outputfolder, "out", "findgenomicpos.fq.tmp"), "|", file.path(ericscriptfolder, "lib", "perl", "xa2multi.pl"), "-", "|","grep -v -e \'^\\@\' -",">", file.path(outputfolder, "out", "findgenomicpos.fq.sam")))
  xx.pos <- read.delim(file.path(outputfolder, "out", "findgenomicpos.fq.sam"), sep = "\t", header= F)
  id.pos <- as.character(xx.pos[[1]])
  flag.pos <- as.character(xx.pos[[2]])
  chr.pos <- as.character(xx.pos[[3]])
  if (length(grep("chr", chr.pos)) > 0) {
    chr.pos <- gsub("chr", "", chr.pos)
  }
  pos.pos <- as.numeric(as.character(xx.pos[[4]]))
  mapq.pos <- as.numeric(as.character(xx.pos[[5]]))
  for (i in 1: length(ensgenename1)) {
    if (i %in% ix.na.pos_1) {
      ## for 5' gene
      ix.mypos <- which(id.pos == paste(i, "_1", sep = ""))
      chr.pos.ix <- chr.pos[ix.mypos]
      flag.pos.ix <- flag.pos[ix.mypos]
      pos.pos.ix <- pos.pos[ix.mypos]
      mapq.pos.ix <- mapq.pos[ix.mypos] 
      ix.okpos <- which(chr.pos.ix == chr1[i] & pos.pos.ix >= as.numeric(genestart1[i]) & pos.pos.ix <= as.numeric(geneend1[i]))
      if (length(ix.okpos) > 1) {
        ix.okpos <- ix.okpos[which.max(mapq.pos.ix[ix.okpos])]
      }
      if (length(ix.okpos) > 0) {
        if (flag.pos.ix[ix.okpos] == 16) {
          genpos_1[i] <- pos.pos.ix[ix.okpos] 
        } else {
          genpos_1[i] <- pos.pos.ix[ix.okpos] + 24
        }
      }
    }
    ## for 3' gene
    if (i %in% ix.na.pos_2) {
      ix.mypos <- which(id.pos == paste(i, "_2", sep = ""))
      chr.pos.ix <- chr.pos[ix.mypos]
      flag.pos.ix <- flag.pos[ix.mypos]
      pos.pos.ix <- pos.pos[ix.mypos]
      mapq.pos.ix <- mapq.pos[ix.mypos] 
      ix.okpos <- which(chr.pos.ix == chr2[i] & pos.pos.ix >= as.numeric(genestart2[i]) & pos.pos.ix <= as.numeric(geneend2[i]))
      if (length(ix.okpos) > 1) {
        ix.okpos <- ix.okpos[which.max(mapq.pos.ix[ix.okpos])]
      }
      if (length(ix.okpos) > 0) {
        if (flag.pos.ix[ix.okpos] == 16) {
          genpos_2[i] <- pos.pos.ix[ix.okpos] + 24
        } else {
          genpos_2[i] <- pos.pos.ix[ix.okpos] 
        }
      }
    }
  }
  # NEW Find GenomicPosition (25nt_2)
  ix.na.pos_1 <- which(genpos_1 == 0)
  ix.na.pos_2 <- which(genpos_2 == 0)
  left_junction.trim <- substr(left_junction, 1, 25)
  right_junction.trim <- substr(right_junction, 26, 50)
  for (i in 1: length(ensgenename1)) {
    if (i == 1) {
      cat(paste("@", i, "_", 1, "\n", left_junction.trim[i], "\n+\n", gsub(", ", "", toString(rep("I", nchar(left_junction.trim[i])))), "\n", "@", i, "_", 2, "\n", right_junction.trim[i],"\n+\n", gsub(", ", "", toString(rep("I", nchar(right_junction.trim[i])))), sep = ""), sep = "\n", file = file.path(outputfolder, "out", "findgenomicpos.fq"), append = F)    
    } else {
      cat(paste("@", i, "_", 1, "\n", left_junction.trim[i], "\n+\n", gsub(", ", "", toString(rep("I", nchar(left_junction.trim[i])))), "\n", "@", i, "_", 2, "\n", right_junction.trim[i],"\n+\n", gsub(", ", "", toString(rep("I", nchar(right_junction.trim[i])))), sep = ""), sep = "\n", file = file.path(outputfolder, "out", "findgenomicpos.fq"), append = T)  
    }
  }
  system(paste("bwa aln", "-R 50", genomeref, file.path(outputfolder, "out", "findgenomicpos.fq"), ">", file.path(outputfolder, "out", "findgenomicpos.fq.sai"), "2>>", file.path(outputfolder, "out", ".ericscript.log")))
  system(paste("bwa samse", "-n 50", genomeref, file.path(outputfolder, "out", "findgenomicpos.fq.sai"), file.path(outputfolder, "out", "findgenomicpos.fq"), ">", file.path(outputfolder, "out", "findgenomicpos.fq.tmp"), "2>>", file.path(outputfolder, "out", ".ericscript.log")))
  system(paste("cat", file.path(outputfolder, "out", "findgenomicpos.fq.tmp"), "|", file.path(ericscriptfolder, "lib", "perl", "xa2multi.pl"), "-", "|","grep -v -e \'^\\@\' -",">", file.path(outputfolder, "out", "findgenomicpos.fq.sam")))
  
  
  xx.pos <- read.delim(file.path(outputfolder, "out", "findgenomicpos.fq.sam"), sep = "\t", header= F)
  id.pos <- as.character(xx.pos[[1]])
  flag.pos <- as.character(xx.pos[[2]])
  chr.pos <- as.character(xx.pos[[3]])
  if (length(grep("chr", chr.pos)) > 0) {
    chr.pos <- gsub("chr", "", chr.pos)
  }
  pos.pos <- as.numeric(as.character(xx.pos[[4]]))
  mapq.pos <- as.numeric(as.character(xx.pos[[5]]))
  for (i in 1: length(ensgenename1)) {
    if (i %in% ix.na.pos_1) {
      ## for 5' gene
      ix.mypos <- which(id.pos == paste(i, "_1", sep = ""))
      chr.pos.ix <- chr.pos[ix.mypos]
      flag.pos.ix <- flag.pos[ix.mypos]
      pos.pos.ix <- pos.pos[ix.mypos]
      mapq.pos.ix <- mapq.pos[ix.mypos] 
      ix.okpos <- which(chr.pos.ix == chr1[i] & pos.pos.ix >= as.numeric(genestart1[i]) & pos.pos.ix <= as.numeric(geneend1[i]))
      if (length(ix.okpos) > 1) {
        ix.okpos <- ix.okpos[which.max(mapq.pos.ix[ix.okpos])]
      }
      if (length(ix.okpos) > 0) {
        if (flag.pos.ix[ix.okpos] == 16) {
          genpos_1[i] <- pos.pos.ix[ix.okpos] 
        } else {
          genpos_1[i] <- pos.pos.ix[ix.okpos] + 24
        }
      }
    }
    ## for 3' gene
    if (i %in% ix.na.pos_2) {
      ix.mypos <- which(id.pos == paste(i, "_2", sep = ""))
      chr.pos.ix <- chr.pos[ix.mypos]
      flag.pos.ix <- flag.pos[ix.mypos]
      pos.pos.ix <- pos.pos[ix.mypos]
      mapq.pos.ix <- mapq.pos[ix.mypos] 
      ix.okpos <- which(chr.pos.ix == chr2[i] & pos.pos.ix >= as.numeric(genestart2[i]) & pos.pos.ix <= as.numeric(geneend2[i]))
      if (length(ix.okpos) > 1) {
        ix.okpos <- ix.okpos[which.max(mapq.pos.ix[ix.okpos])]
      }
      if (length(ix.okpos) > 0) {
        if (flag.pos.ix[ix.okpos] == 16) {
          genpos_2[i] <- pos.pos.ix[ix.okpos] + 24
        } else {
          genpos_2[i] <- pos.pos.ix[ix.okpos] 
        }
      }
    }
  }
  }
  
#   # refine genomic coordinates
#   genpos_1.recal <- genpos_1
#   genpos_2.recal <- genpos_2
#   for (i in 1: length(ensgenename1)) {
#     ix.ref <- which(generef == ensgenename1[i])
#     if (strand1[i] == "+") {
#       exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonEnd[ix.ref]), ",")))
#     } else {
#       exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonStart[ix.ref]), ",")))
#     }
#     ix.exon <- which.min(abs(genpos_1[i] - exonpos))
#     mydiff <- abs(genpos_1[i] - exonpos[ix.exon])
#     if (mydiff <= 3) {
#       genpos_1.recal[i] <- exonpos[ix.exon]
#     }
#     
#     ix.ref <- which(generef == ensgenename1[i])
#     if (strand2[i] == "+") {
#       exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonStart[ix.ref]), ",")))
#     } else {
#       exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonEnd[ix.ref]), ",")))
#     }
#     ix.exon <- which.min(abs(genpos_2[i] - exonpos))
#     mydiff <- abs(genpos_2[i] - exonpos[ix.exon])
#     if (mydiff <= 3) {
#       genpos_2.recal[i] <- exonpos[ix.exon]
#     }
#   }
#   mynames <- names(SummaryMat)
#   SummaryMat$Breakpoint1 <- genpos_1.recal
#   SummaryMat$Breakpoint2 <- genpos_2.recal
  
  SummaryMat$Breakpoint1 <- genpos_1
  SummaryMat$Breakpoint2 <- genpos_2  
  SummaryMat$EricScore <- ericscore
  SummaryMat$Blacklist <- myblacklist
  save(SummaryMat, file=file.path(outputfolder, paste(samplename, ".Summary.recalc.RData", sep = "")))
  n.spanning <- as.numeric(as.character(SummaryMat$spanningreads))
  n.crossing <- as.numeric(as.character(SummaryMat$crossingreads))
  oddity.spanningreads <- rep(0, length(genpos_1))
  oddity.spanningreads[which(n.spanning == 1 & n.crossing >= 10)] <- 1    
  
  if (dim(SummaryMat)[1] > 0) {
    write.table(SummaryMat, file = file.path(outputfolder,paste(samplename,".results.recalc.total.tsv", sep = "")), sep = "\t", row.names = F, quote = F)
    ix.sorting.score <- sort(ericscore, decreasing = T, index.return = T)$ix
    ericscore.sorted <- ericscore[ix.sorting.score]
    myblacklist.sorted <- myblacklist[ix.sorting.score]
    oddity.spanningreads.sorted <- oddity.spanningreads[ix.sorting.score]
    SummaryMat.sorted <- SummaryMat[ix.sorting.score, ]
    SummaryMat.Filtered <-  SummaryMat.sorted[which(ericscore.sorted > 0.5 & myblacklist.sorted == "" & oddity.spanningreads.sorted == 0), ]
    write.table(SummaryMat.Filtered[, -16], file = file.path(outputfolder,paste(samplename,".results.recalc.filtered.tsv", sep = "")), sep = "\t", row.names = F, quote = F)
  }
  cat("done. \n")
  cat(paste("[EricScript] Breakpoint position corrected files are in ", outputfolder, ".\n", sep = ""))
  
} else {  
  cat("[EricScript] No files of results found in ", outputfolder, ". Nothing to be done, Exit!\n", sep =" ")
}

