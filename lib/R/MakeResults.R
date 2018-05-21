## MakeResults.R v0.5
## different read count-based method for gene expression level estimation
## added machine-learning based algorithm as summarization score
## new method to retrieve genomic coordinates
## genome reference from db data
## edited gene fusions separator 

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]
ericscriptfolder <- split.vars[3]
readlength <- as.numeric(split.vars[4])
verbose <- as.numeric(split.vars[5])
refid <- as.character(split.vars [6])
dbfolder <- as.character(split.vars[7])
#genomeref <- as.character(split.vars[8])

flag.ada <- require(ada, quietly = T)
if (flag.ada == F) {
  require(kernlab, quietly = T)  
}
load(file.path(outputfolder,"out",paste(samplename,".chimeric.RData", sep = "")))
load(file.path(outputfolder,"out",paste(samplename,".DataMatrix.RData", sep = "")))
load(file.path(outputfolder,"out", paste(samplename,".ids_homology.RData", sep  = "")))
load(file.path(dbfolder, "data", refid, "EnsemblGene.GeneInfo.RData"))
load(file.path(dbfolder, "data", refid, "EnsemblGene.Structures.RData"))
load(file.path(dbfolder, "data", refid, "EnsemblGene.GeneNames.RData"))
load(file.path(dbfolder, "data", refid, "EnsemblGene.Sequences.RData"))
load(file.path(ericscriptfolder, "lib","data", "_resources", "BlackList.RData"))
load(file.path(ericscriptfolder, "lib","data", "_resources", "DataModel.RData"))
genomeref <- file.path(dbfolder, "data", refid, "allseq.fa")
  
checkselfhomology.fa <- scan(file = file.path(outputfolder,"out", paste(samplename, ".checkselfhomology.fa", sep = "")), sep= "\n", what = "", quiet = T)
seq.length <- nchar(checkselfhomology.fa)[2]
ix.id.fa <- which(checkselfhomology.fa %in% paste(">",info.id.and.homology[,1], sep = ""))
ix.fa <- ix.id.fa + 1
left_junction <- substr(checkselfhomology.fa[ix.fa], 1, (seq.length/2))
right_junction <- substr(checkselfhomology.fa[ix.fa], (seq.length/2 + 1), seq.length)
tmp_id <- unlist(strsplit(checkselfhomology.fa[ix.id.fa], "----", fixed = T))
left_id <- tmp_id[seq(1, length(ix.id.fa)*2, by = 2)]
right_id <- paste(">",tmp_id[seq(2, length(ix.id.fa)*2, by = 2)], sep = "")
exone.tmp <- unlist(strsplit(as.character(DataMatrix[,1]), "----", fixed = T))
exone.tmp1 <- exone.tmp[seq(1, dim(DataMatrix)[1]*2, by = 2)]
exone.tmp2 <- exone.tmp[seq(2, dim(DataMatrix)[1]*2, by = 2)]
ensgenename1 <- unlist(strsplit(exone.tmp1, "_"))[seq(1, dim(DataMatrix)[1]*2, by = 2)]
ensgenename2 <- unlist(strsplit(exone.tmp2, "_"))[seq(1, dim(DataMatrix)[1]*2, by = 2)]
ix.idd <- which(DataMatrix[,1] %in% info.id.and.homology[,1])
DataMatrixF <- DataMatrix[ix.idd ,]
homology <- info.id.and.homology[,2]
flag.dup.a <- info.id.and.homology[,3]
flag.dup.b <- info.id.and.homology[,4]
geneinfo.id <- as.character(EnsemblGene.Structures$EnsemblGene)
geneinfo.chr <- as.character(EnsemblGene.Structures$Chromosome)
geneinfo.start <- as.character(EnsemblGene.Structures$geneStart)
geneinfo.end <- as.character(EnsemblGene.Structures$geneEnd)
geneinfo.description <- as.character(EnsemblGene.GeneInfo$Description)
geneinfo.genename <- as.character(EnsemblGene.GeneInfo$GeneName)
geneinfo.strand <- as.character(EnsemblGene.Structures$Strand)
ensgenename12 <- as.character(DataMatrixF[,1])
exone.tmp <- unlist(strsplit(as.character(DataMatrixF[,1]), "----", fixed = T))
exone.tmp1 <- exone.tmp[seq(1, dim(DataMatrixF)[1]*2, by = 2)]
exone.tmp2 <- exone.tmp[seq(2, dim(DataMatrixF)[1]*2, by = 2)]
ensgenename1 <- unlist(strsplit(exone.tmp1, "_"))[seq(1, dim(DataMatrixF)[1]*2, by = 2)]
ensgenename2 <- unlist(strsplit(exone.tmp2, "_"))[seq(1, dim(DataMatrixF)[1]*2, by = 2)]
geneid1 <- rep("",length(ensgenename1))
geneid2 <- rep("",length(ensgenename1))
description1 <- rep("",length(ensgenename1))
description2 <- rep("",length(ensgenename1))
status1 <- rep("",length(ensgenename1))
status2 <- rep("",length(ensgenename1))
chr1 <- rep("", length(ensgenename1))
chr2 <- rep("", length(ensgenename1))
genestart1 <- rep("", length(ensgenename1))
genestart2 <- rep("", length(ensgenename1))
geneend1 <- rep("", length(ensgenename1))
geneend2 <- rep("", length(ensgenename1))
biotype1 <- rep("", length(ensgenename1))
biotype2 <- rep("", length(ensgenename1))
strand1 <- rep("", length(ensgenename1))
strand2 <- rep("", length(ensgenename1))

for (i in 1:length(ensgenename1)) {
  ix.gene1 <- which(geneinfo.id == ensgenename1[i])
  ix.gene1.info <- which(EnsemblGene.GeneInfo$EnsemblGene == ensgenename1[i])
  if (length(ix.gene1)!=0) {
    geneid1[i] <- geneinfo.genename[ix.gene1.info]
    description1[i] <- geneinfo.description[ix.gene1.info]
    chr1[i] <- geneinfo.chr[ix.gene1]
    strand1[i] <- geneinfo.strand[ix.gene1]
    genestart1[i] <- geneinfo.start[ix.gene1]
    geneend1[i] <- geneinfo.end[ix.gene1]
  }
  ix.gene2 <- which(geneinfo.id == ensgenename2[i])
  ix.gene2.info <- which(EnsemblGene.GeneInfo$EnsemblGene == ensgenename2[i])
  if (length(ix.gene2)!=0) {
    geneid2[i] <- geneinfo.genename[ix.gene2.info]
    description2[i] <- geneinfo.description[ix.gene2.info]
    chr2[i] <- geneinfo.chr[ix.gene2]
    strand2[i] <- geneinfo.strand[ix.gene2]
    genestart2[i] <- geneinfo.start[ix.gene2]
    geneend2[i] <- geneinfo.end[ix.gene2]
    
  }
}

exone1 <- unlist(strsplit(exone.tmp1, "_"))[seq(2, dim(DataMatrixF)[1]*2, by = 2)]
exone2 <- unlist(strsplit(exone.tmp2, "_"))[seq(2, dim(DataMatrixF)[1]*2, by = 2)]
genename1 <- geneid1
genename2 <- geneid2
diff.adj.tot <- rep(NA, length(ensgenename1))
n.crossing <- rep(NA, length(ensgenename1))
for (hh in 1: length(ensgenename1)) {
  hhix <- which(MyGF$id1 == ensgenename1[hh] & MyGF$id2 == ensgenename2[hh])
  diff.adj.tot[hh] <- MyGF$diffpos[hhix]
  n.crossing[hh] <- MyGF$nreads[hhix]
}
n.spanning <- as.numeric(as.character(DataMatrixF[,4])) 
edge.score <- round(as.numeric(as.character(DataMatrixF[,6])), digits = 4)
gjs.score <- round(as.numeric(as.character(DataMatrixF[,7])), digits = 4)
unique.score <- round(as.numeric(as.character(DataMatrixF[,8])), digits = 4)
isize.score <- round(as.numeric(as.character(DataMatrixF[,9])), digits = 4)
ins.size <- round(as.numeric(as.character(DataMatrixF[,3])), digits = 2)
adj <- rep("Not Adjacent", length(ensgenename1))
adj[which(diff.adj.tot == 1)] <- "Adjacent"
if (refid == "homo_sapiens") {
  
  acceptable.chrs <- c(seq.int(1,22), "X","Y")
  ix.chr1 <- which(chr1 %in% acceptable.chrs)
  ix.chr2 <- which(chr2 %in% acceptable.chrs)
  ix.chr.tmp1 <- intersect(ix.chr1, ix.chr2)
  ix.chr.tmp2 <- intersect(grep("NNN", right_junction, invert = T), grep("NNN", left_junction, invert = T))
  ix.chr <- intersect(ix.chr.tmp1, ix.chr.tmp2)
  
  
  if (length(ix.chr) > 0) {
    genename1 <- genename1[ix.chr]
    genename2 <- genename2[ix.chr]
    chr1 <- chr1[ix.chr]
    chr2 <- chr2[ix.chr]
    genestart1 <- genestart1[ix.chr]
    genestart2 <- genestart2[ix.chr]
    geneend1 <- geneend1[ix.chr]
    geneend2 <- geneend2[ix.chr]
    exone1 <- exone1[ix.chr]
    exone2 <- exone2[ix.chr]
    ensgenename1 <- ensgenename1[ix.chr]
    ensgenename12 <- ensgenename12[ix.chr]
    ensgenename2 <- ensgenename2[ix.chr]
    DataMatrixF <- rbind(DataMatrixF[ix.chr,])
    homology <- homology[ix.chr]
    flag.dup.a <- flag.dup.a[ix.chr]
    flag.dup.b <- flag.dup.b[ix.chr]
    adj <- adj[ix.chr]
    n.spanning <- n.spanning[ix.chr]
    n.crossing <- n.crossing[ix.chr]
    edge.score <- edge.score[ix.chr]
    unique.score <- unique.score[ix.chr]
    gjs.score <- gjs.score[ix.chr]
    isize.score <- isize.score[ix.chr]
    ins.size <- ins.size[ix.chr]
    biotype1 <- biotype1[ix.chr]
    biotype2 <- biotype2[ix.chr]
    status1 <- status1[ix.chr]
    status2 <- status2[ix.chr]
    description1 <- description1[ix.chr]
    description2 <- description2[ix.chr]
    left_id <- left_id[ix.chr]
    left_junction <- left_junction[ix.chr]
    right_id <- right_id[ix.chr]
    right_junction <- right_junction[ix.chr]
    strand1 <- strand1[ix.chr]
    strand2 <- strand2[ix.chr]
  }
}
stats <- read.delim(file.path(outputfolder, "out", paste(samplename, ".stats", sep = "")), sep = "\t", header = F)
gene.stats <- as.character(stats[,1])
length.stats <- as.numeric(stats[,2])
nreads.stats <- as.numeric(stats[,3])
gene1.exp <- rep(0, length(ensgenename1))
gene2.exp <- rep(0, length(ensgenename1))
gene12.exp <- round(((as.numeric(DataMatrixF[,2]) + as.numeric(DataMatrixF[,4]))*readlength/as.numeric(DataMatrixF[,5])), digits = 2)
#gene12.exp <- rep(0, length(ensgenename1))
for (iexpr in 1: length(ensgenename1)) {
  ix1.stats <- which(gene.stats == ensgenename1[iexpr])
  ix2.stats <- which(gene.stats == ensgenename2[iexpr])
  #    ix12.stats <- which(gene.stats == ensgenename12[iexpr])
  gene1.exp[iexpr] <- round(nreads.stats[ix1.stats]*readlength/length.stats[ix1.stats], digits = 2)
  gene2.exp[iexpr] <- round(nreads.stats[ix2.stats]*readlength/length.stats[ix2.stats], digits = 2)
  #    gene12.exp[iexpr] <- round(nreads.stats[ix12.stats]*readlength/length.stats[ix12.stats], digits = 2)
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
        genpos_1[i] <- pos.pos.ix[ix.okpos] - 1
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
        genpos_2[i] <- pos.pos.ix[ix.okpos] - 1 
      }
    }
  }
}
}

genpos_1.recal <- genpos_1
genpos_2.recal <- genpos_2

# # refine genomic coordinates
# for (i in 1: length(ensgenename1)) {
#   ix.ref <- grep(ensgenename1[i], GeneNames)
#   ix.ref.table <- which(GeneNames[ix.ref] == geneinfo.id)
#   if (strand1[i] == "+") {
#     exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonEnd[ix.ref.table]), ",")))
#   } else {
#     exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonStart[ix.ref.table]), ",")))
#   }
#   ix.exon <- which.min(abs(genpos_1[i] - exonpos))
#   mydiff <- abs(genpos_1[i] - exonpos[ix.exon])
#   if (mydiff <= 3) {
#     genpos_1.recal[i] <- exonpos[ix.exon]
#   }
#   
#   ix.ref <- grep(ensgenename2[i], GeneNames)
#   ix.ref.table <- which(GeneNames[ix.ref] == geneinfo.id)
#   if (strand2[i] == "+") {
#     exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonStart[ix.ref.table]), ",")))
#   } else {
#     exonpos <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonEnd[ix.ref.table]), ",")))
#   }
#   ix.exon <- which.min(abs(genpos_2[i] - exonpos))
#   mydiff <- abs(genpos_2[i] - exonpos[ix.exon])
#   if (mydiff <= 3) {
#     genpos_2.recal[i] <- exonpos[ix.exon]
#   }
# }


nreads.score <- pmin( n.crossing, n.spanning)/pmax(n.crossing, n.spanning)
myscores <- cbind(gjs.score, edge.score, nreads.score, gene12.exp)
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
ix.repeated <- unique(c(which(is.na(genpos_1)), which(is.na(genpos_2))))
ix.norepeated <- intersect(which(is.na(genpos_1) == F), which(is.na(genpos_2) == F))
genpos_1.recal[which(genpos_1.recal == 0)] <- "Unable to predict breakpoint position"
genpos_2.recal[which(genpos_2.recal == 0)] <- "Unable to predict breakpoint position"
myblacklist <- rep("", length(genename1))
ix.bl <- which((genename1 %in% gene.bl1 & genename2 %in% gene.bl2) | (genename1 %in% gene.bl2 & genename2 %in% gene.bl1))
if (length(ix.bl) > 0) {
  for (bli in 1: length(ix.bl)) {
    ix.bli <- which((gene.bl1 == genename1[ix.bl[bli]] & gene.bl2 == genename2[ix.bl[bli]]) | (gene.bl2 == genename1[ix.bl[bli]] & gene.bl1 == genename2[ix.bl[bli]]))
    myblacklist[ix.bl[bli]] <- paste("Frequency:", sum(freq.bl[ix.bli]))
  }
}
oddity.spanningreads <- rep(0, length(genpos_1))
oddity.spanningreads[which(n.spanning == 1 & n.crossing >= 10)] <- 1    
fusion.type <- rep("inter-chromosomal",length(genpos_1))
fusion.type[which(chr1==chr2)] <- "intra-chromosomal"
fusion.type[intersect(intersect(which(chr1==chr2) , which(adj == "Adjacent")), which(strand1 == strand2))] <- "Read-Through"
fusion.type[intersect(intersect(which(chr1==chr2) , which(adj == "Adjacent")), which(strand1 != strand2))] <- "Cis"
junctionsequence <- paste(tolower(left_junction), right_junction, sep = "")
SummaryMat <- cbind(genename1, genename2, chr1, genpos_1.recal, strand1, chr2, genpos_2.recal , strand2, ensgenename1, ensgenename2, n.crossing, n.spanning,ins.size, homology, fusion.type, myblacklist,description1, description2, junctionsequence, gene1.exp, gene2.exp, gene12.exp, edge.score, gjs.score, nreads.score, ericscore)
colnames(SummaryMat) <- c("GeneName1", "GeneName2", "chr1","Breakpoint1", "strand1", "chr2" ,"Breakpoint2", "strand2","EnsemblGene1", "EnsemblGene2", "crossingreads","spanningreads","mean insertsize","homology","fusiontype", "Blacklist","InfoGene1", "InfoGene2", "JunctionSequence", "GeneExpr1", "GeneExpr2", "GeneExpr_Fused", "ES", "GJS", "US","EricScore")
SummaryMat <- data.frame(SummaryMat)
save(SummaryMat, file=file.path(outputfolder, paste(samplename, ".Summary.RData", sep = "")))
if (dim(SummaryMat)[1] > 0) {
  write.table(SummaryMat, file = file.path(outputfolder,paste(samplename,".results.total.tsv", sep = "")), sep = "\t", row.names = F, quote = F)
  ix.sorting.score <- sort(ericscore, decreasing = T, index.return = T)$ix
  ericscore.sorted <- ericscore[ix.sorting.score]
  myblacklist.sorted <- myblacklist[ix.sorting.score]
  oddity.spanningreads.sorted <- oddity.spanningreads[ix.sorting.score]
  SummaryMat.sorted <- SummaryMat[ix.sorting.score, ]
  SummaryMat.Filtered <-  SummaryMat.sorted[which(ericscore.sorted > 0.5 & myblacklist.sorted == "" & oddity.spanningreads.sorted == 0), ]
  write.table(SummaryMat.Filtered[, -16], file = file.path(outputfolder,paste(samplename,".results.filtered.tsv", sep = "")), sep = "\t", row.names = F, quote = F)
} else
{
  Results <- "No Chimeric Transcript found!"
  write.table(Results, file = file.path(outputfolder,paste(samplename,".results.total.tsv", sep = "")), sep = "\t", row.names = F, col.names = F, quote = F)
  
}














