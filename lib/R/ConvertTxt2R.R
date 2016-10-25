vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
ericscriptfolder <- split.vars [1]
refid <- split.vars[2]
dbfolder <- split.vars[3]
tmpfolder <- split.vars[4]

xx <- read.delim(file.path(tmpfolder, "genepos.txt"), sep = "\t", header = F)
chrs <- as.character(xx[[2]])
unique.chrs <- unique(chrs)
geneid <- as.character(xx[[1]])
genepos <- as.numeric(as.character(xx[[3]]))
xx.strand <- read.delim(file.path(tmpfolder, "strand.txt"), sep = "\t", header = F)
strand <- as.character(xx.strand[[2]])
## sorting genes by genomic pos
ix.srt <- rep(NA, dim(xx)[1])
count <- 0
for ( i in 1: length(unique.chrs)) {
  ix.chr <- which(chrs == unique.chrs[i])
  tmp <- sort(genepos[ix.chr], index.return = T)
  ix.srt[(count + 1):(count + length(ix.chr))] <- ix.chr[tmp$ix]
  count <- count + length(ix.chr)
}

geneid.srt <- geneid[ix.srt]
genepos.str <- genepos[ix.srt]
strand.srt <- strand[ix.srt]

EnsemblGene.GenePosition <- xx[ix.srt, ]
names(EnsemblGene.GenePosition) <- c("EnsemblGene", "Chromosome", "Position")
save(EnsemblGene.GenePosition, file = file.path(dbfolder, "data", refid, "EnsemblGene.GenePosition.RData"))

xx <- read.delim(file.path(tmpfolder, "exonstartend.mrg.txt"), sep = "\t", header = F)

exgeneid <- as.character(xx[[4]])
exchr <- as.character(xx[[1]])
exstart.tmp <- as.numeric(as.character(xx[[2]])) + 1
exend.tmp <- as.character(xx[[3]])
ix.dup <- which(!duplicated(exgeneid))
exstart <- rep("", length(ix.dup))
exend <- rep("", length(ix.dup))
excount <- rep(NA, length(ix.dup))
start <- rep("", length(ix.dup))
end <- rep("", length(ix.dup))
strand1 <- rep("", length(ix.dup))
for (i in 1: (length(ix.dup) - 1)) {
  if (strand[i] == "-1") {
    strand1[i] <- "-"
  } else {
    strand1[i] <- "+"    
  }
  ix.tmp <- ix.dup[i]:(ix.dup[i+1] - 1)
  start[i] <- exstart.tmp[ix.tmp][1]
  end[i] <- exend.tmp[ix.tmp][length(ix.tmp)]
  exstart[i] <- toString(exstart.tmp[ix.tmp])
  exend[i] <- toString(exend.tmp[ix.tmp])
  excount[i] <- length(ix.tmp)
}



EnsemblGene.Structures <- cbind(exgeneid[ix.dup], exchr[ix.dup], strand1, start, end, excount, exstart, exend)[ix.srt, ]
EnsemblGene.Structures <- data.frame(EnsemblGene.Structures)
names(EnsemblGene.Structures) <- c("EnsemblGene", "Chromosome", "Strand", "geneStart", "geneEnd", "exonCount", "exonStart", "exonEnd")
save(EnsemblGene.Structures, file = file.path(dbfolder, "data", refid, "EnsemblGene.Structures.RData"))

xx <- read.delim(file.path(tmpfolder, "geneinfo.txt"), sep = "\t", header = F)
EnsemblGene.GeneInfo <- xx[ix.srt, ]
names(EnsemblGene.GeneInfo) <- c("EnsemblGene", "GeneName", "Description")
save(EnsemblGene.GeneInfo, file = file.path(dbfolder, "data", refid, "EnsemblGene.GeneInfo.RData"))

if (refid == "homo_sapiens") {
xx <- read.delim(file.path(tmpfolder, "paralogs.txt"), sep = "\t", header = F)
EnsemblGene.Paralogs <- xx
names(EnsemblGene.Paralogs) <- c("EnsemblGene", "Paralogs")
save(EnsemblGene.Paralogs, file = file.path(dbfolder, "data", refid, "EnsemblGene.Paralogs.RData"))
}
