vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]
ericscriptfolder <- split.vars[3]
minreads <- as.numeric(split.vars[4])
MAPQ <-  as.numeric(split.vars[5])
refid <- as.character(split.vars[6])
dbfolder <- as.character(split.vars[7])

filein <- file.path(outputfolder, "out", paste(samplename, ".filtered.out", sep = ""))
xx <- readLines(filein, n = 1)
if (length(xx) == 0) {
  myflag <- 0
  cat(myflag, file = file.path(outputfolder, "out", ".ericscript.flag"))
  stop("No lines available in ",filein,". No discordant reads found with MAPQ set to ", MAPQ, ". Try to decrease MAPQ parameter and run again EricScript. Exit!")
} else {
  myflag <- 1
  cat(myflag, file = file.path(outputfolder, "out", ".ericscript.flag"))
}
x <- read.delim(filein, sep = "\t", header = F)
load(file.path(dbfolder,"data", refid, "EnsemblGene.GenePosition.RData"))
flag <- x[,1]
#ix.flag <- which((flag > 63 & flag < 70) | (flag > 95 & flag < 118) | flag == 161 | flag == 181) 
ix.flag <- which((flag > 63 & flag < 70) | (flag > 95 & flag < 112) | flag == 161) 
id_1 <- as.character(x[ix.flag,2])
id_2 <- as.character(x[ix.flag,5])
pos_1 <- as.numeric(as.character(x[ix.flag,3]))
pos_2 <- as.numeric(as.character(x[ix.flag,6]))
rm(x)
genename <- as.character(EnsemblGene.GenePosition$EnsemblGene)
id1 <- c()
id2 <- c()
nreads <- c()
diffpos <- c()
generef <- unique(id_1)
for (i in 1: length(generef)) {
  ix.generef <- which(genename == generef[i])	
  ix.gene <- which(id_1 == generef[i])
  tmp <- sort(summary(as.factor(id_2[ix.gene]), maxsum = length(unique(id_2[ix.gene]))),  decreasing = T)
  tmp.genename <- names(tmp)
  tmp.weight <- as.numeric(tmp)
  if ((max(tmp.weight) >= minreads) & (length(which(tmp.weight >= minreads)) <= 10)) {
    ix.maxnodes <- which(tmp.weight >= minreads)
    tmp.genename <- tmp.genename[ix.maxnodes]
    tmp.weight <- tmp.weight[ix.maxnodes]
    for (j in 1:length(tmp.weight)) {
      ix.genelink <- which(genename == tmp.genename[j])
      if (length(ix.genelink)!=0) {
        id1 <- c(id1, generef[i])
        id2 <- c(id2, tmp.genename[j])
        nreads <- c(nreads, tmp.weight[j])
        diffpos <- c(diffpos, abs(ix.generef-ix.genelink))
      }
    }
  } 
  
}
## filter paralogs if paralogs exist

if (file.exists(file.path(dbfolder,"data", refid, "EnsemblGene.Paralogs.RData"))) {
  
  load(file.path(dbfolder, "data", refid, "EnsemblGene.Paralogs.RData"))  
  paralogs.flag <- rep(0, length(id1))
  
  if (length(id1) == 0) {
    myflag <- 0
    cat(myflag, file = file.path(outputfolder, "out", ".ericscript.flag"))
    stop("No discordant reads found with minimum reads set to ", minreads, ". Exit!")
  }
  for (i in 1: length(id1)) {
    ix.paralogs <- which(EnsemblGene.Paralogs$EnsemblGene == id1[i])
    paralogs <- as.character(EnsemblGene.Paralogs$Paralogs[ix.paralogs])
    if (length(grep(id2[i], paralogs)) > 0) {
      paralogs.flag[i] <- 1
    }
  }
  ##
  paralogs.filter <- which(paralogs.flag == 0)
  id1f <- id1[paralogs.filter]
  id2f <- id2[paralogs.filter]
  nreadsf <- nreads[paralogs.filter]
  diffposf <- diffpos[paralogs.filter]
  if (length(id1) == 0) {
    myflag <- 0
    cat(myflag, file = file.path(outputfolder, "out", ".ericscript.flag"))
    stop("No discordant reads found with minimum reads set to ",minreads,". Exit!")
  } 
  nfus <- length(id1f)
  MyGF <- vector("list", 6)
  names(MyGF) <- c("id1", "id2", "nreads", "pos1", "pos2", "diffpos")
  MyGF$id1 <- id1f
  MyGF$id2 <- id2f
  MyGF$nreads <- nreadsf
  MyGF$diffpos <- diffposf
  MyGF$pos1 <- vector("list", nfus)
  MyGF$pos2 <- vector("list", nfus)
  for (i in 1: (nfus)) {
    ix.pos <- which((id_1 == id1f[i]) & (id_2 ==id2f[i]))
    MyGF$pos1[[i]] <- pos_1[ix.pos]
    MyGF$pos2[[i]] <- pos_2[ix.pos]
  }
} else {
  
  nfus <- length(id1)
  MyGF <- vector("list", 6)
  names(MyGF) <- c("id1", "id2", "nreads", "pos1", "pos2", "diffpos")
  MyGF$id1 <- id1
  MyGF$id2 <- id2
  MyGF$nreads <- nreads
  MyGF$diffpos <- diffpos
  MyGF$pos1 <- pos_1
  MyGF$pos2 <- pos_2
  
}

save(MyGF, file = file.path(outputfolder, "out", paste(samplename, ".chimeric.RData", sep = "")))

