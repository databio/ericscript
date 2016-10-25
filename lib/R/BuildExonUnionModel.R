vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
ericscriptfolder <- split.vars [1]
refid <- split.vars[2]
dbfolder <- split.vars [3]
tmpfolder <- split.vars [4]


formatfasta <- function(myfasta, step = 50) {  
  totalchar <- nchar(myfasta)
  if (totalchar > step) {
    steps <- seq(1, totalchar, by = step)
    newfasta <- rep("", (length(steps) - 1))
    for (j in 1: (length(steps) - 1)) {
      aa <- substr(myfasta, steps[j], (steps[j] + (step - 1)))
      newfasta[j] <- aa 
    }
    if ((totalchar - tail(steps, n = 1)) > 0) {
      newfasta <-  c(newfasta, substr(myfasta, steps[j+1], totalchar))
    }
  } else
  {
    newfasta <- substr(myfasta, 1, totalchar)
  }
  return(newfasta)
}

convertToComplement <- function(x) {
  
  bases <- c("A", "C", "G", "T")
  #xx <- unlist(strsplit(toupper(x), NULL))
  xx <- rev(unlist(strsplit(toupper(x), NULL)))
  paste(unlist(lapply(xx, function(bbb) {
    if (bbb=="A") compString <- "T"
    if (bbb=="C") compString <- "G"
    if (bbb=="G") compString <- "C"
    if (bbb=="T") compString <- "A"
    if (!bbb %in% bases) compString <- "N"
    return(compString)
  })), collapse="")
  
}

refid.folder <- file.path(dbfolder, "data", refid)
if (file.exists(refid.folder) == F) {
  dir.create(refid.folder)
}
x <- scan(file.path(tmpfolder, "subseq.fa"), what = "", quiet = T)
x.bed <- read.delim(file.path(tmpfolder, "exonstartend.mrg.txt"), sep = "\t", header = F)
refid.bed <- paste(as.character(x.bed[[1]]), paste((as.numeric(as.character(x.bed[[2]])) + 1), as.character(x.bed[[3]]), sep = "-"), sep = ":")
tmp <- grep(">", x)
genomicreg <- substr(x[tmp], 2, nchar(x[tmp]))
sequences.tmp <- rep("", length(tmp))
for (i in 1: (length(tmp) - 1)) { 
  sequences.tmp[i] <- gsub(", ", "", toString(x[(tmp[i] + 1):(tmp[i+1] - 1)]))
}
sequences.tmp[length(tmp)] <- gsub(", ", "", toString(x[(tmp[length(tmp)] + 1): length(x)]))
genenames.tmp <- as.character(x.bed[[4]])
unique.genenames <- unique(genenames.tmp)
strand.tmp <- read.delim(file.path(tmpfolder, "strand.txt"), sep = "\t", header = F)
strand <- strand.tmp[[2]]
sequences <- rep("", length(unique.genenames))
for (i in 1: length(unique.genenames)) {
    genenames1 <- paste(">", unique.genenames[i], sep = "")
    ix.gene <- which(genenames.tmp == unique.genenames[i])
    ix.refid <- which(genomicreg %in% refid.bed[ix.gene])
    if (strand[i] == "-1") {
      seqtmp0 <- gsub(", ", "", toString(sequences.tmp[ix.refid]))
      sequences[i] <- convertToComplement(seqtmp0)
    } else {
      sequences[i] <- gsub(", ", "", toString(sequences.tmp[ix.refid]))
    }
    if (nchar(sequences[i])  == 0) {
      sequences[i] <- "NNNNNNN"
    }
      
    if (i == 1) {
      cat(genenames1, file = file.path(refid.folder, "EnsemblGene.Reference.fa"), append = F, sep = "\n")
    } else {
      cat(genenames1, file = file.path(refid.folder, "EnsemblGene.Reference.fa"), append = T, sep = "\n")
    }
    cat(formatfasta(sequences[i]), file = file.path(refid.folder, "EnsemblGene.Reference.fa"), append = T, sep = "\n")  
}
ix.emptyseq <- which(nchar(sequences) == 0) 
GeneNames <- unique.genenames
if (length(ix.emptyseq) > 0) {
GeneNames <- GeneNames[-ix.emptyseq]
sequences <- sequences[-ix.emptyseq]
}
save(GeneNames, file = file.path(refid.folder, "EnsemblGene.GeneNames.RData"))
save(sequences, file = file.path(refid.folder, "EnsemblGene.Sequences.RData"))

