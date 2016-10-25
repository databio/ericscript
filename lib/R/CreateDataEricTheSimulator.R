options(stringsAsFactors=F)
vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
refid <- split.vars[1]
dbfolder <- split.vars[2]
tmpfolder <- split.vars[3]

## read transcript genomic info
xx <-  read.delim(file.path(tmpfolder, "transcripts.txt"), sep = "\t", header = F)
geneid <- xx[[1]]
transcriptid <- xx[[2]]
exonstart <- xx[[3]]
exonend <- xx[[4]]
chr <- xx[[5]]
strandtmp <- xx[[6]]
strand <- rep("+", length(strandtmp))
strand[strandtmp == "-1"] <- "-"
## read transcript cdna
xxseq <- scan(file.path(tmpfolder, "transcripts.fa"), what = list(seq="", id=""), sep = "\t", quiet = T)
seqtmp <- xxseq[[1]]
transcriptid.seqtmp <- xxseq[[2]]
rm (xx, xxseq)
unique.transcriptid <- unique(transcriptid)
EnsemblGene.Structures <- c()
GeneNames <- rep("", length(unique.transcriptid))
sequences <- rep("", length(unique.transcriptid))
for (i in 1: length(unique.transcriptid)) {
  ix <- which(transcriptid == unique.transcriptid[i])
  ixseq <- which(transcriptid.seqtmp == unique.transcriptid[i])
  if (length(ixseq) > 0) {
    sequences[i] <-  seqtmp[ixseq]    
  } else {
    sequences[i] <-  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
  }
  ixsrt <- sort(exonstart[ix], index.return = T)$ix
  genestart <- min(c(exonstart[ix], exonend[ix]))
  geneend <- max(c(exonstart[ix], exonend[ix]))
  exonStart <- toString(exonstart[ix[ixsrt]])
  exonEnd <- toString(exonend[ix[ixsrt]])
  exoncount <- length(ix)
  mychr <- unique(chr[ix])
  mystrand <- unique(strand[ix])
  GeneNames[i] <- unique(geneid[ix])
  EnsemblGene.Structures <- rbind(EnsemblGene.Structures, c(unique.transcriptid[i], mychr, mystrand, genestart, geneend, exoncount, exonStart, exonEnd))
}
colnames(EnsemblGene.Structures) <- c("EnsemblGene", "Chromosome", "Strand", "geneStart", "geneEnd", "exonCount", "exonStart", "exonEnd")
EnsemblGene.Structures <- data.frame(EnsemblGene.Structures)
save(EnsemblGene.Structures, GeneNames, sequences, file = file.path(dbfolder, "data", refid, "EnsemblGene.Transcripts.RData"))




