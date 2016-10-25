vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]
ericscriptfolder <- split.vars[3]
readlength <- max(as.numeric(split.vars[4]))
refid <- as.character(split.vars[5])
dbfolder <- as.character(split.vars[6])


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


load(file.path(outputfolder,"out", paste(samplename,".chimeric.RData", sep = "")))
load(file.path(dbfolder, "data", refid, "EnsemblGene.GeneNames.RData"))
load(file.path(dbfolder, "data", refid, "EnsemblGene.Sequences.RData"))
load(file.path(dbfolder, "data", refid, "EnsemblGene.Structures.RData"))
load(file.path(outputfolder, "out", paste(samplename,".chimeric.RData", sep = "")))
id1 <- MyGF$id1
id2 <- MyGF$id2
junctions <- rep(NA, length(id1))
ids_fasta <- rep("", length(id1))
sequences.fasta <- rep("", length(id1))
fasta.file <- c()
maxgap <- 300
for (i in 1: length(id1)) {
	ix.genetable1 <- which(EnsemblGene.Structures$EnsemblGene == id1[i])
	ix.genetable2 <- which(EnsemblGene.Structures$EnsemblGene == id2[i])
	ix.gene1 <- which(GeneNames == id1[i])
	ix.gene2 <- which(GeneNames == id2[i])
	min.pos1 <- min(MyGF$pos1[[i]]) - 2*readlength
	max.pos1 <- max(MyGF$pos1[[i]]) + readlength - 1
	if (min.pos1 < 1) {min.pos1 <- 1}
	min.pos2 <- min(MyGF$pos2[[i]])
	max.pos2 <- max(MyGF$pos2[[i]]) + 2*readlength
	a <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonStart[ix.genetable1]), ",")))
	b <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonEnd[ix.genetable1]), ",")))
	strand1 <- as.character(EnsemblGene.Structures$Strand[ix.genetable1])
	if (strand1 == "+") {
		tmp.sum1 <- cumsum((b - a ))
	} else {
		tmp.sum1 <- cumsum(rev(b - a))		
	}
	exonenumber1 <- which(tmp.sum1 >= max.pos1)[1]
	if (is.na(exonenumber1)) {exonenumber1 <- length(tmp.sum1)}
	a2 <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonStart[ix.genetable2]), ",")))
	b2 <- as.numeric(unlist(strsplit(as.character(EnsemblGene.Structures$exonEnd[ix.genetable2]), ",")))
	strand2 <- as.character(EnsemblGene.Structures$Strand[ix.genetable2])
	if (strand2 == "+") {
		tmp.sum2 <- cumsum((b2 - a2))
	} else {
		tmp.sum2 <- cumsum(rev(b2 - a2))
	}
	exonenumber2 <- which(tmp.sum2 >= min.pos2)[1]
	if (is.na(exonenumber2)) {exonenumber2 <- length(tmp.sum2)}
	id.gf1 <- paste(id1[i], exonenumber1, sep = "_")
	fasta.gf1.tmp0 <- sequences[ix.gene1]
	start.end.exons <- c(0,tmp.sum1)
	fasta.gf1 <- substr(fasta.gf1.tmp0, min.pos1, (max.pos1 + maxgap - 1))	
	id.gf2 <- paste(id2[i], exonenumber2, sep = "_")
	fasta.gf2.tmp0 <- sequences[ix.gene2]
	if (max.pos2 > nchar(fasta.gf2.tmp0)) {max.pos2 <- nchar(fasta.gf2.tmp0)}
	start.end.exons <- c(0,tmp.sum2)
	fasta.gf2 <- substr(fasta.gf2.tmp0, (min.pos2 - maxgap), max.pos2)
	id.fastaGF <- paste(">",id.gf1,"----",id.gf2," junction@",nchar(fasta.gf1),sep = "")
	sequences.fasta[i] <- paste(fasta.gf1, fasta.gf2, sep = "")
	fasta.gf12 <- formatfasta(sequences.fasta[i])
	ids_fasta[i] <- paste(id.gf1,id.gf2, sep = "----")
	junctions[i] <- nchar(fasta.gf1)
	fastaGF <- c(id.fastaGF, fasta.gf12)
	fasta.file <- c(fasta.file, fastaGF)
}
save(junctions, file = file.path(outputfolder, "out", paste(samplename,".junctions.RData", sep = "")))
save(sequences.fasta, file = file.path(outputfolder, "out", paste(samplename,".sequences_fasta.RData", sep = "")))
save(ids_fasta, file = file.path(outputfolder, "out",  paste(samplename, ".ids_fasta.RData", sep = "")))
cat(fasta.file, file = file.path(outputfolder,"out", paste(samplename,".EricScript.junctions.fa",sep = "")), sep = "\n")






