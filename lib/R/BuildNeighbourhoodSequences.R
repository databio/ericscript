vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]
z <- read.delim(file.path(outputfolder,"out",paste(samplename,".intervals.pileup", sep = "")), sep = "\t", header = F)
load(file.path(outputfolder,"out",paste(samplename,".ids_filtered.RData", sep = "")))
load(file.path(outputfolder,"out",paste(samplename,".junctions.recalibrated.RData", sep = "")))
load(file.path(outputfolder, "out",  paste(samplename, ".ids_fasta.RData", sep = "")))
id.pileup <- as.character(z[,1])
pos.pileup <-  as.numeric(as.character(z[,2]))
sequence.pileup <- as.character(z[,3])
unique.ids.pileup <- unique(id.pileup)
width <- 100
fasta.file <-  c()
for (i in 1:length(id.filtered)) {
	ix.id <- which(id.pileup == id.filtered[i])
	ix.ref <- which(ids_fasta == id.filtered[i])
	junction <- junctions.recalibrated[ix.ref]
	ix.id.pileup <- which(id.pileup == id.filtered[i])
	ix.junction1 <- which(pos.pileup[ix.id.pileup] == junction)
	ix.junction2 <- which(pos.pileup[ix.id.pileup] == (junction + 1))
	seq.vec <- rep("N", width)
	pos.seq <- seq.int((junction-(width/2-1)), (junction + (width/2)))
	ix.pos.tmp <- which(pos.seq %in% pos.pileup[ix.id.pileup])
	seq.vec[ix.pos.tmp] <- sequence.pileup[ix.id.pileup]
	if ((length(ix.junction1)!=0) & (length(ix.junction2)!=0)) {
		query.sequence <- character(length = 1)
		for (ii in 1:length(seq.vec)) {
			query.sequence <- paste(query.sequence, seq.vec[ii], sep = "")
		}
		ids_fasta_query <-  paste(">", id.filtered[i],sep = "")
		fasta.file <-  c(fasta.file, c(ids_fasta_query, query.sequence))
	}
}
cat(fasta.file, sep = "\n", file = file.path(outputfolder, "out",paste(samplename,".checkselfhomology.fa", sep = "")))
cat(file.path(outputfolder, "out", paste(samplename,".checkselfhomology.fa", sep = "")), file = file.path(outputfolder, "out", ".link"))
