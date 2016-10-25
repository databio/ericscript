vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]
readlength <-  as.numeric(split.vars[3])
verbose <- as.numeric(split.vars[4])
grep.readlength <- c("grep")
for (i in 1: length(readlength)) {
	grep.readlength <- paste(grep.readlength, " -v -e MD:Z:", readlength[i], sep = "")
}

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

TryRecalibration <- function(outputfolder, verbose) {
    
	if (verbose == 0) {
		x <- system(paste("blat -tileSize=8 -fine", file.path(outputfolder,"out",".tmp.ref.fa"), file.path(outputfolder, "out", ".link"), file.path(outputfolder, "out",".recalibrated.junctions.blat"), " 1>> ", file.path(outputfolder, "out",".ericscript.log")))
	} else {
	  	x <- system(paste("blat -tileSize=8 -fine", file.path(outputfolder,"out",".tmp.ref.fa"), file.path(outputfolder, "out", ".link"), file.path(outputfolder, "out",".recalibrated.junctions.blat")))
	}
	yy <- readLines(file.path(outputfolder, "out", ".recalibrated.junctions.blat"), n = 6)
	if (length(yy) > 5) {
		xx <- read.delim(file.path(outputfolder, "out", ".recalibrated.junctions.blat"), sep = "\t", skip = 5, header = F)
		gapsize <- xx[,8]
	}
	if (all(gapsize <= 3) | (length(yy) <= 5)) {
		if (verbose == 0) {
			x <- system(paste("blat -tileSize=8", file.path(outputfolder,"out",".tmp.ref.fa"), file.path(outputfolder, "out", ".link"), file.path(outputfolder, "out",".recalibrated.junctions.blat"), " 1>> ", file.path(outputfolder, "out",".ericscript.log")))
		} else {
			x <- system(paste("blat -tileSize=8", file.path(outputfolder,"out",".tmp.ref.fa"), file.path(outputfolder, "out", ".link"), file.path(outputfolder, "out",".recalibrated.junctions.blat")))
		}
		
	}
	
}

load(file.path(outputfolder,"out",paste(samplename,".junctions.RData", sep = "")))
load(file.path(outputfolder,"out",paste(samplename,".ids_fasta.RData", sep = "")))
load(file.path(outputfolder,"out",paste(samplename,".sequences_fasta.RData", sep = "")))
cat(file.path(outputfolder,"out",".tmp.query.fa"), file = file.path(outputfolder,"out",".link"))
sequences <- sequences.fasta
recal.left <- rep(0, length(ids_fasta))
recal.right <- rep(0, length(ids_fasta))
count.recal <- rep(0, length(ids_fasta))
count.total <- rep(0, length(ids_fasta))
sequences.recal <- sequences
junctions.recalibrated <- junctions
for (i in 1:length(ids_fasta)) {
	junction.tmp <- junctions[i]
	x <- system(paste("samtools view ", file.path(outputfolder, "aln", paste(samplename,".remap.sorted.bam", sep = ""))," ", ids_fasta[i], " | ", grep.readlength, " | awk '((($5==0) && ($6==\"*\")) || ($5>=0))' | awk '{ print \">\" $1\"_\"$2,\"\\n\"$10}' > ", file.path(outputfolder,"out",".tmp.query.fa"),sep = ""))
	xx <- readLines(file.path(outputfolder,"out",".tmp.query.fa"), n = 1)
	if (length(xx) != 0) {
		x <- cat(paste(">",ids_fasta[i],sep=""), sequences[i],sep = "\n", file = file.path(outputfolder,"out",".tmp.ref.fa"))
		try.recal <- TryRecalibration(outputfolder, verbose)
		rm(xx)
		yy <- readLines(file.path(outputfolder, "out", ".recalibrated.junctions.blat"), n = 6)
		if (length(yy) > 5) {
			xx <- read.delim(file.path(outputfolder, "out", ".recalibrated.junctions.blat"), sep = "\t", skip = 5, header = F)
			gapsize <- xx[,8]
			gapstarts <- strsplit(as.character(xx[,21]), ",")
			blocksize <- strsplit(as.character(xx[,19]), ",")
			if (any(gapsize > 3)) {
				ix.0 <- which(gapsize > 3)
				if (length(ix.0) > 0) {
					gapstarts1 <- gapstarts[ix.0]
					gapstarts.tmp1 <- c()
					gapstarts.tmp2 <- c()
					for (jgap in 1:length(gapstarts1)) {
						ccc <- as.numeric(gapstarts1[[jgap]])
						gapstarts.tmp1 <- c(gapstarts.tmp1, ccc[1])
						gapstarts.tmp2 <- c(gapstarts.tmp2, ccc[2])
					}
					ix.gap.in.junct <- ix.0[which((gapstarts.tmp1 <= junction.tmp) & (gapstarts.tmp2 >= (junction.tmp - 10 + 1)))]
					gaps <- gapsize[ix.gap.in.junct]
					unique.gaps <- unique(gaps)
					rr <- tabulate(gaps)
					max.rr <- max(rr)
					if( max.rr >= 1) {
						gap.length <- which.max(rr)
						ix.gaps <- which(gapsize == gap.length) 
						a <- rep(0, length(ix.gaps))
						b <- rep(0, length(ix.gaps))
						aa <- rep(0, length(ix.gaps))
						for (jj in 1:length(ix.gaps)) {
							gapstarts.tmp <- as.numeric(gapstarts[[ix.gaps[jj]]])
							blocksize.tmp <- as.numeric(blocksize[[ix.gaps[jj]]])
							a[jj] <- gapstarts.tmp[1] + blocksize.tmp[1] - 1
							b[jj] <- gapstarts.tmp[2]
							aa[jj] <- gapstarts.tmp[1]
						}
						max.a <- max(tabulate(a))
						max.b <- max(tabulate(b))
						my.a <- which.max(tabulate(a))
						my.b <- which.max(tabulate(b))
						if ((abs(max.a-max.b)/max.a) < 0.31) {
							count.total[i] <- length(gaps)
							count.recal[i] <- max.rr
							recal.left[i] <- my.a + 1
							recal.right[i] <- my.b + 1 
							sequences.recal[i] <- paste(substr(sequences[i], 1, recal.left[i]), substr(sequences[i], recal.right[i], nchar(sequences[i])), sep = "")
							junctions.recalibrated[i] <- recal.left[i]
						}
					}
					
				}
			}
		}
		
	} 
}
ids_fasta.recalibrated <- paste(">", ids_fasta, " junction@", junctions.recalibrated, sep = "")
ref.recalibrated <- c()
for (i in 1:length(ids_fasta.recalibrated)) {
	ref.recalibrated <- c(ref.recalibrated, c(ids_fasta.recalibrated[i], formatfasta(sequences.recal[i])))
}
write(ref.recalibrated, file = file.path(outputfolder,"out",paste(samplename,".EricScript.junctions.recalibrated.fa", sep = "")), ncolumns = 1, sep = "")
Recalibrated.Data <- cbind(recal.left, recal.right, junctions, count.recal, count.total)
colnames(Recalibrated.Data) <- c("Left_Junction", "Right_Junction", "Junction", "Recal_Count", "Total_Count")
save(sequences.recal, file = file.path(outputfolder,"out",paste(samplename,".sequences.recalibrated.RData", sep = "")))
save(Recalibrated.Data, file = file.path(outputfolder,"out",paste(samplename,".Recalibrated.Data.RData", sep = "")))
save(junctions.recalibrated, file = file.path(outputfolder,"out",paste(samplename,".junctions.recalibrated.RData", sep = "")))


