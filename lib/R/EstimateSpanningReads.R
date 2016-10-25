## EstimateSpanningReads v2: exclude soft clipped reads
vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]
readlength <- max(as.numeric(split.vars[3]))
load(file.path(outputfolder,"out",paste(samplename,".junctions.recalibrated.RData", sep = "")))
load(file.path(outputfolder, "out", paste(samplename, ".ids_fasta.RData", sep = "")))
load(file.path(outputfolder, "out", "isize.RData"))
mysigma <- readlength/4
myref <- sort(dnorm(seq(1, readlength, by = 1), readlength/2, mysigma), decreasing = T)
DataMatrix <- matrix(NA, nrow = length(ids_fasta), ncol = 9)
for (i in 1: length(ids_fasta)) {
	junction.tmp <- junctions.recalibrated[i]
	x <- system(paste("samtools view ", file.path(outputfolder, "aln", paste(samplename,".remap.recal.sorted.rmdup.bam", sep = ""))," ", ids_fasta[i],  ":", junction.tmp, "-", junction.tmp + 1, " | awk '($5>0)'  | cut -f 1,4,6 > ", file.path(outputfolder,"out",".spanningreads.sam"),sep = ""))
	x <- system(paste("samtools view ", file.path(outputfolder, "aln", paste(samplename,".remap.recal.sorted.rmdup.bam", sep = ""))," ", ids_fasta[i]," | awk '($9>0)'  | cut -f 9 - > ", file.path(outputfolder,"out",".insertsize.sam"),sep = ""))
	x <- system(paste("samtools view ", file.path(outputfolder, "aln", paste(samplename,".remap.recal.sorted.rmdup.bam", sep = ""))," ", ids_fasta[i]," | awk '($4<", junction.tmp, ") && ($8>",junction.tmp+1,")' | cut -f 1,4 -  > ", file.path(outputfolder,"out",".crossingreads.sam"),sep = ""))
	
	spanningreads.tmp <- scan(file.path(outputfolder,"out",".spanningreads.sam"), sep = "\t", what = list("", 1, ""), quiet = T)
	ix.nosc <- sort(intersect(grep("^[0-9]*S", spanningreads.tmp[[3]], perl = T, invert = T), grep("[0-9]*S$", spanningreads.tmp[[3]], perl = T, invert = T)))
	spanningreads <- vector("list", length = 3)
  if (length(ix.nosc) > 0) {
    spanningreads[[1]] <-  spanningreads.tmp[[1]][ix.nosc]   
    spanningreads[[2]] <-  spanningreads.tmp[[2]][ix.nosc]   
    spanningreads[[3]] <-  spanningreads.tmp[[3]][ix.nosc]   
  } else {
    spanningreads <- spanningreads.tmp
  }
	crossingreads <- scan(file.path(outputfolder,"out",".crossingreads.sam"), sep = "\t", what = list("", 1), quiet = T)
	insert.size <- mean(abs(as.numeric(readLines(file.path(outputfolder,"out",".insertsize.sam"))))) - readlength
	id.spanningreads <- spanningreads[[1]]
	pos.spanningreads <- spanningreads[[2]]
	id.crossingreads <- crossingreads[[1]]
	pos.crossingreads <- crossingreads[[2]]
	spanning.score <- 0
	edge.score <- 0
	range.pos.crossingreads <- 0
	if (length(pos.crossingreads) > 0) {
		range.pos.crossingreads <- junction.tmp - min(pos.crossingreads)
	}
	n.crossingreads <- length(which(id.crossingreads %in% id.spanningreads == F))
	n.spanningreads <- 0
	gjs <- 0
	unique.score <- 0
	us.prob <- 0
	insertsize.score <- 0
	if (length(pos.spanningreads) > 0) {
		pos <- pos.spanningreads - junction.tmp + readlength
		us.pos <- tabulate(pos, nbins = readlength)
		
		if (sum(us.pos) > 0) {
			us.mult <- floor(sum(us.pos)/readlength)
			us.residuals <- sum(us.pos)/readlength - floor(sum(us.pos)/readlength) 
			us.refdistr <- rep(us.mult, readlength)
			if (us.residuals > 0) {
				for (kk in 1: (sum(us.pos) - us.mult*readlength)) {
					us.refdistr[kk] <- us.refdistr[kk] + 1  
				}
			} else {
				us.refdistr[1] <- us.refdistr[1] + 1
				us.refdistr[2] <- us.refdistr[2] - 1
			}
			
			ff <- which(sort(us.pos!=0))
			us.prob <- 1- sum(abs(sort(us.pos)[ff]-sort(us.refdistr)[ff]))/sum(us.pos)
		}

		mynorm <- sum(myref[1:length(unique(pos))])
		prob <- sum(dnorm(unique(pos),  readlength/2, mysigma))
		gjs <- prob/mynorm
		
		
		left.spanningreads <- pos.spanningreads[(pos.spanningreads <= (junction.tmp - round(readlength/3)))]
		right.spanningreads <- pos.spanningreads[(pos.spanningreads > (junction.tmp - round(readlength/3)))]
		n.left.spanningreads <- length(left.spanningreads )
		n.right.spanningreads <- length(right.spanningreads)
		spanning.score <- 1- abs(n.left.spanningreads - n.right.spanningreads)/(n.left.spanningreads + n.right.spanningreads)
		
		if (length(left.spanningreads) > 0) {
			left.score <- mean((junction.tmp - readlength)-left.spanningreads)
		} else {
			left.score <- 0
		}
		if (length(right.spanningreads) > 0) {
			right.score <- mean(right.spanningreads - junction.tmp)
		} else {
			right.score <- 0
		}
		edge.score <- 1 - 1.1^(mean(c(left.score,right.score)))
		insertsize.score <- dnorm(insert.size, isize.mean, isize.sd)/dnorm(isize.mean, isize.mean, isize.sd)
		n.spanningreads <- length(pos.spanningreads)
		
	}
	
	DataMatrix[i,] <- c(ids_fasta[i], n.crossingreads , insert.size, n.spanningreads, range.pos.crossingreads, edge.score, gjs, us.prob, insertsize.score)
	
}
colnames(DataMatrix) <- c("id", "nreads","mean_ins_size","nreads_junc", "rangepos", "edgescore", "gjs", "uniformity.score", "isize.score")
save(DataMatrix, file = file.path(outputfolder,"out",paste(samplename, ".DataMatrix.RData", sep = "")))

filecon <- file(file.path(outputfolder,"out", paste(samplename, ".intervals", sep = "")), open = "w")
ix.filter <- sort(unique(intersect(which(DataMatrix[,4] > 0), which(DataMatrix[,3] != "NaN"))))
if (length(ix.filter) > 0) {
  width <- 100
  id.filtered <- ids_fasta[ix.filter]
  save(id.filtered, file = file.path(outputfolder, "out",paste(samplename,".ids_filtered.RData", sep = "")))
  for (i in 1:length(ix.filter)) {
	  ix.ref <- ix.filter[i]
	  junction <- junctions.recalibrated[ix.ref]
	  pileup.interval <- seq.int((junction - (width/2 - 1)), (junction + (width/2)))
	  pileup.interval[which(pileup.interval < 1)] <- 1
	  pileup.interval <- unique(pileup.interval)
	  cat(paste(rep(id.filtered[i], length(pileup.interval)), pileup.interval, sep = " "), file = filecon, sep = "\n", append = T)
  }
  close(filecon)
  myflag <- 1
	cat(myflag, file = file.path(outputfolder, "out", ".ericscript.flag"))
} else {
  myflag <- 0
  cat(myflag, file = file.path(outputfolder, "out", ".ericscript.flag"))
	stop("No chimeric transcripts found. Exit!")
} 

