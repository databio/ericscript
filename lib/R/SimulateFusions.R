### simulate data [revised].

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))

readlength <- as.numeric(split.vars[1])
outputfolder <- split.vars[2]
ericscriptfolder <- split.vars[3]
verbose <- as.numeric(split.vars[4])
ins.size <- as.numeric(split.vars[5])
sd.inssize <- as.numeric(split.vars[6])
ngenefusion <- as.numeric(split.vars[7])
min.coverage <- as.numeric(split.vars[8])
max.coverage <- as.numeric(split.vars[9])
nsims <- as.numeric(split.vars[10])
BE.data <- as.numeric(split.vars[11])
IE.data <- as.numeric(split.vars[12])
background.data_1 <- as.character(split.vars[13])
background.data_2 <- as.character(split.vars[14])
nreads.background <- as.numeric(split.vars[15])
dbfolder <- as.character(split.vars[16])
refid <- as.character(split.vars[17])
  
mysyndata <- file.exists(file.path(dbfolder, "data", refid, "EnsemblGene.Transcripts.RData"))
if (mysyndata == T) {
  cat("[EricScript simulator] Load genes data ...")
  load(file.path(dbfolder, "data", refid, "EnsemblGene.Transcripts.RData"))
} else {
  cat( paste("[EricScript simulator] You need to download", refid, "data before running EricScript Simulator. Exit.\n"))
  system(paste("rm -r", outputfolder))
  quit()
}

# myurl <- "http://dl.dropbox.com/u/3629305/EnsemblGene.Transcripts.RData"
# if (mysyndata == T) {
# 	cat("[EricScript simulator] Load genes data ...")
# 	load(file.path(dbfolder, "data", "EnsemblGene.Transcripts.RData"))
# } else {
# 	cat("[EricScript simulator] Retrieving genes data ...")
# 	download.file(myurl, destfile = file.path(dbfolder, "data", "EnsemblGene.Transcripts.RData"), quiet = T)
# 	load(file.path(dbfolder, "data", "EnsemblGene.Transcripts.RData"))
# 	cat(" done.\n")
# 	cat("[EricScript simulator] Load genes data ...")
# 	
# }

flag.background <- 0
if (nchar(background.data_1) > 2 & nchar(background.data_2) > 2) {
	flag.background <- 1
}

dataset <- c()
if (BE.data == 1) {
	dataset <- c(dataset, "BE")
}
if (IE.data == 1) {
	dataset <- c(dataset, "IE")
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

## evaluate n.backgound reads


TranscriptNames <- as.character(EnsemblGene.Structures$EnsemblGene)
acceptable.chrs <-  c(seq(1,22), "X", "Y")
mycoverage <- seq(min.coverage, max.coverage, length.out = ngenefusion)
minlength <- ins.size + 2*sd.inssize
if (refid == "homo_sapiens") {
  ix.geneok <- which((EnsemblGene.Structures$Chromosome %in% acceptable.chrs))
} else {
  ix.geneok <- seq(1, length(EnsemblGene.Structures$Chromosome))
}
genenameok <- as.character(EnsemblGene.Structures$EnsemblGene)[ix.geneok]
strandok <- as.character(EnsemblGene.Structures$Strand)[ix.geneok]
ix.goodseq <- which((nchar(sequences) > 2*minlength) & (TranscriptNames %in% genenameok) & is.na(GeneNames) == F)
sequences <- sequences[ix.goodseq]
GeneNames <- GeneNames[ix.goodseq]
TranscriptNames <- TranscriptNames[ix.goodseq] 


formatted.count.tmp <-paste("00000", seq(1, nsims), sep = "")
formatted.count <- substr(formatted.count.tmp, nchar(formatted.count.tmp) - 4, nchar(formatted.count.tmp))
formatted.count.tmp <-paste("00000", seq(1, ngenefusion), sep = "")
formatted.count.fusions <- substr(formatted.count.tmp, nchar(formatted.count.tmp) - 4, nchar(formatted.count.tmp))


for (tt in 1: length(dataset)) {
    dir.create(file.path(outputfolder, dataset[tt]))  
	dir.create(file.path(outputfolder, dataset[tt], "data"))  
	dir.create(file.path(outputfolder, dataset[tt], "reads"))  
    
    for (jj in 1: nsims) {
		
		dir.create(file.path(outputfolder, dataset[tt], "data", paste("sim", formatted.count[jj], sep = "_"))) 
		dir.create(file.path(outputfolder, dataset[tt], "reads", paste("sim", formatted.count[jj], sep = "_"))) 
	}
}
cat(" done.\n")

for (jj in 1: nsims) {
	cat("[EricScript simulator] Generating synthetic dataset", formatted.count[jj], "...")
	
	if (flag.background == 1) {

    myrandomseed <- round(runif(1, 1, 1000))
	  system(paste("seqtk sample -s", myrandomseed, " background.data_1 ", nreads.background, " > ", file.path(outputfolder, "background.reads.1.fq") , sep = ""))
	  system(paste("seqtk sample -s", myrandomseed, " background.data_2 ", nreads.background, " > ", file.path(outputfolder, "background.reads.2.fq") , sep = ""))
	}
	
	ix.gene1 <- rep(0,ngenefusion)
	ix.gene2 <- rep(0,ngenefusion)
	strand1 <- rep(0,ngenefusion)
	strand2 <- rep(0,ngenefusion)
	flag <- 1
	mycount <- 0
	while (flag == 1) {
		trans1 <- sample(TranscriptNames, ngenefusion)
		for (ii in 1: ngenefusion) {
			ix.gene1[ii] <- which(TranscriptNames == trans1[ii])
			strand1[ii] <- strandok[which(genenameok == trans1[ii])]
		}
		gene1 <- GeneNames[ix.gene1]
		trans2 <- sample(TranscriptNames, ngenefusion)
		for (ii in 1: ngenefusion) {
			ix.gene2[ii] <- which(TranscriptNames == trans2[ii])
			strand2[ii] <- strandok[which(genenameok == trans2[ii])]
		}
		gene2 <- GeneNames[ix.gene2]
		ix.gene12 <- c(ix.gene1, ix.gene2)
		if (length(unique(ix.gene12)) == 2*ngenefusion & length(unique(GeneNames[ix.gene12])) == 2*ngenefusion) {
			flag <- 0
		} else
		{flag <- 1}
	}
	
	sequence1 <- sequences[ix.gene1]
	sequence2 <- sequences[ix.gene2]
	
	if ("BE" %in% dataset) {
		
		myref <- c()
		junction1.tot <- c()
		junction2.tot <- c()
		id.fusions <- rep(0, length(sequence1))
		sequence.fusions <- rep(0, length(sequence1))
		sequence.fusions.50bp <- rep(0, length(sequence1))
		
		for (i in 1: length(sequence1)) {
			
			tmp <- seq.int(100,(nchar(sequence1[i]) - 100))
			junction1 <- sample(tmp,1)
			junction1.tot <- c(junction1.tot, junction1)
			tmp <- seq.int(100,(nchar(sequence2[i]) - 100))
			junction2 <- sample(tmp,1)
			junction2.tot <- c(junction2.tot, junction2)
			sequence.fusions[i] <- paste(substr(sequence1[i], 1, junction1), substr(sequence2[i], junction2 + 1, nchar(sequence2[i])), sep = "")
			sequence.fusions.50bp[i] <- paste(substr(sequence1[i], (junction1 - 49), junction1), substr(sequence2[i], junction2 + 1, (junction2 + 50)), sep = "")
			id.fusions[i] <- paste(">", paste(gene1[i], gene2[i], sep = "----"), sep = "")
			myref <- c(myref, c(id.fusions[i], sequence.fusions[i]))
			myref.single <- c(id.fusions[i], formatfasta(sequence.fusions[i]))
			cat(myref.single, file = file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), paste("myref", formatted.count.fusions[i], ".fa", sep = "")), sep = "\n")
		}
		
		GeneFusions <- list()
		GeneFusions[[1]] <-  gene1
		GeneFusions[[2]] <-  gene2
		GeneFusions[[3]] <- junction1.tot
		GeneFusions[[4]] <- junction2.tot
		GeneFusions[[5]] <- sequence.fusions.50bp
		GeneFusions[[6]] <- mycoverage
		GeneFusions[[7]] <- trans1
		GeneFusions[[8]] <- trans2
		names(GeneFusions) <- c("gene1", "gene2", "junction1", "junction2", "junctionseq", "coverage", "trans1", "trans2")
		save(GeneFusions, file = file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "GeneFusions.RData"))
		
		system(paste(">", file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.1.fq")))
		system(paste(">", file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.2.fq")))
		
		for (i in 1: ngenefusion) {
			mynreads <- round(mycoverage[i]*nchar(sequence.fusions[i])/(2*readlength))
			if (verbose == 0) {
				system(paste("wgsim -d ", ins.size, " -r 0.0001 -R 0.001 -s ", sd.inssize, " -N ", mynreads, " -1 ", readlength, " -2 ", readlength," ", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), paste("myref", formatted.count.fusions[i], ".fa", sep = "")), " " ,file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq")," ", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"), " 2>> ", file.path(outputfolder, "wgsim.log"), " 1>> ", file.path(outputfolder, "wgsim.log"), sep = ""))
			} else {
				system(paste("wgsim -d ", ins.size, " -r 0.0001 -R 0.001 -s ", sd.inssize, " -N ", mynreads, " -1 ", readlength, " -2 ", readlength," ", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), paste("myref", formatted.count.fusions[i], ".fa", sep = "")), " " ,file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq")," ", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"), sep = ""))
				
			}
			
			system(paste("cat", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq"), ">>",  file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.1.fq")))
			system(paste("cat", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"), ">>",  file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.2.fq")))
			
		}
		
		if (flag.background == 1) {
			system(paste("cat ", file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.1.fq"), " ",  file.path(outputfolder, "background.reads.1.fq"), " > ", file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "total.reads.1.fq"), sep = ""))
			system(paste("cat ", file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.2.fq"), " ",  file.path(outputfolder, "background.reads.2.fq"), " > ", file.path(outputfolder, "BE", "reads", paste("sim", formatted.count[jj], sep = "_"), "total.reads.2.fq"), sep = ""))
		}
		
		system(paste("rm", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq"))) 
		system(paste("rm", file.path(outputfolder, "BE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"))) 
		
		
	}
	
	if ("IE" %in% dataset) {
		
		
		myref <- c()
		
		Gene.Table <- EnsemblGene.Structures
		junction1.tot <- c()
		junction2.tot <- c()
		id.fusions <- rep(0, length(sequence1))
		sequence.fusions <- rep(0, length(sequence1))
		sequence.fusions.50bp <- rep(0, length(sequence1))
		genename.table <- as.character(Gene.Table[,1])
		
		
		for (i in 1: length(sequence1)) {
			ix.genename.table <- which(genename.table == trans1[i])
			start.exons <- as.numeric(unlist(strsplit(as.character(Gene.Table[ix.genename.table, 7]), ",")))
			end.exons <- as.numeric(unlist(strsplit(as.character(Gene.Table[ix.genename.table, 8]), ",")))
			strand <- as.character(Gene.Table[ix.genename.table, 3])
			if (strand == "+") {
				tmp <- cumsum((end.exons - start.exons))
			} else {
				tmp <- cumsum(rev(end.exons - start.exons))
			}
			if (length(tmp) > 1) {
				junction1 <- sample(tmp,1)
			} else {
				junction1 <- tmp
			}
			junction1.tot <- c(junction1.tot, junction1)
			ix.genename.table <- which(genename.table == trans2[i])
			start.exons <- as.numeric(unlist(strsplit(as.character(Gene.Table[ix.genename.table, 7]), ",")))
			end.exons <- as.numeric(unlist(strsplit(as.character(Gene.Table[ix.genename.table, 8]), ",")))
			strand <- as.character(Gene.Table[ix.genename.table, 3])
			if (strand == "+") {
				tmp <- cumsum((end.exons - start.exons)) 
				tmp <- tmp[-length(tmp)]
			} else {
				tmp <- cumsum(rev(end.exons - start.exons)) 
				tmp <- tmp[-length(tmp)]
			}
			if (length(tmp) > 1) {
				junction2 <- sample(tmp,1)
			} else {
				junction2 <- 1
			}
			junction2.tot <- c(junction2.tot, junction2)
			sequence.fusions[i] <- paste(substr(sequence1[i], 1, junction1), substr(sequence2[i], junction2 + 1, nchar(sequence2[i])), sep = "")
			sequence.fusions.50bp[i] <- paste(substr(sequence1[i], (junction1 - 49), junction1), substr(sequence2[i], junction2 + 1, (junction2 + 50)), sep = "")
			id.fusions[i] <- paste(">", paste(gene1[i], gene2[i], sep = "----"), sep = "")
			myref <- c(myref, c(id.fusions[i], sequence.fusions[i]))
			myref.single <- c(id.fusions[i], formatfasta(sequence.fusions[i]))
			cat(myref.single, file = file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), paste("myref", formatted.count.fusions[i], ".fa", sep = "")), sep = "\n")
			
		}
		
		GeneFusions <- list()
		GeneFusions[[1]] <-  gene1
		GeneFusions[[2]] <-  gene2
		GeneFusions[[3]] <- junction1.tot
		GeneFusions[[4]] <- junction2.tot
		GeneFusions[[5]] <- sequence.fusions.50bp
		GeneFusions[[6]] <- mycoverage
		GeneFusions[[7]] <- trans1
		GeneFusions[[8]] <- trans2
		names(GeneFusions) <- c("gene1", "gene2", "junction1", "junction2", "junctionseq", "coverage", "trans1", "trans2")
		save(GeneFusions, file = file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "GeneFusions.RData"))
		
		
		system(paste(">", file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.1.fq")))
		system(paste(">", file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.2.fq")))
		for (i in 1: ngenefusion) {
			mynreads <- round(mycoverage[i]*nchar(sequence.fusions[i])/(2*readlength))
			if (verbose == 0) {
				system(paste("wgsim -d ", ins.size, " -r 0.0001 -R 0.001 -s ", sd.inssize, " -N ", mynreads, " -1 ", readlength, " -2 ", readlength," ", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), paste("myref", formatted.count.fusions[i], ".fa", sep = "")), " " ,file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq")," ", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"), " 2>> ", file.path(outputfolder, "wgsim.log"), " 1>> ", file.path(outputfolder, "wgsim.log"), sep = ""))
			} else {
				system(paste("wgsim -d ", ins.size, " -r 0.0001 -R 0.001 -s ", sd.inssize, " -N ", mynreads, " -1 ", readlength, " -2 ", readlength," ", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), paste("myref", formatted.count.fusions[i], ".fa", sep = "")), " " ,file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq")," ", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"), sep = ""))
				
			}
			system(paste("cat", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq"), ">>",  file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.1.fq")))
			system(paste("cat", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"), ">>",  file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.2.fq")))
		}
		
		if (flag.background == 1) {
			
			system(paste("cat ", file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.1.fq"), " ",  file.path(outputfolder, "background.reads.1.fq"), " > ", file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "total.reads.1.fq"), sep = ""))
			system(paste("cat ", file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "fusions.reads.2.fq"), " ",  file.path(outputfolder, "background.reads.2.fq"), " > ", file.path(outputfolder, "IE", "reads", paste("sim", formatted.count[jj], sep = "_"), "total.reads.2.fq"), sep = ""))
		}
		system(paste("rm", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.1.fq"))) 
		system(paste("rm", file.path(outputfolder, "IE", "data", paste("sim", formatted.count[jj], sep = "_"), "out.reads.2.fq"))) 
	}
	if (flag.background == 1) {
		system(paste("rm", file.path(outputfolder, "background.reads.1.fq")))
		system(paste("rm", file.path(outputfolder, "background.reads.2.fq")))
	}
	system(paste("rm", file.path(outputfolder, "wgsim.log")))
	cat(" done. \n")
}

