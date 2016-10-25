### calculate statistics

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
resultsfolder <- split.vars[1]
outputfolder <- split.vars[2]
datafolder <- split.vars[3]
algoname <- split.vars[4]
dataset <- split.vars[5]
readlength <- as.numeric(split.vars[6])
normroc <- as.numeric(split.vars[7])
ericscriptfolder <- as.character(split.vars[8])

source(file.path(ericscriptfolder, "lib", "R", "ImportResults.R"))

trapezint <- function (x, y, a, b)  {
## function of the ROC package (http://bioconductor.org) 
    if (length(x) != length(y)) 
	stop("length x must equal length y")
    y <- y[x >= a & x <= b]
    x <- x[x >= a & x <= b]
    if (length(unique(x)) < 2) 
	return(NA)
    ya <- approx(x, y, a, ties = max, rule = 2)$y
    yb <- approx(x, y, b, ties = max, rule = 2)$y
    x <- c(a, x, b)
    y <- c(ya, y, yb)
    h <- diff(x)
    lx <- length(x)
    0.5 * sum(h * (y[-1] + y[-lx]))
}

algonamelist <- c("ericscript", "chimerascan", "defuse", "fusionmap", "shortfuse")
algoname <- tolower(algoname)

if (any(algonamelist %in% algoname) == F) {
	algoid <- "unknown"
} else {
	algoid <- algoname
}

xx <- list.files(resultsfolder, pattern = "sim_")
nsims <- length(xx)
cat("[EricScript calcstats] Found ", nsims, " synthetic data analysis for algorithm ", algoname,". \n", sep = "")

if (nsims > 0) {
	tpr <- rep(NA, nsims)
	fpr <- rep(NA, nsims)
	tpr.5 <- rep(NA, nsims)
	fpr.5 <- rep(NA, nsims)
	tpr.seq <- rep(NA, nsims)
	refpath <- file.path(datafolder, dataset, "data")
	rocs.tpr <- rep(0, 1000)
	rocs.fpr <- rep(0, 1000)
	nosims <- 0
	
	refpath <- file.path(datafolder, dataset, "data")
	
	for (i in 1: nsims) {
		if (i < 10) {
			dataresults <- get(paste("Import", algoid, sep = "_"))(file.path(resultsfolder, paste("sim_", "0000", i, sep = "")))
			if (is.list(dataresults)) {
				load(file.path(refpath, paste("sim_", "0000", i, sep = ""), "GeneFusions.RData"))
				cat ("[EricScript calcstats] Analysing ", paste("sim_", "0000", i, sep = "")," ... ")
			} else if (dataresults == 0) {
				cat("[EricScript calcstats] No results file found for ", paste("sim_", "0000", i, sep = ""), ".\n", sep = "")
				nosims <- nosims + 1
				next
			} else if (dataresults > 1) {
				cat("[EricScript calcstats] Error: ", dataresults, " results file found for ", paste("sim_", "0000", i, sep = ""), ". Only 1 file of results is required. \n", sep = "")
				nosims <- nosims + 1
				next
			}
		} else if (i >= 10 & i < 100) {
			dataresults <- get(paste("Import", algoid, sep = "_"))(file.path(resultsfolder, paste("sim_", "000", i, sep = "")))
			if (is.list(dataresults)) {
				load(file.path(refpath, paste("sim_", "000", i, sep = ""), "GeneFusions.RData"))
				cat ("[EricScript calcstats] Analysing ", paste("sim_", "000", i, sep = "")," ... ")
			} else if (dataresults == 0) {
				cat("[EricScript calcstats] No results file found for ", paste("sim_", "000", i, sep = ""), ".\n", sep = "")
				nosims <- nosims + 1
				next
			} else if (dataresults > 1) {
				cat("[EricScript calcstats] Error: ", dataresults, " results file found for ", paste("sim_", "000", i, sep = ""), ". Only 1 file of results is required. \n", sep = "")
				nosims <- nosims + 1
				next
			}
		} else if (i >= 100 & i < 1000) {
			dataresults <- get(paste("Import", algoid, sep = "_"))(file.path(resultsfolder, paste("sim_", "00", i, sep = "")))
			if (is.list(dataresults)) {
				load(file.path(refpath, paste("sim_", "00", i, sep = ""), "GeneFusions.RData"))
				cat ("[EricScript calcstats] Analysing ", paste("sim_", "00", i, sep = "")," ... ")
			} else if (dataresults == 0) {
				cat("[EricScript calcstats] No results file found for ", paste("sim_", "00", i, sep = ""), ".\n", sep = "")
				nosims <- nosims + 1
				next
			} else if (dataresults > 1) {
				cat("[EricScript calcstats] Error: ", dataresults, " results file found for ", paste("sim_", "00", i, sep = ""), ". Only 1 file of results is required. \n", sep = "")
				nosims <- nosims + 1
				next
			}
		} else if (i >= 1000 & i < 10000) {
			dataresults <- get(paste("Import", algoid, sep = "_"))(file.path(resultsfolder, paste("sim_", "0", i, sep = "")))
			if (is.list(dataresults)) {
				load(file.path(refpath, paste("sim_", "0", i, sep = ""), "GeneFusions.RData"))
				cat ("[EricScript calcstats] Analysing ", paste("sim_", "0000", i, sep = "")," ... ")
			} else if (dataresults == 0) {
				cat("[EricScript calcstats] No results file found for ", paste("sim_", "0", i, sep = ""), ".\n", sep = "")
				nosims <- nosims + 1
				next
			} else if (dataresults > 1) {
				cat("[EricScript calcstats] Error: ", dataresults, " results file found for ", paste("sim_", "0", i, sep = ""), ". Only 1 file of results is required. \n", sep = "")
				nosims <- nosims + 1
				next
			}
		} else if (i >= 10000 & i < 100000) {
			dataresults <- get(paste("Import", algoid, sep = "_"))(file.path(resultsfolder, paste("sim_", i, sep = "")))
			if (is.list(dataresults)) {
				load(file.path(refpath, paste("sim_", i, sep = ""), "GeneFusions.RData"))
				cat ("[EricScript calcstats] Analysing ", paste("sim_",  i, sep = "")," ... ")
			} else if (dataresults == 0) {
				cat("[EricScript calcstats] No results file found for ", paste("sim_", i, sep = ""), ".\n", sep = "")
				nosims <- nosims + 1
				next
			} else if (dataresults > 1) {
				cat("[EricScript calcstats] Error: ", dataresults, " results file found for ", paste("sim_", i, sep = ""), ". Only 1 file of results is required. \n", sep = "")
				nosims <- nosims + 1
				next
			}
		}
		
		id1.simul <- GeneFusions[[1]]
		id2.simul <- GeneFusions[[2]]
		seq.simul <- GeneFusions[[5]]
		ngenefusions <- length(id1.simul)
		if (!exists("cov.tpr")) {
			cov.tpr <- rep(0, ngenefusions)
		}
		gene1ens <- as.character(dataresults$gene5)
		gene2ens <- as.character(dataresults$gene3)
		nreads <- as.numeric(as.character(dataresults$nreads))
		score <- as.numeric(as.character(dataresults$score))
		if (normroc > 1) {
			score <- score/normroc
		}
		seq <- as.character(dataresults$seq)
		
		tpr[i] <- length(which(gene1ens %in% id1.simul & gene2ens %in% id2.simul))/ngenefusions
		ix.tpr <- which(id1.simul %in% gene1ens & id2.simul %in%gene2ens)
		ix.fpr <- which((gene1ens %in% id1.simul & gene2ens %in% id2.simul) == F)
		cov.tpr[ix.tpr] <- cov.tpr[ix.tpr] + 1
		
		tpr.5[i] <- length(which(gene1ens %in% id1.simul & gene2ens %in% id2.simul & nreads > 5))/ngenefusions
		fpr.5[i] <- length(which((gene1ens %in% id1.simul & gene2ens %in% id2.simul) == F & nreads > 5))/length(gene1ens)
		
		fpr[i] <- (length(gene1ens) - tpr[i]*ngenefusions)/length(gene1ens)
		
		
		rocs.tpr <- colSums(rbind(rocs.tpr, tabulate(score[which(gene1ens %in% id1.simul & gene2ens %in% id2.simul)]*1000, nbins = 1000)))
		rocs.fpr <- colSums(rbind(rocs.fpr, tabulate(score[which((gene1ens %in% id1.simul & gene2ens %in% id2.simul) == F)]*1000, nbins = 1000)))
		
		ix.correctseq <- c()
		for (ii in 1: length(ix.tpr)) {
			ix.tmp <- which(gene1ens == id1.simul[ix.tpr[ii]] & gene2ens == id2.simul[ix.tpr[ii]])
			if (length(ix.tmp) > 1) { 
				ix.tmp <- ix.tmp[1]
			}
			if (length(agrep(toupper(seq[ix.tmp]), seq.simul[ix.tpr], max = 5)) > 0 & is.na(seq[ix.tmp]) == F) {
				ix.correctseq <- c(ix.correctseq, ix.tpr[ii])
			}
		}
		tpr.seq[i] <- length(ix.correctseq)/length(ix.tpr)
		
		
		cat ("done.\n")
		
	}
	
	
	roc.total <- rocs.tpr + rocs.fpr
	ntot <- sum(roc.total)
	ntot.tpr <- sum(rocs.tpr)
	ntot.fpr <- sum(rocs.fpr)
	sens <- rep(0, 1000)
	spec <- rep(0, 1000)
	for (i in 1: 1000) {
		if (i == 1) {
			sens[i] <- sum(rocs.tpr)/ntot.tpr
			spec[i] <- 1 - sum(rocs.fpr)/ntot.fpr 
		} else {
			sens[i] <- sum(rocs.tpr[-c(1:i)])/ntot.tpr
			spec[i] <- 1 - sum(rocs.fpr[-c(1:i)])/ntot.fpr
		}
	}
	
	stats <- list()
	stats$algorithm <- algoname
	stats$dataset <- dataset
	stats$readlength <- readlength
	stats$totalsims <- nsims
	stats$nsims <- nsims - nosims
	stats$meantpr <- mean(tpr, na.rm = T)
	stats$meanfpr <- mean(fpr, na.rm = T)
	stats$meantpsr <- mean(tpr.seq, na.rm = T)
	stats$auc <- trapezint(sens, 1 - spec, 0, 1)
	stats$meantpr5 <- mean(tpr.5, na.rm = T)
	stats$meanfpr5 <- mean(fpr.5, na.rm = T)
	stats$tpr <- tpr
	stats$fpr <- fpr
	stats$tpsr <- tpr.seq
	stats$tpr5 <- tpr.5
	stats$fpr5 <- fpr.5
	stats$scoring_sensitivity <- sens
	stats$scoring_specificity <- spec
	stats$covtpr <- cov.tpr/(nsims-nosims)
	
	save(stats, file = file.path(outputfolder, paste(algoname, dataset, readlength, "stats","RData", sep = ".")))
} else {
	cat ("[EricScript calcstats] Error: no directories containing results on synthetic data have been found in ", resultsfolder, ". Exit.\n", sep = "")
	
}

