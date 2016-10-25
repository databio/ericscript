## MakeResults.R v0.2
## different read count-based method for gene expression level estimation
## added machine-learning based algorithm as summarization score
vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]

Results <- "No Chimeric Transcript found!"
write.table(Results, file = file.path(outputfolder,paste(samplename,".results.total.tsv", sep = "")), sep = "\t", row.names = F, col.names = F, quote = F)
	















