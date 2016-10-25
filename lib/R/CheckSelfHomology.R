vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
samplename <- split.vars [1]
outputfolder <- split.vars[2]
x <- read.delim(file.path(outputfolder,"out",paste(samplename, ".checkselfhomology.blat", sep = "")), sep = "\t", header = F)
query.fa <- readLines(file.path(outputfolder,"out",paste(samplename, ".checkselfhomology.fa", sep = "")))
id.fa <- substr(query.fa[seq(1, length(query.fa), by = 2)], 2, nchar(query.fa[seq(1, length(query.fa), by = 2)]))
id.query.nolabel <- as.character(x[,1])
unique.ids.nolabel <- unique(id.query.nolabel)
ix.fa <- which(id.fa %in% unique.ids.nolabel)
seq.fa <- query.fa[seq(2, length(query.fa), by = 2)][ix.fa]
id.target.f <-  as.character(x[,2])
id.match.f <- as.numeric(as.character(x[,4]))
start.match <- as.numeric(as.character(x[,7]))
end.match <-  as.numeric(as.character(x[,8]))
diff.match <- end.match - start.match + 1
ix.junct.nohomology <- rep(0, length(unique.ids.nolabel))
ix.junct.homology <- rep(0, length(unique.ids.nolabel))
flag.dup.a <- rep(0, length(unique.ids.nolabel))
flag.dup.b <- rep(0, length(unique.ids.nolabel))
homology.list <- vector("list", length(unique.ids.nolabel))
for (i in 1:length(unique.ids.nolabel)) {
	ix.id <- which(id.query.nolabel == unique.ids.nolabel[i])
	query <- unique.ids.nolabel[i]
	target <- id.target.f[ix.id]
	match <- as.numeric(id.match.f[ix.id])
	query.tmp <- unlist(strsplit(query, "----", fixed = T))
	query_a <- unlist(strsplit(query.tmp[1], "_"))[1]
	query_b <- unlist(strsplit(query.tmp[2], "_"))[1]
	diff.match.tmp <- diff.match[ix.id]
	start.match.tmp <- start.match[ix.id]
	end.match.tmp <- end.match[ix.id]
	width <- 100 - length(grep("N", unlist(strsplit(seq.fa[i], ""))))
	ix.c <- seq.int(1,length(target))
	ix.a <- which(target %in% query_a)
	ix.b <- which(target %in% query_b)
	ix.ab <- c(ix.a, ix.b)
	if (((length(ix.a) > 0) & (length(ix.b) > 0)) | (length(ix.a) > 1 & length(ix.b) == 1) | (length(ix.a) == 1 & length(ix.b) > 1)) {
		if ((length(ix.a) > 1 & length(ix.b) == 1)) {
			myflag <- length(which(start.match.tmp[ix.a] %in% c((start.match.tmp[ix.b]-3):(start.match.tmp[ix.b]+3)))) + length(which(end.match.tmp[ix.a] %in% c((end.match.tmp[ix.b]-3):(end.match.tmp[ix.b]+3))))
		} else if ((length(ix.a) == 1 & length(ix.b) > 1)) {
			myflag <- length(which(start.match.tmp[ix.b] %in% c((start.match.tmp[ix.a]-3):(start.match.tmp[ix.a]+3)))) + length(which(end.match.tmp[ix.b] %in% c((end.match.tmp[ix.a]-3):(end.match.tmp[ix.a]+3))))
		} else {
			myflag <- 0
		}
		if (myflag == 0) {
			if (max(diff.match.tmp[ix.ab]) < round(0.8*width) & ((length(ix.ab) > 2) & any(diff.match.tmp[ix.ab] < 30) | (length(ix.ab) == 2))) {
				if (length(ix.ab) != 0) {
					ix.c <- ix.c[-ix.ab]
				}
				
				if(length(ix.c) != 0) {
					unique.ids.homology <- unique(target[ix.c])
					homology.list[[i]] <- vector("list", length(unique.ids.homology))
					for (j in 1:length(unique.ids.homology)) {
						ix.id.homology <- which(target[ix.c] == unique.ids.homology[j])
						max.match <- max(match[ix.c][ix.id.homology])
						homology.list[[i]][[j]] <- cbind(unique.ids.homology[j], max.match)
					}
					ix.junct.homology[i] <- 1
				}
				if(length(ix.c) == 0) {
					ix.junct.nohomology[i] <- 1
					
				}
				if (length(ix.a) > 1) {
					flag.dup.a[i] <- 1
				}
				if (length(ix.b) > 1) {
					flag.dup.b[i] <- 1
				}
				
			}
			
		}
	}
}
ix.junct <- sort(c(which(ix.junct.nohomology == 1), which(ix.junct.homology == 1)))
if (length(ix.junct) == 0) {
	myflag <- 0
	cat(myflag, file = file.path(outputfolder, "out", ".ericscript.flag"))
	stop("No putative gene fusions pass the self-homology filter. Exit!")		
}

info.homology <- rep("", length(ix.junct))
for (i in 1: length(ix.junct)) {
	list.tmp <- homology.list[[ix.junct[i]]]
	if (is.null(list.tmp) == F) {
		info.homology.tmp <- c()
		n.homo <- length(list.tmp)
		if (n.homo > 30) {
			n.homo <- 30
			info.homology.tmp <- "More than 30 homologies found: "
		}
		for (j in 1: n.homo) {
			info.homology.tmp <- paste(info.homology.tmp, paste(list.tmp[[j]][1,1]," (",list.tmp[[j]][1,2],"%)", sep = "" ), sep = "")
			if (n.homo > 1 & j < n.homo) {
				info.homology.tmp <- paste(info.homology.tmp, ", ", sep = "")
			} else 
			{
				info.homology[i] <- info.homology.tmp
			}
		}
		
	}
}

info.id.and.homology <- cbind(unique.ids.nolabel[ix.junct], info.homology, flag.dup.a[ix.junct], flag.dup.b[ix.junct])
save(info.id.and.homology, file = file.path(outputfolder,"out",paste(samplename,".ids_homology.RData", sep  = "")))
