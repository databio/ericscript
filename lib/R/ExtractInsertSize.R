vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
outputfolder <- split.vars[1]
bwa_aln <- as.numeric(as.character(split.vars[2]))
if (bwa_aln == 0) {
  xx <- readLines(file.path(outputfolder, "out", ".ericscript.log"))
  ix.isize <- grep("mean and std.dev:", xx)
  if (length(ix.isize) > 0) {
    isize.tmp <- strsplit(tail(xx[ix.isize], n = 1), ": ")[[1]][2]
    isize.tmp1 <- unlist(strsplit(isize.tmp, ",", fixed = T))
    isize.mean <- as.numeric(substr(isize.tmp1[1], 2, nchar(isize.tmp1[1])))
    isize.sd <- as.numeric(substr(isize.tmp1[2], 1, nchar(isize.tmp1[2])-1))
  } else {
    isize.mean <- 200
    isize.sd <- 40
  }
  save(isize.mean, isize.sd, file = file.path(outputfolder, "out", "isize.RData"))
} else {
  xx <- readLines(file.path(outputfolder, "out", ".ericscript.log"))
  ix.isize <- grep("inferred external isize from", xx)
  if (length(ix.isize) > 0) {
    isize.tmp <- strsplit(tail(xx[ix.isize], n = 1), ": ")[[1]][2]
    isize.tmp1 <- unlist(strsplit(isize.tmp, "+/-", fixed = T))
    isize.mean <- as.numeric(isize.tmp1[1])
    isize.sd <- as.numeric(isize.tmp1[2])
  } else {
    isize.mean <- 200
    isize.sd <- 40
  }
  save(isize.mean, isize.sd, file = file.path(outputfolder, "out", "isize.RData"))
  
}
