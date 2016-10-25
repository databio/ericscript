### retrieverefid ver 3: changed ensembl ftp filesystem ver 2
### added user selectable ensembl version
vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
ericscriptfolder <- split.vars [1]
flagprint <- split.vars [2]
dbfolder <-  split.vars [3]
ensversion <- split.vars [4]

mynetrc <- "machine ftp.ensembl.org login anonymous password -"
if (file.exists("~/.netrc")) {
  netrc.data <- scan("~/.netrc", what = "", sep = "\n", quiet = T)
  if (any(netrc.data %in% mynetrc) == F) {
    cat(mynetrc, file = "~/.netrc", append = T, sep = "\n")
  }  
} else {
  cat(mynetrc, file = "~/.netrc", sep = "\n")
}
system(paste("sh", file.path(ericscriptfolder, "lib", "bash", "Ftp2Ensembl.sh"), ericscriptfolder, ensversion))
xx.tmp <- readLines(file.path(ericscriptfolder, "lib", "data", "_resources", ".ftplist1"))
xx.tmp1 <- strsplit(xx.tmp, " ")
xx <- rep("", length(xx.tmp))
for (i in 1: length(xx.tmp1)) {
  if (length(xx.tmp1[[i]]) > 0) {
    xx[i] <- xx.tmp1[[i]][length(xx.tmp1[[i]])]
  } else {
    xx[i] <- "" 
  }
}

mybreaks <- which(xx == "")
ensrefid.tmp <- xx[1: mybreaks[1]]
ensrefid <- c()
ensrefid.path <- c()
for ( i in 1: length(ensrefid.tmp)) {
  #  ix.start <- which(xx == paste("./", ensrefid.tmp[i], "/dna:", sep ="")) + 1
  ix.start <- which(xx == paste(ensrefid.tmp[i], "/dna:", sep ="")) + 1
  if (length(ix.start) != 0) {
    ix.end <- grep("./", xx[ix.start: length(xx)])[1] + ix.start - 1
    ensrefid <- c(ensrefid, ensrefid.tmp[i])
    ensrefid.path <- c(ensrefid.path, grep("dna.toplevel", xx[ix.start: ix.end], value = T))
  }
}
ensversion0 <- ensversion
if (ensversion == 0) {
  xx.tmp <- readLines(file.path(ericscriptfolder, "lib", "data", "_resources", ".ftplist0"))
  xx.tmp1 <- xx.tmp[grep("[0-9] release-", xx.tmp, perl = T)]
  xx.tmp2 <- strsplit(xx.tmp1, " release-")
  xx.tmp3 <- rep(NA, length(xx.tmp2))
  for (i in 1: length(xx.tmp2)) {
    xx.tmp3[i] <- as.numeric(unlist(strsplit(xx.tmp2[[i]][2], " "))[1])
  }
  ensversion <- max(xx.tmp3)
}
if (flagprint != 0) {
if (ensversion0 != 0) {
cat("Selected Ensembl version:", ensversion, "\n")
} else {
cat("Current Ensembl version:", ensversion, "\n")
}
if (file.exists(file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData")) & any(file.exists(file.path(dbfolder, "data", ensrefid)))) {
  load(file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData"))
  cat("Installed Ensembl version:", version, "\n")
} else {
  cat("Installed Ensembl version:", "No database installed", "\n")
}
  cat("Available reference IDs:\n", paste("\t", ensrefid, "\n"))
}
if (file.exists(file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData"))) {
  load(file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData"))
  if (ensversion > version | ensversion0 != 0) {
    refid <- ensrefid
    refid.path <- ensrefid.path
    version <- ensversion
    save(refid, refid.path, version, file = file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData"))
    cat(version, file = file.path(ericscriptfolder, "lib", "data", "_resources", "Ensembl.version"))
    cat(1, file = file.path(ericscriptfolder, "lib", "data", "_resources", ".flag.updatedb"))
  } else {
    cat(0, file = file.path(ericscriptfolder, "lib", "data", "_resources", ".flag.updatedb"))
  }
} else {
  refid <- ensrefid
  refid.path <- ensrefid.path
  version <- ensversion
  save(refid, refid.path, version, file = file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData"))
  cat(version, file = file.path(ericscriptfolder, "lib", "data", "_resources", "Ensembl.version"))
  cat(1, file = file.path(ericscriptfolder, "lib", "data", "_resources", ".flag.updatedb"))  
}



