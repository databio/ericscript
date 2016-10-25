vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
ericscriptfolder <- split.vars[1]
refid <- split.vars[2]
dbfolder <- as.character(split.vars[3])

flag.dbexists <- 1

mydbdata.homo <- c("EnsemblGene.Reference.fa", "EnsemblGene.Sequences.RData", "EnsemblGene.GenePosition.RData", "EnsemblGene.Structures.RData", "EnsemblGene.GeneInfo.RData", "EnsemblGene.Paralogs.RData", "EnsemblGene.GeneNames.RData")
mydbdata <- c("EnsemblGene.Reference.fa", "EnsemblGene.Sequences.RData", "EnsemblGene.GenePosition.RData", "EnsemblGene.Structures.RData", "EnsemblGene.GeneInfo.RData","EnsemblGene.GeneNames.RData")

xx <- file.exists(file.path(dbfolder, "data", refid))
  if (xx) {
    xx.files <- list.files(file.path(dbfolder, "data", refid))
    if (refid == "homo_sapiens") {
      xx1 <- all( mydbdata.homo %in% xx.files )
    } else {
      xx1 <- all( mydbdata %in% xx.files )
    }
    if (!xx1) {
      flag.dbexists <- 0
      cat("[EricScript] Some required db files were not found for", refid, "genome. Please run \"ericscript.pl --downdb --refid", refid, "\" to solve this.\n")
    }
  } else {
    flag.dbexists <- 0
    cat("[EricScript] DB data for", refid, "genome do not exist. Set correct -db option or run \" ericscript.pl --downdb --refid", refid, "\" to solve this.\n")
  }

## check bwa version
yy <- file.exists(file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version"))
if (yy) {
  prev.version.bwa <- scan(file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version"), what = "", quiet = T, sep = "\n")
  system(paste("bwa", "2>&1", "|", "grep ersion", ">", file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp")))
  curr.version.bwa <- scan(file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp"), what = "", quiet = T, sep = "\n")
  if (curr.version.bwa != prev.version.bwa) {
    cat("[EricScript] Updating BWA indexes for", refid,"... ")    
    system(paste("bwa index", file.path(file.path(dbfolder, "data", refid, "EnsemblGene.Reference.fa")), "1>>", file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp"), "2>>", file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp")))    
    cat("done.\n")
    cat(curr.version.bwa, file = file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version"))
  }
} else {
  system(paste("bwa", "2>&1", "|", "grep ersion", ">", file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp")))
  curr.version.bwa <- scan(file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp"), what = "", quiet = T, sep = "\n")
  version.a <- gsub("Version: ", "", strsplit(curr.version.bwa, ".", fixed = T)[[1]][1])
  version.b <- strsplit(curr.version.bwa, ".", fixed = T)[[1]][2]
  version.c <- gsub("[a-z]", "", strsplit(strsplit(curr.version.bwa, ".", fixed = T)[[1]][3], "-")[[1]][1])
  version.tot <- as.numeric(paste(version.a, version.b, version.c, sep = ""))
  if (version.tot >= 74) {
    cat("[EricScript] Updating BWA indexes for", refid, "...")    
    system(paste("bwa index", file.path(file.path(dbfolder, "data", refid, "EnsemblGene.Reference.fa")), "1>>", file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp"), "2>>", file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version.tmp")))    
    cat("done.\n")
    system(paste("bwa", "2>&1", "|", "grep ersion", ">", file.path(ericscriptfolder, "lib", "data", "_resources", ".bwa.version")))
   } else {
    flag.dbexists <- 0
    cat("[EricScript] BWA version >= 0.7.4 is required. Exit.\n") 
  }
}

mydbdata.bwa <- c("EnsemblGene.Reference.fa.bwt", "EnsemblGene.Reference.fa.pac", "EnsemblGene.Reference.fa.ann", "EnsemblGene.Reference.fa.amb", "EnsemblGene.Reference.fa.sa")

if (xx) {
  xx.files.bwa <- list.files(file.path(dbfolder, "data", refid))
  xx1 <- all( mydbdata.bwa %in% xx.files.bwa )
  if (!xx1) {
    flag.dbexists <- 0
    cat("[EricScript] Some required files (bwa indexes) were not found for", refid, "genome. Please run \"ericscript.pl --downdb --refid", refid, "\" to solve this.\n")
  }  
}

cat(flag.dbexists, file = file.path(ericscriptfolder, "lib", "data", "_resources", ".flag.dbexists"))


