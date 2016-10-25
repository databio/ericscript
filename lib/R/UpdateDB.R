vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
ericscriptfolder <- split.vars [1]
dbfolder <- split.vars[2]

mydblist.tmp <- list.files(file.path(dbfolder, "data")) 
if (length(mydblist.tmp) > 1) {
  mydblist <- mydblist.tmp[-which(mydblist.tmp == "_resources")]
  flag <- scan(file.path(ericscriptfolder, "lib", "data", "_resources", ".flag.updatedb"), what = "numeric", quiet = T)
  if (flag == 0) {
    cat("[EricScript] Nothing to update. Exit.\n", sep = "")    
  } else
  {  
    cat("[EricScript] Found a new release of Ensembl Gene. Updating database for ", toString(mydblist),".\n", sep = "")    
    for (i in 1: length(mydblist)) {
      system(paste("sh", file.path(ericscriptfolder, "lib", "bash", "BuildSeq.sh"), ericscriptfolder, mydblist[i]))
    }
  }
} else {
  cat("[EricScript] No database was found in ", file.path(dbfolder, "data"), ". Please run ericscript.pl --downdb to download your databases.\n", sep = "")
}


