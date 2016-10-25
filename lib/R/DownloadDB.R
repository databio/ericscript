vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
ericscriptfolder <- split.vars [1]
user.refid <- split.vars [2]
tmpfolder <- split.vars [3]
ensversion <- as.numeric(split.vars [4])
####

load(file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData"))
ix.refid <- which(refid == user.refid)
if (length(ix.refid) == 0) {
  cat("\n[EricScript] Error: No data available for genome ", user.refid, ". Run ericscript.pl --printdb to view the available genomes.\n", sep = "")
  cat(0, file = file.path(tmpfolder, ".refid.flag"))
  quit( save = "no")
}
cat(1, file = file.path(tmpfolder, ".refid.flag"))
myrefid <- refid[ix.refid]
myrefid.path <- refid.path[ix.refid]

if (ensversion == 0) {
  ensversion <- version
}
## create xml queries

tmp0 <- unlist(strsplit(myrefid, "_"))
myrefid.xml <- paste(gsub(", ", "", toString(c(substr(tmp0[1], 1, 1), tmp0[2]))), "gene", "ensembl", sep = "_")

## genepos file
fileout <- file.path(tmpfolder, "genepos.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >", file = fileout, sep = "\n", append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", sep = ""), file = fileout, sep = "\n", append = T)
cat("<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
cat("<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
  cat("<Filter name = \"biotype\" value = \"processed_transcript,protein_coding\"/>", file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, sep = "\n", append = T)
} else {
  cat("<Filter name = \"biotype\" value = \"protein_coding\"/>", file = fileout, sep = "\n", append = T)
}
#### start attributes
cat("<Attribute name = \"ensembl_gene_id\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"chromosome_name\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"start_position\" />", file = fileout, sep = "\n", append = T)
#### end attributes
cat("</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

## geneinfo file
fileout <- file.path(tmpfolder, "geneinfo.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >", file = fileout, sep = "\n", append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", sep = ""), file = fileout, sep = "\n", append = T)
cat("<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
cat("<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
  cat("<Filter name = \"biotype\" value = \"processed_transcript,protein_coding\"/>", file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, sep = "\n", append = T)
} else {
  cat("<Filter name = \"biotype\" value = \"protein_coding\"/>", file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("<Attribute name = \"ensembl_gene_id\" />", file = fileout, sep = "\n", append = T)
if (ensversion > 0 & ensversion <= 75) {
  cat("<Attribute name = \"external_gene_id\" />", file = fileout, sep = "\n", append = T)  
} else {
  cat("<Attribute name = \"external_gene_name\" />", file = fileout, sep = "\n", append = T)
}
cat("<Attribute name = \"description\" />", file = fileout, sep = "\n", append = T)
#### end attributes
cat("</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

## exonstartend file
fileout <- file.path(tmpfolder, "exonstartend.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >", file = fileout, sep = "\n", append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", sep = ""), file = fileout, sep = "\n", append = T)
cat("<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
cat("<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
  cat("<Filter name = \"biotype\" value = \"processed_transcript,protein_coding\"/>", file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, sep = "\n", append = T)
} else {
  cat("<Filter name = \"biotype\" value = \"protein_coding\"/>", file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("<Attribute name = \"ensembl_gene_id\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"exon_chrom_start\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"exon_chrom_end\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"chromosome_name\" />", file = fileout, sep = "\n", append = T)
#### end attributes
cat("</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

## strand file
fileout <- file.path(tmpfolder, "strand.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >", file = fileout, sep = "\n", append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", sep = ""), file = fileout, sep = "\n", append = T)
cat("<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
cat("<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
  cat("<Filter name = \"biotype\" value = \"processed_transcript,protein_coding\"/>", file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, sep = "\n", append = T)
} else {
  cat("<Filter name = \"biotype\" value = \"protein_coding\"/>", file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("<Attribute name = \"ensembl_gene_id\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"strand\" />", file = fileout, sep = "\n", append = T)
#### end attributes
cat("</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)


## paralogs file (if it exists)
if (myrefid == "homo_sapiens") {
  fileout <- file.path(tmpfolder, "paralogs.xml")
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
  cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >", file = fileout, sep = "\n", append = T)
  cat("", file = fileout, sep = "\n", append = T)
  cat(paste("<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", sep = ""), file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
  if (myrefid == "homo_sapiens") {
    cat("<Filter name = \"biotype\" value = \"processed_transcript,protein_coding\"/>", file = fileout, sep = "\n", append = T)
    cat("<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, sep = "\n", append = T)
  } else {
    cat("<Filter name = \"biotype\" value = \"protein_coding\"/>", file = fileout, sep = "\n", append = T)
  } 
  #### start attributes
  cat("<Attribute name = \"ensembl_gene_id\" />", file = fileout, sep = "\n", append = T)
  cat("<Attribute name = \"hsapiens_paralog_ensembl_gene\" />", file = fileout, sep = "\n", append = T)
  #### end attributes
  cat("</Dataset>", file = fileout, sep = "\n", append = T)
  cat("</Query>", file = fileout, sep = "\n", append = T)
}

##  transcripts (eric the simulator)
fileout <- file.path(tmpfolder, "transcripts.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >", file = fileout, sep = "\n", append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", sep = ""), file = fileout, sep = "\n", append = T)
cat("<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
cat("<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
  cat("<Filter name = \"biotype\" value = \"processed_transcript,protein_coding\"/>", file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, sep = "\n", append = T)
} else {
  cat("<Filter name = \"biotype\" value = \"protein_coding\"/>", file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("<Attribute name = \"ensembl_gene_id\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"ensembl_transcript_id\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"exon_chrom_start\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"exon_chrom_end\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"chromosome_name\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"strand\" />", file = fileout, sep = "\n", append = T)
#### end attributes
cat("</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

fileout <- file.path(tmpfolder, "transcripts_cdna.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >", file = fileout, sep = "\n", append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", sep = ""), file = fileout, sep = "\n", append = T)
cat("<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
cat("<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
  cat("<Filter name = \"biotype\" value = \"processed_transcript,protein_coding\"/>", file = fileout, sep = "\n", append = T)
  cat("<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, sep = "\n", append = T)
} else {
  cat("<Filter name = \"biotype\" value = \"protein_coding\"/>", file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("<Attribute name = \"ensembl_transcript_id\" />", file = fileout, sep = "\n", append = T)
cat("<Attribute name = \"cdna\" />", file = fileout, sep = "\n", append = T)
#### end attributes
cat("</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

## download gene data
system(paste("perl", file.path(ericscriptfolder, "lib", "perl", "retrievefrombiomart.pl"), file.path(tmpfolder, "genepos.xml"), ensversion, "| sort -u - >", file.path(tmpfolder, "genepos.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl", "retrievefrombiomart.pl"), file.path(tmpfolder, "geneinfo.xml"), ensversion, "| sort -u - >", file.path(tmpfolder, "geneinfo.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl", "retrievefrombiomart.pl"), file.path(tmpfolder, "exonstartend.xml"), ensversion, "| sort -u - >", file.path(tmpfolder, "exonstartend.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl", "retrievefrombiomart.pl"), file.path(tmpfolder, "strand.xml"), ensversion, "| sort -u - >", file.path(tmpfolder, "strand.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl", "retrievefrombiomart.pl"), file.path(tmpfolder, "transcripts.xml"), ensversion, "| sort -u - >", file.path(tmpfolder, "transcripts.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl", "retrievefrombiomart.pl"), file.path(tmpfolder, "transcripts_cdna.xml"), ensversion, ">", file.path(tmpfolder, "transcripts.fa")))
if (myrefid == "homo_sapiens") {
  system(paste("perl", file.path(ericscriptfolder, "lib", "perl", "retrievefrombiomart.pl"), file.path(tmpfolder, "paralogs.xml"), ensversion, "| sort -u - >", file.path(tmpfolder, "paralogs.txt")))
  acc.chrs <- c(1:22, "X", "Y")
  cat (acc.chrs, file = file.path(tmpfolder, "chrlist"), sep = "\n")
}
## download seq data
download.file(file.path("ftp://ftp.ensembl.org/pub", paste("release-", ensversion, sep = ""), "fasta", myrefid, "dna", myrefid.path), destfile = file.path(tmpfolder, "seq.fa.gz"), quiet = T)



