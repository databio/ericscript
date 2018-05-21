vars.tmp <- commandArgs()
#cat("\nvars.tmp: ", vars.tmp, "\n", sep="")
vars <- vars.tmp[length(vars.tmp)]
#cat("\nvars: ", vars, "\n", sep="")
split.vars <- unlist(strsplit(vars, ","))
#cat("\nsplit.vars: ", split.vars, "\n", sep="")
ericscriptfolder <- split.vars [1]
#cat("\nsplit.vars [1]: ", ericscriptfolder, "\n", sep="")
user.refid <- split.vars [2]
#cat("\nsplit.vars [2]: ", user.refid, "\n", sep="")
tmpfolder <- split.vars [3]
#cat("\nsplit.vars [3]: ", tmpfolder, "\n", sep="")
ensversion <- as.numeric(split.vars [4])

###############################################################################

load(file.path(ericscriptfolder, "lib", "data", "_resources", "RefID.RData"))
ix.refid <- which(refid == user.refid)
#cat("\nThe provided refid is: ", ix.refid, "\n", sep="")
if (length(ix.refid) == 0) {
    cat("\n[EricScript] Error: No data available for genome ", user.refid, 
        ". Run ericscript.pl --printdb to view the available genomes.\n", 
        sep = "")
    cat(0, file = file.path(tmpfolder, ".refid.flag"))
    quit( save = "no")
}
cat(1, file = file.path(tmpfolder, ".refid.flag"))
myrefid <- refid[ix.refid]
myrefid.path <- refid.path[ix.refid]

if (ensversion == 0) {
    ensversion <- version
}

###############################################################################
####                                                                       ####
####                        create xml queries                             ####
####                                                                       ####
###############################################################################

tmp0   <- unlist(strsplit(myrefid, "_"))
tmpfil <- tempfile()
cat(myrefid, file = tmpfil, sep = "\n")
fields <- count.fields(tmpfil, sep = "_")
#cat("\nfields: ", fields, "\n")
if(fields == 2) {
    myrefid.xml <- paste(gsub(", ", "", 
                              toString(c(substr(tmp0[1], 1, 1), tmp0[2]))),
                         "gene", "ensembl", sep = "_")
} else {
    myrefid.xml <- paste(gsub(", ", "", 
                              toString(c(substr(tmp0[1], 1, 1), tmp0[3]))),
                         "gene", "ensembl", sep = "_")
}
unlink(tmpfil)
#cat("\n", myrefid.xml, "\n")

###############################################################################
## genepos file
###############################################################################

fileout <- file.path(tmpfolder, "genepos.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<!DOCTYPE Query>", file = fileout, sep = "\n", append = T)
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"",
    "0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n", 
    file = fileout, append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("\t<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >",
          sep = ""), file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"status\" value = \"KNOWN\"/>", 
    # file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"transcript_status\" value = \"KNOWN\"/>", 
    # file = fileout, sep = "\n", append = T)

if (myrefid == "homo_sapiens") {
    cat("\t\t<Filter name = \"biotype\" value = \"processed_transcript,",
        "protein_coding\"/>\n", file = fileout, append = T)
    cat("\t\t<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout,
        sep = "\n", append = T)
} else {
    cat("\t\t<Filter name = \"biotype\" value = \"protein_coding\"/>",
        file = fileout, sep = "\n", append = T)
}
#### start attributes
cat("\t\t<Attribute name = \"ensembl_gene_id\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"chromosome_name\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"start_position\" />", file = fileout,
    sep = "\n", append = T)
#### end attributes
cat("\t</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

###############################################################################
## geneinfo file
###############################################################################

fileout <- file.path(tmpfolder, "geneinfo.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<!DOCTYPE Query>", file = fileout, sep = "\n", append = T)
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = ",
    "\"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n",
    file = fileout, append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("\t<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >",
          sep = ""), file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout,
    # sep = "\n", append = T)
# cat("\t\t<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout,
    # sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
    cat("\t\t<Filter name = \"biotype\" value = \"",
        "processed_transcript,protein_coding\"/>\n", 
        file = fileout, append = T)
    cat("\t\t<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout, 
        sep = "\n", append = T)
} else {
    cat("\t\t<Filter name = \"biotype\" value = \"protein_coding\"/>", 
        file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("\t\t<Attribute name = \"ensembl_gene_id\" />", file = fileout, sep = "\n",
    append = T)
if (ensversion > 0 & ensversion <= 75) {
    cat("\t\t<Attribute name = \"external_gene_id\" />", 
    file = fileout, sep = "\n", append = T)  
} else {
    cat("\t\t<Attribute name = \"external_gene_name\" />", file = fileout,
        sep = "\n", append = T)
}
cat("\t\t<Attribute name = \"description\" />", file = fileout, sep = "\n",
    append = T)
#### end attributes
cat("\t</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

###############################################################################
## exonstartend file
###############################################################################

fileout <- file.path(tmpfolder, "exonstartend.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<!DOCTYPE Query>", file = fileout, sep = "\n", append = T)
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header =",
    "\"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n",
    file = fileout, append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("\t<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", 
          sep = ""), file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"status\" value = \"KNOWN\"/>", 
    # file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout,
    # sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
    cat("\t\t<Filter name = \"biotype\" value = \"processed_transcript,",
        "protein_coding\"/>\n", file = fileout, append = T)
    cat("\t\t<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout,
        sep = "\n", append = T)
} else {
    cat("\t\t<Filter name = \"biotype\" value = \"protein_coding\"/>",
        file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("\t\t<Attribute name = \"ensembl_gene_id\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"exon_chrom_start\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"exon_chrom_end\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"chromosome_name\" />", file = fileout,
    sep = "\n", append = T)
#### end attributes
cat("\t</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

###############################################################################
## strand file
###############################################################################

fileout <- file.path(tmpfolder, "strand.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<!DOCTYPE Query>", file = fileout, sep = "\n", append = T)
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = ",
    "\"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n",
    file = fileout, append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("\t<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >", 
          sep = ""), file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout,
    # sep = "\n", append = T)
# cat("\t\t<Filter name = \"transcript_status\" value = \"KNOWN\"/>", file = fileout,
    # sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
    cat("\t\t<Filter name = \"biotype\" value = \"processed_transcript,",
        "protein_coding\"/>\n", file = fileout, append = T)
    cat("\t\t<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout,
        sep = "\n", append = T)
} else {
    cat("\t\t<Filter name = \"biotype\" value = \"protein_coding\"/>",
        file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("\t\t<Attribute name = \"ensembl_gene_id\" />", 
    file = fileout, sep = "\n", append = T)
cat("\t\t<Attribute name = \"strand\" />", 
    file = fileout, sep = "\n", append = T)
#### end attributes
cat("\t</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

###############################################################################
## paralogs file (if it exists)
###############################################################################

if (myrefid == "homo_sapiens") {
    fileout <- file.path(tmpfolder, "paralogs.xml")
    cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
        file = fileout, sep = "\n")
    cat("<!DOCTYPE Query>", file = fileout, sep = "\n", append = T)
    cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" ",
        "header = \"0\" uniqueRows = \"0\" count = \"\" ",
        "datasetConfigVersion = \"0.6\" >\n", file = fileout, append = T)
    cat("", file = fileout, sep = "\n", append = T)
    cat(paste("\t<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >",
              sep = ""), file = fileout, sep = "\n", append = T)
    # cat("\t\t<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout,
        # sep = "\n", append = T)
    # cat("\t\t<Filter name = \"transcript_status\" value = \"KNOWN\"/>",
        # file = fileout, sep = "\n", append = T)
    if (myrefid == "homo_sapiens") {
        cat("\t\t<Filter name = \"biotype\" value = \"processed_transcript,",
            "protein_coding\"/>", file = fileout, sep = "\n", append = T)
        cat("\t\t<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout,
            sep = "\n", append = T)
    } else {
        cat("\t\t<Filter name = \"biotype\" value = \"protein_coding\"/>",
            file = fileout, sep = "\n", append = T)
    } 
    #### start attributes
    cat("\t\t<Attribute name = \"ensembl_gene_id\" />", file = fileout,
        sep = "\n", append = T)
    cat("\t\t<Attribute name = \"hsapiens_paralog_ensembl_gene\" />",
        file = fileout, sep = "\n", append = T)
    #### end attributes
    cat("\t</Dataset>", file = fileout, sep = "\n", append = T)
    cat("</Query>", file = fileout, sep = "\n", append = T)
}

###############################################################################
##  transcripts (eric the simulator)
###############################################################################

fileout <- file.path(tmpfolder, "transcripts.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<!DOCTYPE Query>", file = fileout, sep = "\n", append = T)
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = ",
    "\"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n",
    file = fileout, sep = "\n", append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("\t<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >",
          sep = ""), file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout,
    # sep = "\n", append = T)
# cat("\t\t<Filter name = \"transcript_status\" value = \"KNOWN\"/>",
    # file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
    cat("\t\t<Filter name = \"biotype\" value = \"processed_transcript,",
        "protein_coding\"/>\n", file = fileout, append = T)
    cat("\t\t<Filter name = \"with_hgnc\" excluded = \"0\"/>",
        file = fileout, sep = "\n", append = T)
} else {
    cat("\t\t<Filter name = \"biotype\" value = \"protein_coding\"/>", \
        file = fileout, sep = "\n", append = T)
} 
#### start attributes
cat("\t\t<Attribute name = \"ensembl_gene_id\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"ensembl_transcript_id\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"exon_chrom_start\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"exon_chrom_end\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"chromosome_name\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"strand\" />", file = fileout, sep = "\n", append = T)
#### end attributes
cat("\t</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

###############################################################################
####                        transcripts_cdna                               ####
###############################################################################

fileout <- file.path(tmpfolder, "transcripts_cdna.xml")
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = fileout, sep = "\n")
cat("<!DOCTYPE Query>", file = fileout, sep = "\n", append = T)
cat("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = ",
    "\"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n",
    file = fileout, append = T)
cat("", file = fileout, sep = "\n", append = T)
cat(paste("\t<Dataset name = \"", myrefid.xml,"\" interface = \"default\" >",
          sep = ""), file = fileout, sep = "\n", append = T)
# cat("\t\t<Filter name = \"status\" value = \"KNOWN\"/>", file = fileout,
    # sep = "\n", append = T)
# cat("\t\t<Filter name = \"transcript_status\" value = \"KNOWN\"/>",
    # file = fileout, sep = "\n", append = T)
if (myrefid == "homo_sapiens") {
    cat("\t\t<Filter name = \"biotype\" value = \"processed_transcript,",
        "protein_coding\"/>\n", file = fileout, append = T)
    cat("\t\t<Filter name = \"with_hgnc\" excluded = \"0\"/>", file = fileout,
        sep = "\n", append = T)
} else {
    cat("\t\t<Filter name = \"biotype\" value = \"protein_coding\"/>",
        file = fileout, sep = "\n", append = T)
}
#### start attributes
cat("\t\t<Attribute name = \"ensembl_transcript_id\" />", file = fileout,
    sep = "\n", append = T)
cat("\t\t<Attribute name = \"cdna\" />", file = fileout,
    sep = "\n", append = T)
#### end attributes
cat("\t</Dataset>", file = fileout, sep = "\n", append = T)
cat("</Query>", file = fileout, sep = "\n", append = T)

###############################################################################
####                          download gene data                           ####
###############################################################################

cat("\nDownloading files...\n")
system(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                               "retrievefrombiomart.pl"),
             file.path(tmpfolder, "genepos.xml"),
             ensversion, "| sort -u - >",
             file.path(tmpfolder, "genepos.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                               "retrievefrombiomart.pl"),
             file.path(tmpfolder, "geneinfo.xml"),
             ensversion, "| sort -u - >",
             file.path(tmpfolder, "geneinfo.txt")))
cat("\nDownloading exonstartend.txt file\n")
cat(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                               "retrievefrombiomart.pl"),
             file.path(tmpfolder, "exonstartend.xml"),
             ensversion, "| sort -u - >",
             file.path(tmpfolder, "exonstartend.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                               "retrievefrombiomart.pl"),
             file.path(tmpfolder, "exonstartend.xml"),
             ensversion, "| sort -u - >",
             file.path(tmpfolder, "exonstartend.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                               "retrievefrombiomart.pl"),
             file.path(tmpfolder, "strand.xml"),
             ensversion, "| sort -u - >",
             file.path(tmpfolder, "strand.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                               "retrievefrombiomart.pl"),
             file.path(tmpfolder, "transcripts.xml"),
             ensversion, "| sort -u - >",
             file.path(tmpfolder, "transcripts.txt")))
system(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                               "retrievefrombiomart.pl"),
             file.path(tmpfolder, "transcripts_cdna.xml"),
             ensversion, ">",
             file.path(tmpfolder, "transcripts.fa")))
if (myrefid == "homo_sapiens") {
    system(paste("perl", file.path(ericscriptfolder, "lib", "perl",
                                   "retrievefrombiomart.pl"),
                 file.path(tmpfolder, "paralogs.xml"),
                 ensversion, "| sort -u - >",
                 file.path(tmpfolder, "paralogs.txt")))
    acc.chrs <- c(1:22, "X", "Y")
    cat (acc.chrs, file = file.path(tmpfolder, "chrlist"), sep = "\n")
}
## download seq data
#download.file(file.path("ftp://ftp.ensembl.org/pub",
#                        paste("release-", ensversion, sep = ""),
#                        "fasta", myrefid, "dna", myrefid.path),
#              destfile = file.path(tmpfolder, "seq.fa.gz"), quiet = T)

download.file(file.path("ftp://ftp.ensembl.org/pub", 
                        paste("release-", ensversion, sep = ""),
                        "fasta", myrefid.path),
              destfile = file.path(tmpfolder, "seq.fa.gz"), quiet = T)