#!/bin/bash

ericscriptfolder=$1
ensversion=$2
if [ $ensversion -eq 0 ]; then
   rsync --list-only rsync://ftp.ensembl.org/ensembl/pub/ . > .ftplist0
   mv .ftplist0 $ericscriptfolder/lib/data/_resources/
   rsync --list-only -av rsync://ftp.ensembl.org/ensembl/pub/current_fasta/ . > .ftplist1
   mv .ftplist1 $ericscriptfolder/lib/data/_resources/
else 
   fasta_path="release-"$ensversion"/fasta/"
   rsync --list-only rsync://ftp.ensembl.org/ensembl/pub/ . > .ftplist0
   mv .ftplist0 $ericscriptfolder/lib/data/_resources/
   rsync --list-only -av rsync://ftp.ensembl.org/ensembl/pub/${fasta_path}/ . > .ftplist1
   mv .ftplist1 $ericscriptfolder/lib/data/_resources/
fi