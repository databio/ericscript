#!/bin/bash
ericscriptfolder=$1
refid=$2
dbfolder=$3
ensversion=$4
myrandomn=$RANDOM
tmpfolder=$dbfolder/".tmp_"$myrandomn
mkdir $tmpfolder
printf "[EricScript] Downloading $refid data. This process may take from few minutes to few hours depending on the selected genome ..."
R --slave --args $ericscriptfolder,$refid,$tmpfolder,$ensversion < $ericscriptfolder/lib/R/DownloadDB.R
flagrefid=`cat $tmpfolder/.refid.flag` 
if [ $flagrefid -eq 1 ]
then
bedtools sort -i $tmpfolder/exonstartend.txt | bedtools merge -c 4 -o collapse -i - | cut -d ',' -f1 - | awk '{print $4"\t"($2-1)"\t"$3"\t"$1}' - > $tmpfolder/exonstartend.mrg.txt
seqtk subseq $tmpfolder/seq.fa.gz $tmpfolder/exonstartend.mrg.txt > $tmpfolder/subseq.fa
printf "done.\n"
printf "[EricScript] Creating database for $refid ..."
R --slave --args $ericscriptfolder,$refid,$dbfolder,$tmpfolder < $ericscriptfolder/lib/R/BuildExonUnionModel.R
R --slave --args $ericscriptfolder,$refid,$dbfolder,$tmpfolder < $ericscriptfolder/lib/R/ConvertTxt2R.R
R --slave --args $refid,$dbfolder,$tmpfolder < $ericscriptfolder/lib/R/CreateDataEricTheSimulator.R 
if [ $refid == "homo_sapiens" ]
then
seqtk subseq -l 50 $tmpfolder/seq.fa.gz $tmpfolder/chrlist > $dbfolder/data/$refid/allseq.fa
else
gunzip -c -d $tmpfolder/seq.fa.gz > $dbfolder/data/$refid/allseq.fa
fi
printf "done.\n"
printf "[EricScript] Building reference indexes with BWA for transcriptome and genome ..."
bwa index $dbfolder/data/$refid/allseq.fa 1>> $tmpfolder/.tmp.log 2>> $tmpfolder/.tmp.log
bwa index $dbfolder/data/$refid/EnsemblGene.Reference.fa 1>> $tmpfolder/.tmp.log 2>> $tmpfolder/.tmp.log
printf "done.\n"
fi
printf "[EricScript] Removing temporary files ..."
rm -r $tmpfolder
printf "done.\n"
