#!/bin/bash
## bash pipeline
rnum=$1
. ~/.ericscript.$rnum.vars 
rm ~/.ericscript.$rnum.vars
if [ $verbose -eq 1 ]; then
touch $outputfolder/out/.ericscript.log
fi
printf "\n[EricScript] Starting EricScript analysis for sample $samplename.\n"
if [ $flagbin -eq 0 ]; then
readlength=$((`head -n 2 $reads_1 | tail -n 1 | wc -m` - 1))
else
readlength=$((`gunzip -c $reads_1 | head -n 2 | tail -n 1 | wc -m` - 1))
fi
if [ $ntrim -lt $readlength -a $ntrim -ne -1 ]; then
if [ $ntrim -lt 36 -a $ntrim -ge 0 ]; then
printf "[EricScript] Minimum trimming value is 36. Reads will be trimmed to 36 nt.\n"
ntrim=36
fi
if [ $ntrim -eq -1 -a $readlength -ge 70 ]; then
ntrim=50
fi
printf "[EricScript] Trimming PE reads to $ntrim nt ..."
if [ $flagbin -eq 0 ]; then
perl $ericscriptfolder/lib/perl/trimfq.pl $reads_1 $ntrim $outputfolder/aln/$samplename.1.fq.trimmed
perl $ericscriptfolder/lib/perl/trimfq.pl $reads_2 $ntrim $outputfolder/aln/$samplename.2.fq.trimmed
printf "done. \n"
else
gunzip -c $reads_1 | perl $ericscriptfolder/lib/perl/trimfq.pl - $ntrim $outputfolder/aln/$samplename.1.fq.trimmed
gunzip -c $reads_2 | perl $ericscriptfolder/lib/perl/trimfq.pl - $ntrim $outputfolder/aln/$samplename.2.fq.trimmed
fi
if [ $verbose -eq 0 ]; then
printf "[EricScript] Aligning with bwa ..."
if [ $bwa_aln -eq 1 ]; then
bwa aln -R 5 -t $nthreads $myref $outputfolder/aln/$samplename.1.fq.trimmed > $outputfolder/aln/"$samplename"_1.sai 2>> $outputfolder/out/.ericscript.log
bwa aln -R 5 -t $nthreads $myref $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/"$samplename"_2.sai 2>> $outputfolder/out/.ericscript.log
fi
if [ $MAPQ -gt 0 ]; then
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/"$samplename".sam 2>> $outputfolder/out/.ericscript.log
else
bwa mem -t $nthreads $myref $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/"$samplename".sam 2>> $outputfolder/out/.ericscript.log
fi
else
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
else
bwa mem -Y -t $nthreads $myref $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
fi
cat $outputfolder/aln/tmp.sam | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/"$samplename".sam
fi
else
printf "[EricScript] Aligning with bwa ...\n"
if [ $bwa_aln -eq 1 ]; then
bwa aln -t $nthreads $myref $outputfolder/aln/$samplename.1.fq.trimmed > $outputfolder/aln/"$samplename"_1.sai
bwa aln -t $nthreads $myref $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/"$samplename"_2.sai
fi
if [ $MAPQ -gt 0 ]; then
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/"$samplename".sam
else
bwa mem -t $nthreads $myref $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed > $outputfolder/aln/"$samplename".sam
fi
else
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/"$samplename".sam
else
bwa mem -Y -t $nthreads $myref $outputfolder/aln/$samplename.1.fq.trimmed $outputfolder/aln/$samplename.2.fq.trimmed | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/"$samplename".sam
fi
fi
fi
else
if [ $ntrim -ge $readlength ]; then
printf "[EricScript] Selected trimming value is greater equal to read length. Reads will not be trimmed.\n"
fi
if [ $verbose -eq 0 ]; then
printf "[EricScript] Aligning with bwa ..."
if [ $bwa_aln -eq 1 ]; then
bwa aln -R 5 -t $nthreads $myref $reads_1 > $outputfolder/aln/"$samplename"_1.sai 2>> $outputfolder/out/.ericscript.log
bwa aln -R 5 -t $nthreads $myref $reads_2 > $outputfolder/aln/"$samplename"_2.sai 2>> $outputfolder/out/.ericscript.log
fi
if [ $MAPQ -gt 0 ]; then
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $reads_1 $reads_2 > $outputfolder/aln/"$samplename".sam 2>> $outputfolder/out/.ericscript.log
else
bwa mem -t $nthreads $myref $reads_1 $reads_2 > $outputfolder/aln/"$samplename".sam 2>> $outputfolder/out/.ericscript.log
fi
else
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $reads_1 $reads_2 > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
else
bwa mem -Y -t $nthreads $myref $reads_1 $reads_2 > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
fi
cat $outputfolder/aln/tmp.sam | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/"$samplename".sam
fi
else
printf "[EricScript] Aligning with bwa ...\n"
if [ $bwa_aln -eq 1 ]; then
bwa aln -t $nthreads $myref $reads_1 > $outputfolder/aln/"$samplename"_1.sai 
bwa aln -t $nthreads $myref $reads_2 > $outputfolder/aln/"$samplename"_2.sai 
fi
if [ $MAPQ -gt 0 ]; then
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $reads_1 $reads_2 > $outputfolder/aln/"$samplename".sam 
else
bwa mem -t $nthreads $myref $reads_1 $reads_2 > $outputfolder/aln/"$samplename".sam
fi
else
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P -c 0.001 $myref $outputfolder/aln/"$samplename"_1.sai $outputfolder/aln/"$samplename"_2.sai $reads_1 $reads_2 | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/"$samplename".sam
else
bwa mem -Y -t $nthreads $myref $reads_1 $reads_2 | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/"$samplename".sam
fi
fi
fi
fi
printf "done. \n"
R --slave --args $outputfolder,$bwa_aln <  $ericscriptfolder/lib/R/ExtractInsertSize.R
printf "[EricScript] Extracting discordant alignments ... "
grep -v '^@' $outputfolder/aln/"$samplename".sam | awk -v mapq="$MAPQ" '(($7!="=") && ($7!="*") && ($5>=mapq)) { print }' | cut -f2,3,4,5,7,8 > $outputfolder/out/"$samplename".filtered.out
R --slave --args $samplename,$outputfolder,$ericscriptfolder,$minreads,$MAPQ,$refid,$dbfolder < $ericscriptfolder/lib/R/MakeAdjacencyMatrix.R
myflag=`cat $outputfolder/out/.ericscript.flag`
if [ $myflag -eq 0 ]; then
printf "done. \n"
printf "[EricScript] No chimeric transcripts found! Writing results ..."
R --slave --args $samplename,$outputfolder <  $ericscriptfolder/lib/R/MakeEmptyResults.R
printf "done. \n"
exit 1
fi
printf "done. \n"
printf "[EricScript] Building exon junction reference ... "
R --slave --args $samplename,$outputfolder,$ericscriptfolder,$readlength,$refid,$dbfolder < $ericscriptfolder/lib/R/BuildFasta.R
printf "done. \n"
## Aligning to putative junction reference
if [ $verbose -eq 0 ]; then
printf "[EricScript] Aligning to exon junction reference ... "
bwa index $mynewref 1>> $outputfolder/out/.ericscript.log 2>> $outputfolder/out/.ericscript.log
if [ $bwa_aln -eq 1 ]; then
bwa aln -t $nthreads $mynewref $reads_1 > $outputfolder/aln/"$samplename"_1.remap.sai 2>> $outputfolder/out/.ericscript.log
bwa aln -t $nthreads $mynewref $reads_2 > $outputfolder/aln/"$samplename"_2.remap.sai 2>> $outputfolder/out/.ericscript.log
fi
if [ $MAPQ -gt 0 ]; then
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P $mynewref $outputfolder/aln/"$samplename"_1.remap.sai $outputfolder/aln/"$samplename"_2.remap.sai $reads_1 $reads_2 > $outputfolder/aln/$samplename.remap.sam 2>> $outputfolder/out/.ericscript.log
else
bwa mem -t $nthreads $mynewref $reads_1 $reads_2 > $outputfolder/aln/$samplename.remap.sam 2>> $outputfolder/out/.ericscript.log
fi
else
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P $mynewref $outputfolder/aln/"$samplename"_1.remap.sai $outputfolder/aln/"$samplename"_2.remap.sai $reads_1 $reads_2 > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
else
bwa mem -Y -t $nthreads $mynewref $reads_1 $reads_2 > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
fi
cat $outputfolder/aln/tmp.sam | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/$samplename.remap.sam
fi
samtools view -@ $nthreads -bS -o $outputfolder/aln/$samplename.remap.bam $outputfolder/aln/$samplename.remap.sam 1>> $outputfolder/out/.ericscript.log 2>> $outputfolder/out/.ericscript.log
samtools sort -@ $nthreads $outputfolder/aln/$samplename.remap.bam $outputfolder/aln/$samplename.remap.sorted 1>> $outputfolder/out/.ericscript.log 2>> $outputfolder/out/.ericscript.log
samtools index $outputfolder/aln/$samplename.remap.sorted.bam 1>> $outputfolder/out/.ericscript.log
else
printf "[EricScript] Aligning to exon junction reference ... \n"
bwa index $mynewref
if [ $bwa_aln -eq 1 ]; then
bwa aln -t $nthreads $mynewref $reads_1 > $outputfolder/aln/"$samplename"_1.remap.sai 
bwa aln -t $nthreads $mynewref $reads_2 > $outputfolder/aln/"$samplename"_2.remap.sai 
fi
if [ $MAPQ -gt 0 ]; then
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P $mynewref $outputfolder/aln/"$samplename"_1.remap.sai $outputfolder/aln/"$samplename"_2.remap.sai $reads_1 $reads_2 > $outputfolder/aln/$samplename.remap.sam 
else
bwa mem -t $nthreads $mynewref $reads_1 $reads_2 > $outputfolder/aln/$samplename.remap.sam
fi
else
if [ $bwa_aln -eq 1 ]; then
bwa sampe -P $mynewref $outputfolder/aln/"$samplename"_1.remap.sai $outputfolder/aln/"$samplename"_2.remap.sai $reads_1 $reads_2 | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/$samplename.remap.sam 
else
bwa mem -Y -t $nthreads $mynewref $reads_1 $reads_2 | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/$samplename.remap.sam
fi
fi
samtools view -@ $nthreads -bS -o $outputfolder/aln/$samplename.remap.bam $outputfolder/aln/$samplename.remap.sam
samtools sort -@ $nthreads $outputfolder/aln/$samplename.remap.bam $outputfolder/aln/$samplename.remap.sorted
samtools index $outputfolder/aln/$samplename.remap.sorted.bam
fi
printf "done. \n"
## Recalibrating junctions
printf "[EricScript] Recalibrating junctions ... "
R --slave --args $samplename,$outputfolder,$readlength,$verbose < $ericscriptfolder/lib/R/RecalibrateJunctions.R
cat $outputfolder/out/$samplename.EricScript.junctions.recalibrated.fa $myref > $mynewref_recal
printf "done. \n"
## Aligning not properly mapped reads
if [ $verbose -eq 0 ]; then
printf "[EricScript] Aligning to recalibrated junction reference ... "
bwa index $mynewref_recal 1>> $outputfolder/out/.ericscript.log 2>> $outputfolder/out/.ericscript.log
if [ $bwa_aln -eq 1 ]; then
bwa aln -R 5 -t $nthreads $mynewref_recal $reads_1 > $outputfolder/aln/"$samplename"_1.remap.recal.sai 2>> $outputfolder/out/.ericscript.log
bwa aln -R 5 -t $nthreads $mynewref_recal $reads_2 > $outputfolder/aln/"$samplename"_2.remap.recal.sai 2>> $outputfolder/out/.ericscript.log
bwa sampe -P $mynewref_recal $outputfolder/aln/"$samplename"_1.remap.recal.sai $outputfolder/aln/"$samplename"_2.remap.recal.sai $reads_1 $reads_2 > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
else
bwa mem -Y -t $nthreads $mynewref_recal $reads_1 $reads_2 > $outputfolder/aln/tmp.sam 2>> $outputfolder/out/.ericscript.log
fi
cat $outputfolder/aln/tmp.sam | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/$samplename.remap.recal.sam
samtools view -@ $nthreads -bt $mynewref_recal -o $outputfolder/aln/$samplename.remap.recal.bam $outputfolder/aln/$samplename.remap.recal.sam 1>> $outputfolder/out/.ericscript.log 2>> $outputfolder/out/.ericscript.log
samtools sort -@ $nthreads $outputfolder/aln/$samplename.remap.recal.bam $outputfolder/aln/$samplename.remap.recal.sorted 1>> $outputfolder/out/.ericscript.log 2>> $outputfolder/out/.ericscript.log
samtools rmdup $outputfolder/aln/$samplename.remap.recal.sorted.bam $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam 1>> $outputfolder/out/.ericscript.log 2>> $outputfolder/out/.ericscript.log
samtools index $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam 1>> $outputfolder/out/.ericscript.log 
samtools view -@ $nthreads -b -h -q 1 $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam > $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.q1.bam
samtools index $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.q1.bam
else
printf "[EricScript] Aligning to recalibrated junction reference ... \n"
bwa index  $mynewref_recal 
if [ $bwa_aln -eq 1 ]; then
bwa aln -R 5 -t $nthreads $mynewref_recal $reads_1 > $outputfolder/aln/"$samplename"_1.remap.recal.sai 
bwa aln -R 5 -t $nthreads $mynewref_recal $reads_2 > $outputfolder/aln/"$samplename"_2.remap.recal.sai 
bwa sampe -P $mynewref_recal $outputfolder/aln/"$samplename"_1.remap.recal.sai $outputfolder/aln/"$samplename"_2.remap.recal.sai $reads_1 $reads_2 | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/$samplename.remap.recal.sam 
else
bwa mem -Y -t $nthreads $mynewref_recal $reads_1 $reads_2 | $ericscriptfolder/lib/perl/xa2multi.pl > $outputfolder/aln/$samplename.remap.recal.sam
fi
samtools view -@ $nthreads -bt $mynewref_recal -o $outputfolder/aln/$samplename.remap.recal.bam $outputfolder/aln/$samplename.remap.recal.sam 
samtools sort -@ $nthreads $outputfolder/aln/$samplename.remap.recal.bam $outputfolder/aln/$samplename.remap.recal.sorted
samtools rmdup $outputfolder/aln/$samplename.remap.recal.sorted.bam $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam 
samtools index $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam 
samtools view -@ $nthreads -b -h -q 1 $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam > $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.q1.bam
samtools index $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.q1.bam
fi
printf "done. \n"
samtools idxstats $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.q1.bam > $outputfolder/out/$samplename.stats
rm $outputfolder/aln/*.sam
## Estimating spanning reads
printf "[EricScript] Scoring candidate fusions ..."
R --slave --args $samplename,$outputfolder,$readlength <  $ericscriptfolder/lib/R/EstimateSpanningReads.R
myflag=`cat $outputfolder/out/.ericscript.flag`
if [ $myflag -eq 0 ]; then
printf "done. \n"
printf "[EricScript] No chimeric transcripts found! Writing results ..."
R --slave --args $samplename,$outputfolder <  $ericscriptfolder/lib/R/MakeEmptyResults.R
printf "done. \n"
exit 1
fi
printf "done. \n"
printf "[EricScript] Filtering candidate fusions ..."
if [ $verbose -eq 0 ]; then
samtools mpileup -A -f $mynewref_recal -l $outputfolder/out/$samplename.intervals $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam > $outputfolder/out/$samplename.remap.recal.sorted.rmdup.pileup 2>> $outputfolder/out/.ericscript.log
else
samtools mpileup -A -f $mynewref_recal -l $outputfolder/out/$samplename.intervals $outputfolder/aln/$samplename.remap.recal.sorted.rmdup.bam > $outputfolder/out/$samplename.remap.recal.sorted.rmdup.pileup
fi
cut -f1,2,3 $outputfolder/out/$samplename.remap.recal.sorted.rmdup.pileup | grep -e '[0-9]----[a-z | A-Z]' - >  $outputfolder/out/$samplename.intervals.pileup
R --slave --args $samplename,$outputfolder <  $ericscriptfolder/lib/R/BuildNeighbourhoodSequences.R
if [ $verbose -eq 0 ]; then
blat $myref $outputfolder/out/.link $outputfolder/out/$samplename.checkselfhomology.blat -out=blast8 1>> $outputfolder/out/.ericscript.log
else
blat $myref $outputfolder/out/.link $outputfolder/out/$samplename.checkselfhomology.blat -out=blast8
fi
R --slave --args $samplename,$outputfolder <  $ericscriptfolder/lib/R/CheckSelfHomology.R
myflag=`cat $outputfolder/out/.ericscript.flag`
if [ $myflag -eq 0 ]; then
printf "done. \n"
printf "[EricScript] No chimeric transcripts found! Writing results ..."
R --slave --args $samplename,$outputfolder <  $ericscriptfolder/lib/R/MakeEmptyResults.R
printf "done. \n"
exit 1
fi
printf "done. \n"
## Writing results
if [ $verbose -eq 0 ]; then
printf "[EricScript] Writing results ... "
else
printf "[EricScript] Writing results ... \n"
fi
R --slave --args $samplename,$outputfolder,$ericscriptfolder,$readlength,$verbose,$refid,$dbfolder <  $ericscriptfolder/lib/R/MakeResults.R
printf "done. \n"
rm $outputfolder/out/*fai
if [ $removetemp -eq 1 ]; then
printf "[EricScript] Removing temporary files ... "
rm -r $outputfolder/aln
rm -r $outputfolder/out
printf "done. \n"
fi 
printf "[EricScript] Open $outputfolder/$samplename.results* to view the results of EricScript analysis.\n"
