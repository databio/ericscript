#!/usr/bin/perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use Cwd 'abs_path';

our ($verbose, $help, $man);
our ($samplename, $reads_1, $reads_2, $outputfolder, $minreads, $removetemp, $nthreads, $MAPQ, $checkdb, $demo, $refid, $printdb, $downdb, $dbfolder, $recalc, $bwa_aln, $genomeref);
our ($simulator, $readlength, $ntrim, $insize, $sd_insize, $ngenefusion, $min_cov, $max_cov, $nsims, $be, $ie, $background_1, $background_2, $nreads_background, $ensversion);
our($calcstats, $resultsfolder, $datafolder, $algoname, $dataset, $normroc);
my @command_line = @ARGV;
GetOptions('verbose|v'=>\$verbose, 
'help|h'=>\$help, 
'man|m'=>\$man, 
'samplename|name=s'=>\$samplename,
'outputfolder|o=s'=>\$outputfolder, 
'dbfolder|db=s'=>\$dbfolder, 
'background_1=s'=>\$background_1, 
'background_2=s'=>\$background_2, 
'refid=s'=>\$refid, 
'minreads|minr=i'=>\$minreads,
'remove'=>\$removetemp,
'nthreads|p=i'=>\$nthreads,
'readlength|rl=i'=>\$readlength,
'ntrim=i'=>\$ntrim,
'insize=f'=>\$insize,
'sd_insize=f'=>\$sd_insize,
'ngenefusion=i'=>\$ngenefusion,
'min_cov=i'=>\$min_cov,
'max_cov=i'=>\$max_cov,
'nsims=i'=>\$nsims,
'nreads_background=i'=>\$nreads_background,
'checkdb'=>\$checkdb,
'ie'=>\$ie,
'be'=>\$be,
'demo'=>\$demo,
'simulator'=>\$simulator,
'recalc'=>\$recalc,
'printdb'=>\$printdb,
'checkdb'=>\$checkdb,
'downdb'=>\$downdb,
'bwa_aln'=>\$bwa_aln,
'calcstats'=>\$calcstats,
'resultsfolder=s'=>\$resultsfolder, 
'datafolder=s'=>\$datafolder, 
'algoname=s'=>\$algoname, 
'dataset=s'=>\$dataset, 
'normroc=i'=>\$normroc,
'ensversion=i'=>\$ensversion,
'MAPQ=f'=>\$MAPQ) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);

my $sysmsg;
my $sysflag=0;
$sysmsg = qx/samtools 2>&1/ || '';		
if ($sysmsg !~ m/bamshuf/) {
	print STDERR "[EricScript] Error: SAMtools >= 0.1.19 not found! Please install and add it to your PATH.\n";
	$sysflag=1;
}
if ($sysmsg =~ m/htslib/) {
    print STDERR "[EricScript] Error: SAMtools >= 1.0 detected! EricScript is not yet compatible with it. Please use samtools 0.1.19 to run EricScript.\n";
    $sysflag=1;
}
$sysmsg = qx/bwa 2>&1/ || '';		
if ($sysmsg !~ m/sampe/) {
	print STDERR "[EricScript] Error: BWA not found! Please install and add it to your PATH.\n";
	$sysflag=1;
} else {
	$sysmsg = qx/bwa mem 2>&1/ || '';
	if ($sysmsg !~ m/-Y/) {
	    print STDERR "[EricScript] Error: BWA >= 0.7.12 not found! Please install and add it to your PATH.\n";
	    $sysflag=1;
	}
}
#} elsif ($sysmsg !~ m/mem/) {
#    print STDERR "[EricScript] Error: BWA >= 0.7.4 not found! Please install and add it to your PATH.\n";
#    $sysflag=1;
#}
$sysmsg = qx/blat 2>&1/ || '';		
if ($sysmsg !~ "tileSize") {
	print STDERR "[EricScript] Error: BLAT not found! Please install and add it to your PATH.\n";
	$sysflag=1;
}
$sysmsg = qx/R --help 2>&1/ || '';		
if ($sysmsg !~ "Rdiff") {
	print STDERR "[EricScript] Error: R not found! Please install and add it to your PATH.\n";
	$sysflag=1;
}
$sysmsg = qx/seqtk 2>&1/ || '';
if ($sysmsg !~ "subseq") {
    print STDERR "[EricScript] Error: Seqtk not found! Please install and add it to your PATH.\n";
    $sysflag=1;
}
$sysmsg = qx/bedtools 2>&1/ || '';
if ($sysmsg !~ "sample") {
    print STDERR "[EricScript] Error: bedtools (>=2.18) not found! Please install and add it to your PATH.\n";
    $sysflag=1;
}

if ($sysflag == 1) {
	exit(100);
}

my $myscriptname = "ericscript.pl";
my $myscript = abs_path($0);
my $ericscriptfolder=substr $myscript, 0, (length($myscript) - length($myscriptname) - 1) ;

$calcstats ||= 0;
$simulator ||= 0;
$checkdb ||= 0;
$demo ||= 0;
$printdb ||=0;
$downdb ||=0;
$recalc ||=0;
$bwa_aln ||= 0;
$ensversion ||=0;
#if ($bwa_aln != 0) {
#	$bwa_aln = 1;
#}
$dbfolder ||= "$ericscriptfolder/lib";

if ( ($ensversion > 0) && ($ensversion < 70) ) {
	pod2usage ("[EricScript] Error in selecting Ensembl release. Minimum supported version is 70.\n");
}

if ($downdb != 0) {
	$refid or pod2usage ("[EricScript] Syntax error. \"refid\" parameter needs to be specified. Run ericscript.pl --printdb to see the list of available refid.\n");
} else {
	$refid ||= "homo_sapiens";
}

if ($printdb != 0) {
    system("R --slave --args $ericscriptfolder,$printdb,$dbfolder,$ensversion < $ericscriptfolder/lib/R/RetrieveRefId.R");
	exit(100);
}

if ($recalc != 0) {
    system("R --slave --args $ericscriptfolder,$refid,$dbfolder < $ericscriptfolder/lib/R/CheckDB.R");
    my $flagdb;
    open FILE, "< $ericscriptfolder/lib/data/_resources/.flag.dbexists";
    $flagdb = <FILE>;
    if ($flagdb == 0) {
        pod2usage ();
    }

     	## check inputs !! ONLY FOR DEBUG purposes
    my $file3 = File::Spec->catfile ($genomeref);
    -f $file3 or pod2usage ("[EricScript] Error: please specify a genome reference file.\n");
    my $file4 = File::Spec->catfile ("$genomeref.pac");
    -f $file4 or pod2usage ("[EricScript] Error: BWA indexes for $genomeref not found. Create BWA indexes then run EricScript.\n");
	$outputfolder or pod2usage ("[EricScript] Error: Please specify where past analysis is stored by using --outputfolder");
#	$outputfolder = abs_path($outputfolder);
	if (-d $outputfolder) { 
		my $abs_outputfolder = abs_path($outputfolder);
		system("R --slave --args $ericscriptfolder,$abs_outputfolder,$dbfolder,$refid,$genomeref < $ericscriptfolder/lib/R/CalcBreakpointPositions.R");	
    } else {
        die "[EricScript] Error: output folder $outputfolder does not exist.\n";
	}

}

if ($checkdb == 0 && $demo == 0 && $simulator == 0 && $calcstats == 0 && $downdb == 0 && $recalc == 0) {
	
	@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
	@ARGV == 2 or pod2usage ("[EricScript] Syntax error.\n");
	($reads_1, $reads_2) = @ARGV;

		## check db existence
	system("R --slave --args $ericscriptfolder,$refid,$dbfolder < $ericscriptfolder/lib/R/CheckDB.R");
	my $flagdb;
	open FILE, "< $ericscriptfolder/lib/data/_resources/.flag.dbexists";
	$flagdb = <FILE>;
	if ($flagdb == 0) {
		exit(100);
	}	

		## check inputs
	my $file1 = File::Spec->catfile ($reads_1);
	-f $file1 or pod2usage ("[EricScript] Error: the required file $reads_1 does not exist.\n");
	my $file2 = File::Spec->catfile ($reads_2);
    -f $file2 or pod2usage ("[EricScript] Error: the required file $reads_2 does not exist.\n");
#	my $file3 = File::Spec->catfile ($genomeref);
#    -f $file3 or pod2usage ("[EricScript] Error: please specify a valid genome reference file.\n");
#	my $file4 = File::Spec->catfile ("$genomeref.pac");
#    -f $file4 or pod2usage ("[EricScript] Error: BWA indexes for $genomeref not found. Create BWA indexes then run EricScript.\n");


	my $userhome = $ENV{HOME}; 
	$samplename ||= 'MyEric';
	$outputfolder ||= "$userhome/$samplename";
#	$outputfolder = abs_path($outputfolder);
	if (-d $outputfolder) {
		die "[EricScript] Error: output folder $outputfolder already exists.\n";
	}
	mkdir($outputfolder) || die "[EricScript] Error: the directory $outputfolder is not writable by the current user. \n";
	mkdir("$outputfolder/aln");
	mkdir("$outputfolder/out");
	$verbose ||= 0;			
	$minreads ||= 3;
	$ntrim ||= -1;
	$removetemp ||= 0;
	$MAPQ ||= 20;
	$nthreads ||= 4;
	my $myref="$dbfolder/data/$refid/EnsemblGene.Reference.fa";
	my $mynewref="$outputfolder/out/$samplename.EricScript.junctions.fa";
	my $mynewref_recal="$outputfolder/out/$samplename.EricScript.fa";
	my $flagbin= 1;
	if (-T $reads_1) {
		$flagbin = 0;
	}
	## write vars
	my $abs_outputfolder = abs_path($outputfolder);
	my $range = 10000;
	my $rnum = int(rand($range));
	my $varfile = "$userhome/.ericscript.$rnum.vars";
	open(FILE, ">", "$varfile") or die "Couldn't open: $!";
	print FILE "samplename=\"$samplename\"\n";
	print FILE "outputfolder=\"$abs_outputfolder\"\n";
	print FILE "dbfolder=\"$dbfolder\"\n";
	print FILE "reads_1=\"$reads_1\"\n";
	print FILE "reads_2=\"$reads_2\"\n";
	print FILE "minreads=$minreads\n";
	print FILE "ntrim=$ntrim\n";
	print FILE "MAPQ=$MAPQ\n";
	print FILE "removetemp=$removetemp\n";
	print FILE "flagbin=$flagbin\n";
	print FILE "nthreads=$nthreads\n";
	print FILE "ericscriptfolder=\"$ericscriptfolder\"\n";
	print FILE "myref=\"$myref\"\n";
	print FILE "mynewref=\"$mynewref\"\n";
	print FILE "mynewref_recal=\"$mynewref_recal\"\n";
	print FILE "refid=\"$refid\"\n";
	print FILE "verbose=$verbose\n";
	print FILE "bwa_aln=$bwa_aln\n";
	close (FILE); 
	
	system("cp", "$varfile", "$outputfolder/out/.ericscript.vars");
	system("bash", "$ericscriptfolder/lib/bash/RunEric.sh", "$rnum");
} elsif ($checkdb != 0) {
	if (-d $dbfolder) {
		print STDERR "[EricScript] Checking installed Database.\n";
		system("R --slave --args $ericscriptfolder,$printdb,$dbfolder,$ensversion < $ericscriptfolder/lib/R/RetrieveRefId.R");
		system("R --slave --args $ericscriptfolder,$dbfolder < $ericscriptfolder/lib/R/UpdateDB.R");
		exit(100);
	} else {
		die "[EricScript] Error: the directory $dbfolder does not exist.\n";
	}
} elsif ($downdb != 0) {
	if ( -d $dbfolder ) {
		if ( !-d "$dbfolder/data" ) {
			mkdir ("$dbfolder/data");
		}
	system("R --slave --args $ericscriptfolder,$printdb,$dbfolder,$ensversion < $ericscriptfolder/lib/R/RetrieveRefId.R");
	system("bash $ericscriptfolder/lib/bash/BuildSeq.sh $ericscriptfolder $refid $dbfolder $ensversion");	
	} else {
		die "[EricScript] Error: the directory $dbfolder does not exist.\n";
	}

} elsif ($demo != 0) {
        ## check db
#    my $file3 = File::Spec->catfile ($genomeref);
#    -f $file3 or die ("[EricScript] Error: please specify a genome reference file.\n");
    system("R --slave --args $ericscriptfolder,$refid,$dbfolder < $ericscriptfolder/lib/R/CheckDB.R");
    my $flagdb;
    open FILE, "< $ericscriptfolder/lib/data/_resources/.flag.dbexists";
    $flagdb = <FILE>;
    if ($flagdb == 0) {
        pod2usage ();
    }
	my $userhome = $ENV{HOME}; 
	$samplename ='demo';
	$outputfolder ||= "$userhome/ericscript_demo";
#    $outputfolder = abs_path($outputfolder);
	if (-d $outputfolder) {
		die "[EricScript] Error: output folder $outputfolder already exists.\n";
	}	
	mkdir($outputfolder) || die "[EricScript] Error: the directory $outputfolder is not writable by the current user. \n";
	mkdir("$outputfolder/aln");
	mkdir("$outputfolder/out");	
	$verbose=0;			
	$minreads ||= 3;
	$ntrim ||= -1;
	$removetemp ||= 0;
	$MAPQ ||= 20;
	$nthreads ||= 4;
	$refid ||= "homo_sapiens";
	my $myref="$dbfolder/data/$refid/EnsemblGene.Reference.fa";
	my $mynewref="$outputfolder/out/$samplename.EricScript.junctions.fa";
	my $mynewref_recal="$outputfolder/out/$samplename.EricScript.fa";
	my $reads_1="$ericscriptfolder/lib/demo/myreads_1.fq.gz";
	my $reads_2="$ericscriptfolder/lib/demo/myreads_2.fq.gz";	
	my $flagbin= 1;
	if (-T $reads_1) {
		$flagbin = 0;
	}
	my $abs_outputfolder = abs_path($outputfolder);
    my $range = 10000;
    my $rnum = int(rand($range));
    my $varfile = "$userhome/.ericscript.$rnum.vars";
	open(FILE, ">", "$varfile") or die "Couldn't open: $!";
	print FILE "samplename=\"$samplename\"\n";
	print FILE "outputfolder=\"$abs_outputfolder\"\n";
	print FILE "dbfolder=\"$dbfolder\"\n";
	print FILE "reads_1=\"$reads_1\"\n";
	print FILE "reads_2=\"$reads_2\"\n";
	print FILE "flagbin=$flagbin\n";
	print FILE "minreads=$minreads\n";
	print FILE "ntrim=$ntrim\n";
	print FILE "MAPQ=$MAPQ\n";
	print FILE "removetemp=$removetemp\n";
	print FILE "nthreads=$nthreads\n";
	print FILE "ericscriptfolder=\"$ericscriptfolder\"\n";
	print FILE "myref=\"$myref\"\n";
	print FILE "mynewref=\"$mynewref\"\n";
	print FILE "mynewref_recal=\"$mynewref_recal\"\n";
	print FILE "refid=\"$refid\"\n";
	print FILE "verbose=$verbose\n";
	print FILE "bwa_aln=$bwa_aln\n";
	close (FILE); 
	
	system("cp", "$varfile", "$outputfolder/out/.ericscript.vars");
	system("bash", "$ericscriptfolder/lib/bash/RunEric.sh", "$rnum");
	
} 
elsif ($simulator != 0) {
	$sysmsg = qx/wgsim --help 2>&1/ || '';	
	my $userhome = $ENV{HOME}; 
	if ($sysmsg !~ "outer") {
    	print STDERR "[EricScript] Error: wgsim not found! Please install and add it to your PATH.\n";
    	$sysflag=1;
	}
	if ($sysflag == 1) {
    	exit(100);
	}
    $outputfolder ||= "$userhome/ericscript_simulator";
#    $outputfolder = abs_path($outputfolder);
	if (-d $outputfolder) {  
        die "[EricScript] Error: output folder $outputfolder already exists.\n";
   }
    mkdir($outputfolder) || die "[EricScript] Error: the directory $outputfolder is not writable by the current user. \n";
    $verbose ||= 0;
    $readlength ||= 75;
	$insize ||= 200;
	$sd_insize ||= 50;
	$ngenefusion ||= 50;
	$min_cov ||= 1;
	$max_cov ||= 50;
	$nsims ||= 10;
	$be ||= 0;
	$ie ||= 0;
	$refid ||= "homo_sapiens";
	if ($be == 0 && $ie == 0) {
		$ie=1;
	}
	$background_1 ||= 0;
	$background_2 ||= 0;
	$nreads_background ||= 200000;
	my $abs_outputfolder = abs_path($outputfolder);
	my $simcommand="R --slave --args $readlength,$abs_outputfolder,$ericscriptfolder,$verbose,$insize,$sd_insize,$ngenefusion,$min_cov,$max_cov,$nsims,$be,$ie,$background_1,$background_2,$nreads_background,$dbfolder,$refid < $ericscriptfolder/lib/R/SimulateFusions.R";
	system($simcommand);	

}
elsif ($calcstats != 0) {

	my $userhome = $ENV{HOME}; 
	$outputfolder ||= "$userhome/ericscript_stats";
#    $outputfolder = abs_path($outputfolder);
  	if (-d $outputfolder) {  
#  	    die "[EricScript] Error: output folder $outputfolder already exists.\n";
 	} else {
	    mkdir($outputfolder) || die "[EricScript] Error: the directory $outputfolder is not writable by the current user. \n";	
	}
	$resultsfolder || pod2usage ("[EricScript] Error: Argument resultsfolder is not specified! \n");
	$datafolder || pod2usage ("[EricScript] Error: Argument datafolder is not specified! \n");
	$algoname || pod2usage ("[EricScript] Error: Argument algoname is not specified! \n");
	$dataset || pod2usage ("[EricScript] Error: Argument dataset is not specified! \n");
	$readlength || pod2usage ("[EricScript] Error: Argument readlength is not specified! \n");
	$normroc ||= 1;
    my $abs_outputfolder = abs_path($outputfolder);
	-e $resultsfolder ||  die "[EricScript] Error: the folder $resultsfolder does not exist. \n";
	-e $datafolder || die "[EricScript] Error: the folder $datafolder does not exist. \n";
	my $datafolder1="$datafolder/$dataset";
	-e $datafolder1 || die "[EricScript] Error: no $dataset synthetic data have been found in $datafolder. \n";

	my $calcstatscommand="R --slave --args $resultsfolder,$abs_outputfolder,$datafolder,$algoname,$dataset,$readlength,$normroc,$ericscriptfolder < $ericscriptfolder/lib/R/CalcStats.R";
	system($calcstatscommand);

}


=head1 SYNOPSIS
 
 ericscript.pl [arguments] <reads_1.fq(.gz)> <reads_2.fq(.gz)>
 
	Optional arguments:
	-h, --help                      print help message
	-m, --man                       print complete documentation
	-v, --verbose                   use verbose output
	-name, --samplename <string>	what's the name of your sample?
	-o, --outputfolder <string>	where the results will be stored
	-db, --dbfolder <string>	where database is stored. Default is ERICSCRIPT_FOLDER/lib/
	-minr, --minreads <int>		minimum reads to consider discordant alignments [3]
	-p, --nthreads <int>		number of threads for the bwa aln process [4]
	-ntrim <int>			trim PE reads from 1st base to $ntrim. Default is no trimming. Set ntrim=0 to don't trim reads.
	--MAPQ <int>			minimum value of mapping quality to consider discordant reads. For MAPQ 0 use a negative value [20]
	--remove			remove all temporary files.
	--demo				Run a demonstration of EricScript on simulated reads.
 	--refid				Genome reference identification. Run ericscript.pl --printdb to see available refid [homo_sapiens].
	--bwa_aln			Use BWA ALN instead of BWA MEM to search for discordant reads.

	Subcommands:
	--checkdb			Check if your database is up-to-date, based on the latest Ensembl release.
	--downdb			Download, build database. refid parameter need to be specified.
	--simulator			Generate synthetic gene fusions with the same recipe of the ericscript's paper
	--calcstats			Calculate the statistics that we used in our paper to evaluate the performance of the algorithms.
	
	--------
	arguments for databases subcommands (downdb, checkdb):

 		-db, --dbfolder <string>	where database is stored. Default is ERICSCRIPT_FOLDER/lib/
 		--refid				Genome reference identification. Run ericscript.pl --printdb to see available refid [homo_sapiens].
 		--printdb			Print a list of available genomes and exit.
		--ensversion		Download data of a specific Ensembl version (>= 70). Default is the latest one.
 
	-------
	arguments for simulator:
 		-o, --outputfolder <string>	where synthetic datasets will be stored [HOME/ericscript_simulator]
		-rl, --readlength <int>		length of synthetic reads [75]
	 	--refid				Genome reference identification. Run ericscript.pl --printdb to see available refid [homo_sapiens].
		-v, --verbose               	use verbose output
		--insize			parameter of wgsym. Outer distance between the two ends [200]
		--sd_insize			parameter of wgsym. Standard deviation [50]
		--ngenefusion			The number of synthetic gene fusions per dataset? [50]
		--min_cov			Minimum coverage to simulate [1]
		--max_cov			Maximum coverage to simulate [50]
		--nsims				The number of synthetic datasets to simulate [10]
		--be				Use --be to generate Broken Exons (BE) data [no]
		--ie 				Use --ie to generate Intact Exons (IE) data [yes]
		-db, --dbfolder			where database is stored. Default is ERICSCRIPT_FOLDER/lib/ 
		--background_1			Fastq file (forward)  for generating background reads. 
		--background_2			Fastq file (reverse) for generating background reads. 
		--nreads_background		The number of reads to extract from background data [200e3].

	-------
	arguments for calcstats:
 		-o, --outputfolder <string>	where statistics file will be stored [HOME/ericscript_calcstats]
		--resultsfolder <string>	path to folder containing algorithm results.
		--datafolder <string>		path to folder containing synthetic data generated by ericscript simulator.
		--algoname <string>		name of the algorithm that generated results. 
		--dataset <string>		type of synthetic data to considered for calculating statistics. IE or BE? 
		-rl, --readlength <int>		length of synthetic reads 
		--normroc <int>		 	factor to normalize the score given by the algorithm.
			
 ericscript.pl automatically runs a pipeline to detect chimeric transcripts in
 paired-end RNA-seq samples. It is also able to generate datasets with synthetic gene fusions.
 More information about running EricScript Simulator and EricScript CalcStats can be 
 found at http://ericscript.sourceforge.net
 
 Version: 0.5.5b
 
=head1 OPTIONS
 
=over 8

=item B<--help>
 
 print a brief usage message and detailed explanation of options.
 
=item B<--man>
 
 print the complete manual of the program.
 
=item B<--verbose>
 
 use verbose output.
 
 
=item B<--samplename>
 
 Choose a name for your sample. Default is "MyEric"
 
=item B<--outputfolder>
 
 Folder that will contain all the results of the analysis. Default is YOUR_HOME/SAMPLENAME
 
=item B<--dbfolder>
 
 Folder that contains transcriptome sequences and information of the downloaded species. Default is
 ERICSCRIPT_FOLDER/lib

=item B<--minreads>
 
 Minimum reads to consider discordant alignments. Default is 3 reads with minimum MAPQ.

=item B<-ntrim>
 
 trim PE reads from 1st base to $ntrim. Trimmed reads will be used only for the first alignment (identification 
 of discordant reads). Setting ntrim to values lower than orginal read length allows EricScript to 
 increase its sensitivity, especially when the length of reads is 75nt or 100 nt. 
 Default is no trimming. Set ntrim=0 to don't trim reads.
 
=item B<--nthreads >
 
 Number of threads for the bwa aln process.
 
=item B<--MAPQ >
 
 minimum value of mapping quality to consider discordant reads. For MAPQ 0 use a negative value. Default is 20.
 
=item B<--remove>
 
 remove all temporary files. By default, all temporary files will be kept for 
 user inspection, but this will easily clutter the directory.
 
=item B<--checkdb>
 
 Check if your database is up-to-date, based on the latest Ensembl release.

=item B<--downdb>
 
 Download, build database. refid parameter need to be specified.

=item B<--refid>

 Genome reference identification. Run ericscript.pl --printdb to see available refid.[homo_sapiens]

=item B<--ensversion>
    
 Download data of a specific version of Ensembl. Default is downloading the latest version of Ensembl.
 Minimum supported version is 70.

=item B<--printdb>

 Print a list of available genomes and exit.

=item B<--demo>
 
 Run a demonstration of EricScript on simulated reads.
 
=back
 
=head1 DESCRIPTION
 
 EricScript (chimEric tranScript Detection Algorithm) is a computational framework for the discovery of 
 gene fusion products in paired end RNA-seq data.
 EricScript is freely available to the community for non-commercial use under GPLv3 license. For 
 questions or comments, please contact matteo.benelli@gmail.com.
 
 
=cut
