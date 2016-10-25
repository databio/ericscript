#!/usr/bin/perl
use warnings;
use strict;

my($listin) = $ARGV[0];
my($ntrim) = $ARGV[1];
my($outfile) = $ARGV[2];

open LIST, "${listin}" or die $!;
open OUT, ">$outfile";
my $a;
my $count = 1;
while (<LIST>)

   {
	if ($count++ % 2 == 0) {
   	$a = substr($_, 0, $ntrim);
   	print OUT "$a\n";
   	} else 
   	{
   	print OUT "$_"; 
   	}      
   
   }

close LIST;
close OUT;
