# script by BioMart 
# feb16:  mod with selection of genome release

use strict;
use LWP::UserAgent;

my $ensversion = $ARGV[1];

open (FH,$ARGV[0]) || die ("\nUsage: perl webExample.pl Query.xml\n\n");

my $xml;
while (<FH>){
    $xml .= $_;
}
close(FH);

my $path = "";

if ($ensversion == 70) {
        $path = "http://jan2013.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 71) {
        $path = "http://apr2013.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 72) {
        $path = "http://jun2013.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 73) {
        $path = "http://sep2013.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 74) {
        $path = "http://dec2013.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 75) {
        $path = "http://feb2014.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 76) {
        $path = "http://aug2014.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 77) {
        $path = "http://oct2014.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 78) {
        $path = "http://dec2014.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 79) {
        $path = "http://mar2015.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 80) {
        $path = "http://may2015.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 81) {
        $path = "http://jul2015.archive.ensembl.org/biomart/martservice?";
}
elsif ($ensversion == 82) {
        $path = "http://sep2015.archive.ensembl.org/biomart/martservice?";
} 
elsif ($ensversion > 82) {
	$path = "http://www.ensembl.org/biomart/martservice?";
}
my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
my $ua = LWP::UserAgent->new;
my $response;

$ua->request($request, 
	     sub{   
		 my($data, $response) = @_;
		 if ($response->is_success) {
		     print "$data";
		 }
		 else {
		     warn ("Problems with the web server: ".$response->status_line);
		 }
	     },1000);

