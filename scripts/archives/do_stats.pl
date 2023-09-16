#!/usr/bin/perl -w

BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use Packages::General;

sub ProcessFile();
sub CalculateStatistics();

if (! @ARGV || $#ARGV < 2) {
   die "usage: $0 datafile start end\n";
}

my ($datfile, $start_time, $end_time) = @ARGV;
my (%Eng);

die "Invalid start time: Expected decimal, Got $start_time\n"
    if ( (! IsDecimal($start_time) ) && (! IsInteger($start_time) ) );

die "Invalid end time: Expected decimal, Got $end_time\n"
    if ( (! IsDecimal($end_time) ) && (! IsInteger($end_time) ) );

die "Cannot locate $datfile: $!\n"
    if (! -e $datfile);

ProcessFile();
CalculateStatistics();

sub ProcessFile() {
     my ($curr_line, $isValid);

     $isValid = 0;

     open INFILE, $datfile or die "Cannot open $datfile: $!\n";
     while (<INFILE>) {
	chomp;
	$curr_line = $_;
	if ($curr_line =~ /^\s+(\-?\d+\.\d+)\s+(\-\d+\.\d+)/) {
		if ($1 >= $start_time && $1 <= $end_time) {
			$Eng{$1} = $2;
			$isValid = 1;
		}
	}
     }
    
    close INFILE;
 	
    die "ERROR: $datfile is invalid\n" if (! $isValid);

}

sub CalculateStatistics() {

    my ($datString);
    for (keys %Eng) {
	$datString .= $Eng{$_} . " ";
    }

    chop $datString;

    my ($avg, $stdev, $total) = STDev($datString);

    $total = $total / $avg;
    
    printf "%6s%13s%13s\n", CenterText("Total", 6), 
			  CenterText("Average", 13),
			  CenterText("Std. Dev", 13); 
    printf "%6d%13.4f%13.4f\n", $total, $avg, $stdev;
}
