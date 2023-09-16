#!/usr/bin/perl -w

# This script will calculate the appropriate bin size for a series of data

use strict;
sub GetDataVals(@);
sub GetBinSize(@);

die "usage: $0 datafile\n"
    if (! @ARGV);

my ($datavals, $result, $range, $total_pts);
my ($datafile) = @ARGV;

die "Error cannot access file $datafile: $!\n"
    if (! -e $datafile or ! -r $datafile or ! -T $datafile);

$datavals = GetDataVals($datafile);
$total_pts = ($#{ $datavals } + 1);
print "Obtained $total_pts data vals\n";
($range, $result) = GetBinSize($datavals);

print "TOTAL PTS: $total_pts\n";
print "BIN SIZE: $result\n";
print "TOTAL BINS: " . int($range/$result) . "\n";

sub GetDataVals(@) {
    my ($infile) = $_[0];
    my (@DATA);

    open INFILE, $infile or die "Cannot open datafile $infile: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s+\d+\s+(\-?\d+\.\d+)/) {
	    push @DATA, $1;
	}
    }
    close INFILE;

    die "Error: $infile is invalid\n"
	if ($#DATA == -1);

    return (\@DATA);

}

sub GetBinSize(@) {
# This will implement the nonparametric density estimation method too calculate
# the bin size.
# This result was also obtained by Freedman and Diaconis sumarize in: 
# Izenman, A. J. 1991. "Recent developments in nonparametric density estimation."
#                       Journal of the American Statistical Association, 86(413):205-224. 
#
# W = 2 (IQR) N**-1/3 
#    where  IQR is the interquartile range 
#    (the 75th percentile minus the 25th percentile). 

    my (@datavals) = sort { $a <=> $b } @{ $_[0] };

    my ($IQR, $first_percentile, $third_percentile, $max_val, $min_val, $range);

    $first_percentile = sprintf("%.0f", .25 * ($#datavals + 1));
    $first_percentile = ($datavals[$first_percentile + 1] + $datavals[$first_percentile])/2;

    $third_percentile = sprintf("%.0f", .75 * ($#datavals + 1));
    $third_percentile = ($datavals[$third_percentile + 1] + $datavals[$third_percentile])/2;

    $IQR = $third_percentile - $first_percentile;
    my ($result) = 2 * ($IQR) / (($#datavals +1) ** (1/3));

    $max_val = $datavals[$#{ $datavals }];
    $min_val = $datavals[0];

    $range = $max_val - $min_val;

    print "MIN VAL: $min_val\n";
    print "MAX VAL: $max_val\n";

    return ($range, $result);


}
