#!/usr/bin/perl -w
#
use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub parseDataFile;
sub getEquilPoint;
sub numerically { ($a<=>$b); }

my ($dataFile, $avg, $prec, $col);
my ($DATA, $equilPoint, @AVG, $i);

$|++;
&init;
print "Creating reverse running averages...";
$DATA = parseDataFile($dataFile, $col);
print "Done\nObtaining equil point...";
$equilPoint = getEquilPoint($DATA, $avg, $prec);
print "Done\nEquil Pt: $equilPoint\n";
open OUTFILE, "> avg.dat";
for $i (@AVG) {
	print OUTFILE "$i\n";
}
close OUTFILE;
sub getEquilPoint {
	my ($DATA, $avg, $prec) = @_;
	my (@tsteps, $i, $sum, $runningAvg, $count, $equilPt);

	if ($prec) {
		$prec = "%.${prec}f";
	} else {
		$prec = "%d";
	}
	@tsteps = sort numerically keys %{ $DATA };
	$runningAvg = $count = $sum = 0;
	$count = 1;
	$equilPt = $tsteps[$#tsteps];
	for $i (reverse @tsteps) {
		$sum += $DATA->{$i};
		push @AVG, ($sum/$count);
		$runningAvg = sprintf("$prec", ($sum/$count));
		$equilPt = $i if ($runningAvg == $avg) ;
		$count++;
	}

	return $equilPt;
}

sub parseDataFile {
	my ($datFile, $scol) = @_;
	my (%DATA, $tmp);

	$scol--;
	open DATAFILE, $datFile or die "ERROR: Cannot open datafile $datFile: $!\n";
	while (<DATAFILE>) {
		chomp;
		if ($_ =~ /^\s*(\d+\.?\d*)\s+(-?\d+\.?\d*.*)/) {
			@{ $tmp } = split /\s+/,$2;
			next if ($scol > scalar(@{ $tmp }));
			$DATA{$1}=$tmp->[$scol-1];
		}
	}
	close DATAFILE;
	die "ERROR: datafile $datFile is not valid!\n" if (! %DATA);
	return \%DATA;
}

sub init {
	my (%OPTS, $tmp, $count);

	getopt('dac',\%OPTS);
	die "usage: $0 -d datafile -a average -c (column=2)\n" 
		if (! exists($OPTS{d})) or ! exists($OPTS{a});
	($dataFile, $avg, $col) = ($OPTS{d}, $OPTS{a}, $OPTS{c});

	print "Initializing...";
	die "ERROR: Expected integer or decimal for avg. Got \"$avg\"\n"
		if ($avg !~ /^(\d+\.?\d*)$/);

	die "ERROR: datafile \"$dataFile\" is not valid: $!\n"
		if (! -e $dataFile or ! -r $dataFile or ! -T $dataFile);
	$prec = 0;
	if ($avg =~ /\.(\d+)/) {
		$tmp = $1;
		while ($tmp =~ /(\d)/g) {
			$count++;
			$prec = $count if ($1 > 0);
		}
	}
	$col = 2 if (!defined($col) or $col !~ /^\d+$/);
	print "Done\n";
}
