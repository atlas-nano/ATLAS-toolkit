#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(H2others);
use General qw(FileTester);
use Getopt::Std qw(getopt);

sub usage;

my ($h_file);
my ($HMAT, $cell);

$|++;
&init;
print "Getting H matrix from $h_file...";
$HMAT = getHMatrix($h_file);
$cell = H2others($HMAT);
print "$cell\n";

sub getHMatrix {
	my ($inFile) = $_[0];
	my ($hmat, $i, $rec);

	open INFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
	while (<INFILE>) {
		chomp;
		if ($_ =~ /^\s*(-?\d+\.\d*)\s+(-?\d+\.\d*)\s+(-?\d+\.\d*)/) {
			$rec = ([$1,$2,$3]);
			push @{ $hmat }, $rec;
		}
	}
	close INFILE;
	
	die "ERROR: Expected 3 columns while reading $inFile. Got " . scalar(@{ $hmat }) . "\n" if ($#{ $hmat } < 2);
	print "Done\nThe H Matrix:\n";
	for $i (0 .. $#{ $hmat }) {
		printf "%12.7f %12.7f %12.7f\n", $hmat->[$i][0], $hmat->[$i][1], $hmat->[$i][2];
	}

	return $hmat;
}

sub init {
	my (%OPTS, $atomSel);
	getopt('h',\%OPTS);
	$h_file = $OPTS{h};
	&usage if (! defined($h_file));

	print "Initializing...";
	FileTester($h_file);
	print "Done\n";
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -h hmatrix"
Arguments:
  -h hmatrix: the Hmatrix
  h00  h01 h02
  h10  h11 h12
  h20  h21 h22
DATA
die "\n";

}
