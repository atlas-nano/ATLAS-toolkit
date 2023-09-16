#!/usr/bin/perl -w
use strict;

if (! @ARGV) {
	die "usage: $0 sequence_file [output fasta file]\n";
}

my ($in_file, $out_file) = @ARGV;
my ($outString);

if (! $out_file) {
	$out_file = $in_file;
}

die "Cannot locate $in_file: $!\n"
	if (! -e $in_file);

open INFILE, $in_file || die "Cannot open $in_file: $!\n";
while (<INFILE>) {
	chomp;
	$outString .= $_ . "\n";
}

close INFILE;

$outString =~ s/\s+//g;

open OUTFILE, "> $out_file" || die "Cannot write to $out_file: $!\n";
print OUTFILE $outString;
close OUTFILE;
