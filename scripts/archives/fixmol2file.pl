#!/usr/bin/perl -w
use strict;

if (! @ARGV) {
    die "usage: $0 input_mol_file [output_file]\n";
}

my ($infile, $outfile) = @ARGV;
my ($i, $inline, @outarray);

if (! $outfile) {
    $outfile = $infile;
}

-e $infile or die "Cannot locate $infile: $!\n";

print "Reading $infile....";

open INFILE, $infile or die "Cannot open $infile: $!\n";
while (<INFILE>) {
    chomp;
    $inline = $_;
    if ($_ =~ /^(\s+\d+\s+\w+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+)N\.3(\s+\d+\s+)(\w+)(.+)/) {
	if (! ($3 =~ /LYS/i) ) {
	    $inline = $1 . "N.4" . $2 . $3 . $4;
	}
    }
    push @outarray, $inline;
}

print "Sucess\n";
close INFILE;

print "Writing changes to $outfile....";
open OUTFILE, "> $outfile" or die "Cannot write to $outfile: $!\n";
for $i (@outarray) {
    print OUTFILE "$i\n";
}

close OUTFILE;
print "Sucess\n";
