#!/usr/bin/perl -w

my $in_fle = $ARGV[0];

if (! $in_fle) {
    die "usage: flip.pl pdbfile\n";
}

-e $in_fle or die "Cannot find $in_fle: $!\n";

open INFILE, $in_fle or die "Cannot open $in_fle: $!\n";
while (<INFILE>) {
    $instring = $_;
    chomp($instring);
    $instring = lc($instring);

    $instring =~ s/at/ta/g;
    #$instring =~ s/ta/at/g;
    $instring =~ s/gc/cg/g;
    #$instring =~ s/cg/gc/g;

    $instring = uc($instring);
    push @outarray, $instring;
}

close INFILE;

open OUTFILE, "> " . $in_fle or die "Cannot write to $in_fle: $!\n";
for ($i=0; $i<=$#outarray; $i++) {
    print OUTFILE "$outarray[$i]\n";
}

close OUTFILE;
