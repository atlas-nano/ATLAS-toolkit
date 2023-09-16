#!/usr/bin/perl -w
use strict;

sub ExecuteLeapCmd();
sub RenumberPdbFile();
sub FixLine($);

my ($infile, $outfile, $start_nm, $end_nm) = @ARGV;
my ($outline, $is_multiple, $startfile, $startoutfile);

if (!@ARGV) {
    die "usage: $0 infile [outfile] [start] [end]\n";
}

$is_multiple = 0;
if ($start_nm and $end_nm) {
    if ($start_nm =~ /^\d+$/ and $end_nm =~ /^\d+$/) {
	$is_multiple = 1;
	if ($infile =~ /^(.+)_\d+\.pdb/) {
	    $infile = $1;
	}
	print "Ready to convert " . ($end_nm - $start_nm + 1) . " file(s).\n";
	print "File template: $infile\n";
    } else {
	$start_nm = 1;
	$end_nm =1;
    }
} else {
    $start_nm = 1;
    $end_nm =1;
}

$startfile = $infile;
$startoutfile = $outfile;

while ($start_nm <=$end_nm) {
    if ($is_multiple) {
	$infile = $startfile . "_";
	$infile .= sprintf("%0" . length($end_nm) . "d", $start_nm);
	$infile .= ".pdb";

	$outfile = $startoutfile . "_";
	$outfile .= sprintf("%0" . length($end_nm) . "d", $start_nm);
	$outfile .= ".pdb";
	print "Converting $infile -> $outfile...";
    } else {
	if (! $outfile) {
	    $outfile = $infile; 
	    print "Converting $infile...";
	} else {
	    print "Converting $infile -> $outfile...";
	}
    }
    if (! -e $infile) { 
	print  "ERROR: Cannot locate $infile: $!\n";
    } else {
    
	$outline = "";
    
    
	open INFILE, $infile or die "Cannot open $infile: $!\n";
	while (<INFILE>) {
	    chomp;
	    $outline .= FixLine($_);
	}
	
	close INFILE;
	
	open OUTFILE, "> $outfile" or die "Cannot write too $outfile: $!\n";
	print OUTFILE "$outline";
	close OUTFILE;
	
	ExecuteLeapCmd();
#	RenumberPdbFile();
	print "Done\n";
	
    }
    $start_nm++;
}

sub FixLine($) {
    my ($inline) = $_[0];

    if ($inline =~ /^ATOM/) {
	$inline =~ s/  H2\* / H2\'1 /g;
	$inline =~ s/ (\d)HN(\d) /  H$2$1 /g;
	$inline =~ s/ (\d)H5M /  H7$1 /g;
    }
    $inline .= "\n";
    return $inline;
}

sub ExecuteLeapCmd() {

my ($my_cmd) = "cp /ul/tpascal/scripts/newleaprc leaprc";
system $my_cmd;

open MYLEAPRC, ">> leaprc" or die "Cannot write to leaprc: $!\n";
print MYLEAPRC "model = loadpdb $outfile\n";
print MYLEAPRC "savepdb model $outfile\n";
print MYLEAPRC "quit\n";

close MYLEAPRC;

system "tleap > junk";

#system "rm -f junk leaprc leap.log";

}

sub RenumberPdbFile() {

    my ($outline, $counter, $tmpdata);

    $counter = 1;
    open INFILE, $outfile or die "Cannot open $outfile: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^ATOM\s+\d+(.+)$/) {
	    if ($_ !~ /H5T|H3T/)  {
		$outline .= sprintf("ATOM%7d", $counter);
		$outline .= "$1\n";
		$counter++;
	    }
	} else {
	    $outline .= "$_\n";
	}
    }

    close INFILE;
    
    open OUTFILE, "> $outfile" or die "Cannot write too $outfile: $!\n";
    print OUTFILE "$outline";
    close OUTFILE;
    
}
