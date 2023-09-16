#!/usr/bin/perl -w

use strict;

if (! @ARGV or $#ARGV < 1) {
    die "usage: $0 pdbfile residue_name [type] [save name]\n";

}
sub ReadPDBFile();
sub WritePDBFile();
sub Numercially;

my ($infile, $name, $mytype, $outfile) = @ARGV;
my (%PFile);

die "Cannot locate $infile: $!\n"
    if (! -e $infile);

$outfile = $infile
    if (! $outfile);

$mytype = "GR"
    if (! $mytype);

if (length($name) > 3) {
    $name = substr($name, 0, 3);
}

$name = uc($name);

ReadPDBFile();
WritePDBFile();

sub ReadPDBFile() {

    my ($rec, $counter);
    open INFILE, $infile || die "Cannot open $infile: $!\n";
    $counter = 1;
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\w+\s+\d+\s+\w+\s+\w+\s+\d+\s+(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*/) {
	    $rec = (
		    {
			"X" => $1,
			"Y" => $2,
			"Z" => $3,
		    }
		    );
	    $PFile{$counter} = $rec;
	    $counter++;
	}
    }
    close INFILE;

    die "Invalid pdbfile $infile\n"
	if ($counter == 1);

}

sub WritePDBFile() {
    my ($my_index, $my_res, $outline, $res_nm, $leap);

    $my_res = 1;

    $leap = "$name = loadpdb " . $outfile . "\n";
    $leap .= "bondbydistance $name\n";
    for $my_index (sort Numerically keys %PFile) {
	$res_nm = "C_";
	if ( ($my_index % 6) == 0) {
	    $res_nm .= "6";
	} else {
	    $res_nm .= ($my_index % 6);
	}	    

	$outline .= sprintf("ATOM   %4d%5s $name%6d", 
			    $my_index, $res_nm, $my_res);
	$outline .= sprintf(" %11.3f%8.3f%8.3f\n",
			    $PFile{$my_index}->{"X"},
			    $PFile{$my_index}->{"Y"},
			    $PFile{$my_index}->{"Z"});

	$leap .= "set $name." . $my_res . "." . $res_nm . " type $mytype\n";
	$leap .= "set $name." . $my_res . "." . $res_nm . " charge 0.0\n";

	if ( ($my_index % 6) == 0) {
#	    $outline .= "TER\n";
	    $my_res++;
	}
    }

    open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
    print OUTFILE $outline;
    close OUTFILE;
    open OUTFILE, "> leapadd" || die "Cannot create temporary files\n";
    print OUTFILE $leap;
    close OUTFILE;
}

sub Numerically {
    ($a <=> $b);
}
