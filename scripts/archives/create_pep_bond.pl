#!/usr/bin/perl -w
use strict;

sub GetTerminusInfo();
sub CreateCerius2File();
sub ExecuteCmd();

if (! @ARGV) {
    die "usage: $0 bgffile [savefile]\n";
}

my ($infile, $outfile) = @ARGV;
my ($nterminus, $cterminus, $cterm_con, $nterm_con) = (0, 0, 0, 0);

if (! $outfile) {
    $outfile = $infile;
}

-e $infile or die "Cannot locate $infile: $!\m";

GetTerminusInfo();
print "Running cerius2 for $infile. Saving as $outfile...";
CreateCerius2File();
ExecuteCmd();
print "Done\n";

sub CreateCerius2File() {
    my ($i);

    open OUTFILE, "> tmp_file" or die "Cannot create tmp_file: $!\n";
    print OUTFILE "FILES/LOAD_FORMAT  BGF\nFILES/LOAD  \"$infile\"\n";
    print OUTFILE "SKETCHER/SKETCHER_DIALOG 1\n";
    print OUTFILE "SKETCHER/TEMPLATE_FILE  \"./Cerius2-Models/templates/organic/acetyl\"\n";
    print OUTFILE "SKETCHER/ADD_TEMPLATE AT   -1.616   -1.738    0.000\n";
    print OUTFILE "SKETCHER/TEMPLATE_FILE  \"./Cerius2-Models/templates/organic/amino\"\n";
    print OUTFILE "SKETCHER/ADD_TEMPLATE AT    5.371    2.939    0.000\n";
    print OUTFILE "SKETCHER/TEMPLATE_FILE  \"./Cerius2-Models/templates/organic/methyl\"\n";
    print OUTFILE "SKETCHER/ADD_TEMPLATE AT    6.396    3.098    0.000\n";
    print OUTFILE "SKETCHER/DELETE_ATOM ATOM Atom(" . $nterm_con . ")\n";
    print OUTFILE "SKETCHER/DELETE_ATOM ATOM Atom(" . $cterminus . ")\n";
    print OUTFILE "SKETCHER/FUSE ATOMS Atom(" . ($cterminus + 14) . ") Atom(";
    print OUTFILE ($cterminus + 8) . ")\n";
    print OUTFILE "SKETCHER/FUSE ATOMS Atom(" . ($cterminus + 10) . ") Atom(";
    print OUTFILE $cterm_con . ")\n";
    print OUTFILE "SKETCHER/FUSE ATOMS Atom(" . ($cterminus + 2) . ") Atom(";
    print OUTFILE $nterminus . ")\n";
    print OUTFILE "SKETCHER/CLEAN  START\n";

    for $i (0 .. 50) {
	print OUTFILE "SKETCHER/CLEAN  CONTINUE\n";
    } 

    print OUTFILE "SKETCHER/CLEAN  END\n";
    print OUTFILE "FORCE-FIELD/LOAD_FORCE_FIELD  \"././Cerius2-Resources/FORCE-FIELD/DREIDING2.21\"\n";
    print OUTFILE "FORCE-FIELD/CALCULATE_TYPING\n";
    print OUTFILE "CHARGE/CALCULATION_METHOD  \"Charge-Equilibration\"\n";
    print OUTFILE "CHARGE/CALCULATE\nFORCE-FIELD/SETUP_EXPRESSION\n";
    print OUTFILE "FILES/SAVE_FORMAT  BGF\n";
    print OUTFILE "FILES/SAVE  \"" . $outfile . "\"\n";

    close OUTFILE;
}

sub ExecuteCmd() {
    if (! open OUTFILE, "cerius2 -n tmp_file |") {
	die "Cannot run cerius2: $!\n";
    }
    while (<OUTFILE>) {
	chomp;
    }
    close OUTFILE;

}

sub GetTerminusInfo() {
    
    open INFILE, $infile or die "Cannot open $infile: $!\n";
    while (<INFILE>) {
	if ($_ =~ /^ATOM\s+(\d+)\s+N\s+\w+\s+1/) {
	    $nterminus = $1;
	}elsif ($_ =~ /^ATOM\s+(\d+)\s+OXT/) {
	    $cterminus = $1;
	}elsif ($cterminus > 0 && $_ =~ /^CONECT\s+$cterminus\s+(\d+)/) {
	    $cterm_con = $1;
	}elsif ($nterminus > 0 && $_ =~ /^CONECT\s+$nterminus\s+\d+\s+\d+\s+(\d+)/) {
	    $nterm_con = $1;
	}
    }
    close INFILE;

    if ($nterminus && $cterminus && $cterm_con && $nterm_con) {
	print "SUCESS: nt: $nterminus ct: $cterminus ct_con: $cterm_con nt_con: $nterm_con\n";
    } else {
	die "ERROR: Unable too parse bgf file\nnt: $nterminus ct: $cterminus ct_con: $cterm_con nt_con: $nterm_con\n";
    }
}
