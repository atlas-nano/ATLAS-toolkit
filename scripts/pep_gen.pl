#!/usr/bin/perl -w
use strict;

sub FindFiles();
sub ValidateSequence(@);
sub CreatePepetide(@);
sub GetPName(@);
sub FixMol2Name(@);
sub FixBGFName(@);
sub RecordData(@);
sub GetTerminusInfo(@);
sub CreateCerius2File(@);
sub ExecuteCmd();
sub FixTermini(@);
if (! @ARGV) {
    die "usage: $0 sequence [param file]\n";
}

my ($pepsequence) = $ARGV[0];

my (%PeptideList, @peptides, $i, $parfile, $template); 
my (@filenames, @pepname, $outcmd, $bgf_string, $mol2_string, $fsm_string);
my ($nterminus, $cterminus, $cterm_con, $nterm_con) = (0, 0, 0, 0);


system "clear";
FindFiles();
ValidateSequence($pepsequence);

mkdir "backups", 0777 if (! -d "backups");


print "\n\n---==== Creating BGF FILES ====---\n";
print "----------------------------------\n";

for $i (0..$#peptides) {
    if ($peptides[$i] =~ /\w+/) {
	print "Creating $filenames[$i].bgf..";
	CreatePepetide($peptides[$i], $filenames[$i], $pepname[$i]);

	print "Fixing termini of $filenames[$i]...";
	FixTermini($filenames[$i], $pepname[$i]);
	print "Sucess\n";

	FixBGFName($filenames[$i], $pepname[$i]);
	system "cp $filenames[$i]" . ".bgf backups/$filenames[$i]" . ".bgf";
	
	$outcmd = "Score_single_file.pl $filenames[$i]" . ".bgf min";
	system "$outcmd >> outdata";
	
	system "mv $filenames[$i]" . ".bgf backups/$filenames[$i]" . ".back.bgf";
	system "mv $filenames[$i]" . ".fin.bgf $filenames[$i]" . ".bgf";
	
	print "Creating $filenames[$i]" . ".mol2...";
	system "bgf2mol2  $filenames[$i]" . ".bgf lig >> outdata";
	system "/ul/tpascal/scripts/fixmol2file.pl $filenames[$i]" . ".mol2 >> outdata";
	
	FixMol2Name($filenames[$i], $pepname[$i]);
	print "Sucess\n";
	
	print "Creating $filenames[$i].v400fsm.bgf...";
	system "bgf2bgf4fsm.pl $filenames[$i]" . ".bgf cmp";
	print "Sucess\n\n";
	
	system "rm -f $filenames[$i]*.table dreidii322-mpsim.par";
	$bgf_string .= "$filenames[$i]" . ".bgf ";
	$mol2_string .= "$filenames[$i]" . ".mol2 ";
	$fsm_string .= "$filenames[$i].v400fsm.bgf ";
    }
}

system "cat $bgf_string > Ligand.bgf";
system "cat $mol2_string > Ligand.mol2";
system "cat $fsm_string > Ligand.v400fsm.bgf";

sub FindFiles() {
    my ($tmp);

    $parfile = "/ul/tpascal/ff/dreidii322-quanta";
    $tmp = $ARGV[1];

    if ($tmp) {
	if (-e $tmp) {
	    if ($tmp =~ s/\.par//) {
		$parfile = $tmp;
	    }
	}
    }
    $template = "/ul/tpascal/scripts/pep_template.txt";
    -e $template or die "Cannot locate required template $template: $!\n";
    print "Using param file: $parfile\n";
}

sub ValidateSequence(@) {
    my (@tmplist, @peptidelist, $start_index);
    my ($increment, $curr_frag, $pepline, $counter, $j, $peptidename);
    my ($insequence) = $_[0];

    $PeptideList{"A"} = "ala";
    $PeptideList{"R"} = "arg";
    $PeptideList{"N"} = "asn";
    $PeptideList{"D"} = "asp";
    $PeptideList{"C"} = "cys";
    $PeptideList{"Q"} = "gln";
    $PeptideList{"E"} = "glu";
    $PeptideList{"G"} = "gly";
    $PeptideList{"H"} = "his";
    $PeptideList{"I"} = "ile";
    $PeptideList{"L"} = "leu";
    $PeptideList{"K"} = "lys";
    $PeptideList{"M"} = "met";
    $PeptideList{"F"} = "phe";
    $PeptideList{"P"} = "pro";
    $PeptideList{"S"} = "ser";
    $PeptideList{"T"} = "thr";
    $PeptideList{"W"} = "trp";
    $PeptideList{"Y"} = "tyr";
    $PeptideList{"V"} = "val";

    @tmplist = split //, $insequence;
    for $i (@tmplist) {
	$i = uc($i);
	if ($PeptideList{$i}) {
	    push @peptidelist, uc($PeptideList{$i});
	}
    }

    if ($#peptidelist >1) {
	for $i (0 .. ($#peptidelist -1) ) {
	    for $increment ($i .. $#peptidelist - 1) {
		$curr_frag = "";
		$pepline = "";
		$peptidename = "";
		$counter = $increment;
		while ( ($counter) <= $#peptidelist ) {
		    $j = ($counter - $increment) + $i;
		    $curr_frag .= $peptidelist[$j] . "_";
		    $pepline .= "      Selections/" . $peptidelist[$j] . "_l\n";
		    $pepline .= "%       /biograf/biogv330/peptide_l/";
		    $pepline .= $peptidelist[$j] . "_l.bgf\n";
		    $peptidename .= GetPName($peptidelist[$j]);
		    $counter++;
		}
		if ($pepline =~ /\w+/) {
		    RecordData $curr_frag, $pepline, $peptidename;
#		    while ($curr_frag =~ /^(.*)HIS(.*)$/) {
#			$curr_frag =~ s/HIS/HSD/;
#			$peptidename =~ s/H/H1/;
#			RecordData $curr_frag, " ", $peptidename;

#			$peptidename =~ s/H1/H2/;
#			$curr_frag =~ s/HSD/HSP/;
#			RecordData $curr_frag, " ", $peptidename;
#		    }
		}
	    }
	    $increment = $i;
	}
    } else {
	print "The sequence entered: $pepsequence, is invalid or not long enough\n";
    }
}

sub CreatePepetide(@) {
    my ($intext, $flename, $p_name) = @_;
    my ($inline, @outfile);

    if (system ("cp $template temp.macro")) {
	die "There was an error creating the temporary files.\nEnsure that you have write access to the current directory\n";
    }

    open INFILE, "temp.macro" or die "Cannot open temp.macro: $!\n";
    while (<INFILE>) {
	chomp;

	$inline = $_;
	$inline =~ s/pepname/$p_name/;
	$inline =~ s/filenms/$intext/;
	$inline =~ s/outname/$flename/;

	push @outfile, $inline;
    }
    close INFILE;

    open OUTFILE, "> newtemp.macro" or die "Cannot write too temp.macro: $!\n";
    for $i (@outfile) {
	print OUTFILE "$i\n";
    }
    close OUTFILE;

    if (! system "/biograf/bio_msc batbio $parfile newtemp.macro > out") {
	print "Sucess\n";
	system "rm -f out temp.macro newtemp.macro logfile*";
    } else {
	print "Failure\n";
    }

}

sub GetPName(@) {
    my ($instr) = $_[0];
    my ($pkeys, $returnval);

    if ($instr eq "HSD") {
	$returnval = "H1";
    } elsif ($instr eq "HSP") {
	$returnval = "H2";
    } else {
	for $pkeys (keys %PeptideList) {
	    if (uc($PeptideList{$pkeys}) eq uc($instr)) {
		$returnval = $pkeys;
		last;
	    }
	}
    }
    return $returnval;
}

sub FixMol2Name(@) {

    my ($filename, $ligandcode) = @_;
    my ($inline, @writearray, $i, $wrongname);

    $wrongname = $filename . ".bgf";
    $filename .= ".mol2";

    -e $filename or die "Cannot locate $filename: $!\n";
    open INFILE, $filename or die "Cannot open $filename: $!\n";
    while (<INFILE>) {
	chomp;
	$inline = $_;
	$inline =~ s/$wrongname/$ligandcode/;
	push @writearray, $inline;
    }
    close INFILE;

    open OUTARRAY, "> $filename" or die "Cannot write too $filename: $!\n";
    for $i (@writearray) {
	print OUTARRAY "$i\n";
    }
    close OUTARRAY;
}

sub FixBGFName(@) {
    my ($filename, $ligandcode) = @_;
    my ($inline, @writearray, $i);

    $filename .= ".bgf";

    -e $filename or die "Cannot locate $filename: $!\n";
    open INFILE, $filename or die "Cannot open $filename: $!\n";
    while (<INFILE>) {
	chomp;
	$inline = $_;
	if ($_ =~ /^DESCRP/) {
	    $inline = "DESCRP $ligandcode";
	}
	push @writearray, $inline;
    }
    close INFILE;

    if ($#writearray > 1) {
	open OUTARRAY, "> $filename" or die "Cannot write too $filename: $!\n";
	for $i (@writearray) {
	    print OUTARRAY "$i\n";
	}
	close OUTARRAY;
    } else {
	die "Error writing too $filename\n";
    }

}

sub RecordData(@) {
    my ($curr_frag, $pepline, $peptidename) = @_;
    my ($i, $ispresent);

    chop $curr_frag;
    $ispresent = 0;
    for $i (0 .. $#peptides) {
	if ($peptides[$i] eq $pepline) {
	    $ispresent = 1;
	}
    }

    if (! $ispresent) {
	push @filenames, $curr_frag;
	push @peptides, $pepline;
	if (length($peptidename) < 10) {
	    push @pepname , $peptidename;
	} else {
	    push @pepname, substr ($peptidename, 0, 9) . "_";
	}
    }
}

sub FixTermini(@) {

    my ($infile) = $_[0] . ".bgf";
    my ($pepnm) = $_[1];

    if (GetTerminusInfo($infile)) {
	CreateCerius2File($infile, $pepnm);
	ExecuteCmd();
	system "rm -f tmp_file";
    }
}

sub CreateCerius2File(@) {
    my ($i, $outfile);
    my ($infile) = $_[0];
    my ($pepnm) = $_[1];
    my ($oldname) = $infile;

    $oldname =~ s/\.bgf//g;

    $outfile = $infile;

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
#    print OUTFILE "MODEL/RENAME Model (" . $oldname . ") \"" . $pepnm . "\"\n";
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

sub GetTerminusInfo(@) {
    
    my ($i, $outfile);
    my ($infile) = $_[0];
    $outfile = $infile;

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
	return 1;
    } else {
	return 0;
    }
}
