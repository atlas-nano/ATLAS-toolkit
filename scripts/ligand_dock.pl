#!/usr/bin/perl -w
use strict;
use File::Basename;

sub GetLigands();
sub GetFiles(@);

if (! @ARGV || $#ARGV < 2) {
    die "usage: $0 protein_prefix ligand_dir save_dir_nm\n";
}

my ($protein_prefix, $ligand_dir, $save_dir_nm) = @ARGV;
my ($curr_file, $protein_files, @Ligands, $protein_name);
my ($counter, $l_counter, $isvalid, $out_cmd);

$protein_name = basename($protein_prefix);

if ($ligand_dir =~ /(.+)\/$/) {
    $ligand_dir = $1;
}

die "Cannot locate $ligand_dir" 
    if (! -d $ligand_dir);

for (".bgf", ".mol2", ".pdb", ".v400fsm.bgf") {
    $curr_file = $protein_prefix . $_;
    die "Cannot locate $curr_file $!\n"
	if (! -e $curr_file );
    $protein_files .= "$curr_file ";
}

$save_dir_nm = "/temp1/" . $ENV{USER} . "/" . $save_dir_nm;
system "mkdir -p $save_dir_nm";
chdir $save_dir_nm;

GetLigands();

for $l_counter (@Ligands) {
    $isvalid = 1;
    GetFiles($l_counter);
    print "LIGAND: $l_counter\n--=== SCANBINDSITE ===--\n";
    $out_cmd = "ScanBindSite.pl -l " . $l_counter . ".mol2 -p ";
    $out_cmd .= $protein_name . " -m pass -a 20 -e 90 >& out";
   
    if (! system $out_cmd) {

	for $counter ("-gridbox.pdb", "-vls.sph", "-vlssph.pdb") {
	    $curr_file = $protein_name . $counter;
		if (! -e $curr_file ) {
		    print "WARNING: Cannot locate required file $curr_file :$! \n";
		    print "Cannot run HierDock for ligand: " . $l_counter . "\n";
		    $isvalid = 0;
		    last;
		} else {
		    system "cp $curr_file ./PTN/";
		}
	}
	    
	if ($isvalid) {
	    print "--=== HIERDOCK ===--\n";
	    if ( system "HierDock $protein_name $l_counter >& HD.out") {
		print "HIERDOCK Error\n";
	    };

	    chdir "../";
	    print "\n";
	}
    } else {
	print "ERROR: ScanBindSite Error\n";
    }
}

sub GetFiles(@) {
    my ($in_file) = $_[0];
    system "mkdir -p $in_file";
    chdir $in_file;
    system "mkdir -p PTN";
    system "mkdir -p MOL2";
    system "mkdir -p BGF";

    system "cp $ligand_dir" . "/" . $in_file . ".mol2 .";
    system "cp $ligand_dir" . "/" . $in_file . ".mol2 ./MOL2";
    system "cp $ligand_dir" . "/" . $in_file . ".v400fsm.bgf ./BGF";
    system "cp $protein_files ./PTN";
    system "cp $protein_files .";
}

sub GetLigands() {

    opendir(LIGS, $ligand_dir) || die "Cannot read from $ligand_dir: $!\n";
    my (@lig_files) = grep { /^.+\.mol2$/ } readdir LIGS;    

    close LIGS;

    for $counter (@lig_files) {
	if ($counter =~ /^(.+)\.mol2$/) {
	    $counter = $1;
	
	    $isvalid = 1;
	    for $l_counter (".bgf", ".v400fsm.bgf") {
		$curr_file = $counter . $l_counter;
		if (! -e ($ligand_dir . "/" . $curr_file) ) {
		    $isvalid = 0;
		    last;
		}
	    }
	    if ($isvalid) {
		push @Ligands, $counter;
	    }
	}
    }
    die "ERROR: Cannot find any valid ligand files\n"
	if ($#Ligands == -1);
}
