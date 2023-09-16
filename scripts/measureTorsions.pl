#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo);
use General qw(FileTester GetTorsion);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString GetAtmData);

sub usage;
sub measureTorsions;
sub writeData;
sub numerically { ($a<=>$b); }
sub getCons;

my ($bgfFile, $saveName, $selection);
my ($ATOMS, $BONDS, $SELECTIONS, $DATA);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, undef) = GetBGFFileInfo($bgfFile,1);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTIONS });
print "Done\nCalculating torsions...";
$DATA = measureTorsions($ATOMS, $BONDS, $SELECTIONS);
print "Done\nSaving data to $saveName...";
&writeData($DATA, $saveName);
print "Done\n";

sub writeData {
    my ($data, $save) = @_;
    my ($i);

    open OUTDATA, "> $save" or die "ERROR: Cannot write to $save: $!\n";
    print OUTDATA "#atom1 atom2 atom3 atom4 torsion\n";
    for $i (@{ $data }) {
	printf OUTDATA "$i->{KEY} %.2f\n",$i->{VAL};
    }
    close OUTDATA or die "ERROR: Cannot finalize $save: $!\n";
}
    
sub measureTorsions {
    my ($atoms, $bonds, $list) = @_;
    my ($i, $j, $k, $cons, $atomList, $data, $rec, $tlist, $used, $tmp, $adata);

    $used = ();
    $cons = ();
    for $i (sort numerically keys %{ $list }) {
	$atomList = ();
	$tmp = ();
	$tmp->{$i}  = 1;
	$used->{$i} = 1;
	&getCons($i, \%{ $cons }, 3, $atoms, $bonds, $list, \%{ $used });
	delete $cons->{$i} if (! keys %{ $cons->{$i} });
     }
     die "ERROR: No valid torsions found for atoms selection\n" if (! keys %{ $cons }); 
     $tlist = getAtomList($cons);
     for $j (@{ $tlist }) {
	for $k (0 .. $#{ $j }) {
	    $tmp = ();
	    $tmp->{$j->[$k]} = 1;
	    $adata = GetAtmData($atoms, $tmp);
	    %{ $atomList->[$k] } = %{ $adata->{$j->[$k]} };
	}
	$rec = ();
	$rec->{VAL} = GetTorsion(@{ $atomList });
	$rec->{KEY} = "$i @{ $j }";
	push @{ $data }, $rec;
    }
    return $data;
}

sub getAtomList {
    my ($cons) = $_[0];
    my ($i, $j, $k, $l, $alist, $rec);

    for $i (keys %{ $cons }) {
	for $j (keys %{ $cons->{$i} }) {
	    for $k (keys %{ $cons->{$i}{$j} }) {
		for $l (keys %{ $cons->{$i}{$j}{$k} }) {
		   $rec = ();
		    @{ $rec } = ($i,$j,$k,$l);
		    push @{ $alist }, $rec;
		}
	    }
	}
    }

    return $alist;
}

sub getCons {
    my ($curr, $cons, $depth, $atoms, $bonds, $list, $used) = @_;
    my ($i);

    $depth--;
    for $i (@{ $bonds->{$curr} }) {
	if (exists($list->{$i}) and ! exists($used->{$i})) {
	    $used->{$i} = 1;
	    if ($depth > 0) {
		&getCons($i, \%{ $cons->{$curr} }, $depth, $atoms, $bonds, $list, $used);
	    } else {
		$cons->{$curr}{$i} = 1;
	    }
	}
    }
}

sub init {
    my (%OPTS, $atomSel);
    getopt('bas',\%OPTS);
    ($bgfFile, $saveName, $atomSel) = ($OPTS{b},$OPTS{s},$OPTS{a});
    for ($bgfFile, $atomSel) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    $selection = BuildAtomSelectionString($atomSel);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\w+$//;
	$saveName .= "torsions.dat";
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -a atom_selection
Arguments:
  bgf_file: name of bgf_file
  save_name: name of data file to save
  atom selection options:
    any valid bgf field expression. E.g. resname eq 'WAT' will select
    all the "WAT" residues while index > 10 will select all indices > 10.
    combine multiple expressions to make complicated selections: e.g.
    (xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
