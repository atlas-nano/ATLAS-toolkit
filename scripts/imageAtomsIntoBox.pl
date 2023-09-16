#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester CoM);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols GetAtmData SelectAtoms BuildAtomSelectionString 
				  RemoveMolLongBonds ReimageMols AddMolsToSelection);
use BOX qw(GetCellFromHeader GetBox);
use REPLICATE qw(GetBoxDisplacementTensor);

sub usage;
sub getdist;
sub makeAtomsMols;


my ($bgfFile, $saveName, $selection, $molOpt, $reaxOpt);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $MOLS, $SELECT);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
#$BOX = GetCellFromHeader($HEADERS);
$BOX = GetBox($ATOMS, undef, $HEADERS);
&GetBoxDisplacementTensor($BOX);
print "Done\nSelecting relevant atoms...";
$SELECT = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection!\n" if (! keys %{ $SELECT });
&AddMolsToSelection($SELECT, $ATOMS) if ($molOpt); 
$MOLS = GetMols($ATOMS, $BONDS, $SELECT);
$MOLS = &makeAtomsMols($ATOMS) if ($reaxOpt);
print "Done\nReimaging atoms into box...";
&RemoveMolLongBonds($ATOMS, $BONDS, $MOLS, $BOX);
&ReimageMols($ATOMS, $MOLS, $BOX);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub makeAtomsMols {
	my ($atoms) = $_[0];
	my ($i, $mols);

	for $i (keys %{ $atoms }) {
		$mols->{$i}{INDEX} = $i;
		$mols->{$i}{MEMBERS}{$i} = $i;
		$mols->{$i}{MOLSIZE} = 1;
	}

	return $mols;
}

sub init {
	my (%OPTS, $atomSel);
	getopt('bomsr',\%OPTS);
	($bgfFile, $saveName, $atomSel, $molOpt, $reaxOpt) = ($OPTS{b},$OPTS{s}, $OPTS{o}, $OPTS{m}, $OPTS{r});
	
	&usage if (!defined($bgfFile));
	
	print "Initializing...";
	FileTester($bgfFile);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
	$atomSel = "index>0" if (! defined($atomSel));
	$selection = BuildAtomSelectionString($atomSel);
	$molOpt = 0 if (! defined($molOpt) or $molOpt !~ /yes|1/i);
	$molOpt = 1 if ($molOpt =~ /yes|1/i);
	$reaxOpt = 0 if (! defined($reaxOpt) or $reaxOpt !~ /yes|1/i);
	$reaxOpt = 1 if ($reaxOpt =~ /yes|1/i);
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -o (atom selection=all) -m (molopt=no) -r (reaxff_bgf=no)
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  molopt: select entire molecule = no
  reaxff_bgf: BGF file generated from reaxff trajectory = no
  atom selection:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
