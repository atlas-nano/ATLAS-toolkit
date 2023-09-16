#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms insertHeaderRemark);
use General qw(FileTester HasCell LoadElements AddElementField);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);
use BOX qw(GetBox);

sub usage;

my ($bgfFile, $saveName, $selection, $OM);
my ($ATOMS, $BONDS, $SELECTIONS, $sMOLS, $HEADERS, $ELEMENTS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&AddElementField($ATOMS, $ELEMENTS);
$BOX = &GetBox($ATOMS, undef, $HEADERS) if (HasCell($HEADERS));
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
&AddMolsToSelection($SELECTIONS, $ATOMS);
$sMOLS = &GetMols($ATOMS, $BONDS, $SELECTIONS);
print "Done\nAdding OM particle to selected molecules using O-M length " . ${OM}->{dist} . " and H-O-M angle " . ${OM}->{angle} . " ...";
&addOMParticles($ATOMS, $BONDS, $sMOLS, $OM);
print "Done\nCreating BGF file $saveName...";
&insertHeaderRemark($HEADERS, "REMARK " . localtime() . " added OM particle at OM-dist " . ${OM}->{dist} . " and angle " . ${OM}->{angle});
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub addOMParticles {
	my ($atoms, $bonds, $solvM, $om_opts) = @_;
	my ($i);

	for $i (values %{ $solvM }) {
		next if(! isWat($atoms, $i));
		print "";
	}
}

sub isWat {
	my ($atoms, $molecule) = @_;
	my ($i, $eleStr);

	for $i (keys %{ $molecule->{MEMBERS} }) {
		$eleStr .= $atoms->{$i}{ELEMENT}{SYMBOL};
	}
	return 0 if ($eleStr !~ /^(HHO|HOH|OHH)$/);
	return 1;
}

sub init {
	my (%OPTS, $om_opt_str, $atomSel);
	getopt('bosa',\%OPTS);
	($bgfFile, $saveName, $om_opt_str, $atomSel) = ($OPTS{b},$OPTS{s},$OPTS{o},$OPTS{a});
	for ($bgfFile, $om_opt_str) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);

	$atomSel = "resname eq 'WAT'" if (! defined($atomSel));
	$selection = BuildAtomSelectionString($atomSel);

	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/\.om\.bgf/;
	}
	if ($om_opt_str =~ /(\d+\.?\d*)\s+(\d+\.?\d*)/) {
		$OM->{dist} = $1;
		$OM->{angle} = $2;
	} elsif ($om_opt_str =~ /SWM4-DP/i) {
		$OM->{dist} = 0.23808;
		$OM->{angle} = 52.26;
	} elsif ($om_opt_str =~ /SWM4-NDP/i) {
		$OM->{dist} = 0.24034;
		$OM->{angle} = 52.26;
	} else {
		die "ERROR: Expected 'om_dist om_angle' or swm4-dp or swm4-ndp for om_particle_options. Got '$om_opt_str'. Aborting...";
	}
	$ELEMENTS = LoadElements();
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -o om_particle_options -a (atom_selection) -s (save_name) 
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  om_particle_options:
	expected 'om_dist hom_angle' or swm4-dp or swm4-ndp
  atom_selection:
  	any valid bgf field expression. default is "resname eq 'WAT'"
DATA
die "\n";

}
