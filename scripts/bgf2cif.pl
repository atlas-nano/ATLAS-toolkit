#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader GetBGFAtoms);
use General qw(FileTester HasCell LoadElements);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);
use BOX qw(GetBox);

sub usage;
sub CreateCIF;
sub numerically {($a<=>$b)}
sub findElement;

my ($bgfFile, $saveName, $selection, $molOpt);
my ($ATOMS, $BONDS, $SELECTIONS, $BGF, $CONS, $tmp, $HEADERS, $ELEMENTS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
&AddMolsToSelection($SELECTIONS, $ATOMS) if ($molOpt);
($BGF, $CONS, $tmp) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $CONS });
print "Done\nCreating BGF file $saveName...";
&CreateCIF($BGF, $CONS, $saveName, $HEADERS);
print "Done\n";

sub CreateCIF {
	my ($atoms, $bonds, $fName, $headers) = @_;
	my ($tmp, $i, $j, $cell_info, $fields, $box, $xyzlabel);

	for $i (keys %{ $atoms }) {
		$atoms->{$i}{ELEMENT} = findElement($ELEMENTS,$atoms->{$i}{FFTYPE});
	}

	@{ $fields } = ("XCOORD","YCOORD","ZCOORD");
	$cell_info = "";
	$xyzlabel = "_atom_site_coord_x\n_atom_site_coord_y\n_atom_site_coord_z";
	if (HasCell($HEADERS)) {
		$box = GetBox($atoms, undef, $headers);
		@{ $fields } = ("FA","FB","FC");
		$xyzlabel = "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z";
		$cell_info=sprintf("%-17s%14.4f\n%-17s%14.4f\n%-17s%14.4f\n%-17s%14.4f\n%-17s%14.4f\n%-17s%14.4f\n",
			"_cell_length_a",$box->{X}{len},"_cell_length_b",$box->{Y}{len},"_cell_length_c",$box->{Z}{len},
			"_cell_angle_alpha",$box->{X}{angle},"_cell_angle_beta",$box->{Y}{angle},"_cell_angle_gamma",$box->{Z}{angle});
	}
	open CIFFILE, "> $fName" or die "ERROR: Cannot write to: $fName: $!\n";
	$tmp = localtime;
	my ($sname) = basename($0);
	print CIFFILE <<DATA;
data_block_1
_audit_creation_date	'$tmp'
_audit_creation_method	'generate by $sname'

$cell_info

_symmetry_space_group_name_H-M   'P 1'
_symmetry_Int_Tables_number       1


loop_
_atom_site_label
$xyzlabel
DATA
	for $i (sort numerically keys %{ $atoms }) {
		printf CIFFILE "%2s",$atoms->{$i}{ELEMENT};
		for $j (@{ $fields }) {
			printf CIFFILE "%11.6f",$atoms->{$i}{$j}
		}
		printf CIFFILE "\n";
	}
	close CIFFILE or die "ERROR: Cannot close $fName: $!\n";
}

sub findElement {
	my ($elements, $sfield) = @_;
	my ($i, $eleNum, $eleHash);

	$sfield =~ /([A-Z]+)/i;
	$sfield = $1;
	for $i (keys %{ $elements }) {
		$eleHash->{lc $elements->{$i}{SYMBOL} } = $i;
	}
	if (exists($eleHash->{lc $sfield })) {
		$eleNum = $eleHash->{lc $sfield};
	} else {
		for $i (keys %{ $elements }) {
			if ($sfield =~ /^$elements->{$i}{SYMBOL}/) {
				$eleNum = $i;
				last;
			}
		}
	}
	die "ERROR: Cannot find element type of fftype $sfield\n"
		if (!defined($eleNum));
	return $elements->{$eleNum}{SYMBOL};
}

sub init {
	my (%OPTS, $atomSel);
	getopt('boms',\%OPTS);
	($bgfFile, $saveName, $atomSel, $molOpt) = ($OPTS{b},$OPTS{s},$OPTS{o}, $OPTS{m});
	&usage if (! defined($bgfFile));
	print "Initializing...";
	FileTester($bgfFile);
	$atomSel = "index>0" if (!defined($atomSel));
	$selection = BuildAtomSelectionString($atomSel);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$//;
		$saveName .= ".cif";
	}
	$molOpt = 0 if (! defined($molOpt) or $molOpt !~ /yes|1/i);
	$molOpt = 1 if ($molOpt =~ /yes|1/i);
	$ELEMENTS = LoadElements();
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -o options -m (entire mol = no)
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  entire mol: select all atoms belonging to selected molecule even if not selected.
		default = no
  atom selection options:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
	Default is "index>0"
DATA
die "\n";

}
