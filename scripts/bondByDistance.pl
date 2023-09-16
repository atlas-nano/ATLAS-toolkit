#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester GetBondLength LoadFFs LoadElements AddElementField ReadFFs);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString CreateBondsByDistance);
use BOX qw(GetBox);
use CERIUS2 qw(GenerateUFFParms);

sub usage;

my ($bgfFile, $saveName, $selection, $max_bonds, $ffList);
my ($ATOMS, $BONDS, $SELECTIONS, $HEADERS, $BOX, $FF, $FFILES, $ELEMENTS, $tBonds);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$BOX = GetBox($ATOMS, undef, $HEADERS);
&AddElementField($ATOMS, $ELEMENTS);
print "Done\n";
($FFILES, undef) = ReadFFs($ffList);
$FF = LoadFFs($FFILES);
&saveUFFParms($ATOMS, $BONDS, $FF, \%{ $tBonds }) if (exists($FF->{PARMS}{UFF}) and keys %{ $FF->{PARMS}{UFF} });
&GenerateUFFParms($ATOMS, $BONDS, $FF) if (exists($FF->{PARMS}{UFF}) and keys %{ $FF->{PARMS}{UFF} });
&removeTbonds($BONDS,$tBonds) if (defined($tBonds));
print "Selecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTIONS });
&checkBonds($ATOMS, $SELECTIONS, $BONDS, $FF);
print "Done\nCreating bonds...";
&CreateBondsByDistance($ATOMS,$BONDS,$BOX,$SELECTIONS,$FF,$max_bonds);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub removeTbonds {
	my ($bonds, $tmp) = @_;
	my ($i, $j, $k);

	for $i (keys %{ $tmp }) {
		for $j (keys %{ $tmp->{$i} }) {
			$k = 0;
			while ($k <= $#{ $bonds->{$i} }) {
				if ($bonds->{$i}[$k] == $j) {
					splice @{ $bonds->{$i} }, $k, 1;
				} else {
					$k++;
				}
			}
		}
	}
}

sub saveUFFParms {
	my ($atoms, $bonds, $ff, $tmpBonds) = @_;
	my ($i, $j, $a1, $a2, $mass, $fftype, $ffMass, $valid, $alist);

	for $i (keys %{ $atoms }) {
		$fftype =  $atoms->{$i}{FFTYPE};
		if (exists($ff->{ATOMTYPES}{$fftype})) {
			if(!exists($alist->{$fftype}{1})) {
				$alist->{$fftype}{1}=  $i;
			} else {
				$alist->{$fftype}{2} = $i;
			}
			$ff->{ATOMTYPES}{$fftype}{USED} = 1;
			$ff->{ATOMTYPES}{$fftype}{ELENUM} = $atoms->{$i}{ELEMENT}{NUM};
			next;
		}
		$mass = $atoms->{$i}{ELEMENT}{MASS};
		$valid = 0;
		for $j (keys %{ $ff->{ATOMTYPES} }) {
			$ffMass = $ff->{ATOMTYPES}{$j}{MASS};
			if(($mass-0.1)<$ffMass and $ffMass<($mass+.1)) {
				%{ $ff->{ATOMTYPES}{$fftype} } = %{$ff->{ATOMTYPES}{$j} };
				$ff->{ATOMTYPES}{$fftype}{USED} = 1;
				$ff->{ATOMTYPES}{$fftype}{ELENUM} = $atoms->{$i}{ELEMENT}{NUM};
				%{ $ff->{PARMS}{UFF}{$fftype} } = %{$ff->{PARMS}{UFF}{$j} };
				$valid = 1;
				last;
			}
		}
		die "ERROR: Cannot locate any forcefield parameter for atom $i ($fftype)!\n"
			if (! $valid);
		$alist->{$fftype}{1} = $i;
	}	

	#create temporary bonds for UFF generation
	for $i (keys %{ $alist }) {
		$a1 = $alist->{$i}{1};
		for $j (keys %{ $alist }) {
			$a2 = $alist->{$j}{1};
			$a2 = $alist->{$j}{2} if ($i eq $j);
			next if (bonded($bonds, $a1,$a2));
			push @{ $bonds->{$a1} }, $a2;
			push @{ $bonds->{$a2} }, $a1;
			$tmpBonds->{$a1}{$a2} = 1;
			$tmpBonds->{$a2}{$a1} = 1;
		}
	}
}

sub bonded {
	my ($bonds, $a1, $a2) = @_;
	my ($i, $isbond);

	$isbond = 0;
	for $i (@{ $bonds->{$a1} }) {
		if($i == $a2) {
			$isbond = 1;
			last;
		}
	}
	return $isbond if ($isbond);
	for $i (@{ $bonds->{$a2} }) {
		if ($i == $a1) {
			$isbond = 1;
			last
		}
	}
	return $isbond;
}

sub checkBonds {
	my ($atoms, $select, $bonds, $ff) = @_;
	my ($i, $j, $ffi, $ffj, $checked);

	for $i (keys %{ $select }) {
		$ffi = $atoms->{$i}{FFTYPE};
		for $j (@{ $bonds->{$i} }) {
			next if (! exists($select->{$j}));
			$ffj = $atoms->{$j}{FFTYPE};
			next if ((exists($checked->{$i}) and exists($checked->{$i}{$j})) or 
				     (exists($checked->{$j}) and exists($checked->{$j}{$i})));
			if ((! exists($ff->{BONDS}{$ffi}) or ! exists($ff->{BONDS}{$ffi}{$ffj})) and
				(! exists($ff->{BONDS}{$ffj}) or ! exists($ff->{BONDS}{$ffj}{$ffi}))) {
				die "ERROR: No valid ff entry for bond $ffi - $ffj found!\n";
			}
			$checked->{$ffi}{$ffj} = 1;
		}
	}
}

sub init {
	my (%OPTS, $atomSel);
	getopt('bsafm',\%OPTS);
	($bgfFile, $saveName, $atomSel, $ffList, $max_bonds) = ($OPTS{b},$OPTS{s},$OPTS{a}, $OPTS{f}, $OPTS{m});
	&usage if (! defined($bgfFile));

	print "Initializing...";
	FileTester($bgfFile);
	$atomSel = "index>0" if (!defined($atomSel));
	$selection = BuildAtomSelectionString($atomSel);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
	$ffList="UFF" if (!defined($ffList));
	$max_bonds = 99999 if(!defined($max_bonds) or (defined($max_bonds) and $max_bonds !~ /\d+/));
	$ELEMENTS = &LoadElements;
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s (save_name) -a (atom selections) -f (forcefield(s)) -m (max_bonds)
Arguments:
	bgf_file: name of bgf_file
	save_name: name of file to save
	forcefield: 1 or more Cerius2/MPSIM forcefields for bond information. Default is the UFF forcefield
	max_bonds: maximum number of bonds. Default is unlimited
	atom selection options:
		any valid bgf field expression. E.g. resname eq 'WAT' will select
		all the "WAT" residues while index > 10 will select all indices > 10.
		combine multiple expressions to make complicated selections: e.g.
		(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
