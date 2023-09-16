#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester);
use FileFormats qw(createBGF addHeader GetBGFFileInfo);
use ManipAtoms qw(GetMols);

sub init;
sub updateAtomTypes;
sub updateMolData;
sub updateResData;
sub updateFieldData;
sub numerically { $a<=>$b; }

my ($resFile, $saveFile, $bgfFile, $FIELDS, $molOpt);
my ($ATOMS, $BONDS, $MOLS, $HEADERS, $resATOMS, $resBONDS, $resMOL);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\nParsing residue library file $resFile...";
($resATOMS, $resBONDS, undef) = GetBGFFileInfo($resFile);
$resMOL = GetMols($resATOMS, $resBONDS);
$resATOMS = updateResData($resATOMS) if (!$molOpt);
&updateFieldData($ATOMS, $FIELDS);
print "Done\nUpdating atom info based on library...";
&updateAtomTypes($ATOMS, $resATOMS, $FIELDS) if (!$molOpt);
&updateMolData($ATOMS, $MOLS, $resATOMS, $resMOL, $FIELDS) if ($molOpt);
print "Done\nCreating BGF file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateMolData {
	my ($atoms, $mols, $refAtoms, $refMol, $fields) = @_;
	my ($i, $j, $f, $refMolSize, $molSize, $refMolAtoms, $molAtoms, $refAtom, $molAtom);
	
	$refMolSize = $refMol->{1}{MOLSIZE};
	@{ $refMolAtoms } = sort numerically keys %{ $refMol->{1}{MEMBERS} };
	for $i (keys %{ $mols }) {
		$molSize = $mols->{$i}{MOLSIZE};
		@{ $molAtoms } = sort numerically keys %{ $mols->{$i}{MEMBERS} };
		next if ($molSize != $refMolSize);
		for $j (0 .. ($refMolSize-1)) {
			$refAtom = $refMolAtoms->[$j];
			$molAtom = $molAtoms->[$j];
			for $f (keys %{ $fields }) {
				$atoms->{$molAtom}{$f} = $refAtoms->{$refAtom}{$f};
			}
		}
	}
}

sub updateFieldData {
	my ($atoms, $fields) = @_;
	my ($i, @tmp, $fatom, $fList);

	@tmp = keys %{ $atoms };
	$fatom = $atoms->{ $tmp[0] };

	for $i (keys %{ $fields }) {
		$fList .= "$i ";
		delete $fields->{$i} if (! exists($fatom->{$i}));
	}

	die "ERROR: field \"$fList\" is invalid!\n" if (scalar(keys %{ $fields }) == 0);
}

sub updateAtomTypes {
	my ($atoms, $resInfo, $fields) = @_;
	my ($i, $j, $resName, $atmName);

	for $i (keys %{ $atoms }) {
		($resName, $atmName) = ($atoms->{$i}{RESNAME}, $atoms->{$i}{ATMNAME});
		$resName =~ s/\s//g;
		$atmName =~ s/\s//g;
		if (exists($resInfo->{$resName}) and exists($resInfo->{$resName}{$atmName})) {
			for $j (keys %{ $fields }) {
				$atoms->{$i}{$j} = $resInfo->{$resName}{$atmName}{$j};
			}
		}
	}
}

sub updateResData {
	my ($resInfo) = $_[0];
	my (%RES, $i, $resName, $atmName);

	for $i (keys %{ $resInfo }) {
		$resName = $resInfo->{$i}{RESNAME};
		$atmName = $resInfo->{$i}{ATMNAME};
		$resName =~ s/\s//g;
		$atmName =~ s/\s//g;
		%{ $RES{$resName}{$atmName} } = %{ $resInfo->{$i} };
	}
	
	return \%RES;
}

sub init {
	my (%OPTS, $fList);
	getopt('brsfm',\%OPTS);
	
	for ("b", "r") {
		die "usage: $0 -b bgf structure file -r residue library (bgf) -s (save name (optional) -f (fields = fftype) -m (use mols)\n"
			if (! exists($OPTS{$_}));
	}

	($bgfFile, $saveFile, $resFile, $fList, $molOpt) = ($OPTS{b}, $OPTS{s}, $OPTS{r}, $OPTS{f}, $OPTS{m});
	print "Initializing...";
	FileTester($bgfFile);
	FileTester($resFile);
	if (! defined($saveFile)) {
		$saveFile = basename($bgfFile);
		$saveFile =~ s/\.\w+$//;
		$saveFile .= "_typed.bgf";
	}

	$fList = "FFTYPE" if (! defined($fList));
	while ($fList =~ /(\S+)/g) {
		$FIELDS->{uc($1)} = 1;
	}
	$molOpt = 0 if (! defined($molOpt) or $molOpt !~ /1|yes/i);
	$molOpt = 1 if ($molOpt =~ /1|yes/i);
	print "Done\n";
}
