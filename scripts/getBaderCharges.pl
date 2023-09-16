#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester LoadElements AddElementField);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox);

sub usage;

my ($bgfFile, $saveName, $baderQfile, $eleMods);
my ($ATOMS, $BONDS, $HEADERS, $ELEMENTS, $bCharges);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&AddElementField($ATOMS, $ELEMENTS, undef);
print "Done\nParsing Bader charge file ${baderQfile}...";
$bCharges = parseBaderQfile($baderQfile);
print "Done\nUpdating charges...";
&updateCharges($ATOMS, $bCharges, $eleMods);
print "Done\nWriting update BGF file ${saveName}...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub updateCharges {
	my ($atoms, $bQs, $eMods) = @_;
	my ($i, $nval);

	die "ERROR: Number of atoms in bgf file (" . scalar(keys %{ $atoms }) . 
	") is not equal to number of charges in bader file (" . scalar(keys %{ $bQs }) . 
		")!\n" if(scalar(keys %{ $atoms }) != scalar(keys %{ $bQs }));
	for $i (keys %{ $atoms }) {
		$nval = $atoms->{$i}{ELEMENT}{NVALENCE};
		$nval = $eMods->{ $atoms->{$i}{FFTYPE} } if(defined($eMods) and exists($eMods->{ $atoms->{$i}{FFTYPE} }));
		$atoms->{$i}{CHARGE} = $nval - $bQs->{$i};
	}
}

sub parseBaderQfile {
	my ($inFile) = $_[0];
	my ($charges, $i, $qfield, $tmp);

	$qfield = 5; #volume

	open BADERQFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
	while(<BADERQFILE>) {
		chomp;
		if ($_ =~ /ATOMIC/) {
			$qfield = 0; #atomic
		} elsif ($_ =~ /^\s*(\d+.*)/) {
			@{ $tmp } = split /\s+/,$1;
			next if ($#{ $tmp } != 6);
			$charges->{ $tmp->[$qfield] } = $tmp->[4];
		}
	}
	close BADERQFILE;
	die "ERROR: No valid info found in $inFile...Cannot continue\n" if (! defined($charges));
	return $charges;
}

sub init {
	my (%OPTS, $atomSel, $mOpts);
	getopt('bqsm',\%OPTS);
	($bgfFile, $saveName, $baderQfile, $mOpts) = ($OPTS{b},$OPTS{s},$OPTS{q},$OPTS{m});
	for ($bgfFile, $baderQfile) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	FileTester($baderQfile);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/.baderQ\.bgf/;
	}
	$ELEMENTS = LoadElements();
	if(defined($mOpts)) {
		while($mOpts =~ /(\S+):(\S+)/g) {
			$eleMods->{$1} = $2;
		}
	}
	print "Done\n";
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -q bader_charge_file -s (save_bgf_name) 
Arguments:
  bgf_file: name of bgf_file
  bader_charge_file: location of bader charge file. Either the atomic or volume charges
  save_bgf_name: (Optional) name of file to save
DATA
die "\n";

}
