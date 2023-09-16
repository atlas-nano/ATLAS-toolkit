#!/usr/bin/perl -w
#TODO - add impropers!!
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetPDBFileInfo addHeader createHeaders createBGF);
use General qw(FileTester);
use ManipAtoms qw(CreateBondsByDistance);
use constant PI => 4*atan2(1.0,1.0);

sub init;
sub parseGromacsFF;

my ($ffFile, $saveName, $pdbFile);
my ($ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Parsing PDB file $pdbFile...";
($ATOMS, $BONDS) = GetPDBFileInfo($pdbFile);
&CreateBondsByDistance($ATOMS,$BONDS) if (! defined($BONDS));
print "Done\nParsing GROMACS forcefield file $ffFile...";
&parseGromacsFF($ffFile,$ATOMS);
print "Done\nCreating BGF file $saveName...";
$HEADERS = createHeaders(undef, $saveName);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS,$BONDS,$saveName);
print "Done\n";

sub parseGromacsFF {
	my ($inFile, $atoms) = @_;
	my ($i, $fields, @vals, $counter, $valid, $index);

	$valid = $counter = 0;
	@{ $fields } = ("FFTYPE","RESNUM","RESNAME","ATMNAME","CGNR","CHARGE","MASS");
	open GROMACSFF, $inFile or die "ERROR: Cannot open GROMACS force field $inFile: $!\n";
	while(<GROMACSFF>) {
		chomp;
		if ($_ =~ /^\[\s*(\S+)\s*\]/) {
			$valid = 0;
			$valid = 1 if ($1 eq "atoms");
			next;
		}
		next if (! $valid or $_ !~ /^\s*(\d+)\s*(.*)/);
		$index = $1;
		@vals = split /\s+/,$2;
		next if ($#vals<6);
		for $i (0 .. 6) {
			$atoms->{$index}{$fields->[$i]} = $vals[$i];
		}
		delete $atoms->{$index}{LABEL};
		for $i ("ATMNUM","INDEX") {
			$atoms->{$index}{$i} =~ s/^\s+//;
			$atoms->{$index}{$i} =~ s/\s+$//;
		}
		$counter++
	}
	close GROMACSFF;
	die "ERROR: $inFile is not valid!\n" if (! $counter);
}

sub init {
	my (%OPTS);

	getopt('tps',\%OPTS);
	($ffFile, $saveName, $pdbFile) = ($OPTS{t}, $OPTS{s}, $OPTS{p});
	die "usage: $0 -t gromacs force field (itp) -p pdb file -s (save prefix)\n"
		if (! exists($OPTS{p}) or ! exists($OPTS{t}));

	print "Initializing...";
	FileTester($ffFile);
	FileTester($pdbFile);

	if (! defined($saveName)) {
		$saveName = basename($pdbFile);
		$saveName =~ s/\.\w+$//;
		$saveName .= ".bgf";
	}
	print "Done\n";
}
