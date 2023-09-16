#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(createBGF GetBGFFileInfo createMOL2 GetMOL2FileInfo
			     addHeader);
use Packages::General qw(FileTester);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub convertFile;
sub getAmberAtomTypes;
sub updateAtomTypes;

my ($inFile, $outFile, $fileType);
my ($parseFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS, $atomTypes);
$|++;
&init;
print "Parsing $fileType file $inFile...";
($ATOMS, $BONDS, $HEADERS) = $parseFile->($inFile, 1);
print "Done\n"; 
if ($fileType eq "BGF") {
    print "Converting BGF -> MOL2...";
    $inFile = &convertFile($inFile);
    print "Done\n";
} else {
    delete $ATOMS->{HEADER};
}
print "Typing with AMBER FF atom types...";
$atomTypes = &getAmberAtomTypes($inFile);
print "Done\nUpdating atom types...";
&updateAtomTypes($ATOMS, $atomTypes);
print "Done\nCreating $fileType file $outFile...";
&addHeader($ATOMS, $HEADERS) if ($fileType eq "BGF");
$saveFile->($ATOMS, $BONDS, $outFile);
print "Done\n";

sub updateAtomTypes {
    my ($atoms, $types) = @_;
    my ($i);
    
    for $i (keys %{ $atoms }) {
	die "ERROR: No atom type found for atom $i!\n" if (! defined($types->{$i}));
	$atoms->{$i}{FFTYPE} = $types->{$i};
    }
}

sub getAmberAtomTypes {
    my ($inFile) = $_[0];
    my ($acFile, $typeCmd, %ATOMTYPES);

    $acFile = $inFile;
    $acFile =~ s/\.\w+$//g;
    $acFile .= ".ac";
    $typeCmd = $ENV{AMBERHOME} . "/exe/" if ($ENV{AMBERHOME});
    $typeCmd .= "atomtype -i $inFile -f mol2 -o $acFile -p amber >& _ambertype";
    die "Error while executing $typeCmd! See ATOMTYPE.INF\n" if (system($typeCmd));

    open ACFILE, $acFile or die "ERROR: Cannot open $acFile: $!\n";
    while (<ACFILE>) {
	chomp;
	if ($_ =~ /^ATOM\s*(\d+).*\s+(\S+)\s*$/) {
	    $ATOMTYPES{$1} = $2;
	}
    }
    close ACFILE;
    die "ERROR: $acFile is invalid!\n" if (! %ATOMTYPES);
    system("rm -fr _ambertype ATOMTYPE.INF $acFile");
    return \%ATOMTYPES;
}

sub convertFile {
    my ($bgfFile) = $_[0];
    my ($mol2File, $cnvCmd);

    $mol2File = $bgfFile;
    $mol2File =~ s/\.\w+$//;
    $mol2File .= ".mol2";
    $cnvCmd = "/ul/tpascal/scripts/bgf2mol2.pl -b $bgfFile -m $mol2File >& _bgf2mol2cnv";
    die "Error while executing $cnvCmd. See _bgf2mol2cnv file!\n" if (system($cnvCmd));
    system ("rm -fr _bgf2mol2cnv");
    return $mol2File;
}
       
sub init {
    my (%OPTS);
    getopt('ist',\%OPTS);
    
    for ("i", "t") {
	die "usage: $0 -i input file -t file type (mol2|bgf) -s (save name - optional)\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($inFile, $fileType, $outFile) = ($OPTS{i}, $OPTS{t}, $OPTS{s});
    FileTester($inFile);
    die "ERROR: Expected \"mol2\" or \"bgf\" for file type. Got $fileType!\n"
	if ($fileType !~ /^(mol2|bgf)/i);
    $fileType = uc($fileType);
    if ($fileType eq "MOL2") {
	$parseFile = \&GetMOL2FileInfo;
	$saveFile = \&createMOL2;
    } else {
	$parseFile = \&GetBGFFileInfo;
	$saveFile = \&createBGF;
    }
    
    if (! defined($outFile)) {
	$outFile = $inFile;
	$outFile =~ s/\.\w+$//g;
	$outFile .= "_amber." . lc($fileType);
    }
    print "Done\n";
}
