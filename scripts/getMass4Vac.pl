#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms AddMass);
use General qw(GetSelections FileTester LoadElements LoadFFs);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetAtmList);
use ManipAtoms qw(FindElement);

sub usage;
sub getMass;
sub printMass;

my ($bgfFile, $saveName, $selection, $FFs);
my ($SELECT, $ATOMS, $BONDS, $MASS, $PARMS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile,0);
print "Done\n";
$PARMS = LoadFFs($FFs);
AddMass($ATOMS, $PARMS);
print "Parsing atom/residue selection...";
$SELECT = GetSelections($selection, 0);
$SELECT = GetAtmList($SELECT, $ATOMS);
die "ERROR: No valid atoms selected!\n" if (! keys %{ $SELECT });
print "Done\nGetting Masses...";
$MASS = getMass($SELECT, $ATOMS, $PARMS);
print "Done\nSaving masses to $saveName...";
&writeMass($MASS, $saveName);
print "Done\n";

sub writeMass {
    my ($mass, $savefile) = @_;

    open OUTDATA, "> $savefile" or die "ERROR: Cannot write to $savefile: $!\n";
    print OUTDATA $mass;
    close OUTDATA;
}

sub getMass {
    my ($select, $atoms, $parms) = @_;
    my ($masses, $i, $ffType, $eleNum, $atomMass, $ELEMENTS, $eleList, $USED);
    
    $ELEMENTS = &LoadElements;

    for $i (keys %{ $select }) {
	$ffType = $atoms->{$i}{FFTYPE};
	next if (exists($USED->{$ffType}));
	$USED->{$ffType} = 1;
	$atomMass = $atoms->{$i}{MASS};
	(undef, $eleNum) = FindElement($ffType, $parms->{ATOMTYPES}, $ELEMENTS);
	$masses .= "{\"$ffType\", $eleNum, $atomMass },\n";
    }
    return $masses;
}
 
sub init {
    my (%OPTS, $select, $fflist);
    getopt('bosf',\%OPTS);
    ($bgfFile, $saveName, $select, $fflist) = ($OPTS{b},$OPTS{s},$OPTS{o},$OPTS{f});
    for ($bgfFile, $fflist) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);

    while ($fflist =~ /(\S+)/g) {
	if (-e $1 and -r $1 and -T $1) {
	    push @{ $FFs }, $1;
	}
    }
    die "ERROR: No valid forcefield files found while search \"$fflist\"\n" if (! $FFs);
 
    $select = "*" if (! defined($select));
    if ($select =~ /\s+/) {
	@{ $selection } = split /\s+/, $select;
    } else {
	$selection->[0] = $select;
    }
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -f forcefield(s) -s save_name -o options
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  options:
    [^][:][I|T|N][a|r]
    a   - atom
    r   - residue
    IaX - atom number X
    IrX - residue index X
    TaX - atom type X
    NaX - atom name X
    NrX - residue name X
    Use ":" to specify a range, eg. :Tr1-8 :Ia3-66
    Use "^" to exclude a selection. You can use multiple combinations
    range and exclusion enclosed in quotes, eg, "^:TrIP-IM ^:Ia23-45"
    to exclude residues of type IM and IP and atoms 23-45
DATA

die "\n";

}
