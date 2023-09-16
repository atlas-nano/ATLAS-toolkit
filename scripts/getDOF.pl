#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo);
use General qw(FileTester);
use File::Basename;
use Getopt::Std qw(getopt);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString GetMols);

my ($bgfFile, $selection, $verbose);
my ($BGF, $ATOMS, $BONDS, $dof, $tot, $SELECT);

sub init;
sub getShakeDOF;

$|++;
&init;
print "Parsing BGF file $bgfFile..." if ($verbose);
($BGF, $BONDS) = GetBGFFileInfo($bgfFile, 0, $verbose);
&GetMols($BGF, $BONDS);
print "Done\nGetting Constrained DOFs..." if ($verbose);
$SELECT = SelectAtoms($selection, $BGF);
$dof = getShakeDOF($SELECT, $BGF);
$tot = (scalar(keys %{ $SELECT }) * 3);
print "Done\nTotal DOF: $tot Constrained DOF: $dof Real DOF: " . ($tot - $dof) . "\n" if ($verbose);
print "$dof\n" if (! $verbose);

sub getShakeDOF {
    my ($select, $atoms) = @_;
    my ($i, $dof);
    $dof = 0;
    for $i (keys %{ $select }) {
	if ($atoms->{$i}{FFTYPE} =~ /^H/i) { #if this is a hydrogen
	    $dof++;
	    if (${ $atoms->{$i}{MOLSIZE} } == 3) {
		$dof += 0.5;
	    }
	}
    }
    return $dof;
}

sub init {
    my (%OPTS, $atomSelect);
    
    getopt('bav',\%OPTS);
    ($bgfFile, $atomSelect, $verbose) = ($OPTS{b}, $OPTS{a}, $OPTS{v});
    die "usage:$0 -b bgf file -a (atom selection) -v (verbose = yes)\n"
	if (! defined($bgfFile));

    $verbose = 1 if (! defined($verbose) or $verbose !~ /^\s*(0|no|n)\s*$/);
    $verbose = 0 if ($verbose =~ /^\s*(0|no|n)\s*$/);

    print "Initializing..." if ($verbose);
    FileTester($bgfFile);
    $atomSelect = "index > 0" if (! defined($atomSelect));
    $selection = BuildAtomSelectionString($atomSelect);
    print "Done\n" if ($verbose);
}
