#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use File::Basename;
use Packages::General;
use Packages::FileFormats;

sub main;
sub init;
sub fixFields;
sub getBonds;

die "usage: $0 pdb_file|directory babel_cmd [save directory]\n"
    if (! @ARGV or $#ARGV < 1);

my ($pdb, $babelCmd, $save) = @ARGV;

&main;

sub main {
    my ($FLIST, $file);
    my ($ATOMS, $BONDS, $HEADER, $bgfName, $BOX, $TYPES);
    ($FLIST) = &init;
 
    for $file (@{ $FLIST }) {
	$bgfName = $save . "/" . basename($file);
	$bgfName =~ s/\.\w+$/\.bgf/;
	print "Creating BGF file $bgfName..\r";
	$ATOMS = GetPDBFileInfo($file, 1, " ", " ", 1);
	($TYPES, $BONDS) = getBonds($file, $babelCmd);
	fixFields($ATOMS, $TYPES);
	$HEADER = createHeaders($BOX, $bgfName);
	addHeader($ATOMS, $HEADER);
	createBGF($ATOMS, $BONDS, $bgfName);
	$ATOMS = ();
	$BONDS = ();
	$TYPES = ();
	$HEADER = ();
    }

    print "\nDone\n";
}

sub getBonds {
    my ($file, $cmd) = @_;
    my ($ATOMS, $BONDS);

    $cmd .= " -ipdb $file -obgf _tmp.bgf >& junk";
    die "ERROR: Cannot execute babel command $cmd on file $file\n"
	if (system($cmd));
    -e "_tmp.bgf" or die "ERROR: Babel command produced no output!\n";
    ($ATOMS, $BONDS) = GetBGFFileInfo("_tmp.bgf", 0);
    die "ERROR: Babel command $cmd produced invalid input!\n" if (! $ATOMS);
    return ($ATOMS, $BONDS);
}
    
sub fixFields {
    my ($atoms, $ref) = @_;
    my ($i, $atmName);

    for $i (keys %{ $atoms }) {
	$atmName = Trim($atoms->{$i}{"LABEL"});
	$atoms->{$i}{"LABEL"} = "ATOM";
	$atoms->{$i}{"ATMNAME"} = $atmName;
	$atoms->{$i}{"RESNAME"} = $atoms->{$i}{"RES_NAME"};
	$atoms->{$i}{"RESNUM"} = $atoms->{$i}{"RES_ID"};
	if (! exists($atoms->{$i}{"CHARGE"})) {
	    $atoms->{$i}{"CHARGE"} = $ref->{$i}{"CHARGE"};
	}
	$atoms->{$i}{"LONEPAIRS"} = $ref->{$i}{"LONEPAIRS"};
	$atoms->{$i}{"NUMBONDS"} = $ref->{$i}{"NUMBONDS"};
	$atoms->{$i}{"FFTYPE"} = $ref->{$i}{"FFTYPE"};
    }
}

    
sub init {
    my (@FILES);
    -e $babelCmd or die "ERROR: Cannot access $babelCmd: $!\n";
    
    if (-d $pdb) {
	opendir PDB, $pdb or die "ERROR: Cannot access directory $pdb: $!\n";
	@FILES = grep { /\w+\.pdb$/ } map {"$pdb/$_"} readdir PDB;
	closedir PDB or die "ERROR: Cannot close open directory $pdb: $!\n";
    } elsif (-e $pdb) {
	FileTester($pdb);
	push @FILES, $pdb;
    }

    die "ERROR: No valid PDB files found\n" if (! @FILES);

    if (! $save) {
	$save = "./";
    }

    if (! -d $save) {
	die "ERROR: Cannot access directory $save: $!\n";
    }

    return (\@FILES);
}
