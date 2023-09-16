#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo);
use AMBER qw(CreateAmberTrj);
use BOX qw(GetBox);

sub main;
sub init;
sub numerically;

$|++;

die "usage: $0 bgf_file|directory amber_traj_name\n"
    if (! @ARGV or $#ARGV < 1);
my ($bgfLoc, $saveName) = @ARGV;
my ($scaled) = 0;

&main;

sub main {
    my ($FLIST, $file);
    my ($ATOMS, $BONDS, $BOX, $HEADERS);

    $FLIST = &init;
    print "Creating AMBER trajectory file $saveName\n";
    open OUTFILE, "> $saveName" or die "ERROR: Cannot create Amber Trj file $saveName: $!\n";
    print OUTFILE "TITLE: AMBER Trajectory file created by bgf2ambertrj.pl on " .
	scalar(localtime(time())) . "\n";
    for $file (@{ $FLIST }) {
	print "Reading BGF file $file...";
	($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($file, 1);
	$BOX = GetBox($ATOMS, undef, $HEADERS);
	CreateAmberTrj($ATOMS, $BOX, \*OUTFILE);
        print "Done\r";
    }
    close OUTFILE or die "ERROR: Cannot finalize file $saveName: $!\n";
    print "\nDone\n";
}

sub init {
    my (@FILES);

    if (-d $bgfLoc) {
	opendir BGF, $bgfLoc or die "ERROR: Cannot open directory $bgfLoc: $!\n";
	@FILES = grep { /\.bgf$/ } map {"$bgfLoc/$_"} readdir BGF;
	closedir BGF or die "ERROR: Cannot close directory $bgfLoc: $!\n";
	print "Found " . ($#FILES + 1) . " BGF files\n";
    } elsif (-e $bgfLoc) {
	FileTester($bgfLoc);
	push @FILES, $bgfLoc;
	print "Found BGF file $bgfLoc\n";
    }

    die "ERROR: No bgf file found!\n" if (! @FILES);

    return \@FILES;
}

sub numerically {
    ($a<=>$b);
}
