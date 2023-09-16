#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use General qw(FileTester);
use Getopt::Std qw(getopt);

sub init;
sub getAmberCharges;
sub createBGFfile;

my ($pdbFile, $saveFile);
my ($preFile);

$|++;
&init;
print "Getting AMBER charges for $pdbFile...";
$preFile = getAmberCharges($pdbFile);
print "Done\nCreating bgf file $saveFile...";
&createBGFfile($preFile, $saveFile);
print "Done\n";

sub createBGFfile {
    my ($prefix, $saveName) = @_;
    my ($script) = "/ul/tpascal/scripts/amber2bgf.pl";

    die "ERROR: Cannot locate script $script: $!\n"
	if (! -e $script || ! -r $script || ! -T $script);
    if(system("$script ${prefix}.prmtop ${prefix}.rst7 $saveName >& /dev/null")) {
	die "ERROR while executing \"$script ${prefix}.prmtop ${prefix}.rst7 $saveName\"\n";
    }
}

sub getAmberCharges {
    my ($inFile) = $_[0];
    my ($script, $prefix); 

    $script = "/ul/tpascal/scripts/prep_sander.pl";
    $prefix = basename($inFile);
    $prefix =~ s/\.\w+$//;

    die "ERROR: Cannot find script file $script:$!\n" 
	if (! -e $script || ! -r $script || ! -T $script);
    if(system("$script $inFile $prefix /ul/tpascal/amber/")) {
	die "ERROR while executing \"$script $inFile $prefix /ul/tpascal/amber/\"\n";
    }
    return $prefix;
}

sub init {
    my (%OPTS);
    getopt('ps',\%OPTS);
    die "usage: $0 -p pdbfile -s (savebgf)\n" if (! exists($OPTS{p}));
    print "Initializing...";
    ($pdbFile, $saveFile) = ($OPTS{p}, $OPTS{s});
    FileTester($pdbFile);
    if (! defined($saveFile)) {
	$saveFile = basename($pdbFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_amber.bgf";
    }
    print "Done\n";
}
