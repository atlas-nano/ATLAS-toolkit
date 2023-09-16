#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetPDBFileInfo createBGF addHeader createHeaders);
use General qw(FileTester);

sub init;

my ($pqrFile, $jagOutFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Parsing pqrFile $pqrFile...";
($ATOMS, $BONDS) = GetPDBFileInfo($pqrFile);
print "Done\nCreating BGF file $saveFile...";
$HEADERS = createHeaders(undef,$saveFile);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub init {
    my (%OPTS);
    getopt('ps',\%OPTS);
    
	die "usage: $0 -p pqrfile -s (bgf_file)\n" if (! exists($OPTS{p}));
    
    print "Initializing...";
    ($pqrFile, $saveFile) = ($OPTS{p}, $OPTS{s});
    FileTester($pqrFile);
    
    if (! defined($saveFile)) {
		$saveFile = basename($pqrFile);
		$saveFile =~ s/\.\w+$//;
		$saveFile .= ".bgf";
    }
    print "Done\n";
}
