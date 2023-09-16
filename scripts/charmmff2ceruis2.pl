#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);

use CERIUS2 qw(saveCeriusFF parseCharmmFF);
use General qw(FileTester);
use File::Basename;

my ($charmmTop, $charmmFF, $cerius2FF);;
my ($PARMS);

$|++;
&init;
print "Parsing CHARMM force field $charmmFF...";
&parseCharmmFF($charmmFF, 0, \%{ $PARMS });
print "Done\nCreating CERIUS2 force field $cerius2FF...";
&saveCeriusFF($PARMS, $cerius2FF);
print "Done\n";

sub init {
	my (%OPTS);

	getopt('tps',\%OPTS);
	die "usage: $0 -t charmmTop -p charmmPar -s [saveName]\n"
    	if (!exists($OPTS{t}) or ! exists($OPTS{p}));
	($charmmTop, $charmmFF, $cerius2FF) = ($OPTS{t}, $OPTS{p}, $OPTS{s});
    print "Initializing...";
    FileTester($charmmFF);
    FileTester($charmmTop);
    if (! defined($cerius2FF)) {
		$cerius2FF = basename($charmmFF);
		$cerius2FF =~ s/\.\w+$/\.ff/;
    }
    print "Done\n";
}
