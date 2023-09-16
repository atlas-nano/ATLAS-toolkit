#!/usr/bin/perl -w

use Getopt::Std qw(getopt);
use strict;

$|++;
$FILES = &init;

sub init {
    my ($findCmd, @datFiles, %OPTS);

    getopt('lsf', \%OPTS);
    for ("l", "f") {
	die "ERROR: usage: $0 -l file location -f file name (prefix/suffix) -s savedir -c (burried|rank|energy) -t (total|cuttoff)\n"
	    if (! defined($OPTS{$_}));
    }

    print "Initializing...";
    $findCmd = "find $OPTS{l} -name '$OPTS{f}' -print";
    if (open(FINDCMD, "$findCmd |")) {
        while (<FINDCMD>) {
            chomp;
            push @datFiles, $_;
        }
        close FINDCMD;
    }
    die "ERROR: No valid file found file searching for $OPTS{f} in directory $OPTS{l}\n"
	if (! @datFiles);
    
    if ($OPTS{c} !~ /^(b|r|e)/i) {
	print "Invalid option for selection type. Selecting on burried surface\n";
	$OPTS{c} = "b";
    }

    if ($OPTS{c} eq "b") {
    }
}
