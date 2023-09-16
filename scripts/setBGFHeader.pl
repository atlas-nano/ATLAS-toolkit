#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);

sub init;
sub updateHeader;
sub saveBGFs;

my ($BGFS);

$|++;
&init;
&updateHeader($BGFS);
&saveBGFs($BGFS);

sub saveBGFs {
    my ($bgfs) = $_[0];
    my ($i, $count);
    
    $count = 0;
    for $i (@{ $bgfs }) {
	$count++;
	print "Creating BGF file $i->{FILE}...";
	addHeader($i->{ATOMS}, $i->{HEADERS});
	createBGF($i->{ATOMS}, $i->{BONDS}, $i->{FILE});
	print "Done\r";
    }

    print "Creating BGF files.... ${count} file created                \n";
}

sub updateHeader {
    my ($bgfs) = $_[0];
    my ($i, $sName, $count);
    
    $count = 0;
    for $i (@{ $bgfs }) {
	$count++;
	print "Parsing BGF file $i->{FILE}...";
	($i->{ATOMS}, $i->{BONDS}) = GetBGFFileInfo($i->{FILE}, 0);
	$sName = basename($i->{FILE});
	$sName =~ s/\.\w+$//;
	$sName =~ s/$i->{PREFIX}//;
	@{ $i->{HEADERS} } = ("BIOGRF  321", "DESCRP $i->{PREFIX} $sName", 
			      "REMARK Created by $ENV{USER} \@ $ENV{HOST} at " . scalar(localtime),
			      "FORCEFIELD DREIDING");
	print "Done\r";
    }
    print "Parsing BGF files... $count file found                              \n";
}

sub init {
    my (%OPTS, $findCmd, $rec);
    
    getopt('bp',\%OPTS);
    for ("b", "p") {
	die "usage: $0 -b bgf(s) location -p (save prefix)\n" if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    
    $findCmd = "find $OPTS{b} -name '*.bgf' -print";
    if (-e $OPTS{b}) {
	if ($OPTS{b} =~ /$OPTS{p}/) {
	    $rec->{PREFIX} = $OPTS{p};
	    $rec->{FILE} = $OPTS{b};
	    push @{ $BGFS }, $rec;
	}
    } elsif (open(FINDCMD, "$findCmd |")) {
        while (<FINDCMD>) {
	    $rec = ();
	    $rec->{PREFIX} = $OPTS{p};
	    next if ($_ !~ /$rec->{PREFIX}/);
            chomp;
	    $rec->{FILE} = $_;
            push @{ $BGFS }, $rec;
        }
        close FINDCMD;
    } else {
        die "ERROR: No valid BGF files found while searching $OPTS{b}\n";
    }
    die "ERROR: No valid BFF files in $OPTS{b} found with prefix $OPTS{p}!\n" if (! $BGFS);
    print "Done\n";
}
