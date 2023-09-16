#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);
use File::Basename;

sub init;
sub getEnergies;
sub writeData;

my ($saveName, $FILES, $DATA);

$|++;
$FILES = &init;
$DATA = getEnergies($FILES);
print "Creating $saveName...";
writeData($DATA, $saveName);
print "Done\n";

sub writeData {
    my ($data, $fName) = @_;
    my ($i, $j);

    open OUTDATA, "> $fName" or die "ERROR: Cannot create $fName: $!\n";
    for $i (@{ $data }) {
	print OUTDATA "$i->{FILE}: ";
	for $j (sort {$a cmp $b} keys %{ $i->{GROUP} }) {
	    printf OUTDATA "%11.3f%11.3f", $i->{GROUP}{$j}{ENG}, $i->{GROUP}{$j}{DEV};
	}
	print OUTDATA "\n";
    }
    close OUTDATA;
}

sub getEnergies {
    my ($files) = $_[0];
    my ($i, @DATA, $rec, $strLen);

    for $i (@{ $files }) {
	chomp $i;
	print "Reading $i...\r";
	$strLen = length($i);
	$rec = ();
	open ENGFILE, $i or die "ERROR: Cannot open $i: $!\n";
	while (<ENGFILE>) {
	    if ($_ =~ /^\s+(\d+)\s+(\-?\d+\.\d+)\s+(\d+\.\d+)/) {
		$rec->{GROUP}{$1} = (
				     {
					 "ENG" => $2,
					 "DEV" => $3,
				     }
				     );
	    }
	}
	close ENGFILE;
	$rec->{FILE} = $i;
	push @DATA, $rec if (keys %{ $rec });
    }
    
    die "ERROR: No valid data found in any of the data file!\n" if (! @DATA);
    printf "Reading files...%-${strLen}s\n", "Done";
    return \@DATA;
}
	    
sub init {
    my (%OPTS);
    my ($fileName, $currDir, $findCmd, @FILES);

    getopt('dfs',\%OPTS);
    ($fileName, $currDir, $saveName) = ($OPTS{f},$OPTS{d},$OPTS{s});
    die "usage: $0 -d [directory] -f file name -s [save name]\n" if (! defined($fileName));

    print "Initializing...";
    $currDir = "./" if (! defined($currDir));
    $findCmd = "find $currDir -name '$fileName' -print";
    open FINDCMD, "$findCmd |" or die "ERROR: Cannot execute '$findCmd': $!\n";
    while (<FINDCMD>) {
	push @FILES, $_;
    }
    close FINDCMD;
    die "ERROR: No valid files found while search for $fileName in $currDir!\n" if (! @FILES);

    if (! defined($saveName)) {
	$saveName = basename($fileName);
	$saveName =~ s/\.\w+$/_enthalpies.dat/;
    }

    print "Done\n";
    return \@FILES;
}
