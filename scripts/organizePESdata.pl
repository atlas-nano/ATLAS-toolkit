#!/usr/bin/perl -w

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub parseCSVfile;
sub reorgData;
sub writeCSVdata;
sub numerically { ($a<=>$b); }
sub alphabetically { ($a cmp $b); }

my ($engFile, $refFile, $saveFile);
my ($ENG, $REF, $DATA);

$|++;
&init;
print "Parsing energy file $engFile...";
$ENG = parseCSVfile($engFile);
print "Done\nParsing reference energy file $refFile...";
$REF = parseCSVfile($refFile);
print "Done\nReorganizing data...";
$DATA = reorgData($ENG, $REF);
print "Done\nCreating CSV file $saveFile...";
&writeCSVdata($DATA, $saveFile);
print "Done\n";

sub writeCSVdata {
    my ($data, $saveName) = @_;
    my (@headers, $i, @dist, $j);

    @headers = sort alphabetically keys %{ $data };
    @dist = sort numerically keys %{ $data->{$headers[0]} };
    open CSVDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    print CSVDATA "dist";
    for $i (@headers) {
	print CSVDATA ",$i";
    }
    print CSVDATA "\n";
    
    for $i (@dist) {
	print CSVDATA "$i";
	for $j (@headers) {
	    print CSVDATA ",$data->{$j}{$i}";
	}
	print CSVDATA "\n";
    }
    close CSVDATA;
}

sub reorgData {
    my ($eng, $ref) = @_;
    my ($i, $j, %CSVDATA, $k, $refData, $engData);

    for $i (keys %{ $eng })  {
	die "ERROR: Reference for key $i not found! Aborting\n"
	    if (! exists($ref->{$i}));
	for $j (keys %{ $ref->{$i} }) {
	    $refData = $j;
	}
	for $j (keys %{ $eng->{$i} }) {
	    for $k (keys %{ $eng->{$i}{$j} }) {
		$engData = $k;
	    }
	    $CSVDATA{$i}{$j} = $engData - 2*$refData;
	}
    }

    return \%CSVDATA;
}

sub parseCSVfile {
    my ($csvFile) = $_[0];
    my (%CSVDATA, $data, $key, $curr);

    open CSVFILE, $csvFile or die "ERROR: Cannot open $csvFile: $!\n";
    while (<CSVFILE>) {
	chomp;
	if ($_ =~ /^\s*(\w+)\,(.+)/) {
	    ($key, $data) = ($1, $2);
	    $curr = \%{ $CSVDATA{$key} };
	    while ($data =~ /(\-?\d+\.\d+)/g) {
		$curr->{$1} = ();
		$curr = \%{ $curr->{$1} };
	    }
	}
    }
    close CSVFILE;
    
    die "ERROR: $csvFile is invalid!\n" if (! %CSVDATA);
    return \%CSVDATA;
}

sub init {
    my (%OPTS);
    
    getopt('ers',\%OPTS);
    die "usage: $0 -e energy file csv (from grep) -r reference energy csv file -s (savename)\n"
	if (! exists($OPTS{e}) or ! exists($OPTS{r}));
    print "Initializing...";
    for ("e", "r") {
	die "ERROR: Cannot access csv file $OPTS{$_}: $!\n"
	    if (! -e $OPTS{$_} or ! -r $OPTS{$_} or ! -T $OPTS{$_});
    }
    ($engFile, $refFile, $saveFile) = ($OPTS{e}, $OPTS{r}, $OPTS{s});
    if (! defined($saveFile)) {
	$saveFile = basename($engFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_reorg.csv";
    }
    print "Done\n";
}
