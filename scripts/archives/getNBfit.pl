#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::General qw(FileTester STDev);

sub init;
sub getNBParms;
sub printStats;

my ($FILES, $template, $graceCmd);

$|++;

$FILES = &init;

&getNBParms($FILES, $template, $graceCmd);
&printStats($FILES);

sub printStats {
    my ($data) = $_[0];
    my ($i, %STATS, $j, @tmp);

    for $i (@{ $data }) {
	if ($i->{STATS}{"Correlation coefficient"} < 0.75 or $i->{STATS}{a0} < 0) {
	    $i->{DATAFILE} .= "*";
	    next;
	}
	for $j (keys %{ $i->{STATS} }) {
	    next if ($j !~ /^a\d/);
	    $STATS{$j}{DATA} .= "$i->{STATS}{$j} ";
	}
	@tmp = keys %{ $i->{STATS} };
    }

    for $i (keys %STATS) {
	chomp $STATS{$i};
	($STATS{$i}{AVG},$STATS{$i}{STDEV},$j) = STDev($STATS{$i}{DATA});
    }

    @tmp = sort {$a cmp $b } @tmp;
    printf "%-25s", " ";
    for $i (@tmp) {
	printf "%15s",substr($i,0,15);
    }
    print "\n";

    for $i (@{ $data }) {
	printf "%-25s",$i->{DATAFILE};
	for $j (@tmp) {
	    printf "%15.6f", $i->{STATS}{$j};
	}
	print "\n";
    }
    for $i ("AVG", "STDEV") {
	printf "%-25s",$i;
	for $j (@tmp) {
	    if (! exists($STATS{$j}{$i})) { 
		printf "%15s", " ";
	    } else {
		printf "%15.6f", $STATS{$j}{$i};
	    }
	}
	print "\n";
    }
}

sub getNBParms {
    my ($files, $temp, $grace) = @_;
    my ($i, $numSpaces, $templateData, $templateCopy, $execCmd);

    $execCmd = "$grace -nosafe -hardcopy -hdevice 'PNG' -batch _tmp.dat";
    print "Reading template file $temp...";
    open TEMPLATE, $temp, or die "ERROR: Cannot read template file $temp: $!\n";
    while (<TEMPLATE>) {
	$templateData .= $_;
    }
    close TEMPLATE;
    print "Done\n";

    for $i (@{ $files }) {
	$numSpaces = length($i->{DATAFILE}) + 16;
	print "$i->{DATAFILE}: Getting Parms\r";
	$templateCopy = $templateData;
	$templateCopy =~ s/file\.dat/$i->{DATAFILE}/;
	$templateCopy =~ s/file_avg\.dat/$i->{AVGFILE}/;
	$templateCopy =~ s/file_name/$i->{TITLE}/;
	$templateCopy =~ s/file\.agr/$i->{SAVEFILE}/;
	open TMP, "> _tmp.dat" or die "ERROR: Cannot create tmp file _tmp.dat: $!\n";
	print TMP $templateCopy;
	close TMP;
	open PIPECMD, "$execCmd |" or die "ERROR: Cannot exectute $execCmd: $!\n";
	while (<PIPECMD>) {
	    if ($_ =~ /(a\d) = (\-?\d+\.?\d*)/) {
		$i->{STATS}{$1} = $2;
	    } elsif ($_ =~ /^(.+): (\-?\d+\.?\d*)/) {
		$i->{STATS}{$1} = $2;
	    }
	}
	printf "%${numSpaces}s\r"," ";
    }
    print "Got parms for " . (scalar @{ $files }) . " files\n";
}
sub init {
    my (%OPTS, @FDATA, $findCmd, @datFiles, $rec, $i, $fileTitle, $avgFile, $saveFile);

    getopt('lgst',\%OPTS);

    for ($OPTS{l},$OPTS{g}, $OPTS{t}) {
	die "usage: $0 -t grace template -g grace location -l data file|directory -s [save prefix]\n" 
	    if (! defined($_));
    }
    
    print "Initializing...";
    
    #die "ERROR: Cannot access grace template file $OPTS{t}: $!\n"
#	if (! FileTester($OPTS{t}));
    $template = $OPTS{t};
    
    die "ERROR: Expected file or directory when searching in $OPTS{l}!\n"
	if (! -e $OPTS{l} or ! -d $OPTS{l});
#    die "ERROR: Cannot access grace executable in $OPTS{g}: $!\n"
#	if (! FileTester($OPTS{l}) or ! -x $OPTS{g});
    $graceCmd = $OPTS{g};

    $findCmd = "find $OPTS{l} -name '*.dat' -print";
    if (open(FINDCMD, "$findCmd |")) {
	while (<FINDCMD>) {
	    chomp;
	    push @datFiles, $_;
	}
	close FINDCMD;
    }
 
    die "ERROR: Cannot find any .dat files in $OPTS{l}\n" if (! @datFiles);

    for $i (@datFiles) {
	$avgFile = $i;
	$avgFile =~ s/\.dat$/_avg\.dat/;
	next if (! -e $avgFile);
	$fileTitle = basename($i);
	$saveFile = $i;
	$saveFile =~ s/\.dat$/\.avg/;
	if ($OPTS{s}) {
	    $saveFile = $OPTS{s} . $saveFile;
	}
	if ($fileTitle =~ /^([a-z]+)([A-Z]+[A-Za-z]+)/) {
	    $fileTitle = uc($1) . " " . uc($2);
	}
	
	$rec = (
		{
		    "DATAFILE" => $i,
		    "AVGFILE"  => $avgFile,
		    "TITLE"    => $fileTitle,
		    "SAVEFILE" => $saveFile,
		}
		);
	push @FDATA, $rec;
    }

    die "ERROR: No file with _avg.dat found when searching $OPTS{l}\n" 
	if (! @FDATA);
    
    print "Done\n";
    return \@FDATA;
}

