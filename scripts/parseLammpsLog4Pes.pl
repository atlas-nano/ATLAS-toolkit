#!/usr/bin/perl -w

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub parseLogFile;
sub writeData;
sub numerically { ($a<=>$b); }
sub MIN;

my ($logFile, $saveFile, $offset);
my ($DATA, $STATS);

$|++;
&init;
print "Parsing LAMMPS Log file $logFile...";
$DATA = parseLogFile($logFile);
print "Done\nWriting data to $saveFile...";
&writeData($DATA, \%{ $STATS }, $saveFile, $offset);
print "Done\n==========STATS=========\n";
for my $i (keys %{ $STATS }) {
    printf "%-10s %.2f\n", $i, $STATS->{$i};
}

sub writeData {
    my ($data, $stats, $outFile, $offset) = @_;
    my ($i, $j, @xvals, @yvals, $distfile, $eng, $dist, $header);

    $distfile = $outFile;
    $distfile =~ s/_\w+\.\w+$//;
    $distfile .= "_dist.csv";

    @xvals = sort numerically keys %{ $data };
    shift @xvals;
    @yvals = sort numerically keys %{ $data->{ $xvals[1] } };

    $stats->{min} = 99999999999999999;
    $stats->{max} = -9999999999999999;

    for $i (@xvals) {
	for $j (@yvals) {
	    #next if (! $data->{$i}{$j}); 
	    ($eng->{$i}{$j}, $dist->{$i}{$j}) = MIN($data->{$i}{$j});
            # record stats
            $stats->{X} += $eng->{$i}{$j};
            $stats->{count}++;
            $stats->{min} = $eng->{$i}{$j} if ($eng->{$i}{$j} < $stats->{min});
            $stats->{max} = $eng->{$i}{$j} if ($eng->{$i}{$j} > $stats->{max});
            $stats->{X2} += $eng->{$i}{$j}**2;
	}
    }

    $stats->{avg} = $stats->{X}/$stats->{count};
    $stats->{X2} /= $stats->{count};
    $stats->{stdev} = sqrt($stats->{X2} - $stats->{avg}**2);
    delete $stats->{X};
    delete $stats->{X2};

    open ENGFILE, "> $outFile" || die "ERROR: Cannot write to $outFile: $!\n";
    open DISTFILE, "> $distfile" || die "ERROR: Cannot write to $distfile: $!\n";
    print ENGFILE ",";
    print DISTFILE ",";
    for $i (@xvals) {
	$header = $i;
	$header += $offset->{x} if (exists$offset->{x});
        print ENGFILE "$header,";
        print DISTFILE "$header,";
    }
    print ENGFILE "\n";
    print DISTFILE "\n";

    for $j (@yvals) {
	$header = $j;
	$header += $offset->{y} if (exists($offset->{y}));
        print ENGFILE "$header,";
        print DISTFILE "$header,";
        for $i (@xvals) {
	    #if (! exists($eng->{$i}{$j})) {
		#print ENGFILE ",";
		#print DISTFILE ",";
		#next;
	    #}
	    $eng->{$i}{$j} -= $stats->{$offset->{e}} if (exists($offset->{e}));
	    $dist->{$i}{$j} -= $offset->{d} if (exists($offset->{d}));
            print ENGFILE "$eng->{$i}{$j},";
            print DISTFILE "$dist->{$i}{$j},";
        }
        print ENGFILE "\n";
        print DISTFILE "\n";
    }
    close ENGFILE;
    close DISTFILE;


}

sub MIN {
    my ($vals) = $_[0];
    my ($minE, $minD, @tmp, $i);
    
    @tmp = keys %{ $vals };
    $minE = $vals->{ pop @tmp };
    for $i (@tmp) {
	if ($vals->{$i} < $minE) {
	    $minE = $vals->{$i};
	    $minD = $i;
	}
    }

    return ($minE, $minD);
}

sub parseLogFile {
    my ($inFile) = $_[0];
    my (%VDW, $curr, $stats);

    open LOGFILE, $inFile || die "ERROR: Cannot open $inFile: $!\n";
    while (<LOGFILE>) {
	chomp;
	if ($_ =~ /DATA: (\d+\.?\d*) (\d+\.?\d*) (\d+\.?\d*) (\-?\d+\.\d+)/) {
	    $VDW{$1}{$2}{$3} = $4;
	}
    }
    close LOGFILE;
    die "ERROR: No valid data read!\n" if (! %VDW);
    return \%VDW;
}

sub init {
    my (%OPTS);
    getopt('lsexyd',\%OPTS);
    die "usage: $0 -l log file -s (save file) -e (energy offset: max|min|avg|none(default) -x (x offset) -y (y offset) -d (distance offset)\n" 
	if (! exists($OPTS{l}));
    print "Initializing...";
    ($logFile, $saveFile) = ($OPTS{l}, $OPTS{s});
    die "ERROR: Cannot access log file $logFile: $!\n"
	if (! -e $logFile || ! -r $logFile || ! -T $logFile);
    if (! defined($saveFile)) {
	$saveFile = basename($logFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_pes.csv";
    }
    $offset = ();
    if (exists($OPTS{e})) {
	if ($OPTS{e} =~ /(max|min|avg)/i) {
	    $offset->{e} = lc $1;
	}
    }
    for ("x", "y", "d") {
	if (exists($OPTS{$_}) and $OPTS{$_} =~ /(\-?\d+\.?\d*)/) {
	    $offset->{$_} = $1;
	}
    }

    print "Done\n";
}
