#!/usr/bin/perl -w
use strict;
use File::Basename;

if (! @ARGV or $#ARGV < 2) {
    die "usage: $0 filebase starttime endtime\n";
}

my ($filebase, $starttime, $endtime, $interval) = @ARGV;
my ($temp, $i, $counter, $infile, $outfile, $out_cmd, @fle_list, $file_list);

if ($starttime =~ /(\d+)/) {
    if ($1 > 0) {
	$starttime = $1;
    } else {
	die "Invalid starting time\n";
    }
} else {
    die "Invalid startting time\n";
}

if ($endtime =~ /(\d+)/) {
    if ($1 > 0) {
	$endtime = $1;
    } else {
	die "Invalid ending time\n";
    }
} else {
    die "Invalid ending time\n";
}


if ($starttime > $endtime) {
    $starttime = $temp;
    $starttime = $endtime;
    $endtime = $temp;
}

$counter = $starttime;

$i = $starttime;
while ($i <= $endtime) {
    $counter = sprintf("%0" . length($endtime) . "d", $counter);
    $infile = $filebase . $i . ".gif";
    $outfile = $filebase . $counter . ".gif";
    system "mkdir -p conv";
    if (-e $infile) {
#	if ($infile ne $outfile) {
#   		$out_cmd = "convert $infile $outfile";
#		$outfile = "conv/" . basename($outfile);
#		print "$infile -> $outfile...";
#		if (! system("$out_cmd") ) {
#		    print "Sucess\n";
#		    $counter++;
#		    push @fle_list, $outfile;
#		} else {
#		    print "Failure\n";
#		}
#	} else {
	    $counter++;
	    push @fle_list, $infile;
#	}
    }
    $i += $interval;
}

if ($#fle_list > -1) {

    print "\nCreating movie...";
    $out_cmd = "makemovie -o dna.avi ";
    $out_cmd .= "-f avi -c qt_cvid -r 25 -k 1 ";

    for $i (0 .. $#fle_list) {
	$file_list .= " " . $fle_list[$i];
    }

    $out_cmd .= $file_list;
    system "$out_cmd";
    open OUTFILE, "> movie_cmd" or die "Cannot write to movie_cmd: $!\n";
    print OUTFILE $out_cmd;
    close OUTFILE;
    print "Done\n";

#    system "rm -f $file_list";
} else {
    die "Could not locate any valid files\n";
}









