#!/usr/bin/perl -w
use strict;
use warnings;

sub DetermineLegitInt(@);
sub GetValidFiles();
sub MoveFiles(@);

my ($workingdir, $filebase, $startnum, $endnum, $interval, $totalbases) = @ARGV;
my (@valid_files, $i, $out_cmd);

if (! @ARGV or $#ARGV < 4) {
    die "usage: getmoreparms.pl workingdir filebase startnumber endnumber #bases\n";
}

#-d "./" . $workingdir ."/" or die "Cannot locate $workingdir: $!\n";

chdir "$workingdir";

DetermineLegitInt($startnum, "Starting Number");
DetermineLegitInt($endnum, "Ending Number");

GetValidFiles();

for $i (@valid_files) {
    print "Analysing $i...";
    $out_cmd = "/ul/maiti/Curve_linux/data/gex <<\!\n";
    $out_cmd .= "$i\nn\n0\n\!\n";
    open TMPFILE, "> tmpparmfile" or die "Cannot create temporary file: $!\n";
    print TMPFILE "$out_cmd";
    close TMPFILE;
    system "chmod +x tmpparmfile";
    if (! system "tmpparmfile >& morejunk") {
	if ($i =~ /_(\d+)$/) {
	    print "Sucess\n";
	    system "mv $i" . ".minor $filebase" . "_minor". ".$1" ;
	    system "mv $i" . ".major $filebase" . "_major". ".$1" ;
	    system "mv $i" . ".minor_depth $filebase" . "_minor_depth". ".$1" ;
	    system "mv $i" . ".major_depth $filebase" . "_major_depth". ".$1" ;
	    system "mv $i" . ".radius $filebase" . "_radius". ".$1" ;
	} else { print "Failure: $i\n"; }
    } else { print "WARNING: Failure when processine $i\n"; }
    system "rm -f tmpparmfile morejunk";

}

# move files

MoveFiles("*_minor.*", "../minor_data");
MoveFiles("*_major.*", "../major_data");
MoveFiles("*_minor_depth.*", "../minor_depth");
MoveFiles("*_major_depth.*", "../major_depth");
MoveFiles("*_radius.*", "../radius_data");

# execute program to create graphs
print "\n";
print "Calculating minor groove depth...";

chdir "../minor_depth";
$out_cmd = "/ul/maiti/bgfstr/linux/groove $totalbases 0 $filebase" . "_minor_depth". "." . $startnum . " $interval " . ($endnum - $startnum + 1);

if (! system "$out_cmd >& junk") {
    print "Sucess\n";
} else {
    print "Failure\n";
}

system "mv junk ../minor_groove_depth.out";

chdir "../major_depth";
print "\n";
print "Calculating major groove depth...";

$out_cmd = "/ul/maiti/bgfstr/linux/groove $totalbases 0 $filebase" . "_major_depth". "." . $startnum . " $interval " . ($endnum - $startnum + 1);

if (! system "$out_cmd >& junk") {
    print "Sucess!!\n";
} else {
    print "Failure!!\n";
}

system "mv junk ../major_groove_depth.out";

system "rm -f junk";

sub MoveFiles(@) {

    my ($filespec, $storedir) = @_;

    if (! -d "$storedir" ) {
	mkdir "$storedir";
    }

    system "mv $filespec $storedir";

}

sub GetValidFiles() {
    my $i = 0;
    my $j = 0;
    my $k = 0;
    my $counter = 0;
    my $curr_fle = "";
    my $line_in = "";

    for ($i = $startnum; $i<=$endnum;$i++) {
	$curr_fle = $filebase . "." . $i;
	$curr_fle =~ /(.+)_lis/;
	$line_in = $1 . "_$i";
	if (-e $curr_fle) {
	    system "mv $curr_fle $line_in" . ".lis";
	    push @valid_files, "./$line_in";
	    $counter +=1;
	} elsif (-e "$line_in" . ".lis") {
	    push @valid_files, "./$line_in";
	    $counter +=1;
	} else {
	    print "Warning: Cannot find file $curr_fle\n";
	}
    }

    if ($counter ==0) {
	die "Cannot find any valid files\n";
    }

}

sub DetermineLegitInt(@) {
    my ($inval) = $_[0];
    my ($e_message) = $_[1];

    if (! $inval =~ /^\d+$/) {
	die "Invalid value in $e_message. Expected Integer\n";
    }
}


