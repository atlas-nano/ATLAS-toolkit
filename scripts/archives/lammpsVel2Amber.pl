#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General;
use Packages::FileFormats;
use Packages::LAMMPS;
use Packages::AMBER;

sub parseLammpsData;
sub createAmberTrj;


sub parseLammpsData {
    my ($bgfInfo, $lammpsFile, $resinfo, $outFile, $printStr) = @_;
    my ($inStr, $BOX, $counter, $ATOMS, @dim, $validTStep, $index, $snapC);
    my ($atomC, $HEADERS, $currPos, $progress, $start, $strLen, $tot, $endFrame);
    my ($filesize) = -s $lammpsFile;

    
    $ATOMS = $BOX = $HEADERS = ();
    $counter = $validTStep = $index = $snapC = 0;
    $start = 0;

    $strLen = 0;
    if (defined($SELECT)) {
	@dim = sort numerically keys %{ $SELECT };
	$tot = $#dim + 1;
	$endFrame = pop @dim;
    } else {
	$tot = 0;
    }

    @dim =  ("X", "Y", "Z");
    open LAMMPS, $lammpsFile or die "ERROR: Cannot open LAMMPS trajectory file $lammpsFile: $!\n";
    while (<LAMMPS>) {
	chomp;
	$inStr = $_;
	if ($inStr =~ /^ITEM: TIMESTEP$/) { # either timestep of #atoms
	    if ($ATOMS and $BOX and $HEADERS) {
		#print "SAVING\n";
		createAmberTrj($ATOMS, $BOX, $resinfo, $HEADERS, $outFile);
		$counter = 0;
		$index++;
		if (! $tot) {
		    $currPos = tell(LAMMPS);
		    $strLen = PrintProgress($currPos, $filesize, $start, $printStr);
		} else {
		    $strLen = PrintProgress($index, $tot, $start, $printStr);
		    if ($snapC == $endFrame) {
			$BOX = $ATOMS = $HEADERS = ();
			$validTStep = 0;
			last;
		    }
		}
	    }
	    $BOX = $ATOMS = $HEADERS = ();
	    $validTStep = 0;
	    $counter = 0;
	    #print "$inStr\n";
	    push @{ $HEADERS }, $inStr;
	    $snapC++;
	    if (! $tot || exists($SELECT->{$snapC})) {
		$validTStep = 1;
		if (! $start) {
		    $start = time();
		}
	    }
	} elsif ($validTStep) {
	    if ($inStr =~ /^(\-?\d+\.?\d*) (\-?\d+\.?\d*)$/) { #box info
		push @{ $HEADERS }, $inStr;
		$BOX->{$dim[$counter]} = (
					  {
					      "lo" => $1,
					      "hi" => $2,
					  }
					  );
		$counter++;
	    } elsif ($inStr =~ /^(\d+)(\s\-?\d+.+)/) {
		$atomC = $1;
		$ATOMS->{$atomC}{"bondpointer"} = $bgfInfo->{$atomC}{"bondpointer"};
		$counter = 0;
		$inStr = $2;
		while ($inStr =~ /\s(\-?\d+\.?\d*e?\-?\d*)/g and $counter < 6) {
		    $ATOMS->{$atomC}{$dim[$counter]} = $1;
		    $counter++;
		}
	    } else {
		#print "$inStr\n";
		push @{ $HEADERS }, $inStr;
	    }
	}
    }
    close LAMMPS or die "ERROR: Cannot close $lammpsFile: $!\n";

    if ($ATOMS and $BOX) {
	createAmberTrj($ATOMS, $BOX, $resinfo, $HEADERS, $outFile);
    }
    printf "$printStr%-${strLen}s\n", "Done";
}

sub createAmberTrj {
    my ($ATOMS, $BBOX, $bondList, $HEADERS, $OUTFILE) = @_;
    my ($atomC, $counter, @tmp, $dim, $BOX);

    $counter = 0;
    @tmp = sort numerically keys %{ $ATOMS };
    if (defined($BBOX)) {
	for $dim ("X", "Y", "Z") {
	    $BOX->{$dim} = $BBOX->{$dim}{"hi"} - $BBOX->{$dim}{"lo"};
	}
    }

    for $atomC (@tmp) {
	for $dim ("X", "Y", "Z") {
	    $counter++;
	    if ($scaled and defined($BBOX)) {
		$ATOMS->{$atomC}{$dim} *= $BOX->{$dim}
	    }
	    if (exists($ATOMS->{$atomC}{$dim . "i"}) and defined($BBOX)) {
		printf $OUTFILE "%8.3f", $ATOMS->{$atomC}{$dim} + ($BOX->{$dim} * $ATOMS->{$atomC}{$dim . "i"});
	    } else {
		printf $OUTFILE "%8.3f", $ATOMS->{$atomC}{$dim};
	    }
	    if ($counter == 10) {
		print $OUTFILE "\n";
		$counter = 0;
	    }
	}
    }
    if (defined($BBOX)) {
	for $dim ("X", "Y", "Z") {
	    printf $OUTFILE "%8.3f", $BOX->{$dim};
	}
	print $OUTFILE "\n";
    }
}
