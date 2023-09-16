#!/usr/bin/perl -w
BEGIN {
	unshift @INC, "/home/tpascal/scripts/";
}
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use Packages::Namot2::Namot;
use strict;

sub Initialize;
sub MapChain;
sub GetChain;
sub ReadSpecsFile;
sub AssignVals;
sub Numerically;

die "usage: $0 helix1 helix2 which_helix specsFile\n"
    if (! @ARGV or $#ARGV < 3);

my ($helix1, $helix2, $myHelix, $specsFile) = @ARGV;
my ($is5prime, $isParallel, $total_bases, $periodicity, $half_turn, @crossovers);
my ($c_point, $i, @CROSS, $trans, $counter, @tmp);

my ($opts) = ReadSpecsFile($specsFile);
AssignVals;
Initialize;

#p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("render");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("render");
p5namot::Cmd("set hush WARNING off");
p5namot::Cmd("render");
p5namot::Cmd("set hush REQUESTED on");
p5namot::Cmd("render");

p5namot::Cmd("load pdb na $helix1");
p5namot::Cmd("render"); 
p5namot::Cmd("load pdb na $helix2");
p5namot::Cmd("render");
p5namot::Cmd("trans 2 $trans");
p5namot::Cmd("render");

# open(STDOUT, ">foo.out");
$i = 0;
while ($i < 360) {
    p5namot::Cmd("rotate $myHelix 3 1");
    p5namot::Cmd("render");
    
    for $c_point (@CROSS) {
	next
	    if ($c_point->{"1"} eq "0");
	@tmp = sort Numerically keys %{ $c_point };
	$counter = 0;
	while ( $counter <= ($#tmp - 1)) {
	    p5namot::Cmd("query dist " . $c_point->{$tmp[$counter]} . " " . $c_point->{$tmp[($counter + 1)]});
	      $counter += 2;
	  }
    }
    p5namot::Cmd("echo TER");
    p5namot::Cmd("render");
    $i++;
}

sub Initialize {
    my ($c_point, $rec);

    -e $helix1 or die "Cannot find $helix1\n";
    -e $helix2 or die "Cannot find $helix2\n";
    
    die "Expected integer for total_bases got $total_bases\n"
	if (! $total_bases =~ /^\d+/);

    for $c_point (0 .. $#crossovers) {
	if (! $crossovers[$c_point]{"BASE"} =~ /(\d+)/) {
	    die "Invalid crossover point, expected integer got " . $crossovers[$c_point]{"BASE"} . "\n";
	} else {
	    $crossovers[$c_point]{"BASE"} =~ /(\d+)/;
	    $crossovers[$c_point]{"BASE"} = $1;
	    if ($total_bases < $crossovers[$c_point]{"BASE"}) {
		$rec->{"1"} = 0;
	    } else {
		
		$rec = MapChain($crossovers[$c_point]{"BASE"}, $crossovers[$c_point]{"TYPE"});
	    }

	    push @CROSS, $rec;
	}
    }
    
    die "find_o_angle: Invalid input for 5primein, expected 0 or 1\n"    
	if ($is5prime <0 or $is5prime >1);
    
    die "Invalid helix, expected 1 or 2\n"
	if (! $myHelix =~/\d/);

    $myHelix =~ /(\d)/;
    $myHelix = $1;

    die "Invalid helix, expected 1 or 2\n" 
	if($myHelix > 2 or $myHelix <1);
    

}

sub MapChain {
    my ($currCross, $crossType) = @_;
    my ($ret, $chain1, $chain2, $h1, $h2);

    $chain1 = GetChain($currCross, 1);
    $chain2 = GetChain($currCross, 2);

    $h1 = $h2 = ($currCross+1);

    if ($chain1 == 2) {
	$h1--;
    }

    if ($chain2 == 2) {
	$h2--;
    }

    if ($crossType == 1) {
	$ret = (
		{
		    "1" => "h1:" . $h1 . ":$chain1:P ",
		    "2" => "h2:" . ($h2 + 1) . ":$chain2:P ",
		    "3" => "h1:" . ($h1 + 1) . ":$chain1:P ",
		    "4" => "h2:" . $h2 . ":$chain2:P ",
		}
		);
    } else {
	$ret = (
		{
		    "1" => "h1:" . $h1 . ":$chain1:P ",
		    "2" => "h2:" . $h2 . ":$chain2:P ",
		    "3" => "h1:" . ($h1 + 1) . ":$chain1:P ",
		    "4" => "h2:" . ($h2 + 1) . ":$chain2:P ",
		}
		);
    }

    print $ret->{1} . " " . $ret->{2} . "\n";
    print $ret->{3} . " " . $ret->{4} . "\n";
    return $ret;
}	    

sub GetChain {
# GetChain - gets the chain of a molecule
    my ($curr_group, $whichHelix) = @_;
    my ($returnval);

    $curr_group -= ($half_turn/2);

    while ($curr_group > $periodicity) {
	$curr_group -= $periodicity;
    }

    #print "5Prime " . $is5prime . "\n";

    if ($whichHelix == 1) {
	if ($is5prime and $isParallel) { #Parallel (PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	} elsif (! $is5prime and $isParallel) { # Parallel (anti-PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	} elsif ($is5prime and ! $isParallel) { # AntiParallel (DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	} elsif (! $is5prime and ! $isParallel) { # AntiParallel (anti-DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	}
    } else {
	if ($is5prime and $isParallel) { #Parallel (PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	} elsif (! $is5prime and $isParallel) { # Parallel (anti-PX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	} elsif ($is5prime and ! $isParallel) { # AntiParallel (DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 2;
	    } else {
		$returnval = 1;
	    }
	} elsif (! $is5prime and ! $isParallel) { # AntiParallel (anti-DX)
	    if ($curr_group <= $half_turn) {
		$returnval = 1;
	    } else {
		$returnval = 2;
	    }
	}
    }
	    
    return $returnval;
}

sub ReadSpecsFile {
    my ($cfile) = $_[0];
    my (%specs, $validCounter, $myKey, $isValid, $line_in, $rec);

    $validCounter = $isValid = 0;
    open CROSSFILE, $cfile or die "Error while reading $cfile: $!\n";
    while (<CROSSFILE>) {
	chomp;
	$line_in = $_;
	if ($line_in =~ /^Helix (\d+)\s+(\d+):(\d+)\s+(\d+)\s+bases/) {
	    $myKey = "Helix" . $1;
	    $specs{$myKey} = (
			      {
				  "MajorGroove"   => $2,
				  "MinorGroove"   => $3,
				  "Periodicity"   => ($2 + $3),
				  "TotalBases"    => $4,
			      }
			      );
	    $validCounter++;
	} elsif ($line_in =~ /^Transalation: (\S+)\s+(\S+)\s+(\S+)/) {
	    $specs{"Trans"} = (
			       {
				   "x" => $1,
				   "y" => $2,
				   "z" => $3,
			       }
			       );
	    $validCounter++;
	} elsif ($line_in =~ /^Parallel: ([1|0])/) {
	    $specs{"isParallel"} = $1;
	    $validCounter++;
	} elsif ($line_in =~ /^5PrimeToLeft: ([1|0])/) {
	    $specs{"5prime"} = $1;
	    $validCounter++;
	}elsif ($line_in =~ /CROSSOVERS/) {
            $isValid = 1; # file must have word "CROSSOVERS"
	    $validCounter++;
        }elsif ($isValid and $line_in =~ /^(\d+):(\d+)->(\d+):(\d+) TYPE:([1|2])$/) { # valid crossover spec
	    $rec = (
		    {
			"TYPE" => $5,
			"BASE" => $2,
		    }
		    );
	    push @{ $specs{"crossovers"} }, $rec;
	}
    }
    close CROSSFILE;

    die "ERROR: No crossvers found while reading file $cfile\n"
	if (! $isValid or $validCounter < 5);
    return \%specs;
}

sub AssignVals {
    $is5prime = $opts->{"5prime"};
    $isParallel = $opts->{"isParallel"};
    $total_bases = $opts->{"Helix1"}{"TotalBases"};
    $total_bases = $opts->{"Helix2"}{"TotalBases"}
    if ($opts->{"Helix2"}{"TotalBases"} < $total_bases);
    $half_turn = $opts->{"Helix1"}{"MinorGroove"};
    $periodicity = $opts->{"Helix1"}{"MajorGroove"} + $half_turn;
    @crossovers = @{ $opts->{"crossovers"} };
    $trans = $opts->{"Trans"}{"x"}  . " " . $opts->{"Trans"}{"y"}  . " " . $opts->{"Trans"}{"z"};

}

sub Numerically {
    ($a<=>$b);
}
