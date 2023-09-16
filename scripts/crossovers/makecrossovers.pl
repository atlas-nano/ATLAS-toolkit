#!/usr/bin/perl -w
use strict;
use FindBin ();
use lib "$FindBin::Bin";
use Crossovers;
use Carp;

# This script will generate PX and JX superstructures rather than
# Using Namot2. It will read in the individual helices, create a chain
# topology map and write the corresponding pdb file.

sub Initialize;
sub createTopologyMap;
sub getCrossovers;
sub findTerminal;
sub setCrossovers;
sub createPics;
sub assignVals;
sub updateChainList;
sub input;

sub numerically;
die "usage: $0 helix-1.pdb helix-2.pdb outfile.pdb crossover-spec-file\n"
	if (! @ARGV or $#ARGV < 3);

my ($helix1_loc, $helix2_loc, $output_fle, $crossovers_file) = @ARGV;
my ($periodicity, %ATOMS, $CStructure, @crossovers, $totChain, $is5prime); 
my ($isParallel, $half_turn, $counter, $trans0, %CENTER, $opts);

$|++;

assignVals;
Initialize;

readNamotPDB($helix1_loc, \%ATOMS, 1, $trans0);
readNamotPDB($helix2_loc, \%ATOMS, 2, $opts->{"Trans"});
%{ $CENTER{1} } = %{ $ATOMS{"HELIX1"} };
%{ $CENTER{2} } = %{ $ATOMS{"HELIX2"} };
delete $ATOMS{"HELIX1"};
delete $ATOMS{"HELIX2"};
$ATOMS{1}{2} = reverseRes(\%{ $ATOMS{1}{2} });
$ATOMS{2}{2} = reverseRes(\%{ $ATOMS{2}{2} });

($CStructure, $totChain) = createTopologyMap($is5prime, $isParallel, \@crossovers, \%ATOMS);
writeData($CStructure, $output_fle, $trans0);
createPics($output_fle, $totChain);

sub createTopologyMap {
	my ($bl5, $helicesIsParallel, $cList, $atoms) = @_;
	my ($counter, $curr_chain, $curr_base, $curr_helix, %myCross, $eff_chain, $isDXCross, $isCross, $c2, $index, $cC, @tmp2);
	my (@tmp, $chain_counter, %CrossoverStructure, $base_counter, $isEnd, @chain_lengths, @tmp1, $usedPrevious, $cM, $tC);

	%myCross = setCrossovers($cList);
	$tC = $#{ $cList } + 1;
	$tC *= 2;
	for $index (sort numerically keys %{ $atoms }) {
	@tmp1 = sort numerically keys %{ $atoms->{$index} };
		for $counter (@tmp1) {
			$chain_lengths[($#chain_lengths + 1)] = $opts->{"Helix" . $index}{"TotalBases"} + 1;
			push @tmp, ($counter + (($index - 1) * 2));
		}
	}

	$chain_counter = $base_counter = 1;
	
	$curr_chain = $cC = 0;
	while ($#tmp > -1) {
		$counter = $curr_base = 1;
		$isEnd = $isDXCross = $usedPrevious = $isCross = 0;

		$eff_chain = $tmp[0];
		$curr_helix = int(($eff_chain - 1)/2) + 1;
		$curr_chain = $eff_chain - (($curr_helix - 1) * 2);

		while (! $isEnd ) {
			$isCross = 0;
			if (exists($atoms->{$curr_helix}{$curr_chain}{$curr_base})) {
				if ($curr_base > 1 and ! $usedPrevious) { # DX Crossover in middle of helix so find terminal group
					if (($tC - 1) == $cC) { # if this is the last crossover then join chains
						if (exists($myCross{$curr_helix}{$curr_chain}{($curr_base - 1)})) {
							@tmp1 = split /\,/, $myCross{$curr_helix}{$curr_chain}{($curr_base - 1)};
							$myCross{$curr_helix}{$curr_chain}{$curr_base} = $tmp1[0] . "," . 
							$tmp1[1] . "," . ($tmp1[2] + 1) . "," . $tmp1[3];
							delete $myCross{$tmp1[0]}{$tmp1[1]}{($tmp1[2] + 1)};
							@tmp2 = ($curr_helix,$curr_chain,$curr_base);
						}
					}
					($curr_helix, $curr_chain, $curr_base, $isDXCross) = 
					findTerminal($curr_helix, $curr_chain, $curr_base, \%myCross, $atoms, \@chain_lengths);
#					print "Found Terminal at $curr_helix, $curr_chain, $curr_base\n";
				}
				%{ $CrossoverStructure{$chain_counter}{$base_counter} } = %{ $atoms->{$curr_helix}{$curr_chain}{$curr_base} };
	#			print "$chain_counter : $base_counter - $eff_chain : $curr_base\n";
				delete $atoms->{$curr_helix}{$curr_chain}{$curr_base};
				$base_counter++;
		
				if (exists($myCross{$curr_helix}{$curr_chain}{$curr_base})) { # base in crossover so do crossover
					$cC++;
		#			print "CROSSED from $curr_helix, $curr_chain, $curr_base to ";
					($curr_helix, $curr_chain, $curr_base,$cM) = split /\,/, $myCross{$curr_helix}{$curr_chain}{$curr_base};
#					print "$curr_helix, $curr_chain, $curr_base TYPE: $cM\n";
					if ($cM == 2 and ($tC != $cC)) { 
						$isDXCross = 1;
					} else {
						$isDXCross = 0;
					}
					$isCross = 1;			
				} else {
					$isCross = 0;
				}
				$usedPrevious = 1;
			}

			if (! $isCross) {
				if (! $isDXCross) {
					$curr_base++;
				} else {
					$curr_base--;
				}
			}	

			$eff_chain = ($curr_chain + (($curr_helix - 1) * 2)) - 1;
			if ($curr_base > $chain_lengths[$eff_chain] or $curr_base < 1) {
				$isEnd = 1;
			}
		}
	
		@tmp = updateChainList($tmp[0], \@tmp, \% { $atoms });
	
		$chain_counter++;
		$base_counter = 1;

	}

	$chain_counter--;
	return (\%CrossoverStructure, $chain_counter);
}

sub updateChainList {
	my ($chain_counter, $chainList, $atoms) = @_;
	my (@tmp, @tmp1, $index, $curr_helix, $eff_chain);

	@tmp1 = @{ $chainList};
	
	for $index (1 .. $#tmp1) {
		push @tmp, $tmp1[$index];
	}

	push @tmp, $tmp1[0];

	$index = 0;
	while ($index <= $#tmp) {
		$curr_helix = int(($tmp[$index] - 1)/2) + 1;
		$eff_chain = ($tmp[$index] - (($curr_helix - 1) * 2));
		@tmp1 = keys %{ $atoms->{$curr_helix}{$eff_chain} };
		if ($#tmp1 == -1) {
			delete $atoms->{$curr_helix}{$eff_chain};
			splice(@tmp, $index, 1);
		} else {
			$index++;
		}
	}

	return @tmp;
}

sub setCrossovers {
	my ($cList) = @_;
	my ($b1, $b2, $c1, $c2, %myCross, $counter, $typeCross);

	for $counter (0 .. $#{ $cList }) {
		$c1 = getChain($cList->[$counter]{"1"}{"Base"}, $cList->[$counter]{"1"}{"Helix"}, $opts);
		$c2 = getChain($cList->[$counter]{"2"}{"Base"}, $cList->[$counter]{"2"}{"Helix"}, $opts);
		$b1 = $cList->[$counter]{"1"}{"Base"};
		$b2 = $cList->[$counter]{"2"}{"Base"};
		$typeCross = $cList->[$counter]{"TYPE"};

		if ($typeCross == 1) {
			$myCross{$cList->[$counter]{"1"}{"Helix"}}{$c1}{$b1} = 
			$cList->[$counter]{"2"}{"Helix"} . "," . $c2 . "," . ($b2 + 1) . ",$typeCross";
			$myCross{$cList->[$counter]{"2"}{"Helix"}}{$c2}{$b2} = 
			$cList->[$counter]{"1"}{"Helix"} . "," . $c1 . "," . ($b1 + 1) . ",$typeCross";
		} else {   
			$myCross{$cList->[$counter]{"1"}{"Helix"}}{$c1}{$b1} = 
			$cList->[$counter]{"2"}{"Helix"} . "," . $c2 . "," . $b2 . ",$typeCross";

			$myCross{$cList->[$counter]{"2"}{"Helix"}}{$c2}{($b2 + 1)} = 
			$cList->[$counter]{"1"}{"Helix"} . "," . $c1 . "," . ($b1 + 1) . ",$typeCross";
		}
	}

	return %myCross;
}

sub findTerminal {
	my ($helix, $chain, $resNum, $Cross, $atoms, $chainLen) = @_;
	my ($counter, $isDecreasing, $cM, $isCross);

	$isDecreasing = $isCross = 0;
	$counter = $resNum + 1;
	while (exists($atoms->{$helix}{$chain}{$counter})) {
		if (exists($Cross->{$helix}{$chain}{$counter}) and ! $isCross) {
			($helix, $chain, $counter,$cM) = split /\,/, $Cross->{$helix}{$chain}{$counter};
			
			if ($cM == 2 and ! $isDecreasing) { 
			$isDecreasing = 1;
			} else {
			$isDecreasing = 0;
			}
			$isCross = 1;
		} else {
			$isCross = 0;
		}
		if ($isDecreasing) {
			$counter--;
		} else {
			$counter++;
		}
	}
	
	if (! $isDecreasing) {
		$counter--;
		$isDecreasing = 1;
	} else {
		$counter++;
		$isDecreasing = 0;
	}
	
	return ($helix, $chain, $counter, $isDecreasing);
}

sub Initialize {
	my ($invalid_output_file, $overwrite);
   -e $helix1_loc or die "Error trying to access $helix1_loc\n";
   -e $helix2_loc or die "Error trying to access $helix2_loc\n";
   -e $crossovers_file or die "Error trying to access the crossover specification file $crossovers_file\n";
	
	$invalid_output_file = 0;
	if (-e $output_fle) { 
		$invalid_output_file = 1;
	};

	if ($invalid_output_file) {
		print "The file $output_fle already exists.\n";

		while ($invalid_output_file) {
			$overwrite = input ("Overwrite? [y, n]: ");
			$overwrite =~ s/\s+$//; #remove the final newline
			if ($overwrite =~ /^y|Y$/) {
				last;
			}elsif ($overwrite =~ /^n|N$/) {
				$output_fle = "";
				while (!($output_fle =~ /\w+/)) {
					$output_fle = input ("New file name: ");
					chomp($output_fle);
					$output_fle =~ s/\s+$//;
				}
				$invalid_output_file = 0;
				-e $output_fle or $invalid_output_file = 1;
				if ($invalid_output_file) {
					print "The file $output_fle already exists.\n";
				}
			}
		}	
	}	
	
	$trans0 = (
		   {
		   "XCOORD" => 0,
		   "YCOORD" => 0,
		   "ZCOORD" => 0,
		   }
		   );
}

sub registerCrossover {   
	my ($rec) = (
		 {
			 "1" => (
				 {
				 "Helix" => $_[0],
				 "Base"  => $_[1],
				 }
				 ),
			 "2" => (
				 {
				 "Helix" => $_[2],
				 "Base"  => $_[3],
				 }
				 ),
			 "TYPE" => $_[4],
		 }
		 );
	push @crossovers, $rec;
}

sub createPics {
	my ($output_fle, $chains) = @_;
	my ($pic_fle, @colors, $counter, $script_fle);
	
	@colors = ("blue", "yellow", "red", "green", "cyan", "magenta", "orange", "black");
	$script_fle = $output_fle . ".script";

	$pic_fle = $output_fle;
	$pic_fle =~ s/\.pdb$//;
	
	open SCRIPTFLE, "> $script_fle" or die "Cannot write to script file $script_fle: $!\n";
	print SCRIPTFLE "load pdb na $output_fle\n";
	print SCRIPTFLE "\nwrite pdb $output_fle\n";
	print SCRIPTFLE "write amber " . $pic_fle . "_amber.pdb\n";
	print SCRIPTFLE "# This will color the different strand of the molecule\n";
	print SCRIPTFLE "set depth_cueing on\n";
	print SCRIPTFLE "render\n";
	
	if ($chains > 8) {
		$chains = 8;
	}

	for $counter (0 .. ($chains - 1)) {
		print SCRIPTFLE "set color m1:" . ($counter + 1) . ":*:* "  . $colors[$counter] . "\n";
		print SCRIPTFLE "render\n";
	}

	print SCRIPTFLE "# Create Pictues...\n";
	print SCRIPTFLE "write png " . $pic_fle . ".png\n";
	print SCRIPTFLE "set space *\n";
	print SCRIPTFLE "render\n";
	print SCRIPTFLE "write png " . $pic_fle . "_cpk.png\n";
	print SCRIPTFLE "close\n";
	
	close SCRIPTFLE;

}

sub assignVals {
	my ($counter, $rec);

	$opts = readSpecsFile($crossovers_file);

	$is5prime = $opts->{"5prime"};
	$isParallel = $opts->{"isParallel"};
	$half_turn = $opts->{"Helix1"}{"MinorGroove"};
	$periodicity = $opts->{"Helix1"}{"MajorGroove"} + $half_turn;
	
	for $counter (@{ $opts->{"crossovers"} }) {
		registerCrossover(1, $counter->{"BASE"}, 2, $counter->{"BASE"}, $counter->{"TYPE"});
	}

}


sub numerically {
	($a<=>$b);
}

sub input {
	my ($printstring) = $_[0];
	my ($returnval);

	print "$printstring";
	$returnval = <STDIN>;
	return $returnval;
}


