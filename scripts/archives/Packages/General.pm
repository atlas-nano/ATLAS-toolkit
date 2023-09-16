package Packages::General;

require Exporter;
use strict;
use Cwd;

our (@ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = qw(IsInteger Swap CenterText STDev IsDecimal IsBoolean Trim ParseSanderInfoFile TrjSelections
GetMPSimEnergy FileTester LoadElements val CombineMols CreateGrid GetBondLength GetAngle execCmd CenterOnMol
GetTorsion CrossProduct DotProduct CRadDegrees Rotate CoP PrintProgress GetSelections CoM ParseSelect 
GetSoluteAtoms Permutate);
$VERSION = "1.00";

sub Trim {
    my ($inString) = $_[0];
                                                                                
    for ($inString) {
        s/^\s+//;
        s/\s+$//;
    }
                                                                                
    return $inString;
}

sub FileTester {
    my ($inFile) = $_[0];

    die "Error accessing regular file $inFile: $!\n"
	if (! -e $inFile or ! -r $inFile or ! -T $inFile);
    
}

                                                                             
sub IsBoolean {

    my ($inString) = $_[0];
    
    $inString =~ /^[1|0|true|false]/i ? 
	return 1 :
	return 0;

}
sub IsDecimal {

    my ($inString) = $_[0];
    
    $inString =~ /^\-?\d+\.\d+$/ ? 
	return 1 :
	return 0;

}

sub IsInteger {
    my ($inString) = $_[0];
    
    $inString =~ /^\-?\d+$/ ? 
	return 1 :
	return 0;
}

sub Swap {
    my ($val1, $val2) = @_;
    return ($val2, $val1);
}

sub CenterText {
    my ($in_text, $field_lenght) = @_;
    my ($text_length, $result);

    $text_length = $field_lenght - length($in_text) ;

    if ($text_length > 0) {
	$result = sprintf("%" . ( int($text_length/2) + ($text_length % 2) ) . "s", " ");
	$result .= sprintf("%" . length($in_text) . "s", "$in_text");
	$result .= sprintf("%" . int($text_length/2) . "s", " ");
    }else {
	$result = $in_text;
    }

    return $result;
}

sub STDev {
    my ($dataString) = $_[0];
    my ($x1, $x2, $i, $count, $sigma);

    $count = $sigma = 0;
    while ($dataString =~ /(\-?\d+\.?\d*e?\-?\d*)/g) {
	$count++;
	$x1 += $1;
	$x2 += ($1*$1);	
    }

    if ($count > 0) {
	$x1 /= $count;
	$x2 /= $count;
	$sigma = sqrt($x2 - $x1*$x1);
    }

    return ($x1, $sigma, ($x1 * $count));
}

sub STDev_old {
    my (@datavalues, $n_total, $avg, $result, $i);

    $avg = $result = 0;
    @datavalues = split /\s+/, $_[0];
    $n_total = $#datavalues;
    if ($n_total > -1) {
	$avg = 0.0;
	$result = 0.0;
	
	foreach $i (@datavalues) {
	    $avg += $i;
	}
	
	$avg = $avg/($n_total + 1);
	
	foreach (@datavalues) {
	    $result += ($_ - $avg) **2;
	}
        $n_total = 1 if ($n_total == 0);
	$result = sqrt($result/($n_total - 1));
    }
    return ($avg, $result, ($avg * $n_total));
}

sub GetInfoInput {
    my ($input_string) = $_[0];
    $input_string =~ s/\=//g;
    my (%rec, @valid_data, $counter, $index);
    my (@data_array) = split /\s+/, $input_string;
    
    @valid_data = ("TIME", "TEMP", "PRESS");
                                                                                                                   
    if ($#data_array > -1) {
        for $counter (0 ..$#data_array) {
            for $index (@valid_data) {
                if ($data_array[$counter] =~ /$index/) {
                    if ($data_array[$counter + 1] && $data_array[$counter + 1] =~ /(\d+\.\d+)/) {
                        $rec{$index} = $1;
                    }
                }
            }
        }
    }
                                                                                                                   
    return \%rec;
}
                                                                                                                   

sub ParseSanderInfoFile {
    my ($infile) = $_[0];
    
    my ($rec, %OutData, $valid_info);
    open INFILE, $infile or die "Cannot open $infile: $!\n";
    while (<INFILE>) {
        if ($_ =~ /R M S  F L U C T U A T I O N S/) {
            last;
        }elsif ($_ =~ /NSTEP\s+\=\s*(\d+)\s+(.+)/) {
            $rec = GetInfoInput($2);
            if ($rec) {
                $OutData{$1} = $rec;
                $valid_info = 1;
            }
        }
    }
    close INFILE; 
    
    die "No valid data found in $infile\n"
	if (! $valid_info);

    return \%OutData;
}

sub GetMPSimEnergy {
    my ($curr_fle, $total_bases) = @_;
    my ($data_start, $in_data, $unit_nm, $counter, $eng);
    my ($curr_val, $curr_total, %Eng, @eng_header);

    if (open INPFILE, $curr_fle) {
	while (<INPFILE>) {
	    chomp;
	    $in_data = $_;
	    if ($in_data =~ /^\s*[Atom|GROUP]/ && $#eng_header <= 0) {
		while ($in_data =~ /\s+(\w+)/g) {
		    if ( ($1 !~ /GROUP/i) && ($1 !~ /TOTAL/i) ) {
			push @eng_header, $1;
		    }
		}
	    } elsif ($in_data =~ /^\s*(\d+)\s+/ && $#eng_header > 0) {
		$unit_nm = $1;

		if ($unit_nm == $total_bases +1) {
		    last;
		} else {
		    $counter = 0;
		    $curr_val = $curr_total = 0;
		    while ($in_data =~ /(\-?\d+\.\d+)\s*/g && $counter <= $#eng_header) {
			$eng = $eng_header[$counter];
			$curr_val = $1;
			$curr_total += $curr_val;
			$Eng{$unit_nm}->{$eng} = $curr_val; 
		$counter++;
		    }				
		    $Eng{$unit_nm}->{"TOTAL"} = $curr_total;
		    $curr_val = $curr_total = 0;
		    } 
		}
	    }
	close INPFILE;
	}
	
    die "Error when reading energy file: $curr_fle: $!\n"
	if (! ($#eng_header > 0));
    
    return \%Eng;
}

sub LoadElements {
    my (%ELEMENTS, $indata, $eleNum, $myDir);
    $myDir = "/ul/tpascal/scripts/Packages";
    my ($datFile) = "$myDir/elementList.txt";
    FileTester($datFile);

    open INDATA, $datFile or die "Cannot open file $datFile: $!\n";
    while (<INDATA>) {
	chomp;
	$indata = $_;
	if ($indata =~ /^(\d+)\s+(\*?)\s*\d+/) {
	    $eleNum = $1;
	    if ($2) {
		$ELEMENTS{$1}{"NATURAL"} = 0;
	    } else {
		$ELEMENTS{$1}{"NATURAL"} = 1;
	    }
	    #$indata = $'; 
            if ($indata =~ /^\d+\s+\*?\s*(\d+\.?\d*)\s+(\w+)\s+(\w+)/) {
		$ELEMENTS{$eleNum}{"MASS"} = $1;
		$ELEMENTS{$eleNum}{"NAME"} = $2;
		$ELEMENTS{$eleNum}{"SYMBOL"} = $3;
	    } else {
		delete $ELEMENTS{$eleNum};
	    }
	}
    }

    close INDATA;

    die "ERORR: No valid data found in file $datFile\n"
	if (! %ELEMENTS);

    return \%ELEMENTS;

}

sub val {
    my ($inStr) = $_[0];
    my ($returnStr);

    if (IsInteger($inStr) or IsDecimal($inStr)) {
	$returnStr = $inStr;
    } elsif ($inStr =~ /^(-?\d+\.*\d*)E(\d+)/) {
	$returnStr = $1 * 10**$2;
    } elsif ($inStr =~ /^(-?\d+\.*\d*)E(\+\d+|\-\d+)/) {
	$returnStr = $1 * 10**$2;
    } else {
	$returnStr = $inStr;
    }

    return $returnStr;
}

sub CombineMols {
    my ($mol1, $mol2, $CONN, $connections) = @_;
    my ($atom, $tot_atoms, $tot_res, @tmp, %ATOMS, @BONDS, %CONS, $bond);

    @tmp = sort Numerically keys %{ $mol1 };
    $tot_atoms = $tmp[$#tmp];
    $tot_res = $mol1->{$tot_atoms}{"RESNUM"};
    @tmp = sort Numerically keys %{ $mol2 };

    %ATOMS = %{ $mol1 };
    %CONS = %{ $CONN };
    for $atom (@tmp) {
        %{ $ATOMS{($atom + $tot_atoms)} } = %{ $mol2->{$atom} };
        $ATOMS{($atom + $tot_atoms)}{"RESNUM"} += $tot_res;
	@BONDS = ();
        @BONDS = @{ $connections->{$atom} } if (defined($connections->{$atom}));
        for $bond (@BONDS) {
            push @{ $CONS{($atom + $tot_atoms)} }, ($bond + $tot_atoms);
        }
    }

    return (\%ATOMS, \%CONS);


}

sub GetBondLength(@) {
    my ($bead_1, $bead_2) = @_;
    my ($dist);

    $dist = ($bead_1->{"XCOORD"} - $bead_2->{"XCOORD"})**2;
    $dist += ($bead_1->{"YCOORD"} - $bead_2->{"YCOORD"})**2;
    $dist += ($bead_1->{"ZCOORD"} - $bead_2->{"ZCOORD"})**2;
    $dist = sqrt($dist);

    return $dist;

}

sub GetAngle(@) {
    my ($bead_1, $bead_2, $bead_3, $saveRad) = @_;
    my ($vec_1, $vec_2, $angle, $counter);

    for $counter ("XCOORD", "YCOORD", "ZCOORD") {
        $vec_1->{$counter} = $bead_1->{$counter} - $bead_2->{$counter};
        $vec_2->{$counter} = $bead_3->{$counter} - $bead_2->{$counter};
    }

# Now to calculate the dot product between the two vectors, us the dot product
# A.B = |A||B|cos(O)= sum(AiBi)

    $angle = DotProduct($vec_1, $vec_2);
    $angle = acos($angle);

    if (defined($saveRad) and $saveRad eq "1") {
	return $angle;
    } else {
	return CRadDegrees($angle, 1);
    }
}

sub GetTorsion(@) {
    my ($b1, $b2, $b3, $b4, $saveRad) = @_;
    my ($counter, $v1, $v2, $v3, $v4, $v5);
    my ($normal_a, $normal_b, $answer, $sign);
    my ($pi) = atan2(1,1) *4;

    for $counter ("XCOORD", "YCOORD", "ZCOORD") {
        $v1->{$counter} = $b3->{$counter} - $b2->{$counter};
        $v2->{$counter} = $b4->{$counter} - $b3->{$counter};
        $v3->{$counter} = $b2->{$counter} - $b3->{$counter};
        $v4->{$counter} = $b1->{$counter} - $b2->{$counter};
        $v5->{$counter} = $b3->{$counter} - $b4->{$counter};
    }

    $normal_a = CrossProduct($v2, $v1);
    $normal_b = CrossProduct($v3, $v4);

    $answer = DotProduct($normal_a, $normal_b);
    $sign = DotProduct($normal_b, $v5);
    $answer = acos($answer);
    if ($sign < 0) {
	$answer = (2 * $pi) - $answer;
    }

    if (defined($saveRad) and $saveRad eq "1") {
	return $answer;
    } else {
	$answer = CRadDegrees($answer, 1);

	return $answer;
    }

}

sub CrossProduct(@) {
    my ($n1, $n2) = @_;

    my ($result);

    $result->{"XCOORD"} = ($n1->{"YCOORD"} * $n2->{"ZCOORD"}) - ($n1->{"ZCOORD"} * $n2->{"YCOORD"});
    $result->{"YCOORD"} = ($n1->{"ZCOORD"} * $n2->{"XCOORD"}) - ($n1->{"XCOORD"} * $n2->{"ZCOORD"});
    $result->{"ZCOORD"} = ($n1->{"XCOORD"} * $n2->{"YCOORD"}) - ($n1->{"YCOORD"} * $n2->{"XCOORD"});

    return $result;

}

sub DotProduct(@) {
    my ($n1, $n2) = @_;
    my ($result, $counter, $norm_a, $norm_b);

    for $counter ("XCOORD", "YCOORD", "ZCOORD") {
        $result += ($n1->{$counter} * $n2->{$counter});
        $norm_a += $n1->{$counter}**2;
        $norm_b += $n2->{$counter}**2;
    }

    $result = ($result/(sqrt($norm_a) * sqrt($norm_b)));

    return $result;
}

sub CRadDegrees(@) {

    my ($convertToDegrees) = $_[1];
    my ($inital_angle) = $_[0];
    my ($resulting_angle) = 0.0;

    my ($pi) = atan2(1,1) *4;

    if ($convertToDegrees) { $resulting_angle = $inital_angle * 180 / $pi; }
    else { $resulting_angle = $inital_angle * $pi/180; }

    return $resulting_angle;
}

sub acos(@) {
    my($x) = $_[0];
    if (abs($x) > 1.0) {
        return 0;
    } else {
        return atan2(sqrt(1 - $x * $x), $x);
    }
}


sub Numerically {
    ($a<=>$b);
}

sub Rotate {
   my ($atoms, $rot_array, $coord) = @_;
   my (@rotation_matrix, $counter, $index, $matrix_keys, %F_Data, $i, $offset);
   my (@sinA, @cosA);

   for $counter (keys %{ $atoms }) {
       for $index (keys %{ $atoms->{$counter} }) {
	   $F_Data{$counter}{$index} = $atoms->{$counter}{$index};
       }
   }

   for $counter (0 .. $#{ $rot_array }) {
        $sinA[$counter] = sin($rot_array->[$counter]);
	$cosA[$counter] = cos($rot_array->[$counter]);
   }

   @rotation_matrix = (
                       {
                           "XCOORD" => [1, 0, 0],
                           "YCOORD" => [0, $cosA[0], $sinA[0]],
                           "ZCOORD" => [0, -$sinA[0], $cosA[0]],
                        },

                       {
                           "XCOORD" => [$cosA[1], 0, -$sinA[1]],
                           "YCOORD" => [0, 1, 0],
                           "ZCOORD" => [$sinA[1], 0, $cosA[1]],
                       },
                       {
                           "XCOORD" => [$cosA[2], $sinA[2], 0],
                           "YCOORD" => [-$sinA[2], $cosA[2], 0],
                           "ZCOORD" => [0, 0, 1],
                       },
		       {
			   "XCOORD" => [$cosA[1]*$cosA[2], $sinA[0]*$sinA[1]*$cosA[2]+$cosA[0]*$sinA[2], -$cosA[0]*$sinA[1]*$cosA[2]+$sinA[0]*$sinA[2]],
                           "YCOORD" => [-$cosA[1]*$sinA[2], -$sinA[0]*$sinA[1]*$sinA[2]+$cosA[0]*$cosA[2], $cosA[0]*$sinA[1]*$sinA[2]+$sinA[0]*$cosA[2]],
                           "ZCOORD" => [$sinA[1], -$sinA[0]*$cosA[1], $cosA[0]*$cosA[1]],
		       },
                       );

   for $matrix_keys (keys %{ $atoms }) {
       $counter = \%{ $rotation_matrix[$coord] };
       for $index ("XCOORD", "YCOORD", "ZCOORD") {
	   $atoms->{$matrix_keys}{$index} = ($F_Data{$matrix_keys}{"XCOORD"} * $counter->{$index}[0] +
					    $F_Data{$matrix_keys}{"YCOORD"} * $counter->{$index}[1] +
					    $F_Data{$matrix_keys}{"ZCOORD"} * $counter->{$index}[2]);
       }
   }
   @rotation_matrix = ();

}

sub CenterMol {
   my ($atoms) = @_;
   my ($CENTER, $atomC, $dim);

    $CENTER = CoP($atoms);
    for $atomC (keys %{ $atoms }) {
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $atoms->{$atomC}{$dim} -= $CENTER->{$dim};
	}
   }

   return $CENTER;
}

sub CoP {
    my ($atoms) = @_;
    my (%CENTER, $atomC, $tot, $dim);
    
   $tot = 0;
   for $atomC (keys %{ $atoms }) {
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $CENTER{$dim} += $atoms->{$atomC}{$dim};
        }
        $tot++;
   }

   for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$CENTER{$dim} /= $tot;
   }

    return \%CENTER;
}

sub PrintProgress {
    my ($currPos, $total, $start, $pStr) = @_;
    my ($progress, $str, $end);
    
    $end = time();

    $progress = $currPos/$total;
    
    $str = sprintf("%.2f%% complete %s\r", 
		   100*$progress, getEta(($end - $start), $progress));
    
    print "${pStr}${str}" if (defined($pStr));
    return length($str);
}

sub getEta {
    my ($elapsed, $percentage) = @_;
    my ($totalTime) = $elapsed/$percentage;
    my ($timeLeft) = $totalTime - $elapsed;
    my ($returnStr) = "(";
    

    if ($timeLeft > 60) {
	if ($timeLeft > 3600) {
	    if ($timeLeft > 86400) {
		$returnStr = int($timeLeft/86400) . "d ";
		$timeLeft = $timeLeft % 86400;
	    }
	    $returnStr .= int($timeLeft/3600) . "h ";
	    $timeLeft = $timeLeft % 3600;
	}
	$returnStr .= int($timeLeft/60) . "m ";
	$timeLeft = $timeLeft % 60;
    }

    $returnStr .= int($timeLeft) . "s remaining)          ";
    
    return $returnStr;
}

sub GetSelections {
    my ($opts, $isPrint) = @_;
    my ($index, %SELECT, $saveKey, $counter, $exclude, %EXCLUDED, $interval);

    $isPrint = 1 if (! defined($isPrint));

    for $index (@{ $opts }) {
	$index = Trim($index);
        $exclude = 0;
	if ($index =~ /^\^/) {
	    $exclude = 1;
	    $index =~ s/^\^//;
        }
	if ($index =~ /^([A-Z])([A-Z])(\w+|\*)/i) {
	    $saveKey = getOptName($1, $2);
	    if (defined($saveKey)) {
		print "SELECTED $saveKey $3\n" if ($isPrint);
		$SELECT{$saveKey}{$3} = 1;
		$EXCLUDED{$saveKey}{$3} = 1 if ($exclude);
	    }
	} elsif ($index =~ /^\:([A-Z])([A-Z])(\d+)\-(\d+)\:?(\d*)/i) {
	    $saveKey = getOptName($1, $2);
	    $interval = 1;
	    $interval = $5 if ($5 ne "");
	    if (defined($saveKey)) {
		print "SELECTED $saveKey $3 - $4 every $interval\n" if ($isPrint);
		$counter = $3;
		while ($counter <= $4) {
		    $SELECT{$saveKey}{$counter} = 1;
		    $EXCLUDED{$saveKey}{$counter} = 1 if ($exclude);
		    $counter += $interval;
		}
	    }
	}
    }

    for $saveKey (keys %EXCLUDED) {
	for $counter (keys %{ $EXCLUDED{$saveKey} }) {
	    delete $SELECT{$saveKey}{$counter};
	}
    }
    die "ERROR: No valid options found!\n"
	if (! %SELECT);

    return \%SELECT;
}

sub getOptName {
    my ($optField, $optType) = @_;
    my ($saveKey);

    if (uc($optField) eq "N") {
	if (uc($optType) eq "A") {
	    $saveKey = "ATMNAME";
	} else {
	    $saveKey = "RESNAME";
	} 
    }elsif ($optField eq "T" and uc($optType) eq "A") {
	$saveKey = "FFTYPE";
    } else {
	if (uc($optType) eq "A") {
	    $saveKey = "INDEX";
	} elsif (uc($optType) eq "R") {
	    $saveKey = "RESNUM";
	} elsif (uc($optType) eq "M") {
	    $saveKey = "MOLSIZE";
	} else {
	    $saveKey = "NUMBONDS";
	}
    }
    
    return $saveKey;
}
    
sub CoM {
    my ($atoms) = $_[0];
    my ($i, %COM, $dim, $totalMass, $mass);
    
    $totalMass = 0;
    for $i (keys %{ $atoms }) {
	$mass = 1;
        $mass = $atoms->{$i}{"MASS"} if (exists($atoms->{$i}{MASS}));
        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $COM{$dim} += $atoms->{$i}{$dim} * $mass;
        }
        $totalMass += $mass;
    }
    
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
        $COM{$dim} /= $totalMass;
    }
   
   return \%COM;
}

sub ParseSelect {
    my ($selection) = $_[0];
    my ($i, %SELECT);

    print "SNAPSHOT SELECTION: ";
    if (defined($selection) && Trim($selection) =~ /^(\d+)\s+(\d+)\s+(\d+)/) {
        print "start: $1 end: $2 every: $3...";
        $i = $1;
        while ($i <= $2) {
            $SELECT{$i} = 1;
            $i += $3;
        }
    } else {
        %SELECT = ();
        print "using all";
    }
    print "\n";

    return \%SELECT;
}

sub execCmd {
    my ($cmd, $fileToDelete) = @_;

    die "Error while executing $cmd\n" if (system($cmd));
    system "rm -f $fileToDelete" if (defined($fileToDelete));
}

sub TrjSelections {
    my ($select) = $_[0];
    my (@tmp, $i, %SELECTION, $j, $counter, @tmp2);

    @tmp = split /\s+/, $select;

    $i = 0;
    $counter = 0;
    while ($i <= $#tmp) {
        if ($tmp[$i] =~ /^\:It(\d+)\-(\d+)\:(\d+)$/) {
            $j = $1;
            while ($j <= $2) {
                $SELECTION{$j} = 1;
                $counter++;
            	$j += $3;
            }
        } elsif ($tmp[$i] =~ /^(\d+)$/) {
            $SELECTION{$1} = 1;
            $counter++;
        } elsif (-e $tmp[$i] && -r $tmp[$i] && -T $tmp[$i]) { # if the input is a file
	    if (open INDATA, $tmp[$i]) {
		while (<INDATA>) {
		    if ($_ =~ /^\"(.+)\"/) {
			@tmp2 = split /\s+/, $1;
			for $j (@tmp2) {
			    push @tmp, $j;
			}
		    }
		}
		close INDATA;
	    }
	}
	$i++;
    }

    die "ERROR: No valid selection found. See help file. Got $select\n" if ($select ne "*" && (! %SELECTION || ! keys %SELECTION));

    print "trajectory: ";
    if ($select eq "*") {
        print "selected all frames...";
    } else {
        print "selected $counter frames...";
    }
    return \%SELECTION;
}

sub CenterOnMol {
    my ($atoms, $CENTER) = @_;
    my ($i, @tmp, $j);

    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (keys %{ $atoms }) {
	for $j (@tmp) {
	    $atoms->{$i}{$j} -= $CENTER->{$j};
	}
    }
}

sub GetSoluteAtoms {
    my ($atoms, $moleculeList) = @_;
    my (%SOL, $i, $atm, @tmp);

    for $i (keys %{ $moleculeList }) {
	@tmp = keys %{ $moleculeList->{$i} };
	$atm = pop @tmp;
	next if ($atoms->{$atm}{RESNAME} =~ /WAT|Na|Mg|Cl/i || $atoms->{$atm}{FFTYPE} =~ /OF3C|Cl|Na|Mg|HF3C|OW|HW/i);
	$SOL{$atm} = 1;
	for $atm (@tmp) {
	    $SOL{$atm} = 1;
	}
    }

    return \%SOL;
}

sub Permutate {
    my @head = @{ $_[0] };
    my @tail = @{ $_[1] };
    my (@RET);
    
    @RET = ();
    unless (@head) {
	# stop recursing when there are no elements in the head
	return \@tail;
    } else {
	# for all elements in @head, move one from @head to @tail
	# and call permut() on the new @head and @tail
	my(@newhead,@newtail,$i);
	foreach $i (0 .. $#head) {
	    @newhead = @head;
	    @newtail = @tail;
	    unshift(@newtail, splice(@newhead, $i, 1));
	    push @RET, Permutate([@newhead], [@newtail]);
	}
    }

    return @RET;
}

1;
