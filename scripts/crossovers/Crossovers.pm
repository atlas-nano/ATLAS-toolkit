package Crossovers;
require Exporter;

our (@ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = qw(findAtom getDistance readSpecsFile readNamotPDB writeData getChain Trim cRadDegrees reverseRes isInteger); 

$VERSION = "1.00";

sub findAtom {
    my ($res, $atomName) = @_;
    my ($result, $aC, $atom, $foundAtom);

    $result = $foundAtom = 0;

    for $aC (keys %{ $res }) {
	$atom = \%{ $res->{$aC} };
	if (Trim(lc($atom->{"NAME"})) eq lc($atomName)) {
	    $foundAtom = $aC;
	    last;
	}
    }

    return $foundAtom;
}

sub getDistance {
    my ($atom1, $atom2) = @_;
    my ($dist);

    $dist = 0;
    for ("XCOORD","YCOORD","ZCOORD") {
	$dist += ($atom1->{$_} - $atom2->{$_})**2;
    }

    $dist = sqrt($dist);

    return $dist;
}

sub readSpecsFile {
    my ($cfile) = $_[0];
    my (%specs, $validCounter, $myKey, $isValid, $line_in, $rec);

    $validCounter = $isValid = 0;
    open CROSSFILE, $cfile or die "Error while reading $cfile: $!\n";
    while (<CROSSFILE>) {
	chomp;
	$line_in = $_;
	if ($line_in =~ /^Helix (\d+)\s+(\d+):(\d+)\s+(\d+)\s+bases\s+1st Minor:\s+(\d+)/) {
	    $myKey = "Helix" . $1;
	    $specs{$myKey} = (
			      {
				  "MajorGroove"   => $2,
				  "MinorGroove"   => $3,
				  "Periodicity"   => ($2 + $3),
				  "TotalBases"    => $4,
				  "Minor1"        => $5,
			      }
			      );
	    $validCounter++;
	} elsif ($line_in =~ /^Transalation: (\S+)\s+(\S+)\s+(\S+)/) {
	    $specs{"Trans"} = (
			       {
				   "XCOORD" => $1,
				   "YCOORD" => $2,
				   "ZCOORD" => $3,
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

sub readNamotPDB {
    my ($filenm, $AtomData, $helix, $transVec) = @_;
    my ($inStr, $curr_chain, @tmp, $curr_res, $rec, $counter, $atomCounter);

    $counter = $curr_chain = 1;
    $curr_res = -1;
    $rec = ();
    $atomCounter = 0;
    open NAMOTFILE, $filenm or die "ERROR: Cannot open $filenm: $!\n";
    while (<NAMOTFILE>) {
	chomp;
	$inStr = $_;
	if ($inStr =~ /^ATOM\s+(\d+)(.{5})(.{4})\s+(\w)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $atomCounter++;
	    if ($5 != $curr_res and $curr_res > -1) {
		$AtomData->{$helix}{$curr_chain}{$curr_res} = $rec;
		$rec = ();
		$counter = 1;
		$curr_res = $5;
	    } elsif ($curr_res == -1) {
		$curr_res = $5;
	    }
	    
	    $rec->{$counter} = (
				{
				    "NAME"    => $2,
				    "RESNAME" => $3,
				    "CHAIN"   => $4,
				    "RESNUM"  => $5,
				    "XCOORD"  => $6,
				    "YCOORD"  => $7,
				    "ZCOORD"  => $8,
				    "ATOMNUM" => $1,
				}
				);
	    for ("XCOORD","YCOORD","ZCOORD") {
		$rec->{$counter}{$_} += $transVec->{$_};
		$AtomData->{("HELIX" . $helix)}{$_} += $rec->{$counter}{$_};
	    }

	    $counter++;
	} elsif ($inStr =~ /^TER/) {
	    $AtomData->{$helix}{$curr_chain}{$curr_res} = $rec;
	    $rec = ();
	    $curr_res = -1;
	    $curr_chain++;
	}
    }

    close NAMOTFILE;

    die "ERROR: No atom information found while reading Namot2 PDB file $filenm\n"
	if ($atomCounter == 0);
    for ("XCOORD","YCOORD","ZCOORD") {
	$AtomData->{("HELIX" . $helix)}{$_} /= $atomCounter;
    }
    
    
}

sub writeData {
    my ($HelixData, $save_name, $transVec) = @_;
    my ($atom, $counter, $fmt, $index, $chain, $chainName);
    my ($outData, $header, $fmt2, $res, $resName, $rn, $rc, $atoms);

    $header = "HEADER                                 NUCLEIC ACID             0\n" .
	"TITLE                                            NAMOT Generated MODEL\n" .
	"KEYWDS                                                           MODEL\n";
    $fmt = "%-4s%7d%5s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f\n";
    $fmt2 = "%-4s%7d%9s%2s%4d\n";
    $index = 1;
    
    for $chain (sort Numerically keys %{ $HelixData }) {
	$chainName = lc(chr(64 + $chain ));
	for $res (sort Numerically keys %{ $HelixData->{$chain} }) {
	    for $atoms (sort Numerically keys %{ $HelixData->{$chain}{$res} }) {
		$atom = \%{  $HelixData->{$chain}{$res}{$atoms} };
		for ("XCOORD", "YCOORD", "ZCOORD") {
		    $atom->{$_} -= $transVec->{$_};
		}

		$resName = $atom->{"RESNAME"};
		$outData .= sprintf($fmt, "ATOM", $index, $atom->{"NAME"}, $resName, $chainName, $res,
				    $atom->{"XCOORD"}, $atom->{"YCOORD"}, $atom->{"ZCOORD"}, 0,0);
		$index++;
		$rn = $resName;
	    }
	    $rc = $res;
	}
	$outData .= sprintf($fmt2, "TER", $index, $rn, $chainName, $rc);
	$index++;
    }
    
    if (defined($outData)) {
	$outData = $header . $outData . "ENDMDL\n";

	open OUTDATA, "> $save_name" or die "Cannot write to file $save_name: $!\n";
	print OUTDATA $outData;
	close OUTDATA;
    } else {
	die "ERROR: No data written\n";
    }
}
	    

sub getChain {
# GetChain - gets the chain of a molecule
    my ($curr_group, $whichHelix, $opts) = @_;
    my ($returnval, $minorG, $is5prime, $i, $baseC, $totBases); 
    my ($isParallel, $Minor1, %CrossChain, $majorG, $chain);

    $minorG = $opts->{"Helix1"}{"MinorGroove"};
    $majorG = $opts->{"Helix1"}{"MajorGroove"};
    $is5prime = $opts->{"5prime"};
    $isParallel = $opts->{"isParallel"};
    $Minor1 = $opts->{"Helix" . $whichHelix}{"Minor1"};
    $totBases = $opts->{"Helix" . $whichHelix}{"TotalBases"};
    
    $chain = 1; # specifies which of the chains is facing "inside" i.e. closest to crossover point
    if ((! $is5prime and $whichHelix == 1) or ($whichHelix == 2 and ! $isParallel and $is5prime)) {
	$chain = 2;
    }
	    
    for $i (1 .. $Minor1) {
	$CrossChain{$i} = $chain;
    }
    $baseC = $Minor1;
    if ($chain == 1) {
	$chain = 2;
    } elsif ($chain == 2) {
	$chain = 1;
    }

    while ($baseC < $totBases) {
	for $i (1 .. $minorG) {
	    $CrossChain{($i + $baseC)} = $chain;
	}
	$baseC += $minorG;
	if ($chain == 1) {
	    $chain = 2;
	} elsif ($chain == 2) {
	    $chain = 1;
	}

	for $i (1 .. $majorG) {
	    $CrossChain{($i + $baseC)} = $chain;
	}
	$baseC += $majorG;
	if ($chain == 1) {
	    $chain = 2;
	} elsif ($chain == 2) {
	    $chain = 1;
	}
    }
#    print "BASE: $curr_group Chain $CrossChain{$curr_group}\n";
    return $CrossChain{$curr_group};
}

sub Trim {
    my ($inStr) = $_[0];
 
    $inStr =~ s/^\s+//;
    $inStr =~ s/\s+$//;

    return $inStr;
}

sub cRadDegrees {
                                                                                                                  
    my ($inital_angle, $convertToDegrees) = @_;

    my ($resulting_angle) = 0.0;                                                                      
    my ($pi) = atan2(1,1) *4;
                                                                                                                 
    if ($convertToDegrees) { $resulting_angle = $inital_angle * 180 / $pi; }
    else { $resulting_angle = $inital_angle * $pi/180; }
                                                                                                                  
    return $resulting_angle;
}

sub reverseRes {
    my ($chainData) = $_[0];
    my ($resNum, $atom, $counter, @tmp, %tmpHash);

    @tmp = sort Numerically keys %{ $chainData };

    for $counter (0 .. $#tmp) {
	$resNum = ($tmp[$#tmp] - $tmp[$counter]) + 1;
	%{ $tmpHash{($resNum)} } = %{ $chainData->{$tmp[$counter]} };
    }

    return \%tmpHash;
}


sub Numerically {
    ($a<=>$b);
}

sub isInteger {
    my ($test_str) = $_[0];

    if ($test_str =~ /^\d+/) {
	return 1;
    } else {
	return 0;
    }
}

    
