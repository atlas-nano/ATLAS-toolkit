#!/usr/bin/perl -w
# This script will take a waterbox and insert ions at various points and place a residue dimer water shell apart

BEGIN {
    push (@INC, "/ul/tpascal/scripts");
}
use strict;
use File::Basename;
use Packages::General;
use Packages::FileFormats;
use Packages::BOX;

sub GetParms;
sub ParsePDBFile;
sub GetVals;
sub Numerically;
sub InsertRes;
sub GetResInfo;
sub AddResInfo;
sub CheckOverLap;
sub GetRandomPos;
sub GetMaxDisplace;
sub writePDB;
sub GenerateInsertions;
sub RemoveWaters;
sub RemoveResidue;
sub PlaceResInWatShells;
sub GetSmallestBoxLen;
sub DeterminePlacement;
sub DetermineIonPlacement;
sub GetBoxCenter;
sub CopyHash;
sub EqualizeWaters;
sub SaveFiles;
sub CreateParmScript;
sub DoCmd;

die "usage: $0 paramfile [save_prefix] [is_glucose]\n"
    if (! @ARGV);

my ($paramfile, $savePrefix, $is_glucose) = @ARGV;

my ($PARMS) = GetParms($paramfile);
my ($Atom_Info, $CON, $HEADERS) = GetBGFFileInfo($PARMS->{"BGF_FILE"}, 1);
my ($BOX) = GetBox($Atom_Info, $HEADERS);

$savePrefix = "out"
    if (! defined($savePrefix));

if ($PARMS->{"NumInsertions"} and IsInteger($PARMS->{"NumInsertions"})) {
    print "--====Generating " . $PARMS->{"NumInsertions"} . " structures with different inserions points===--\n";
    GenerateInsertions;
    print "--===Done Insertions===--\n";
}

print "--===Placing " . $PARMS->{"Name"} . " in water shells===--\n";
PlaceResInWatShells($BOX, $PARMS, $Atom_Info);
print "--===Done Water Shells===--\n";

sub PlaceResInWatShells {
    my ($BBOX, $MyParms, $ATOMS) = @_;
    my ($counter, $sLim, $numShells, $index, $tmp1, $tmp2, $insertpoint, $PLACEMENT);
    my ($total_atoms, $total_res, $num_ions, $save_name, $ioninsert, $IONPLACEMENT);
    my ($ionCounter, $numions, @tmp, %FILEDATA, $fileKey, $fdata);
    # This sub will attempt to place the Residue consecutive water shells apart from the center of the box
    # The number of residues is determined the the TotalRes keyword in the input file
    # This should not exceed 8, since the algorithm will attempt to place the residues at the 8 corners
    # of a cubic box, choosing the diagonal placement first to maximize distances
    # Residues will then be placed consecutive watershells apart until 2/3 of the radial box radius is reached
        
    # Determine residue placement in corner of cubic box for each sucessive watershell
    $sLim = (0.66 * (GetSmallestBoxLen($BBOX)/2)) - $MyParms->{"ResRadii"};
    $PLACEMENT = DeterminePlacement($sLim, $MyParms->{"ResRadii"}, $MyParms->{"WaterRadii"}, $BBOX);

    # Determine the placement of the ions. Should be as far from residue as possible to minimize
    # Interaction. Place face-centered on consecutively smaller cubes
    if ($MyParms->{"HasIons"}) {
	$sLim = (0.66 * (GetSmallestBoxLen($BBOX)/2)) - $MyParms->{"IonRadii"};
	if ($MyParms->{"WaterRadii"} > $MyParms->{"IonRadii"}) {
	    $IONPLACEMENT = DetermineIonPlacement($BBOX, $sLim, $MyParms->{"WaterRadii"});
	} else {	    
	    $IONPLACEMENT = DetermineIonPlacement($BBOX, $sLim, $MyParms->{"IonRadii"});
	}
    }

    # Loop over numShells, placing TotalRes residues in corners of a cubic box: length: 2 * numShell + ResRadii
    $total_atoms = $BBOX->{"TOTAL_ATOMS"};
    $total_res = $BBOX->{"TOTAL_RES"};
    $ionCounter = 0;
    @tmp = keys %{ $PLACEMENT };
    $numShells = $#tmp;
    undef @tmp;
    for $counter (0 .. $numShells) {
	print "Water Shell #" . ($counter + 1) . "\n";
	%{ $tmp1 } = %{ $ATOMS };
	for $index (1 .. $MyParms->{"TotalRes"}) {
	    undef $insertpoint;
	    $insertpoint = $PLACEMENT->{$counter}{($index - 1)};
	    $tmp1 = InsertRes($MyParms->{"RESIDUE"}, $insertpoint, $tmp1, 0,
			      ($total_atoms + 1), ($total_res + 1));
	    $total_atoms += $MyParms->{"AtomsInRes"};
	    $total_res++;

	    print "\tInserted " . $MyParms->{"Name"} . " at " . $tmp1->{$total_atoms}{"XCOORD"} . 
		" " . $tmp1->{$total_atoms}{"YCOORD"} . " " . $tmp1->{$total_atoms}{"ZCOORD"} . "\n";
	}
	$save_name = $savePrefix . "_watShell" . ($counter + 1);
	
	if ($MyParms->{"HasIons"} and IsInteger($MyParms->{"HasIons"})) {
	    $numions = $MyParms->{"TotalRes"} * abs($MyParms->{"ResCharge"} / $MyParms->{"IonCharge"});
	    $save_name .= "_wIon";
	    %{ $tmp2 } = %{ $tmp1 };
	    for $counter (1 .. $numions) {		
		undef $ioninsert;
		if (! $IONPLACEMENT->{$ionCounter}) {
		    $ioninsert = GetRandomPos($BBOX, $MyParms->{"IonRadii"});
		} else {
		    $ioninsert = $IONPLACEMENT->{$ionCounter};
		}

		$tmp2 = InsertRes($MyParms->{"IONS"}, $ioninsert, $tmp2, 1, 
				  ($total_atoms + 1), ($total_res + 1));
		if ($tmp2->{"-1"}) {
		    $tmp2->{"-1"} -= 1;
		}
		$ionCounter += 1;
		$total_atoms += 1;
		$total_res += 1;

		print "\tInserted " . $MyParms->{"IonName"} . " at " . $tmp2->{$total_atoms}{"XCOORD"} . 
		    " " . $tmp2->{$total_atoms}{"YCOORD"} . " " . $tmp2->{$total_atoms}{"ZCOORD"} . "\n";
	    }
	    $fileKey = 2 * $counter + 1;
	    $FILEDATA{$fileKey} = (
				   {
				       "DATA" => $tmp2, 
				       "NAME" => $save_name,
				   }
				   );
	    undef $tmp2;
	}

	if (substr($save_name, -5) eq "_wIon") {
	    $save_name = substr($save_name, 0, -5);
	}
	$fileKey = 2 * $counter;
	$FILEDATA{$fileKey} = (
			       {
				   "DATA" => $tmp1, 
				   "NAME" => $save_name,
			       }
			       );
	undef $tmp1;
    }	

    print "\t--==Creating PDB files==--\n";
    $fdata = EqualizeWaters(\%FILEDATA);
    SaveFiles($fdata);
    print "Creating Parameter and Script Files...";
    CreateParmScript($fdata, $MyParms->{"TotalRes"}, $MyParms->{"HasIons"}, $numions);
    print "Done\n";
    undef $fdata;
    undef %FILEDATA;
	
}

sub CreateParmScript {
    my ($fData,  $totRes, $hasIons, $totIons) = @_;
    my ($hashKey, $counter, $savePrefix, $myCmd, $leap_loc);


    for $hashKey (keys %{ $fData }) {
	$savePrefix = $fData->{$hashKey}{"NAME"};

	# Create the .top and .crd files
	if (! defined($is_glucose) or ! IsInteger($is_glucose) or $is_glucose != 1) {
	    $is_glucose = 0;
	    $leap_loc = "";
	} else {
	    $leap_loc = "/ul/tpascal/amber8/leaprc.glucose";
	}

	$myCmd = "/ul/tpascal/scripts/createSanderInp.pl " . $fData->{$hashKey}{"NAME"} . ".pdb " . 
	$fData->{$hashKey}{"NAME"} . " $is_glucose 0 $leap_loc >& _junk";
	DoCmd($myCmd);

	# Create the cluster files
	$counter = $totRes;
	if (($hashKey % 2) == 1 and $hasIons) {
	    $counter += $totIons;
	}
	$myCmd = "/ul/tpascal/scripts/createAmberscript.pl $savePrefix" . ".top " . $savePrefix . ".crd $counter $savePrefix >& _junk";
	DoCmd($myCmd);

    }

}

sub DoCmd {
    my ($sys_cmd) = $_[0];
                                                                                                                   
    if (system($sys_cmd)) {
        die "Error executing command: $sys_cmd\n";
    }
}
                                                                                                                   

sub EqualizeWaters {
    my ($filesData) = $_[0];
    my ($counter, $hashKey, $totRemoved);
    
    $counter = $totRemoved = 0;
    for $hashKey (sort Numerically keys %{ $filesData }) {
	if (exists($filesData->{$hashKey}{"DATA"}{"-1"})) {
	    $totRemoved = $filesData->{$hashKey}{"DATA"}{"-1"}
	    if ($filesData->{$hashKey}{"DATA"}{"-1"} > $totRemoved);
	}
    }
    
 #   print "Max to be removed: $totRemoved\n";
    for $hashKey (keys %{ $filesData }) {
	$counter = $totRemoved;

	if (exists($filesData->{$hashKey}{"DATA"}{"-1"})) {
	    $counter -= $filesData->{$hashKey}{"DATA"}{"-1"};
	}

	if ($counter > 0) {
#	    print "Removing $counter from " .  $filesData->{$hashKey}{"NAME"} . "\n";
	    $filesData->{$hashKey}{"DATA"} =  RemoveWaters($filesData->{$hashKey}{"DATA"}, $counter);
	}
    }

    return $filesData;
}

sub SaveFiles {
    my ($filesData) = $_[0];
    my ($hashKey, $counter) ;

    $counter = 1;
    for $hashKey (sort Numerically keys %{ $filesData }) {
	print "\t$counter" . ". ";
	writePDB($filesData->{$hashKey}{"DATA"}, $filesData->{$hashKey}{"NAME"});
	$counter++;
    }


}

sub GenerateInsertions {
    my ($index, $counter, $insertpoint, $tmp1, $tmp2, $numions, $fileKey);
    my ($ioninsert, $save_name, $total_atoms, $total_res, %FILEDATA, $fdata);
    
    $total_atoms = $BOX->{"TOTAL_ATOMS"};
    $total_res = $BOX->{"TOTAL_RES"};

    for $index (1 .. $PARMS->{"NumInsertions"}) {
	print "Insertion #" . $index . "\n";
	undef $insertpoint;
	$insertpoint = GetRandomPos($BOX, $PARMS->{"ResRadii"});
	$tmp1 = InsertRes($PARMS->{"RESIDUE"}, $insertpoint, $Atom_Info, 0,
			  ($total_atoms + 1), ($total_res + 1));
	
	$total_atoms += $PARMS->{"AtomsInRes"};
	$total_res++;

	print "\tInserted " . $PARMS->{"Name"} . " at " . $tmp1->{$total_atoms}{"XCOORD"} . 
	    " " . $tmp1->{$total_atoms}{"YCOORD"} . " " . $tmp1->{$total_atoms}{"ZCOORD"} . "\n";

	$save_name = $savePrefix . "Insertion" . $index;
	if ($PARMS->{"HasIons"} and IsInteger($PARMS->{"HasIons"})) {
	    $numions = abs($PARMS->{"ResCharge"} / $PARMS->{"IonCharge"});
	    die "ERROR: Residue:Ion Charge Ratio exceeds limit (8). Aborting\n"
		if ($numions > 8); 
	    $save_name .= "_wIon";
	    %{ $tmp2 } =  %{ $tmp1 };
	    for $counter (1 .. $numions) {
		undef $ioninsert;
		$ioninsert = GetMaxDisplace($BOX, $insertpoint, $PARMS->{"IonRadii"}, $counter);
		$tmp2 = InsertRes($PARMS->{"IONS"}, $ioninsert, $tmp2, 1, 
				  ($total_atoms + 1), ($total_res + 1));
		if ($tmp2->{"-1"}) {
		    $tmp2->{"-1"} -= 1;
		}

		$total_atoms += 1;
		$total_res += 1;

		print "\tInserted " . $PARMS->{"IonName"} . " at " . $tmp2->{$total_atoms}{"XCOORD"} . 
		    " " . $tmp2->{$total_atoms}{"YCOORD"} . " " . $tmp2->{$total_atoms}{"ZCOORD"} . "\n";
	    }
	    $fileKey = 2 * $index + 1;
	    $FILEDATA{$fileKey} = (
				   {
				       "DATA" => $tmp2, 
				       "NAME" => $save_name,
				   }
				   );
	    undef $tmp2;
	}
	if (substr($save_name, -5) eq "_wIon") {
	    $save_name = substr($save_name, 0, -5);
	}
	$fileKey = 2 * $index;
	$FILEDATA{$fileKey} = (
				    {
					"DATA" => $tmp1, 
					"NAME" => $save_name,
				    }
				    );
	undef $tmp1;
    }

    print "\t--==Creating PDB files==--\n";
    $fdata = EqualizeWaters(\%FILEDATA);
    SaveFiles($fdata);
    print "Creating Parameter and Script Files...";
    CreateParmScript($fdata, $PARMS->{"TotalRes"}, $PARMS->{"HasIons"}, $numions);
    print "Done\n";
    undef $fdata;
    undef %FILEDATA;
}

sub InsertRes {
    my ($myRes, $insertLoc, $ATOMS, $recordOverlaps, $totalAtoms, $totalRes) = @_;
    my (%CONFIG, $hashKey, $rec, $counter, @overlaps, $oldoverlaps, $localRes, $index);

    
    if (exists($ATOMS->{"-1"})) {
	$oldoverlaps = $ATOMS->{"-1"};
	delete $ATOMS->{"-1"};
    } else {
	$oldoverlaps = 0;
    }

    for $counter (keys %{ $myRes }) {

	CopyHash(\%{ $localRes }, $myRes->{$counter});
	for $index ("XCOORD", "YCOORD", "ZCOORD") {
	    $localRes->{$index} += $insertLoc->{$index};
	}

	$localRes->{"RES_ID"} = $totalRes;
	
	for $hashKey (keys %{ $ATOMS }) {
	    if (IsOverLap($ATOMS->{$hashKey}, $localRes)) {
		if ($ATOMS->{$hashKey}{"RES_NAME"} ne "WAT") {
		    print "OVERLAP of Solute molecules: " . $ATOMS->{$hashKey}{"LABEL"} . " $hashKey" .
			"with " . $localRes->{"LABEL"} . " $counter \n";
		    print "IGNORING\n";
		} else {
		    $rec->{$ATOMS->{$hashKey}{"RES_ID"}} = 1;

		}

	    }
	}

	CopyHash(\%{ $CONFIG{$totalAtoms + $counter} } , $localRes);
    }
    for $hashKey (keys %{ $ATOMS }) {
	if (! $rec->{$ATOMS->{$hashKey}{"RES_ID"}}) {
	    CopyHash(\%{ $CONFIG{$hashKey} } , \%{ $ATOMS->{$hashKey} });
	}
    }
    
    @overlaps = keys %{ $rec };
    $CONFIG{"-1"} = ($#overlaps + 1) + $oldoverlaps;
    
#    print "Recorded " . ($#overlaps + 1) . " overlaps\n";
    return \%CONFIG;

}

sub IsOverLap {
    my ($res1, $res2) = @_;
    my ($hashKey, $result, $dist);

    $result = $dist = 0;
    for $hashKey ("XCOORD", "YCOORD", "ZCOORD") {
	$dist += (($res1->{$hashKey} - $res2->{$hashKey}))**2;
    }
    $dist = sqrt($dist);
    if ($dist <= ($res1->{"RADII"} + $res2->{"RADII"})) {
	
	$result = 1;
	#print "OVERLAP: res1:" . $res1->{"LABEL"} . " " . $res1->{"RADII"} . 
	    #" and res2: " . $res2->{"LABEL"} . " " . $res2->{"RADII"} . "\n"; 
    }

    return $result;

}


sub GetRandomPos {
    my ($BBOX, $radii) = @_;
    my (%POS, $cindex, $clen);

    $clen = 0.9 * ($BBOX->{"xhi"} - $BBOX->{"xlo"});
    $POS{"XCOORD"} = ($BBOX->{"xlo"} + rand($clen));
    if ((($POS{"XCOORD"} - $BBOX->{"xlo"})/$clen) > 0.5) {
	$POS{"XCOORD"} += $radii;
    } else {
	$POS{"XCOORD"} -= $radii;
    }
	
    $clen = 0.9 * ($BBOX->{"yhi"} - $BBOX->{"ylo"});
    $POS{"YCOORD"} = ($BBOX->{"ylo"} + rand($clen));
    if ((($POS{"YCOORD"} - $BBOX->{"ylo"})/$clen) > 0.5) {
	$POS{"YCOORD"} += $radii;
    } else {
	$POS{"YCOORD"} -= $radii;
    }
	
    $clen = 0.9 * ($BBOX->{"zhi"} - $BBOX->{"zlo"});
    $POS{"ZCOORD"} = ($BBOX->{"zlo"} + rand($clen));
    if ((($POS{"XCOORD"} - $BBOX->{"zlo"})/$clen) > 0.5) {
	$POS{"ZCOORD"} += $radii;
    } else {
	$POS{"ZCOORD"} -= $radii;
    }
	
    for $clen ("XCOORD", "YCOORD", "ZCOORD") {
	$POS{$clen} =  sprintf("%.2f", $POS{$clen});
    }
    
    return \%POS;
}

sub GetMaxDisplace {
    my ($BBOX, $resLoc, $radii, $ionCounter) = @_;
    my ($clen, %POS);

    $clen = ($BBOX->{"xhi"} - $BBOX->{"xlo"})/2;
    if (($resLoc->{"XCOORD"} - $BBOX->{"xlo"}) > $clen) {
	$POS{"XCOORD"} = $resLoc->{"XCOORD"} - $clen + $radii;
    } else {
	$POS{"XCOORD"} = $resLoc->{"XCOORD"} + $clen - $radii;
    }
    if ($ionCounter == 4) {
	$POS{"XCOORD"} = -1 * $POS{"XCOORD"};
    }

    $clen = ($BBOX->{"yhi"} - $BBOX->{"ylo"})/2;
    if (($resLoc->{"YCOORD"} - $BBOX->{"ylo"}) > $clen) {
	$POS{"YCOORD"} = $resLoc->{"YCOORD"} - $clen + $radii;
    } else {
	$POS{"YCOORD"} = $resLoc->{"YCOORD"} + $clen - $radii;
    }
    if ($ionCounter == 3 or $ionCounter == 6) {
	$POS{"YCOORD"} = -1 * $POS{"YCOORD"};
    }

    $clen = ($BBOX->{"zhi"} - $BBOX->{"zlo"})/2;
    if (($resLoc->{"ZCOORD"} - $BBOX->{"zlo"}) > $clen) {
	$POS{"ZCOORD"} = $resLoc->{"ZCOORD"} - $clen + $radii;
    } else {
	$POS{"ZCOORD"} = $resLoc->{"ZCOORD"} + $clen - $radii;
    }
    if ($ionCounter == 2 or $ionCounter == 5 or $ionCounter == 8) {
	$POS{"ZCOORD"} = -1 * $POS{"ZCOORD"};
    }

    for $clen ("XCOORD", "YCOORD", "ZCOORD") {
	$POS{$clen} = sprintf("%.2f", $POS{$clen});
    }

    return \%POS;
}

sub GetParms {
    my ($paramfile) = $_[0];
    my ($key_name, $key_val, %PARAMS, $res, $counter);
    die "Error; Cannot access regular file $paramfile: $!\n"
	if (! -e $paramfile or ! -r $paramfile or ! -T $paramfile);

    open PARAMFILE, $paramfile or die "Cannot open $paramfile: $!\n";
    while (<PARAMFILE>) {
	chomp;
	($key_name, $key_val) = GetVals($_);
	if (defined($key_name)) {
	    $PARAMS{$key_name} = $key_val;
	}
    }
    close PARAMFILE;

    die "Error: BGF file keyword not found\n"
	if (! $PARAMS{"BGF_FILE"});
    die "Error: Cannot access regular file " . $PARAMS{"BGF_FILE"} . ": $!\n"
	if (! -e $PARAMS{"BGF_FILE"} or ! -r $PARAMS{"BGF_FILE"} or ! -T $PARAMS{"BGF_FILE"});
    
    for $counter ("Name", "ResName", "ResRadii", "ResCharge", "HasIons", "TotalRes", "WaterRadii") {
	die "ERROR: Keyword not found $counter\n"
	    if (! exists($PARAMS{$counter}));
    }

    if (exists($PARAMS{"ResFile"})) {
	GetResInfo($PARAMS{"ResFile"}, \%PARAMS);
    } else {
	AddResInfo("RESIDUE", $PARAMS{"Name"}, $PARAMS{"ResName"}, 
		   $PARAMS{"ResRadii"}, $PARAMS{"ResCharge"}, \%PARAMS);
    }

    if ($PARAMS{"HasIons"}) {
	AddResInfo("IONS", $PARAMS{"IonName"}, $PARAMS{"IonResName"}, 
		   $PARAMS{"IonRadii"}, $PARAMS{"IonCharge"}, \%PARAMS);
    }

    return \%PARAMS;


}

sub AddResInfo {
    my ($hashkey, $name, $resname, $radii, $charge, $SaveHash) = @_;
    my ($rec);

    $rec = (
	    {
		"HEADER"   => "ATOM",
		"LABEL"    => $name,
		"RES_NAME" => $resname,
		"RES_ID"   => 0,
		"XCOORD"     => "0.0",
		"YCOORD"     => "0.0",
		"ZCOORD"     => "0.0",
		"CHARGE"   => $charge,
		"RADII"    => $radii,
	    }
	    );
    $SaveHash->{$hashkey}{"0"} = $rec;

}

sub GetVals {
    my ($linein) = $_[0];
    my ($return_name, $return_val);

    if ($linein =~ /^PARM_(\S+):\s+(\S+)$/) {
	$return_name = $1;
	$return_val = $2;
    }

    return ($return_name, $return_val);


}

sub GetResInfo {
    my ($resname, $SaveHash) = @_;
    my ($ResAtoms);
    
    die "Cannot locate residue template file $resname: $!\n"
	if (! -e $resname or ! -r $resname or ! -T $resname);

    $ResAtoms = GetPDBFileInfo($resname, 1, " ", " ", 1);
    CopyHash(\%{ $SaveHash->{"RESIDUE"} }, $ResAtoms);
}

sub ParsePDBFile {
    my ($pdb_file) = $_[0];
    my ($Atom_Info) = GetPDBFileInfo($pdb_file, 1, " ", " ", 1);
    my ($counter, $curr_res_name, $curr_res_id, $curr_atom, $BBOX, $tmp, $hashkey, $radii);

    $BBOX = (
	     {
		 "xhi" => -9999,
		 "xlo" => 9999,
		 "yhi" => -9999,
		 "ylo" => 9999,
		 "zhi" => -9999,
		 "zlo" => 9999,
	     }
	     );

    $curr_res_id = 0;
    for $counter (sort Numerically keys %{ $Atom_Info }) {
	$curr_res_id = $Atom_Info->{$counter}{"RES_ID"};
	$curr_res_name = $Atom_Info->{$counter}{"RES_NAME"};
	$curr_atom =  $Atom_Info->{$counter}{"LABEL"};
	if (($curr_res_name eq "WAT" and $curr_atom eq "O") or  
	    (! $curr_res_name eq "WAT")) {
	    $radii = $Atom_Info->{$counter}{"RADII"};
	    $tmp = (
		    {
			"xhi" => $Atom_Info->{$counter}{"XCOORD"} + $radii,
			"xlo" => $Atom_Info->{$counter}{"XCOORD"} - $radii,
			"yhi" => $Atom_Info->{$counter}{"YCOORD"} + $radii,
			"ylo" => $Atom_Info->{$counter}{"YCOORD"} - $radii,
			"zhi" => $Atom_Info->{$counter}{"ZCOORD"} - $radii,
			"zlo" => $Atom_Info->{$counter}{"ZCOORD"} - $radii,
		    }
		    );
	    for $hashkey ("xhi", "yhi", "zhi") {
		if ($BBOX->{$hashkey} < $tmp->{$hashkey}) {
		    $BBOX->{$hashkey} = $tmp->{$hashkey};
		}
	    }
	    
	    for $hashkey ("xlo", "ylo", "zlo") {
		if ($BBOX->{$hashkey} > $tmp->{$hashkey}) {
		    $BBOX->{$hashkey} = $tmp->{$hashkey};
		}
	    }
	}
	$BBOX->{"TOTAL_ATOMS"} = $counter;
	$BBOX->{"TOTAL_RES"} = $curr_res_id;
    }
    
    $BBOX->{"CENTER"} = GetBoxCenter($BBOX);
    return ($Atom_Info, $BBOX);
}

sub Numerically {
    ($a <=> $b);
}

sub writePDB {
    my ($AtomData, $saveName) = @_;
    my ($index, $resName, $resIndex, $atomIndex, $resCounter);
    
    $resName = "";
    $resIndex = $atomIndex = $resCounter = 0;
    
    $saveName .= ".pdb";
    print "Creating $saveName...";
    open OUTDATA, "> $saveName" or die "Cannot write to regular file $saveName: $!\n";
    for $index (sort Numerically keys %{ $AtomData }) {
	next
	    if ($index < 0);
	$atomIndex++;
	if ($AtomData->{$index}{"RES_ID"} != $resCounter) {
	    $resName = $AtomData->{$index}{"RES_NAME"}; 
	    $resIndex++;
	    $resCounter = $AtomData->{$index}{"RES_ID"};
	}
	printf OUTDATA "%-5s%6d%2s%-4s%3s%6d%12.3f%8.3f%8.3f%12.3f%12.3f\n",
	$AtomData->{$index}{"HEADER"}, $atomIndex, " ", $AtomData->{$index}{"LABEL"}, 
	$resName, $resIndex, $AtomData->{$index}{"XCOORD"}, $AtomData->{$index}{"YCOORD"}, 
	$AtomData->{$index}{"ZCOORD"}, $AtomData->{$index}{"CHARGE"}, $AtomData->{$index}{"RADII"};
    }

    close OUTDATA;
    print "\n";
}

sub RemoveWaters {
    my ($AtomData, $numWats) = @_;
    my ($counter, @WATLIST);

    for $counter (keys %{ $AtomData }) {
	next
	    if ($counter == -1);
	if (($AtomData->{$counter}{"RES_NAME"} eq "WAT") and ($AtomData->{$counter}{"LABEL"} !~ /H/i)) {
	    push @WATLIST, $AtomData->{$counter}{"RES_ID"};
	}
    }
    
    for $counter (1 .. $numWats) {
	$AtomData = RemoveResidue($AtomData, int(rand($#WATLIST)));
    }

    return $AtomData;
}

sub RemoveResidue {
    my ($ATOMS, $res_id) = @_;
    my ($hashKey, %NewAtoms);

    for $hashKey (keys %{ $ATOMS }) {
	next
	    if ($hashKey == -1);
	if ($ATOMS->{$hashKey}{"RES_ID"} != $res_id) {
	    $NewAtoms{$hashKey} = $ATOMS->{$hashKey};
	} else {
	    print "";
	}
    }

    return \%NewAtoms;

}

sub GetSmallestBoxLen {
    my ($BOXINFO)  = $_[0];
    my ($returnVal);

    $returnVal = ($BOXINFO->{"xhi"} - $BOXINFO->{"xlo"});
    $returnVal = ($BOXINFO->{"yhi"} - $BOXINFO->{"ylo"})
	if (($BOXINFO->{"yhi"} - $BOXINFO->{"ylo"}) < $returnVal);
    $returnVal = ($BOXINFO->{"zhi"} - $BOXINFO->{"zlo"})
	if (($BOXINFO->{"zhi"} - $BOXINFO->{"zlo"}) < $returnVal);

    return $returnVal;
}
    
sub GetBoxCenter {
    my ($BOXInfo) = $_[0];
    my ($CENTER) = (
		    {
			"XCOORD" => $BOXInfo->{"xlo"} + (($BOXInfo->{"xhi"} - $BOXInfo->{"xlo"})/2),
			"YCOORD" => $BOXInfo->{"ylo"} + ($BOXInfo->{"yhi"} - $BOXInfo->{"ylo"})/2,
			"ZCOORD" => $BOXInfo->{"zlo"} + ($BOXInfo->{"zhi"} - $BOXInfo->{"zlo"})/2,
		    }
		    );

    return $CENTER;
}

sub CopyHash {
    my ($hashRef1, $hashRef2) = @_;
    my ($hashKey);

    for $hashKey (keys %{ $hashRef2 }) {
	$hashRef1->{$hashKey} = $hashRef2->{$hashKey};
    }
}

sub DeterminePlacement {
    my ($sLim, $ResRadii, $WATRadii, $BOXInfo) = @_;
    my ($clen, %CUBE, @PLACE, $displace, $index, $counter);

    $ResRadii = $ResRadii * 2;
    $WATRadii = $WATRadii;

    # @PLACE Array holds the placement information for x, y, z = +/- 1
    @PLACE = (
	      [-1,-1,-1],
	      [ 1, 1, 1],
	      [ 1,-1,-1],
	      [-1, 1, 1],
	      [ 1,-1, 1],
	      [-1, 1,-1],
	      [-1,-1, 1],
	      [ 1, 1,-1],
	      );
	       
    # now using the box center, build consecutive cubic boxes 1 watershell apart
    # the first one should be the radii of the residue
    
    $counter = 0;
    $displace = ($ResRadii + 0.25) - $WATRadii;
    while ($displace < $sLim) {
	$displace += $WATRadii;
	for $index (0 .. $#PLACE) {
	    $CUBE{$counter}{$index} = (
				       {
					   "XCOORD" => ($BOXInfo->{"CENTER"}{"XCOORD"} + ($displace * $PLACE[$index][0])),
					   "YCOORD" => ($BOXInfo->{"CENTER"}{"YCOORD"} + ($displace * $PLACE[$index][1])),
					   "ZCOORD" => ($BOXInfo->{"CENTER"}{"ZCOORD"} + ($displace * $PLACE[$index][2])),
				       }
				       );
	}
	$counter++;
    }

    return \%CUBE;

}

sub DetermineIonPlacement {
    my ($BOXInfo, $boxPos, $boxSpacing) = @_;
    my (%FCCube, $clen, $counter, $index);

    $boxSpacing = $boxSpacing * 2;

    # Now until we can go no further, create face centered poistions for ions on cubes with boxspacing buffer
    $counter = $index = 0;
    while ($boxPos > $boxSpacing) {
	$FCCube{$counter + 1}{"XCOORD"} = ($BOXInfo->{"CENTER"}{"XCOORD"} + $boxPos); # +x
	$FCCube{$counter + 2}{"XCOORD"} = ($BOXInfo->{"CENTER"}{"XCOORD"} - $boxPos); # -x
	$FCCube{$counter + 3}{"YCOORD"} = ($BOXInfo->{"CENTER"}{"YCOORD"} + $boxPos); # +y
	$FCCube{$counter + 4}{"YCOORD"} = ($BOXInfo->{"CENTER"}{"YCOORD"} - $boxPos); # -y
	$FCCube{$counter + 5}{"ZCOORD"} = ($BOXInfo->{"CENTER"}{"ZCOORD"} + $boxPos); # +z
	$FCCube{$counter + 6}{"ZCOORD"} = ($BOXInfo->{"CENTER"}{"ZCOORD"} - $boxPos); # -z
	$counter += 6;
	$boxPos -= $boxSpacing;
    }

    # Finally, make sure that the other dimensions are filled with the default (center) values
    for $index (keys %FCCube) {
	for $counter ("XCOORD", "YCOORD", "ZCOORD") {
	    if (! exists($FCCube{$index}{$counter})) {
		$FCCube{$index}{$counter} = $BOXInfo->{"CENTER"}{$counter};
	    }
	}
    }
	

    return \%FCCube;


}

