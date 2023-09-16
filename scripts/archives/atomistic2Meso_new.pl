#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts");
}

use strict;
use Packages::FileFormats;
use Packages::General;
use File::Basename;
use Packages::Math::Matrix;
use Packages::BOX;

#   atomistic2Meso.pl: This script will open an atomistic bgf file
#   and create the corresponding meso-scale representation

sub main;
sub Initialize;
sub GetParms;
sub Numerically;
sub CreateMesoModel;
sub getStrand;
sub getBeadInfo;
sub isMember;
sub CreateHBondBead;
sub updateCOM;
sub getSaveName;
sub ExtractResidue;
sub GetChain;
sub createBGFHeader;
sub makeBGF;
sub buildBondList;
sub getSearchList;
sub GetAtmMass;
sub getStrandLength;

die "usage: $0 parm_file bgf_file|directory [save_name] [strand length]\n"
    if (! @ARGV or $#ARGV < 1);

my ($parmFile, $fileLoc, $save_name, $sLen) = @ARGV;
my ($FILES, $fC, $ATOMS, $BONDS, $i, $PARMS);
my ($MESO, $bgfFile, $BGF, $CONS, $CHAIN);

&main;

sub main {
    $|++;
    
    $FILES = Initialize($parmFile, $fileLoc);
    print "Obtaining Bead parameters ....";
    $PARMS = GetParms($parmFile);
    print "Done\n\n";
    
    $fC = 0;
    
    for $i(@{ $FILES }) {
        print "$i:....READING...";
	($ATOMS, $BONDS) = GetBGFFileInfo($i, 0);
        $sLen = getStrandLength($ATOMS) if (! defined($sLen));
        next if (! $sLen);
        print "....Creating MESO";
        $MESO = CreateMesoModel($ATOMS, $PARMS);
        $bgfFile = getSaveName($i);
        print "....Writing Meso";
        chdir "structures";
        ($BGF, $CONS) = makeBGF($MESO, $bgfFile, $PARMS, $sLen, $fC);
    }
    
    print "\nAll Tasks Completed\n";

}

sub Initialize {
    my ($pLoc, $fLoc) = @_;
    my (@FILES);

    FileTester($pLoc);

    if (! -T $fLoc) {
	if (-d $fLoc) {
	    opendir FILEDIR, $fLoc or die "ERROR: Cannot access directory $fLoc: $!\n";
	    @FILES = grep { /\.bgf$/ && -f} map { "$fLoc/$_"} readdir FILEDIR;
	    closedir FILEDIR or die "ERROR: Cannot close directory $fLoc: $!\n";
	}
    } else {
	$FILES[0] = $fLoc;
    }

    die "ERROR: No valid pdbfiles found!\n"
	if ($#FILES == -1);

    return \@FILES;
}

sub GetParms {
    my ($parm) = $_[0];
    my ($curr_id, $rec, $mykey, %PARMS_HOLDER);

    $curr_id = -1;
    $rec = ();

    open PARMFILE, $parm or die "Cannot open parameter file $parm: $!\n";
    while(<PARMFILE>) {
	chomp;
	if ($_ =~ /^BEAD_ID: (\d+)/) {
	    if ($curr_id > -1) {
		%{ $PARMS_HOLDER{"BEADS"}{$curr_id} } = %{ $rec };
	    }
	    $curr_id = $1;
	    $rec = ();
	} elsif ($curr_id > -1 and $_ =~ /^BEAD_(\w+): (.+)/) {
	    $rec->{$1} = $2;
	} elsif ($_ =~ /^BOND_\d+: (\w+\-\w+):(\d+)/) {
	    $PARMS_HOLDER{"BONDS"}{$1} = $2;
	} elsif ($_ =~ /^ANGLE_\d+: (\w+\-\w+\-\w+):(\d+)/) {
	    $PARMS_HOLDER{"ANGLES"}{$1} =  $2;
	} elsif ($_ =~ /^TORSION_\d+: (\w+\-\w+\-\w+\-\w+):(\d+)/) {
	    $PARMS_HOLDER{"TORSIONS"}{$1} = $2;
	} elsif ($_ =~ /^MASS_(\w+)\s+(\d+\.*\d*)/) {
	    $PARMS_HOLDER{"MASSES"}{$1} = $2;
	} elsif ($_ =~ /^EQUIV_\d+: (\d+)\s+(\d+)/) {
	    $PARMS_HOLDER{"EQUIVALENCE"}{$1}{$2} = 1;
	    $PARMS_HOLDER{"EQUIVALENCE"}{$2}{$1} = 1;
	}
	
    }
    
    close PARMFILE;

    %{ $PARMS_HOLDER{"BEADS"}{$curr_id} } = %{ $rec };

    die "ERROR: Invalid parameter file $parm\n"
	if ($curr_id == -1);
    return (\%PARMS_HOLDER);
}

sub Numerically {
    ($a<=>$b);
}

sub CreateMesoModel {
    my ($atoms, $parms) = @_;
    my ($atm_counter, $bead_counter, $atm_label, $bead_member, $curr_id, $res_id);
    my (%MESO_MODEL, $atm_res, $base_res, $hbondA, $hbondD, $non_member, $counter);
    my ($result, $isAcceptor, @search_elements, @hblist, $index, $oldResId, $hbType);
    my ($avg, $stdev, $total, $tmp, $hbond_nm, $companion, $strandEnd, $beadName, $ff);

    $res_id = $curr_id = -1;
    $strandEnd = 1;
    $oldResId = 0;

    for $atm_counter (sort Numerically keys %{ $atoms }) {
	$res_id = $atoms->{$atm_counter}{"RESNUM"};
	if (! $oldResId) {
	    $oldResId = $res_id;
	} 

	$atm_label = Trim($atoms->{$atm_counter}{"ATMNAME"});
	$atm_res = Trim($atoms->{$atm_counter}{"RESNAME"});
	$strandEnd = getStrand($atm_res, $strandEnd, $res_id, $oldResId);
	$oldResId = $res_id;
	
	$atm_res =~ s/\d+//g;
	
	for $bead_counter (keys %{ $parms->{"BEADS"} }) {
            $base_res = $parms->{"BEADS"}{$bead_counter}{"RES"};
	    next if (! isMember($base_res,$atm_res));
	    $bead_member = $parms->{"BEADS"}{$bead_counter}{"MEMBERS"};
	    $non_member = $parms->{"BEADS"}{$bead_counter}{"NONMEMBERS"};
	    $hbondA = $parms->{"BEADS"}{$bead_counter}{"HBOND_ACCEPTORS"};
	    $hbondD = $parms->{"BEADS"}{$bead_counter}{"HBOND_DONORS"};
	    $beadName = $parms->{"BEADS"}{$bead_counter}{"NAME"};
	    
	    if ($beadName =~ /PHP|PHO|SUG|SUS/) {
	        ($beadName, $bead_counter) = getBeadInfo($beadName, $bead_counter, $strandEnd, $parms);
	    }
	    $ff = $parms->{"BEADS"}{$bead_counter}{"ELEMENT"};
	    
	    if (isMember($bead_member,$atm_label)) {
		$MESO_MODEL{$res_id}{$bead_counter}{"ATOMLIST"}{$atm_label} = 1;
		$MESO_MODEL{$res_id}{$bead_counter}{"ATOMS"}{$atm_counter} = 1;
		$MESO_MODEL{$res_id}{$bead_counter}{"INDEX"} = $bead_counter;
		updateCOM(\%{ $MESO_MODEL{$res_id}{$bead_counter} },\%{ $atoms->{$atm_counter} });
		$MESO_MODEL{$res_id}{$bead_counter}{"FFTYPE"} = $ff;
		$MESO_MODEL{$res_id}{$bead_counter}{"SOLUTE"} = $atoms->{$atm_counter}{"SOLUTE"};
		$MESO_MODEL{$res_id}{$bead_counter}{"ATMNAME"} = $beadName;
		
		for $counter ("CHARGE", "RADII", "NUMBONDS", "LONEPAIRS", "ELEMENT") {
		    $MESO_MODEL{$res_id}{$bead_counter}{$counter} = $parms->{"BEADS"}{$bead_counter}{$counter};
		}
		
	    } elsif (isMember($non_member,$atm_label)) {
		$MESO_MODEL{$res_id}{$bead_counter}{"ATOMLIST"}{$atm_label} = 1;
	    }

#	    Hydrogen Bonds Acceptors
	    if ($hbondA and isMember($hbondA,$atm_label)) {
		for $counter ("XCOORD", "YCOORD", "ZCOORD") {
		    $MESO_MODEL{$res_id}{$bead_counter}{"HBONDS"}{"ACCEPTORS"}{$atm_label}{$counter} = 
			$atoms->{$atm_counter}{$counter};
		}
	    }
	    
	    #Hydrogen Bond Donors
	    if ($hbondD and isMember($hbondD,$atm_label)) {
		for $counter ("XCOORD", "YCOORD", "ZCOORD") {
		    $MESO_MODEL{$res_id}{$bead_counter}{"HBONDS"}{"DONORS"}{$atm_label}{$counter} = 
			$atoms->{$atm_counter}{$counter};
		}
	    }
	}
    }
    
    for $res_id (keys %MESO_MODEL) {
	for $curr_id (keys %{ $MESO_MODEL{$res_id} }) {
	    $bead_counter = $MESO_MODEL{$res_id}{$curr_id};
            if (! $bead_counter->{"COM"}{"TotalMass"} ) {
                delete $MESO_MODEL{$res_id}{$curr_id};
                next;
            }
	    $bead_counter->{"XCOORD"} = $bead_counter->{"COM"}{"X"}/$bead_counter->{"COM"}{"TotalMass"};
	    $bead_counter->{"YCOORD"} = $bead_counter->{"COM"}{"Y"}/$bead_counter->{"COM"}{"TotalMass"};
	    $bead_counter->{"ZCOORD"} = $bead_counter->{"COM"}{"Z"}/$bead_counter->{"COM"}{"TotalMass"};
	    $index = $total = 1;
	    next
		if (! exists($bead_counter->{"HBONDS"}));

	    for $hbType ("DONORS", "ACCEPTORS") {
                $index = 50 * $total + 1;
		next
		    if (! exists($bead_counter->{"HBONDS"}{$hbType}));
		for $counter (keys %{ $bead_counter->{"HBONDS"}{$hbType} }) {
		    CreateHBondBead(\%{ $MESO_MODEL{$res_id} }, $index, $hbType, $counter, $curr_id);
		    $index++;
		}
                $total++;
	    }
	}
    }
	
    return \%MESO_MODEL;	
}

sub getStrand {
    my ($resName, $strandEnd, $res_id, $oldResId) = @_;
    
    $resName =~ /(\d+)/;
    if (defined($1)) {
	if ($1 == 5) {
	    $strandEnd = 1;
	} else {
	    $strandEnd = -1;
	}
    } else {
	if ($res_id > $oldResId) {
	    $strandEnd *= -1;
	}
    }

    return $strandEnd;
}

sub getBeadInfo {
    my ($beadName, $beadCounter, $strandEnd, $parms) = @_;

    if ($beadName =~ /PHP|PHO/) {
	if ($strandEnd == 1) {
	    $beadName = "PHO";
	    $beadCounter = 1;
	} else {
	    $beadName = "PHP";
	    $beadCounter = 2;
	}
    } elsif ($beadName =~ /SUG|SUS/) {
	if ($strandEnd == 1) {
	    $beadName = "SUG";
	    $beadCounter = 3;
	} else {
	    $beadName = "SUS";
	    $beadCounter = 4;
	}
    }
    
    return ($beadName, $beadCounter);
    
}

sub isMember {
    my ($searchStr, $searchItem) = @_;
    my (@tmp, $returnVal);
    
    $returnVal = 0;
    if (! defined($searchStr)) {
	$returnVal = 0;
    } else {
	
	if ($searchStr !~ /,/) {
	    push @tmp, $searchStr;
	} else {
	    @tmp = split /,/, $searchStr;
	}
    
	for my $counter (@tmp) {
	    if (lc($counter) eq lc($searchItem)) {
		$returnVal =  1;
		last;
	    }
	}
    }

    return $returnVal;
}

sub CreateHBondBead {
    my ($res, $index, $type, $atmName, $parent) = @_;
    my (%Curr_Bead, $counter, $hbond);
    
    $hbond = $res->{$parent}{"HBONDS"}{$type}{$atmName};
    $Curr_Bead{"ATOMLIST"}{$atmName} = 1;
    $Curr_Bead{"INDEX"} = $index;
    
    for $counter ("XCOORD", "YCOORD", "ZCOORD", "ATOMS", "COMPANION") {
	$Curr_Bead{$counter} = $hbond->{$counter};
    }
    $Curr_Bead{"CHARGE"} = 0.0;
    $Curr_Bead{"RADII"} = 1.058;
    $Curr_Bead{"NUMBONDS"} = 1;
    $Curr_Bead{"LONEPAIRS"} = 0;
    $type = substr($type, 0, 1);
    if ($type !~ /D/i) {
        $type = "H";
    }
    $Curr_Bead{"ELEMENT"} = $type;
    $Curr_Bead{"ATMNAME"} = "H" . $index;
    $Curr_Bead{"FFTYPE"} = $type . ($index % 50);
    $Curr_Bead{"PARENT"} = $parent;
    $Curr_Bead{"SOLUTE"} = 1;

    %{ $res->{$index} } = %Curr_Bead;
}

sub updateCOM {
    my ($curr_bead, $atm_counter) = @_;
    my ($atm_mass);

    $atm_mass = GetAtmMass(Trim($atm_counter->{"ATMNAME"}));

    $curr_bead->{"COM"}{"X"} += $atm_counter->{"XCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"Y"} += $atm_counter->{"YCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"Z"} += $atm_counter->{"ZCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"TotalMass"} += $atm_mass;

#    return $curr_bead;
}

sub getSaveName {
    my ($curr_file) = $_[0];

    $curr_file = basename($curr_file);
    if ($save_name) {
	$curr_file = $save_name;
    }

    $curr_file =~ s/\.\w+$/\.bgf/i;

    return $curr_file;

}

sub ExtractResidue {
    my ($curr_res) = $_[0];
    my ($returnval, $tmp);

    for $tmp (keys %{ $curr_res }) {
	if ($curr_res->{$tmp}{"ATMNAME"} !~ /SUS|SUG|PHP|PHO|H\d+/) {
	    $returnval = $curr_res->{$tmp}{"ATMNAME"};
	    last;
	}
    }

    return $returnval;

}

sub GetChain {
    my ($curr_chain, $strand_length) = @_;
    my ($counter, $returnval, $largest_val);

    $counter = 0;
    $counter = int($curr_chain/($strand_length + 1)) + 1;
    $returnval = uc(chr($counter + 64));
    
    return $returnval;
}

sub createBGFHeader {
    my ($HEADER, $muser, $mdate);

    $muser = "tpascal";#$ENV{USER};
    $mdate = (localtime)[3] . "/" .  (localtime)[4] . "/" . ((localtime)[5] + 1900);
    $mdate .= " at " .  sprintf("%2d", (localtime)[2]) . ":" .   sprintf("%2d", (localtime)[1]) . 
	":" . sprintf("%2d", (localtime)[0]);
    push @{ $HEADER }, "BIOGRF  332"; 
    push @{ $HEADER }, "DESCRP MesoDNA"; 
    push @{ $HEADER }, "REMARK Bead Model for Atomistic DNA";
    push @{ $HEADER }, "REMARK Created by $muser on $mdate"; 
    push @{ $HEADER }, "FORCEFIELD DREIDING";
 
    return $HEADER;
}

sub makeBGF {
    my ($Meso, $bgfFile, $parms, $strand_len, $file_no) = @_;
    my ($HEADER, $resC, $beadC, $bead, $index, $BGF, $BONDS);

    $HEADER = createBGFHeader;
    $index = 1;
    print "..atoms";
    for $resC (sort Numerically keys %{ $Meso } ) {
	for $beadC (sort Numerically keys %{ $Meso->{$resC} } ) {
	    $bead = $Meso->{$resC}{$beadC};
	    $bead->{"INDEX"} = $index;
	    $bead->{"LABEL"} = "ATOM";
	    $bead->{"RESNAME"} = Trim(ExtractResidue($Meso->{$resC}));
	    $bead->{"CHAIN"} = GetChain($resC, $strand_len);
	    $bead->{"RESNUM"} = $resC;
	    $bead->{"ID"} = $beadC;
	    $bead->{"ATMNAME"} = Trim($bead->{"ATMNAME"});
	    $bead->{"FFTYPE"} = Trim($bead->{"FFTYPE"});
	    if (exists($bead->{"PARENT"})) {
		$bead->{"PARENTID"} = $Meso->{$resC}{ $bead->{"PARENT"} }{"INDEX"};
	    }
	    
	    %{ $BGF->{$index} } = %{ $bead };
	    $index++;
	}
    }

    print "..bonds";
    $BONDS = buildBondList($BGF, $Meso, $parms, $strand_len);
    print "..saving";
    addHeader($BGF, $HEADER);
    createBGF($BGF, $BONDS, $bgfFile);

    return ($BGF, $BONDS);
}

sub buildBondList {
    my ($BGF, $Meso, $parms, $strand_len) = @_;
    my ($myType, $index, $beadC, $resC, $counter, $bType, $bond);
    my ($searchKey, $searchList, %BONDS, @bondList, %CONS, $parent);

    for $beadC (sort Numerically keys %{ $BGF } ) {
        next if (! $BGF->{$beadC}{"SOLUTE"});
	if (exists($BGF->{$beadC}{"PARENT"})) { #hbond
	    $parent = $BGF->{$beadC}{"PARENTID"};
	    $BONDS{$beadC}{$parent} = 1;
	    $BONDS{$parent}{$beadC} = 1;
	    next;
	}
 	$resC = $BGF->{$beadC}{"RESNUM"};
	$myType = $BGF->{$beadC}{"ID"};
        
        @bondList = ();
	for $counter (0, 1) {
	    next if (! exists($Meso->{$resC + $counter }));
	    for $index (keys %{ $Meso->{$resC + $counter} }) {
		$bond = (
			 {
			     "RES"   => $resC + $counter,
			     "ATOM"  => $Meso->{$resC + $counter}{$index},
			 }
			 );
		push @bondList, $bond;
	    }
	    last if ($resC % $strand_len == 0); # end of strand            
	}

	for $bond (@bondList) {
	    next
		if ($bond->{"ATOM"}{"INDEX"} == $beadC);
	    $bType = $bond->{"ATOM"}{"ID"};
	    $searchList = getSearchList($myType, $bType, $parms);
	    for $searchKey (@{ $searchList }) {
		if (exists($parms->{"BONDS"}{$searchKey}) and
		    ($parms->{"BONDS"}{$searchKey} == abs($resC - $bond->{"RES"}))) {
		    $BONDS{$beadC}{ $bond->{"ATOM"}{"INDEX"} } = 1;
		    $BONDS{ $bond->{"ATOM"}{"INDEX"} }{$beadC} = 1;
		    last;
		}
	    }
	}
    }

    for $bond (keys %BONDS) {
	@{ $CONS{$bond} } = sort Numerically keys %{ $BONDS{$bond} };
    }
    
    return \%CONS;
}

sub getSearchList {
    my ($aType, $bType, $parms) = @_;
    my (@searchList, %keyList, $key1, $key2);

    @{ $keyList{$bType} } = ($aType, "X");
    @{ $keyList{$aType} } = ($bType, "X");
    
    if (exists($parms->{"EQUIVALENCE"}{$bType})) {
	for $key1 (keys %{ $parms->{"EQUIVALENCE"}{$bType} }) {
	    push @{ $keyList{$aType} }, $key1;
	}
    }

    if (exists($parms->{"EQUIVALENCE"}{$aType})) {
	for $key1 (keys %{ $parms->{"EQUIVALENCE"}{$aType} }) {
	    push @{ $keyList{$bType} }, $key1;
	}
    }
    
    for $key1 (@{ $keyList{$bType} }) {
        for $key2 (@{ $keyList{$aType} }) {
            push @searchList, $key1 . "-" . $key2;
        }
    }
    
    return \@searchList;
}

sub GetAtmMass {
    my ($atm_label) = $_[0];
    my ($returnval, $in_label);

    $returnval = 0;
    
    if ($atm_label =~ /^\d*([a-z]+)/i) {
	$in_label = $1;
   
	if ($PARMS->{"MASSES"}{$in_label}) {
	    $returnval = $PARMS->{"MASSES"}{$in_label};
	}
    }
    return $returnval;
}

sub getStrandLength {
    my ($fDATA) = $_[0];
    my ($resNum, $sLen, $sI, $atomC, $atom, $resName);

    $resNum = $sLen = 0;
    $sI = "";

    for $atomC (sort Numerically keys %{ $fDATA }) {
	$atom = $fDATA->{$atomC};
	if ($atom->{"RESNUM"} > $resNum) {
	    $resNum = $atom->{"RESNUM"};
	    $resName = Trim($atom->{"RESNAME"});
	    if ($resName =~ /(\d+)$/) {
		if ($sI eq "") {
		    $sI = $1;
		    $sLen = $resNum;
		} elsif ($sI ne $1) {
		    $sLen = $resNum;
		    last;
		}
	    }
	}
    }

    if (! $sLen) {
	die "ERROR: CANNOT find strand length\n";
    }
    return ($sLen);
}
