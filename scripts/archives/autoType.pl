#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts/";
}

# This script will parse a Cerius2 bgffile, and update the force field types
# to be compatible with amber

use strict;
use File::Basename;
use Packages::General;
use Packages::FileFormats;

sub initialize;
sub parseAmberLib;
sub updateBGF;
sub Numerically;
sub getResList;
sub getCons;
sub determineValidCons;
sub getMaxDepth;
sub getNewAtm;
sub getResData;
sub printFound;
sub findElement;
sub doUpdate;
sub saveNode;
sub syncList;
sub determineRes;

die "usage: $0 bgf_file amber_lib startAtom [save_name]\n"
    if (! @ARGV or $#ARGV < 2);

my ($bgf_file, $amber_lib, $startAtm, $save_name) = @ARGV;
my (%ELEMENTS, $ATMTREE, $RESTREE, %USED, %ERES, $FDATA);
    
initialize;

print "STEP 1. Parsing $bgf_file...";
my ($BGF_INFO, $CONN, $HEADERS) = GetBGFFileInfo($bgf_file, 1);
addCon($BGF_INFO, $CONN);
addElement($BGF_INFO, "ATMNAME");

print "Done\n";

print "STEP 2. Parsing AMBER LIB file $amber_lib...";
my ($RES_INFO) = parseAmberLib($amber_lib);
print "Done\n";

print "STEP 3. Scanning $bgf_file for residue information...\n";
updateBGF($BGF_INFO, $CONN, $RES_INFO, $startAtm);

addHeader($BGF_INFO, $HEADERS);
print "STEP 4. Creating BGF File $save_name...";
createBGF($BGF_INFO, $CONN, $save_name);
print "Done\n";

sub parseAmberLib {
    my ($inFile) = $_[0];
    my (%RESDATA, $i, $start, $resName, $rec, $recC, $connect, $atomC);
    $start = $connect = 0;

    open LIBFILE, $inFile or die "Cannot open AMBER LIB file $inFile: $1\n";
    while (<LIBFILE>) {
	chomp;
	if ($_ =~ /^\!entry\.(\w+)\.unit\.atoms/) {
	    $start = $atomC = 1;
	    $resName = $1;
	    $connect = 0;
	} elsif ($_ =~ /^\!entry\.\w+\.unit\.atomspertinfo/) {
	    $start = 0;
	} elsif ($_ =~ /^\!entry\.(\w+)\.unit\.connectivity/) {
	    $connect = 1;
	    $resName = $1;
	    $start = 0;
	} elsif ($_ =~ /^\!entry\.\w+\.unit\.hierarchy/) {
	    $connect = $start = 0;
	} elsif ($start and $_ =~ /^\s+\"(\S+)\"\s+\"(\S+)\"\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s+(\-?\d+\.\d+)/) {
	    $rec = (
		    {
			"FFTYPE"  => $2,
			"ELEMENT" => substr($2, 0, 1),
			"CHARGE"  => $4,
			"ELENUM"  => $3,
			"ATMNAME" => $1,
			"RESNAME" => $resName,
		    }
		    );
	    $RESDATA{$resName}{$atomC} = $rec;
	    $atomC++;
	} elsif ($connect and exists($RESDATA{$resName}) and
		 $_ =~ /^\s+(\d+)\s+(\d+)\s+(\d+)$/) {
	    if (exists($RESDATA{$resName}{$1}) and exists($RESDATA{$resName}{$2})) {
		$RESDATA{$resName}{$1}{"CONNECT"}{$2} = "";
		$RESDATA{$resName}{$2}{"CONNECT"}{$1} = "";
	    }
	}	
    }
    
    close LIBFILE;

    die "ERROR while parsing AMBER LIB file $inFile\n"
	if (! %RESDATA);
    return \%RESDATA;
}

sub updateBGF {
    my ($atoms, $connections, $resInfo, $startAtm) = @_;
    my ($atomC, $resNum, $myAtm, %ATMS, $isFound, $allRes, $resList);

    $atomC = $startAtm;
    $resNum = $atoms->{$startAtm}{"RESNUM"};
    for $myAtm (keys %{ $atoms }) {
	$ATMS{$myAtm} = " "
	    if ($myAtm < 640);
    }

    while (keys %ATMS and $atomC and exists($BGF_INFO->{$atomC})) {
	$isFound = 0;
	%ERES = (); # Reset the global ERES hash holding the atoms with no match in any of the residues
        $USED{"ATOMS"} = ();
        for $myAtm (keys %{ $atoms }) {
            if (! exists($ATMS{$myAtm})) {
                $USED{"ATOMS"}{$myAtm} = 1;
            }
        }
	$allRes = findElement($atomC, $resInfo); # Search all the residues for current atom's element

#	print "SEARCH STARTED AT $atomC...";
	$BGF_INFO->{$atomC}{"PARENT"} = 0;
	$resList = getResList($atomC, $allRes); # now recursively search the residues matching all the connected atoms
#	print "DONE\n";

	if (keys %{ $resList }) {
	    $resList = getResData($resList); # get the residues with greatest depth
	    if ($resList and keys %{ $resList }) {
		$isFound = 1;
		doUpdate($atomC, $resList, \%ATMS, $resNum);
		$resNum++;
            }
	}

	if (! $isFound) {
	    print "\tSkipping atom $atomC: " . $atoms->{$atomC}{"ATMNAME"} . "\n";
	    delete $ATMS{$atomC};
	    $atoms->{$atomC}{"RESNUM"} = 0;
	}
	$isFound = 0;
	$atomC = getNewAtm(\%ATMS, $atomC);
	$allRes = ();
    }
}

sub doUpdate {
    my ($atomC, $resList, $ATMS, $resNum) = @_;
    my (@RES, $resName, @START, $startP, $j, $currRes, $currAtm, $atmList, $resData);

    @RES = keys %{ $resList };
    $resName = $RES[0];
    if ($#RES > 0) {
	print "\tMultiple residue matches while searching starting at $atomC. Using $resName\n";
    }

    @START = keys %{ $resList->{$resName} };
    $startP = $START[0];
    if ($#START > 0) {
	print "\tMultiple starting points found. Using $startP\n";
    }

    @RES = keys %{ $resList->{$resName}{$startP} };
    $j = 0;
    $BGF_INFO->{$atomC}{"RESNAME"} = $resName;
    $BGF_INFO->{$atomC}{"FFTYPE"} = $RES_INFO->{$resName}{$startP}{"FFTYPE"};
    delete $ATMS->{$atomC};
    $USED{"ATOMS"}{$atomC} = 1;
    $USED{"TMP"}{$atomC} = 1;
    delete $ATMS->{$atomC};

    while ($j <= $#RES) {
        $resData = determineRes($resList->{$resName}{$startP}{$RES[$j]});
        if ( $#{ $resData } > 0) {
	    $j++;
	    next;
	}
        
	$currRes = \%{ $RES_INFO->{$resName}{$resData->[0]} };
	$currAtm = $RES[$j];
	for ("RESNAME", "FFTYPE") {
	    $BGF_INFO->{$currAtm}{$_} = $currRes->{$_};
	}
	$BGF_INFO->{$currAtm}{"RESNUM"} = $resNum;
	$USED{"ATOMS"}{$currAtm} = 1;
	$USED{"TMP"}{$currAtm} = 1;
        delete $ATMS->{$currAtm};
	splice @RES, $j, 1;
    }

    if (@RES) {
	for $j (@RES) {
	    $atmList = determineRes($resList->{$resName}{$startP}{$j});
	    for (@{ $atmList }) {
		next
		    if (exists($USED{"TMP"}{$_}));
		$currRes = \%{ $RES_INFO->{$resName}{$_} };
		$currAtm = $j;
		for ("RESNAME", "FFTYPE") {
		    $BGF_INFO->{$currAtm}{$_} = $currRes->{$_};
		}
                $BGF_INFO->{$currAtm}{"RESNUM"} = $resNum;
		$USED{"ATOMS"}{$currAtm} = 1;
		$USED{"TMP"}{$currAtm} = 1;
		delete $ATMS->{$currAtm};
		last;
	    }
	}
    }

    printFound($resName, \%{ $USED{"TMP"} });
    $USED{"TMP"} = ();
}

sub getResList {
    my ($start, $validRes) = @_;
    my ($currElement, @ATMCONS, $RESLIST, $resIndex, $resName, $resStart, @tmp, $currAtm);
    my ($childList, $resInfo, $maxDepth, $atmC, $i, $parentStr);

    # get the neighbors of this atom
    @ATMCONS = getConns($start);

    if ($#ATMCONS == -1) {
        $RESLIST = saveNode($validRes, $start);
	return $RESLIST; # Exit senario 1
    }

    # find the children of this residue atom that with similar connections
    $currElement = $BGF_INFO->{$start}{"ELEMENT"};
    $childList = determineValidCons($validRes, $currElement, \@ATMCONS, $start);
	
    # Check to see if any of the connections have no match, if so then leaf is terminal, add to ERES
    for $currElement (@ATMCONS) {
	if (! exists($childList->{$currElement->{"ATOM"}})) {
	    $ERES{$currElement->{"ATOM"}} = 1;
	}
    }
    
    # if none of children are similar, return the parent and it's corresponding atom 
    if (! keys %{ $childList }) {
        #$RESLIST = saveNode($validRes, $start);
	#return $RESLIST; # Exit senario 2
        return ();
    }

    for $currElement (@ATMCONS) {
        next
            if (exists($USED{"ATOMS"}{ $currElement->{"ATOM"} }) );
        # update the used atoms/residues
        $USED{"ATOMS"}{$start} = 1;
        
        #remove duplicate atoms from RESLIST and prune childrenlist to exclude residues already used
        ($RESLIST, $childList->{ $currElement->{"ATOM"} }) = syncList($RESLIST, $childList->{ $currElement->{"ATOM"} });
                
        #recursively match the connection list of each present child atom
	($resInfo) = getResList($currElement->{"ATOM"}, \%{ $childList->{$currElement->{"ATOM"}} });
        
        # save the atoms<->residue information of each child
	for $resName (keys %{ $resInfo }) {
	    for $resStart (keys %{ $resInfo->{$resName} }) {
		for $atmC (keys %{ $resInfo->{$resName}{$resStart} }) {
                    $resIndex = $resInfo->{$resName}{$resStart}{$atmC};
                    if (exists($RESLIST->{$resName}{$resStart}{$atmC})) {
                        $RESLIST->{$resName}{$resStart}{$atmC} .= " " . $resIndex;
                    } else {
                        $RESLIST->{$resName}{$resStart}{$atmC} = $resIndex;
                    }                    
		}
	    }
	}
    }
    
    return $RESLIST; #Exit Senario 3

}

sub saveNode {
    my ($validRes, $startAtm) = @_;
    my ($resName, $resIndex, $resStart, $currAtm, $RESLIST, $parentIndex, $parentAtm);

    $RESLIST = ();
    
    for $resName (keys %{ $validRes }) {
	for $resStart (keys %{ $validRes->{$resName} }) {
            for $resIndex( keys %{ $validRes->{$resName}{$resStart} }) {
                if (exists($RESLIST->{$resName}{$resStart}{$startAtm})) {
                    $RESLIST->{$resName}{$resStart}{$startAtm} .= " $resIndex";
                } else {
                    $RESLIST->{$resName}{$resStart}{$startAtm} = $resIndex;
                }
                # save the parent path information
                for $parentIndex (keys %{ $validRes->{$resName}{$resStart}{$resIndex}{"PARENTS"} }) {
                    $parentAtm = $validRes->{$resName}{$resStart}{$resIndex}{"PARENTS"}{$parentIndex};
                    $RESLIST->{$resName}{$resStart}{$parentAtm} = $parentIndex;
                }
	    }
	}
    }
    
    if ($RESLIST) {
        $USED{"ATOMS"}{$startAtm} = 1;
    }
    return $RESLIST;
}

sub determineValidCons {
    my ($resList, $currElement, $ATMCONS, $atmStart) = @_;
    my ($resName, $resIndex, $resStart, $pathStr, $resAtmParent);
    my ($resAtm, $resCon, $currAtom, %RESLIST, $tmp, $parents, $cRes);
    
    for $currAtom (@{ $ATMCONS }) {
        for $resName (keys %{ $resList }) {
            for $resStart (keys %{ $resList->{$resName} }) {
                for $resIndex (keys %{ $resList->{$resName}{$resStart} }) {
                    $pathStr = $resList->{$resName}{$resStart}{$resIndex}{"PATH"};
                    $resAtmParent = \%{ $RES_INFO->{$resName}{$resIndex} };
                    %{ $parents } = %{ $resList->{$resName}{$resStart}{$resIndex}{"PARENTS"} }; 
                    if ($resAtmParent->{"ELEMENT"} eq $currElement and exists($resAtmParent->{"CONNECT"})) {
                        for $resCon (keys %{ $resAtmParent->{"CONNECT"} }) {
                            next
                                if (exists($parents->{$resCon})); # ignore connections to parent
                            $resAtm = \%{ $RES_INFO->{$resName}{$resCon} };
                            if ($resAtm->{"ELEMENT"} eq $currAtom->{"ELEMENT"}) {
                                $cRes = \%{ $RESLIST{$currAtom->{"ATOM"}}{$resName}{$resStart} };
                                $cRes->{$resCon} = (
                                                    {
                                                        "ATOM"    => $currAtom->{"ATOM"},
                                                        "PATH"    => $pathStr . " $resIndex",
                                                    }
                                                   );
                                %{ $cRes->{$resCon}{"PARENTS"} } = %{ $parents };
                                $cRes->{$resCon}{"PARENTS"}{$resIndex} = $atmStart;
                            }
                        }
                    }
                }
            }
        }
    }

    return \%RESLIST;
}

sub getNewAtm {
    my ($AllAtm, $startAtm) = @_;
    my ($aCounter, $retVal);
    
    for $aCounter (sort Numerically keys %ERES) {
	$retVal = $aCounter;
	last;
    }

    if (! defined($retVal)) {
	for $aCounter (sort Numerically keys %{ $AllAtm }) {
	    next
		if ($aCounter < $startAtm);
	    $retVal = $aCounter;
	    last;
	}
   }

    return $retVal;
}
	
sub getConns {
    my ($start) = $_[0];
    my ($i, $element, @ATMCONS);

    for $i (keys %{ $BGF_INFO->{$start}{"CONNECT"} }) {
	next
	    if (exists($USED{"ATOMS"}{$i}));
	$BGF_INFO->{$i}{"FFTYPE"} =~ /(\w+)/;
	$element = $BGF_INFO->{$i}{"ELEMENT"};
	push @ATMCONS, (
			{
			    "ELEMENT" => $element,
			    "ATOM"    => $i,
			}
			);
    }
    
    return @ATMCONS;
}

sub getResData {
    my ($RESLIST) = $_[0];
    my ($maxDepth, $resName, $resStart, $RETLIST, @tmp);
    my ($resSize, %TMP, $resIndex, $i, %RESUSED, $totRes);
    
    $maxDepth = $resSize = 0;

    for $resName (keys %{ $RESLIST }) {
	@tmp = keys %{ $RES_INFO->{$resName} };
	$resSize = $#tmp;
	for $resStart (keys %{ $RESLIST->{$resName} }) {
            %RESUSED = ();
	    @tmp = keys %{ $RESLIST->{$resName}{$resStart} };
            for $i (@tmp) {
                $resIndex = $RESLIST->{$resName}{$resStart}{$i};
                $TMP{$resName}{$resStart}{$i} = $RESLIST->{$resName}{$resStart}{$i};
                $RESUSED{$resIndex} = 1;
                $TMP{$resName}{$resStart}{"COUNT"}++;
            }
            if ($TMP{$resName}{$resStart}{"COUNT"} == $resSize) {
                $TMP{$resName}{$resStart}{"TOTRES"} = 1;
            }
            $maxDepth = $TMP{$resName}{$resStart}{"COUNT"}
                if ($TMP{$resName}{$resStart}{"COUNT"} > $maxDepth);
        }
    }
    
    if ($maxDepth < 3) {
        return ();
    }
   $RESLIST = ();
   for $resName (keys %TMP ) {
        for $resStart (keys %{ $TMP{$resName} }) {
            if ($TMP{$resName}{$resStart}{"COUNT"} == $maxDepth) {
                %{ $RETLIST->{$resName}{$resStart} } = %{ $TMP{$resName}{$resStart} };
                delete $RETLIST->{$resName}{$resStart}{"COUNT"};
                delete $RETLIST->{$resName}{$resStart}{"TOTRES"};
            }
        }
    }
    
    %TMP = ();
    
    return $RETLIST;
}


sub getMaxDepth {
    my ($RESLIST) = $_[0];
    my ($resName, $resStart, $maxDepth, @tmp);

    $maxDepth = 0;
    for $resName (keys %{ $RESLIST }) {
	for $resStart (keys %{ $RESLIST->{$resName} }) {
	    @tmp = keys %{ $RESLIST->{$resName}{$resStart} };
	    $maxDepth = $#tmp
                if ($#tmp > $maxDepth);
	}
    }

    return $maxDepth;
}

sub printFound {
    my ($res, $fList) = @_;
    my ($currAtm, $oldAtm);

    print "\tFOUND $res: ATOMS ";
    for $currAtm (sort Numerically keys %{ $fList }) {
	if (! defined($oldAtm)) { 
	    print "$currAtm";
	    $oldAtm = $currAtm;
	} elsif ( ($currAtm - $oldAtm) == 1) {
	    $oldAtm = $currAtm;
	} else {
	    print " - $oldAtm $currAtm";
	    $oldAtm = $currAtm;
	}
    }

    print " - $oldAtm\n";
}

sub initialize {
    
    for (0 .. 1) {
	FileTester($ARGV[$_]);

    }

    if (! $save_name) {
	$save_name = basename($bgf_file);
	$save_name =~ s/\.bgf$//;
	$save_name .= "_mod.bgf";
    }
    
    die "ERROR: Expected integer for startAtom. Got: $startAtm\n"
	if (! IsInteger(Trim($startAtm)));
    $startAtm = Trim($startAtm);
}

sub Numerically {
    ($a<=>$b);
}

sub findElement {
    my ($atomC, $resList) = @_;
    my ($myElement, $i, $j, $resIndex, $foundRes);

    $myElement = $BGF_INFO->{$atomC}{"ELEMENT"};
    for $i (keys %{ $resList }) {
	for $j (keys %{ $resList->{$i} }) {
	    if (! exists($resList->{$i}{$j}{"ELEMENT"})) {
		delete $resList->{$i}{$j};
		next;
	    }
	    $resIndex = \%{ $resList->{$i}{$j} };
	    if ($resIndex->{"ELEMENT"} eq $myElement) {
		$foundRes->{$i}{$j}{$j} = (
                                           {
                                                "ATOM"    => $atomC,
                                                "PATH"    => "",
                                                "PARENTS" => {},
                                           }
                                          ); 
	    }
	}
    }

    return $foundRes;
}

sub syncList {
    my ($parentList, $childrenList) = @_;
    my ($resName, $resStart, $resIndex, $atmC, $tmp);
    
    for $resName (keys %{ $parentList }) {
        for $resStart (keys %{ $parentList->{$resName} }) {
            for $atmC (keys %{ $parentList->{$resName}{$resStart} }) {
                 $tmp = determineRes($parentList->{$resName}{$resStart}{$atmC});
                $parentList->{$resName}{$resStart}{$atmC} = "";
                for (@{ $tmp }) {
                    $parentList->{$resName}{$resStart}{$atmC} .= Trim($_) . " ";
                }
                chomp $parentList->{$resName}{$resStart}{$atmC};
                if ($#{ $tmp } == 0) {
                    $resIndex = pop @{ $tmp };
                    delete $childrenList->{$resName}{$resStart}{$resIndex};
                }
            }
        }
    }
    
    return ($parentList, $childrenList);
}

sub determineRes {
    my ($resList) = $_[0];
    my (@tmp, %TMP);
    
    @tmp = ();
    %TMP = ();
    if ($resList =~ /\d+\s+\d+/) {
        @tmp = split /\s+/, $resList;
    } else {
        $tmp[0] = Trim($resList);
    }
    for (@tmp) {
        $TMP{$_} = "";
    }
                
    @tmp = keys %TMP;
    
    return \@tmp;
}
