#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use IO::Handle;
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(TrjSelections FileTester GetBondLength CoM LoadFFs CenterText GetStats);
use BOX qw(GetCellFromHeader GetBox);;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use ManipAtoms qw(UnwrapAtoms ImageAtoms GetMols BuildAtomSelectionString 
				SelectAtoms ReimageAtoms GetAtmData);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use constant PI => atan2(1,1);

sub numerically { ($a<=>$b); };
sub alphabetically { ($a cmp $b); };
sub doResidenceTimes;
sub calcResTime;
sub writeDataFileHeaders;
sub groupAtomsByField;
sub writeBINdata;
sub saveCoords;
sub centerOnSolute;
sub getShellID;
sub getSolvShells;
sub getCoM;
sub updateBGF;
sub getChainOffset;
sub init;
sub showUsage;

my ($bgfFile, $shellStr, $saveName, $solvAtomList, $soluAtomList, $reImage, $reCenter, $saveBGF, $numeric);
my ($SELECT, $bin, $trjFile, $field, $pStr, $LAMMPSOPTS, $getByteOffset, $getSnapshots, $trjType, $OUTFILE);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $SOLUTE, $SOLVENT, $WSHELLS, $watShells, $MOLS, $chain_offset);
my ($datFile, $binFile, $BDATA, $nframe, $nshells, $GRPS, $equilPt, $solvTrj, $tstep);

$|++;
$watShells = &init;
print "Getting atom data from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS, $BONDS);
$SOLUTE = SelectAtoms($soluAtomList,$ATOMS);
die "ERROR: No valid solute atoms found!\n" if(! keys%{ $SOLUTE });
$SOLVENT = SelectAtoms($solvAtomList,$ATOMS);
die "ERROR: No valid solvent atoms found!\n" if(! keys%{ $SOLVENT });
$GRPS = &groupAtomsByField($ATOMS,$field,$SOLVENT);
$chain_offset = getChainOffset($ATOMS,$SOLVENT,$numeric);
print "Done\n";
if (! defined($trjFile)) {
    print "Reimaging atoms in box around solute...";
    $SOLUTE = GetAtmData($ATOMS,$SOLUTE);
    $BOX = GetBox($ATOMS, undef, $HEADERS);
    &centerOnSolute($ATOMS, $SOLUTE) if ($reCenter);
    &ReimageAtoms($ATOMS,$BONDS,$MOLS,GetCellFromHeader($HEADERS),$SOLVENT) if ($reImage);
    $SOLVENT = GetMols($ATOMS,$BONDS,$SOLVENT);
    &getCoM($ATOMS,$SOLVENT);
    print "Done\nGetting solvent molecules in shells around solute...";
    $nframe=1;
    $WSHELLS = getSolvShells($SOLUTE,$SOLVENT,$watShells,$BOX,0,1,$bin,\%{ $BDATA });
    print "Done\nCreating $saveName...";
    &updateBGF($ATOMS, $SOLUTE, $SOLVENT, $WSHELLS, $chain_offset, $GRPS, $nshells+1);
    &addHeader($ATOMS, $HEADERS);
    &createBGF($ATOMS, $BONDS, $saveBGF);
    print "Done\n";
} else {
    open $OUTFILE, "> $datFile" or die "ERROR: Cannot write to $datFile: $!\n";
    &writeDataFileHeaders($GRPS, $nshells, $OUTFILE);
    $OUTFILE->autoflush(1);

    $field = scalar keys %{ $ATOMS };
    $getByteOffset->($SELECT, $trjFile, $field);
    if ($trjType == 2) {
        &GetLammpsTrjType($SELECT, $trjFile, "coord", \%{ $LAMMPSOPTS });
        $field = "coord";
    }
    $pStr = "Getting solvent distribution from $trjFile...";
    $getSnapshots->($ATOMS, $trjFile, $SELECT, $field, \&saveCoords, $pStr, $OUTFILE);
    close $OUTFILE;
}
&writeBINdata($BDATA, $bin, $nframe, $binFile) if (defined($bin));
&doResidenceTimes($solvTrj,$tstep,($nshells+1),$saveName) if (defined($tstep));

sub doResidenceTimes {
    my ($trjData,$delT,$nshell,$save) = @_;
    my ($i, $j, $k, $prefix, $resTimes, $tmp, $avg, $stdev);

    $prefix = $save;
    $prefix =~ s/\.\w+$//;
    $prefix .= ".restime.dat";
    print "Writing residence time data to ${prefix}...";
    open RESTIME, "> $prefix" or die "ERROR: Cannot write to $prefix: $!\n";
    $resTimes = calcResTime($trjData,$nshell,$tstep);
    @{ $tmp } = sort alphabetically keys %{ $resTimes };
    printf RESTIME "%-8s","#";
    $j =  $nshell * 12;
    for $i (@{ $tmp }) {
	print RESTIME CenterText($i,$j);
    }
    printf RESTIME "\n%-7s","#";
    for $i (@{ $tmp }) {
	for $j (1 .. $nshell) {
	    printf RESTIME "%12s","Shell${j}";
	}
    }
    print RESTIME "\n";
    for $i (1 .. $nshell) {
	printf RESTIME "%-8s","Shell${i}";
	for $k (@{ $tmp }) {
	    for $j (1 .. $nshell) {
		($avg,$stdev) = (0,0);
		($avg,$stdev) = ($resTimes->{$k}{$i}{$j}{STATS}{AVG},$resTimes->{$k}{$i}{$j}{STATS}{STDEV})
				if(exists($resTimes->{$k}{$i}{$j}{STATS}));
		printf RESTIME "%6.2f,%-5.2f",$avg,$stdev;
	    }
	}
	print RESTIME "\n";
    }
    close RESTIME;
    print "Done\n";
}

sub calcResTime {
    my ($data, $nshell, $delT) = @_;
    my ($i, $j, $l, $sStr, $k, $resTime);
  
    for $i (keys %{ $data }) {
	for $j (1 .. $nshell) {
	    for $l (1 .. $nshell) {
		next if ($j == $l);
		$sStr = "(${j},${j},\[${j},\]\*,${j},${j},)${l}";
		for $k (keys %{ $data->{$i} }) {
		    while($data->{$i}{$k} =~ /$sStr/g) {
			push @{ $resTime->{$i}{$j}{$l}{DATA} }, length($1)*$delT;
		    }
		}
	    }
	}
    }
 
    for $i (keys %{ $resTime }) {
	for $j (keys %{ $resTime->{$i} }) {
	    for $k (keys %{ $resTime->{$i}{$j} }) {
		$resTime->{$i}{$j}{$k}{STATS} = GetStats($resTime->{$i}{$j}{$k}{DATA});
	    }
	}
    }

    return $resTime;
}

sub writeDataFileHeaders {
    my ($grps, $num, $outfile) = @_;
    my ($i, $j, $tmp);

    printf $outfile "%-10s","#";
    @{ $tmp } = sort alphabetically keys %{ $grps };
    for $i (@{ $tmp }) {
	$j = ($num+1)*11;
	print $outfile CenterText($i,$j);
    }
    printf $outfile "\n%10s","#Step";
    for $i (@{ $tmp }) {
	for $j (0 .. $num) {
	    printf $outfile "%11s","Shell" . ($j+1);
	}
    }
    print $outfile "\n";
}

sub groupAtomsByField {
    my ($atoms, $field, $select) = @_;
    my ($i, $grpid, $grps);

    if(defined($field)) {
	for $i (keys %{ $select }) {
	    $grps->{ $atoms->{$i}{$field} }{$i} = 1;
	    $atoms->{$i}{GRPID} = $atoms->{$i}{$field};
	}
    } else {
	for $i (keys %{ $select }) {
	    $grps->{ALL}{$i} = 1;
	    $atoms->{$i}{GRPID} = "ALL";
	}
    }
    return $grps;
}

sub writeBINdata {
    my ($bins, $binSize, $nframe, $bfile) = @_;
    my ($i, $j, $tot, $normalized, $count, $list, $dr, $rlo, $rhi, $nbins, $tmp);

    $tot = (); 
    $nbins = 0;
    print "Writing bin data to $bfile...";
    open BFILE, "> $bfile" or die "ERROR: Cannot write to $bfile: $!\n";
    printf BFILE "#%-11s","";
    @{ $tmp } = sort alphabetically keys %{ $bins };
    for $i (@{ $tmp }) {
	print BFILE CenterText($i,36);
    }
    printf BFILE "\n%-12s","#Distance";
    for $i (@{ $tmp }) {
	printf BFILE "%12s%12s%12s","Probability","Normalized","Integrated";
    }
    print BFILE "\n";
    for $i (keys %{ $bins }) {
	@{ $list } = sort numerically keys %{ $bins->{$i} };
	$nbins = $list->[$#{ $list }] if ($list->[$#{ $list }] > $nbins);
    }
    for $i (0 .. $nbins) {
	$rlo = ($i+1) * $binSize;
	$rhi = $rlo + $binSize;
	$dr = ($rhi**3 - $rlo**3);
	printf BFILE "%-12.3f",$rlo;
	for $j (@{ $tmp }) {
	    $count = 0;
	    $count = $bins->{$j}{$i}/$nframe if(exists($bins->{$j}{$i}));
	     $tot->{$j} += $count;
	     $normalized = $count/(4*PI*$dr*.6023*$nframe)/3;
	     printf BFILE " %11.5g %11.5g %11.5g",$count,$normalized,$tot->{$j};
	}
	print BFILE "\n";
    }
    close BFILE;
    print "Done\n";
}

sub saveCoords {
    my ($atoms, $box, $frameNum, $fileHandle) = @_;
    my ($i, $j, $cMol, $shell_data, $solu, $solv, $count, $CENTER, $tmp);

    $nframe++;
    if (defined($trjFile) and $trjType == 2) { #LAMMPS
        $frameNum = $atoms->{TIMESTEP}[0];
        $box = ConvertLammpsBox($atoms->{"BOX BOUNDS"});
        %{ $box->{X} } = %{ $box->{XCOORD} };
        %{ $box->{Y} } = %{ $box->{YCOORD} };
        %{ $box->{Z} } = %{ $box->{ZCOORD} };
        $atoms = $atoms->{ATOMS};
        if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
            UnwrapAtoms($atoms,  $box, $LAMMPSOPTS->{scaled});
        }
    } elsif (defined($trjFile)) {
        $box = ConvertAmberBox(\%{ $box });
    }
    $solu = GetAtmData($atoms,$SOLUTE);
    $solv = GetMols($atoms,$BONDS,$SOLVENT);
    &getCoM($atoms,$solv);
    $shell_data = getSolvShells($solu,$solv,$watShells,$box,$equilPt,$nframe,$bin,\%{ $BDATA });
    printf $OUTFILE "%-10d",$frameNum;
    @{ $tmp } = sort alphabetically keys %{ $GRPS };
    for $j (@{ $tmp }) {
	for $i (0 .. $nshells) {
	    $count = 0;
	    $count = scalar(@{ $shell_data->{$j}[$i] }) if (exists($shell_data->{$j}) and 
							$i <= $#{ $shell_data->{$j} } and 
							defined($shell_data->{$j}[$i]));
	    printf $OUTFILE " %10d", $count;
	}
    }
    printf $OUTFILE "\n";
    $solu = ();
    $solv = ();
    $atoms = ();
    $shell_data = ();
    $tmp = ();
}

sub centerOnSolute {
    my ($atoms, $solute) = @_;
    my ($com, $i, $j);

    $com = CoM($solute);
    for $i (keys %{ $atoms }) {
	for $j ("XCOORD", "YCOORD", "ZCOORD") {
	    $atoms->{$i}{$j} += $com->{$j};
	}
    }
}

sub getShellID {
    my ($shellDist, $shells) = @_;
    my ($i, $shellID);

    for $i (0 .. $#{ $shells }) {
	if($shellDist<$shells->[$i]) {
	    $shellID = $i;
	    last;
	}
    }
    $shellID = scalar(@{ $shells }) if (!defined($shellID));
    $shellID++;
}

sub getSolvShells {
    my ($solu, $solv, $shells, $box,$startFrame,$currFrame,$bsize,$bdata) = @_;
    my ($index, $d_curr, $dist, $i, $j, $shellID, $shellList, $grpid);

    for $i (keys %{ $solv }) {
	$dist = 9999999999;
	$grpid = $solv->{$i}{GRPID};
	for $j (keys %{ $solu }) {
	    $d_curr = GetBondLength($solv->{$i}{CoM},$solu->{$j},$box);
	    $dist = $d_curr if($d_curr < $dist);
	}
	$shellID = getShellID($dist,$shells);
	if(defined($bsize) and $currFrame >= $startFrame) {
	    $index = int($dist/$bsize);
	    $bdata->{$grpid}{$index}++;
	}
	push @{ $shellList->{$grpid}[$shellID] }, $i;
	$solv->{$i}{SHELL} = $shellID;
	$solvTrj->{$grpid}{$i} .= "$shellID,";
    }

    return $shellList;
}

sub getCoM {
    my ($atoms, $mols) = @_;
    my ($com, $moldata, $i, $ret, $aid, $tmp);

    for $i (keys %{ $mols }) {
	$moldata = GetAtmData($atoms, $mols->{$i}{MEMBERS});
	$mols->{$i}{CoM} = CoM($moldata);
	@{ $tmp } = keys %{ $mols->{$i}{MEMBERS} };
	$aid = pop @{ $tmp };
	$mols->{$i}{GRPID} = $ATOMS->{$aid}{GRPID};
    }
}

sub updateBGF {
    my ($atoms, $solute, $solvent, $shells, $chain_start, $grps, $nshells) = @_;
    my ($i, $j, $k, $chain);

    for $i (sort alphabetically keys %{ $grps }) {
	for $j (keys %{ $solvent }) {
	    next if ($solvent->{$j}{GRPID} ne $i);
	    for $k (keys %{ $solvent->{$j}{MEMBERS} }) {
		$chain = $chain_start;
		for (1 .. $solvent->{$j}{SHELL}) { $chain++; }
		$atoms->{$k}{CHAIN} = $chain;
	    }
	}
	for (1 .. $nshells) { $chain_start++ };
    }
}

sub getChainOffset {
    my ($atoms, $solventAtoms, $isnumeric) = @_;
    my ($i, $tmp, $tmp1, $start);

    if(!$isnumeric) { #alphanumeric
	$tmp = ();
	$start = "A";
	for $i (keys %{ $atoms }) {
	    next if (exists($solventAtoms->{$i}));
	    next if ($atoms->{$i}{CHAIN} =~ /\d/);
	    $tmp->{ $atoms->{$i}{CHAIN} } = 1;
	}
	@{ $tmp1 } = sort alphabetically keys %{ $tmp } if(keys %{ $tmp });
    } else { #numeric
	$tmp = ();
	$start = 1;
	for $i (keys %{ $atoms }) {
	    next if (exists($solventAtoms->{$i}));
	    next if ($atoms->{$i}{CHAIN} !~ /\d/);
	    $tmp->{ $atoms->{$i}{CHAIN} } = 1;
	}
	@{ $tmp1 } = sort numerically keys %{ $tmp } if(keys %{ $tmp });
    }

    if(defined($tmp1)) {
	$start = $tmp1->[$#{ $tmp1 }];
	$start++;
    }
    return $start;
}

sub init {
    my (@tmp, $i, $j, $tmp, $SHELL, %OPTS, $usage, $solvStr, $soluStr, $tSel, $list);

    getopt('bswuvcinltrdgem',\%OPTS);

    ($bgfFile,$shellStr,$saveName,$solvStr,$soluStr,$reImage,$reCenter,
     $numeric,$trjFile,$trjType,$tSel,$bin,$field,$equilPt,$tstep) = 
	($OPTS{b},$OPTS{w},$OPTS{s},$OPTS{v},$OPTS{u},$OPTS{i},$OPTS{c},
	$OPTS{n},$OPTS{l},$OPTS{t},$OPTS{r},$OPTS{d},$OPTS{g},$OPTS{e},$OPTS{m});
    $usage = &showUsage;
    for ("b", "w") {
	die "$usage\n" if (! $OPTS{$_});
    }
    
    print "Initializing...";
    FileTester($bgfFile);
    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$/.solvshell\.grps/;
    }
    undef($field) if (defined($field) and $field !~ /^CHAIN|RESNUM|RESID|RES|MOL|MOLECULE|MOLID|MOLECULEID|MOLSIZE|NUMBONDS|FFTYPE$/i);
    if(defined($field)) {
	$field = uc($field); 
	$field = "RESNUM" if ($field =~ /RESID/); 
	$field = "MOLECULEID" if ($field =~ /MOL/ && $field !~ /SIZE/);
	$field = "MOLSIZE" if ($field =~ /SIZE/);
    }
    $equilPt = 0 if (!defined($equilPt) or $equilPt !~ /^\d+/);
    if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
	undef $tstep if(defined($tstep) and $tstep !~ /^\d+\.?\d*$/);
	$datFile = $saveName;
	$datFile =~ s/\.\w+$//;
	$datFile .= ".trj.dat";
        if (! defined($trjType)) {
            if ($trjFile =~ /\.lammps/) {
                $trjType = "lammps";
            } else {
                $trjType = "amber";
            }
        }

        if (lc($trjType) ne "lammps") {
            $trjType = 1;
            $getSnapshots = \&ParseAmberTrj;
            $getByteOffset = \&GetAmberByteOffset;
        } else {
            $trjType = 2;
            $getSnapshots = \&ParseLAMMPSTrj;
            $getByteOffset = \&GetLammpsByteOffset;
        }
        if (! defined($tSel)) {
            $tSel = "*";
        }
        $list = TrjSelections($tSel);
        for $i (keys %{ $list }) {
            $SELECT->{$i} = $list->{$i};
        }
        die "ERROR: No valid frames selected with selection $tSel!\n"
            if (! keys %{ $SELECT } and $tSel ne "*");
    }
    $bin = 0.1 if (defined($bin) and $bin !~ /^\d+\.?\d*$/); 
    if(defined($bin)) {
	$binFile = $saveName;
	$binFile =~ s/\.\w+$//;
	$binFile .= ".distrib.dat";
    }  
    if ($shellStr =~ /\s+/) {
	@tmp = split /\s+/, $shellStr;
    } else {
	$tmp[0] = $shellStr;
    }

    for $i (@tmp) {
	$tmp->{$i} = 1 if ($i =~ /(\d+\.*\d+)/); 
    }
    die "Expected at least one decimal(integer) value for the shell limts. Got \"$shellStr\"\n" if (! keys %{ $tmp });
    @{ $SHELL } = sort numerically keys %{ $tmp };

    print "total shells: " . ($#tmp + 1) . "...Done\n";
    $nshells = $#tmp +1;
    $solvStr = "resname eq 'WAT'" if(!defined($solvStr));
    $soluStr = "moleculeid==1" if (!defined($soluStr));
    $solvAtomList = BuildAtomSelectionString($solvStr);
    $soluAtomList = BuildAtomSelectionString($soluStr);

    $reImage = 1 if (! defined($reImage) or $reImage !~ /0|no/i);
    $reImage = 0 if ($reImage =~ /0|no/i);
    $reCenter = 0 if (! defined($reCenter) or $reCenter !~ /1|yes/i);
    $reCenter = 1 if ($reCenter =~ /1|yes/i);
    $numeric = 0 if (! defined($numeric) or $numeric !~ /1|yes/i);
    $numeric = 1 if ($numeric =~ /1|yes/i);

    $saveBGF = $saveName;
    $saveBGF =~ s/\.\w+$//;
    $saveBGF .= ".bgf";

    return $SHELL;
}

sub showUsage {

    my ($usage) = <<DATA;
This script will either break up the bgf file into:  solute, ions, shell 1...n, and bulk or plot distribution functions from trajectory
usage: $0 -b bgf file le -w shell distance(s) -s [save name] -c [recenter] -i [reimage] -u [solute] -v [solvent]
	options:
	-b bgf file: the location of the bgf file
	-w shell distances(s): the distance from the solute surface of the 1st solvent shell.
	   you can specify multiple shell distances by enclosing them in "" qoutes
	-s [save name]: (optional) the name to save the group information generated by the script
	   if not specified, will be bgf_file.grps. Each solvent shell will be savename_shellx.bgf
	-c [recenter]: (optional) place the center of mass of the solute at the box center. Default no
	-i [reimage]: (optional) wrap solvent atoms back into simulation cell. Default yes
	-u [solute]: (optional) solute atom selection. Default "moleculeid == 1"
	-v [solvent]: (optional) solvent atom selection. Default "resname eq 'WAT'"
	-n [numeric chain]: (optional). If 0 then chain will be 1,2,3... If 1 then chain will be A,B...
					default: 0
	-l [trajectory file]: (optional). Location of trajectory file
	-t [trajectory type]: (optional). Either Lammps or AMBER
	-r [trajecory selection]: (optional). Selection string for trajectory.
	-d [binsize]: (optional): Size of bins for calculating distributions.
	-g [group field]: (optional). Field to group bined atoms. Default none
	-e [equil frame]: (optional). Equilibration point. Default is frame 1
	-m [trj save timestep]: (optional). Timestep between successive snapshots in trajectory to calc residence time
DATA

    return $usage;
}
