#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use General qw(LoadFFs ReadFFs);
use Getopt::Std qw(getopt);
use AMBER qw(AmberLib);
use File::Basename qw(basename);

sub init;
sub showUsage;
sub getCellOpts;
sub getSolvOpts;
sub embedMols;
sub findLigand;
sub typeLigand;
sub removeHydrogens;
sub runLeap;
sub createMembrane;
sub convertWAT;
sub createLAMMPSfiles;
sub getAtomRange;
sub getNewAtomRange;
sub execCmd;
sub getCellInfo;
sub addIons;
sub removeTmpFiles;
sub removeBadContacts;
sub checkAtomTypes;
sub groupMols;
sub assignChain;

my ($sfile, $solute, $tmp, $tmp1, $ffOpt, $PARMS);
my ($memOpts, $solvOpts, $ffields, $LIBS, $leapPrep, $LIGAND, $soluRange, $tmp2);

$|++;
&init;
print "Gettting data from solute $sfile->{type} file $sfile->{name}...";
$solute->read($sfile->{name}, $sfile->{type});
print "Done\n";
if (! $ffOpt) {
    print "Searching for ligand...";
    $LIGAND = &findLigand($solute, $LIBS, $sfile);
    print "Done\n";
    $leapPrep = &typeLigand($solute, $LIGAND, $sfile) if (defined($LIGAND));
    print "Assigning AMBER atom types to system...";
    $sfile->{pdb} = &removeHydrogens($sfile->{name}, $sfile->{type});
    $sfile->{pdb} = $sfile->{name} if ($sfile->{type} eq "pdb");
    print "Done\n";
    ($tmp,$tmp1) = &runLeap($sfile->{pdb},$leapPrep, $LIGAND);
    $tmp2 = $tmp;
} else {
    $PARMS = LoadFFs($sfile->{ff});
    print "Checking atom types...";
    &checkAtomTypes($solute,$PARMS);
    $tmp = $sfile->{name};
    $tmp2 = $sfile->{save};
    $tmp2 =~ s/\.\w+$/\.bgf/;
    print "Done\n";
}
print "Fixing structure...";
&assignChain($tmp2,$tmp1);
print "Done\nEmbedding system into $memOpts->{size} $memOpts->{type} membrane...";
$sfile->{save} =~ s/\.\w+$/\.bgf/;
$tmp1 = &createMembrane($memOpts);
&embedMols($tmp, $tmp1, $ffields, $sfile->{save});
print "Done\nRemoving bad contacts...";
#&removeBadContacts($sfile->{save}, $ffields);
print "Done\nNeutralizing system...";
&addIons($sfile->{save}, $ffields);
print "Creating $sfile->{save}...";
&convertWAT($sfile->{save}, $solvOpts) if (defined($solvOpts));
&groupMols($sfile->{save});
print "Done\nCreating LAMMPS input and data files...";
&createLAMMPSfiles($sfile->{save}, $ffields);
print "Done\n";
&removeTmpFiles($sfile);

sub assignChain {
    my ($savename, $ligRange) = @_;
    my ($modifyCmd);

    $modifyCmd = "$Bin/modifyAtomData.pl -s $savename -w $savename -a 'index>0' -f 'CHAIN:A'"
;
    &execCmd($modifyCmd);
    if(defined($ligRange)) {
        $modifyCmd = "$Bin/modifyAtomData.pl -s $savename -w $savename -a \"resname eq 'RES'\" -f 'CHAIN:B'";
        &execCmd($modifyCmd);
    }
}

sub groupMols {
    my ($sfile) = $_[0];
    my ($cmdStr);

    $cmdStr = "$Bin/groupAtoms.pl -b $sfile -s $sfile -f chain";
    &execCmd($cmdStr);
}

sub removeTmpFiles {
    my ($sfile) = $_[0];
    my ($prefix);

    print "Removing temporary files...";
    $prefix = $sfile->{name};
    $prefix =~ s/\.\w+$//;
    
    &execCmd("rm -fr ???.frcmod ???.prepi _tmplig.* ${prefix}_noH.* leaprc membrane.bgf _gpcrmembranemd.* leap.log _gpcrmembranecharge.dat ${prefix}_amber.bgf");
    print "Done\n";
}

sub addIons {
    my ($fileName, $ffields) = @_;
    my ($chargeCmd, $neutCmd, $charge, $rCharge, $delCharge, $addChargeCmd);

    $chargeCmd = "$Bin/bgfcharge.pl $fileName";
    $addChargeCmd = "$Bin/modifyAtomData.pl -s $fileName -w $fileName -a 'index == 1' -f \"CHARGE:";
    $neutCmd = "$Bin/addIons.pl -b $fileName -n 0 -f \"$ffields\" -w $fileName -s \"resname eq 'WAT' and chain eq 'W'\" -i ";
    &execCmd($chargeCmd, "_gpcrmembranecharge.dat");
    open CHARGECMD, "_gpcrmembranecharge.dat" or die "ERROR: Cannot execute $chargeCmd: $!\n";
    while (<CHARGECMD>) {
	chomp;
	if ($_ =~ /Total Charge: (\-?\d+\.?\d*)/) {
	    $charge = $1;
	}
    }
    close CHARGECMD;
    die "ERROR: No valid data read while executing cmd $chargeCmd!\n" if (! defined($charge));
    $rCharge = sprintf("%.0f",$charge);
    $delCharge = $rCharge - $charge;
    if (abs($delCharge) > 0.00001) {
	$addChargeCmd .= "+" if ($delCharge > 0);
	$addChargeCmd .= $delCharge;
	&execCmd("${addChargeCmd}\"");
    }
    if ($rCharge > 0) {
	$neutCmd .= "Cl";
    } else {
	$neutCmd .= "Na";
    }
    &execCmd($neutCmd);
}

sub createLAMMPSfiles {
    my ($saveName, $ffields) = @_;
    my ($lammpsStr, $prefix);

    $prefix = basename($saveName);
    $prefix =~ s/\.\w+$//;
    
    $lammpsStr = "$Bin/createLammpsInput.pl -f \"$ffields\" -b $saveName -t gpcr -s $prefix";
    &execCmd($lammpsStr);
}

sub convertWAT {
    my ($filename, $opts) = @_;
    my ($convertStr, $i, $tmp, $j);

    print "Converting water model...";
    $convertStr = "$Bin/modifyAtomData.pl -s $filename -w $filename";
    for $i (keys %{ $opts }) {
	$tmp = "-a \"fftype eq '${i}'\" -f \"";
	for $j (keys %{ $opts->{$i} }) {
	    next if (! exists($opts->{$i}{$j}));
	    $tmp .= "${j}:$opts->{$i}{$j} ";
	}
	&execCmd("${convertStr} ${tmp} \"");
    }
}

sub removeBadContacts {
    my ($bgffile, $ffields) = @_;
    my ($badContactStr);

    $badContactStr = "$Bin/getLammpsAtomEng.pl -b $bgffile -s _tmp_elec_pot.dat -f \"$ffields\"";
    &execCmd($badContactStr);
    $badContactStr = "$Bin/removeBadContact.pl -b $bgffile -s $bgffile -o _tmp_elec_pot.dat";
    &execCmd($badContactStr);
    &execCmd("rm -fr _tmp_elec_pot.dat");
}
   
sub createMembrane {
    my ($options) = $_[0];
    my ($CELL, $i, $DIMS, $replicateStr, $memFile, $shouldReplicate, $trimStr, $multiple);

    $memFile = "$Bin/dat/membranes/" . $options->{type} . "GAFFspc.bgf";
    $replicateStr = "$Bin/replicate.pl -i 0 -c 1 -b $memFile -s ./membrane.bgf -d \"";
    $trimStr = "$Bin/trimCell.pl -b ./membrane.bgf -s _gpcrmembranemd.bgf -m 1 ";    
    $CELL = getCellInfo($options, $memFile);

    $shouldReplicate = 0;
    if ($options->{type} ne "popc" or $options->{size} ne "97x97") {
	for $i ("x", "y") {
	    $DIMS->{$i} = 1;
	    if($options->{$i}>$CELL->{$i}){
	  	$multiple = int($options->{$i}/$CELL->{$i})+1;
	  	$DIMS->{$i} = $multiple;
	  	$shouldReplicate = 1;
	    }
	    $replicateStr .= "$DIMS->{$i} ";
	}
	$replicateStr .= "1\"";
    }    
    if($shouldReplicate) {
	&execCmd($replicateStr);
	$trimStr .= "-c \"$options->{x} $options->{y} $CELL->{z}\"";
	&execCmd($trimStr);
    } else {
	&execCmd("cp $memFile ./_gpcrmembranemd.bgf");
    }
    return "_gpcrmembranemd.bgf";
}

sub runLeap {
    my ($pdbFile, $leapAdd, $ligand) = @_;
    my ($leaprcLoc, $leapCmd, $molname, $convertCmd, $tlist); 
    my ($saveName, $modifyCmd, $i, $atomRange, $tmp, $saltBridgeCmd, $saltBridgeAdd);
    
    $molname = $pdbFile;
    $molname =~ s/\.\w+//;
    $saveName = $molname;
    $saveName =~ s/_noH//;
    $saveName .= "_amber.bgf";

    $modifyCmd = "$Bin/modifyAtomData.pl -s $pdbFile -a \"resname eq 'HSE'\" -f 'RESNAME:HIE' -w $pdbFile -t pdb";
    &execCmd($modifyCmd,0);
    $modifyCmd = "$Bin/modifyAtomData.pl -s $pdbFile -a \"resname eq 'HDD'\" -f 'RESNAME:HID' -w $pdbFile -t pdb";
    &execCmd($modifyCmd,0);
    if (defined($ligand)) {
	for $i (keys %{ $ligand }) {
            ($atomRange, $tmp) = getNewAtomRange($ligand->{$i});
	    $tlist .= "$tmp ";
            $modifyCmd = "$Bin/modifyAtomData.pl -s $pdbFile -a '" . $atomRange . "' -f 'LABEL:HETATM' -w $pdbFile -t pdb";
            &execCmd($modifyCmd);
	}
        $modifyCmd = "$Bin/splitPDBwLigand.pl -p $pdbFile -s $pdbFile -t \"$tlist\"";
        &execCmd($modifyCmd);
    }
    #find salt bridges
    $saltBridgeCmd = "$Bin/saltBridge.pl -f $pdbFile -s ${molname}_saltbridge.leap";
    &execCmd($saltBridgeCmd,0);
    $saltBridgeAdd = "";
    $saltBridgeAdd = "source ${molname}_saltbridge.leap\n" if (-e "${molname}_saltbridge.leap");
    system("sed -i '/CONECT/d' $pdbFile");
    $leaprcLoc = "$ENV{AMBERHOME}/dat/leap/cmd/leaprc_gpcr";
    $leapCmd = "$ENV{AMBERHOME}/exe/tleap";
    $convertCmd = "$Bin/amber2bgf.pl ${molname}.prmtop ${molname}.inpcrd $saveName";
    &execCmd("cp $leaprcLoc ./leaprc");
    open MYLEAPRC, ">> leaprc" or die "Cannot write to leaprc: $!\n";
    print MYLEAPRC "gaff = loadamberparams gaff.dat\n";
    print MYLEAPRC $leapAdd if (defined($leapAdd));
    print MYLEAPRC "prot = loadpdb $pdbFile\n";
    print MYLEAPRC $saltBridgeAdd;
    print MYLEAPRC "saveamberparm prot ${molname}.prmtop ${molname}.inpcrd\n";
    print MYLEAPRC "quit\n";
    close MYLEAPRC;
    &execCmd("$leapCmd");
    print "Converting AMBER files to bgf...";
    &execCmd($convertCmd);
    print "Done\n";
    &execCmd("rm -fr ${molname}_saltbridge.leap",0) if (-e "${molname}_saltbridge.leap");
    return ($saveName,$atomRange);
}

sub removeHydrogens {
    my ($fileName, $fileType) = @_;
    my ($cmd, $saveName); 

    print "Removing all hydrogens and converting to pdb...";
    $saveName = basename($fileName);
    $saveName =~ s/\.\w+$//;
    $saveName .= "_noH.pdb";

    if($fileType eq "bgf") {
	$cmd = "/project/Biogroup/scripts/perl/bgf2pdb_noH.pl";
        &execCmd("${cmd} $fileName", $saveName);
	$fileName = $saveName;
    } 

    $cmd = "$Bin/removePDBh.pl -p $fileName -s $saveName -b 0 >> junk";
    &execCmd($cmd);
    return $saveName;
}

sub typeLigand {
    my ($solu, $ligands, $filedata) = @_;
    my ($antechamberCmd, $i, $extractCmd, $resname, $atomRange, $convertCmd);
    my ($charge, $RES, $topo, $leapAdd, $saveName, $parmchkCmd, $resParm);

    print "Assigning atom types for ligands...";
    $antechamberCmd = "/exec/amber9/exe//antechamber -fo prepi -c am1 -at gaff -pf y";
    $extractCmd = "$Bin/getBGFAtoms.pl -b $filedata->{name} -s _tmplig.bgf -o ";
    $convertCmd = "$Bin/bgf2mol2.pl -v 0 -b _tmplig.bgf -a 1";
    $parmchkCmd = "/exec/amber9/exe//parmchk -f prepi";
    for $resname (keys %{ $ligands }) {
	$charge = 0;
	$RES = uc $resname;
	$RES = substr($resname, 0, 3) if (length($resname) > 3); 
        $resParm = lc $RES . "Parm";
	$topo = "${resname}.res";
	$saveName = "${resname}.prepi";
	($atomRange, undef) = getNewAtomRange($ligands->{$resname});
	for $i (values %{ $ligands->{$resname} }) {
	    $charge += $i->charge;
	}
	$charge = 0 if (abs($charge) < 0.01);
	&execCmd("${extractCmd} '$atomRange'");
	&execCmd($convertCmd);
	&execCmd("${antechamberCmd} -fi mol2 -i _tmplig.mol2 -nc $charge -rn $RES -rf $topo -o $saveName");
        &execCmd("${parmchkCmd} -i ${RES}.prepi -o ${RES}.frcmod");
	$leapAdd .= "loadamberprep $saveName\n";
        $leapAdd .= "$resParm = loadamberparams ${RES}.frcmod\n";
    }
    
    print "Done\n";
    return $leapAdd;
}


sub findLigand {
    my ($soluFile, $amberLibs) = @_;
    my ($resname, $LIGANDS, $count, $i, @atomList, $atomId);

    $count = 0;
    for $i (keys %{ $soluFile->shash->{resid} }) {
        @atomList = keys %{ $soluFile->shash->{"resid"}{$i} };
        next if (! @atomList);
        $atomId = shift @atomList;
        next if (! $atomId or ! $soluFile->shash->{"resid"}{$i}{$atomId});
        $resname = $soluFile->shash->{"resid"}{$i}{$atomId}->resname;	
        $resname = "HIE" if ($resname =~ /his|hse|hdd/i);
        $resname = "CYS" if ($resname =~ /cyx/i);
        $resname = "TYR" if ($resname =~ /ty./i);
	if (! exists($amberLibs->{uc $resname}) and ! exists($LIGANDS->{$resname})) {
	    $LIGANDS->{$resname} = $soluFile->shash->{resid}{$i};
	    $count++;
	}
    }
    print "none found..." if (! $count);
    print "found $count..." if ($count);

    return $LIGANDS;
}

sub embedMols {
    my ($soluFile, $solvFile, $ffields, $saveName) = @_;
    my ($embedStr);

    $embedStr = "$Bin/embedGPCR.pl -m $solvFile -s $soluFile -f \"$ffields\" -c com -w $saveName";
    &execCmd($embedStr);
    &execCmd("rm -fr _out _solv_1cell.bgf _solv_replicate.bgf _solv_trim.bgf");
}

sub createSolventBox {
    my ($solv, $solu, $cell) = @_;
    my ($i, $replicate, $blen, $trim, $box, $bmin, $smin, $offset, $map, $solvBox, $cellScale);

    $map = ({ "a" => "x", "b" => "y", "c" => "z"});
    if (defined($cell->{density})) { #compress/expand the solvent cell to the new density. assume 1 g/cm3
	$cellScale = 1/($cell->{density}**(1/3)); #
	for $i ("x", "y", "z") {
	    $solv->stressCell($i, $cellScale);
	}
    }
    $smin = $solv->getExtrema("min");
    %{ $solvBox } = %{ $solv->vdwbox };
    %{ $solvBox } = %{ $solv->cell } if ($solv->cell->{valid});
    %{ $box } = %{ $solu->vdwbox };
    $replicate = "$Bin/replicate.pl -b ./_solv_1cell.bgf -d \"";
    for $i ("a", "b", "c") {
	$box->{$i}{max} += $cell->{cell}{$i}{max} 
		if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{max}));
        $box->{$i}{min} -= $cell->{cell}{$i}{min} 
		if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{min}));
	$box->{$i}{len} = $box->{$i}{max} - $box->{$i}{min};
	$bmin->{$i} = $box->{$i}{min};
	$blen .= "$box->{$i}{len} ";
	$replicate .= sprintf("%.0f ", (($box->{$i}{len}/$solvBox->{$i}))+1);
    }
    $replicate .= "\" -s _solv_replicate.bgf";
    #move the solvent box to the solute minima
    for $i ("a", "b", "c") {
	$offset->{ $map->{$i }} = $bmin->{$i} - $smin->{$i};
    }
    $solv->moveMol("all", $offset);
    $solv->write("_solv_1cell.bgf", "bgf");
    #replicate the cell by the replication vector calculated above
    die "ERROR while executing \"$replicate\"\n" if (system("${replicate} >& _out.dat"));
    # remove all molecules outside the solute (inflated) cell
    $trim = "$Bin/trimCell.pl -b _solv_replicate.bgf -c \"$blen\" -s _solv_trim.bgf -m 1 -o 2";
    die "ERROR while executing \"$trim\"\n" if (system("${trim} >& _out.dat"));
    undef($solv);
}

sub checkAtomTypes {
    my ($sol, $parms) = @_;
    my ($i);

    for $i (keys %{ $sol->shash->{fftype} }) {
        die "ERROR: FFTYPE \"$i\" not found in forcefield(s)!\n"
	    if (!exists($parms->{ATOMTYPES}{$i}));
    }
}

sub init {
    my (%OPTS, $ffStr, $memTypeStr, $solvTypeStr, $periodStr, $memSizeStr, $solFF);

    $ffOpt = 0;
    getopt('idmswf', \%OPTS);
    die &showUsage . "\n" if (! exists($OPTS{i}));

    print "Initialzing...";
    $ENV{AMBERHOME} = "/home/tpascal/codes/amber11" if (!exists($ENV{AMBERHOME}));
    ($sfile->{name}, $memSizeStr, $memTypeStr, $solvTypeStr, $sfile->{save}, $ffStr) = 
	($OPTS{i}, $OPTS{d}, $OPTS{m}, $OPTS{s}, $OPTS{w}, $OPTS{f});
    $solute =  MolData->new();
    $solute->testFile($sfile->{name});
    if ($sfile->{name} =~ /\.(\w+)$/) {
	$sfile->{type} = lc $1;
    } else {
	die "ERROR: Cannot determine file type from $sfile->{name}!\n";
    }
    $memTypeStr = "popc" if (! defined($memTypeStr) or $memTypeStr !~ /^(popc|pope)$/i);
    $memSizeStr = "97 97" if (! defined($memSizeStr));
    $memOpts = getCellOpts($memTypeStr, $memSizeStr);

    $solFF = "AMBER03.ff";
    if(defined($ffStr) && $sfile->{type} eq "bgf") {
	($sfile->{ff}, undef) = ReadFFs($ffStr,undef);
	$solFF = $ffStr;
	$ffOpt = 1;
    } elsif (defined($ffStr)) {
	print "cannot use custom forcefiled if structure is not bgf...";
    }

    $ffields = "$solFF GAFF.ff spcew.ff";
    ($solvOpts, $ffields) = getSolvOpts($solvTypeStr, $solFF) if(defined($solvTypeStr));
    $sfile->{save} = $solute->getFileName($sfile->{name}) if (! defined($sfile->{save}));
    $LIBS = &AmberLib;
    system("rm -fr gpcrlammpsmd.log");

    print "Done\n";
}

sub getCellOpts {
    my ($memType, $memStr) = @_;
    my ($CELL);
    $CELL->{type} = $memType;
    if ($memStr =~ /^\s*(\d+\.?\d*)\s+(\d+\.?\d*)/) {
	($CELL->{x}, $CELL->{y}) = ($1, $2);
    } else {
	print "invalid membrane size \"$memStr\" encountered.. defaulting to 97 x 97...";
	($CELL->{x}, $CELL->{y}) = (97, 97);
    }
    $CELL->{size} = "$CELL->{x}x$CELL->{y}";
    return $CELL;
}

sub getSolvOpts {
    my ($solventStr, $solFF) = @_;
    my ($SOLVENT, $ffields);
    
    $solventStr = lc $solventStr;
    $ffields = "$solFF GAFF.ff ";
    
    if ($solventStr =~ /tip3/) {
	$SOLVENT->{OW}{CHARGE} = "+0.0068";
	$SOLVENT->{HW}{CHARGE} = "-0.0034";
	if ($solventStr =~ /charmm/) {
	    $ffields .= "tip3_charmm.ff";
	} else {
	    $ffields .= "tip3ew.ff";
	}
    } elsif ($solventStr =~ /tip4/) {
	$SOLVENT->{OW}{CHARGE} = "-0.1394";
	$SOLVENT->{HW}{CHARGE} = "+0.0.697";
	$ffields .= "tip42005.ff";
    } elsif ($solventStr =~ /f3c/) {
	$SOLVENT->{OW}{CHARGE} = "-0.014";
	$SOLVENT->{HW}{CHARGE} = "+0.007";
	$SOLVENT->{OW}{FFTYPE} = "O_3F";
	$SOLVENT->{HW}{FFTYPE} = "H_F";
	$ffields .= "f3cew.ff";
    }
    return ($SOLVENT, $ffields);
}

sub getAtomRange {
    my ($list) = $_[0];
    my (@atomids, $i, $prev, $start, $range);

    @atomids = sort {$a<=>$b} keys %{ $list };
    $i = 1;
    $prev = $start = $atomids[0];;
    while ($i <= $#atomids) {
	if (($atomids[$i] - $prev) > 1) {
	    $range .= ":Ia${start}-${prev} ";
	    $start = $atomids[$i];
	}
	$prev = $atomids[$i];
	$i++;
    }
    $range .= ":Ia${start}-${prev} ";
    return ($range, $atomids[0]);
}

sub getNewAtomRange {
    my ($list) = $_[0];
    my (@atomids, $i, $prev, $start, $range);

    @atomids = sort {$a<=>$b} keys %{ $list };
    $i = 1;
    $prev = $start = $atomids[0];
    while ($i <= $#atomids) {
        if (($atomids[$i] - $prev) > 1) {
            $range .= "(index>=${start} and index<=${prev}) ";
            $start = $atomids[$i];
        }
        $prev = $atomids[$i];
        $i++;
    }
    $range .= "(index>=${start} and index<=${prev}) ";
    return ($range, $atomids[0]);
}

sub getCellInfo {
    my ($options, $memFile) = @_;
    my ($crystxCmd, $CELL);

    #$crystxCmd = "grep crystx $memFile";
    #&execCmd($crystxCmd, "_crystx.dat");
    open CRYSTXCMD, "grep -i crystx $memFile |" or die "ERROR while executing $crystxCmd |: $!\n";
    while (<CRYSTXCMD>) {
	chomp;
	if ($_ = ~/^CRYSTX\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)/) {
	    ($CELL->{x}, $CELL->{y}, $CELL->{z}) = ($1, $2, $3);
	}
    }
    close CRYSTXCMD;
    die "ERROR: No valid cell info found while search $memFile!\n" if (! defined($CELL));
    return $CELL
}

sub execCmd {
    my ($cmdStr, $output) = @_;
    my ($terminal);

    $terminal = 1;
    $terminal = 0 if (defined($output) && ! $output);
    if (! defined($output)) {
	$cmdStr .= ">> gpcrlammpsmd.log";
	$output = "gpcrlammpsmd.log";
    } elsif (defined($output) && ! $output) {
	$cmdStr .= ">& gpcrlammpsmd.log";
	$output = "gpcrlammpsmd.log";
    } else {
	$cmdStr .= ">& $output";
    }
    system("echo \"\ncommand:'$_[0]'\" >> $output");
    if (system($cmdStr)) {
        die "ERROR: Check grcrlammpsmd.log\n" if ($terminal);
    }
}

sub showUsage {
    my $usage = <<DATA;
usage: $0 -i input structure -m (membrane type) -d (membrane xy size) 
	              				-s (solvent type) -f (forcefield) -w (savename)
OPTIONS:
	-i input structure:  REQUIRED. The following file formats are supported: BGF, PDB, MSI, MOL2
	-m membrane type:    OPTIONAL. Either POPC or POPE. Defaults to POPC.
	-d membrane xy size: OPTIONAL. The size of the membrane to construct. Needs to be at least as large as the protein
	                     Will default to 90 x 90.
	-s solvent type:     OPTIONAL. The following predetermined (equilibrated) solvents are available. Default F3C.
			     TIP3: the original Jorgenson TIP3 water model (rigid hydrogens)
			     TIP4: TIP4P with massless pseudo-atom
			     TIP3_CHARMM: TIP3 water model as implemented in CHARMM
			     F3C: F3C water model (no rigid hydrogens)
			     SPC: SPC water model
	-f forcefield:       OPTIONAL. This programs defaults to assigning AMBER charges/atom types for your protein/ligand complex.
			     This is desirable in most instances, however can be turned off using this flag to specify the forcefield(s)
			     commensurate with your strucuture. You should check that your current forcefield is compatible with AMBER
			     which is used for the lipid
	-w savename:         OPTIONAL. Will assume \$prefix_mod.\$suffix if not specified.
DATA

    return $usage;
}
