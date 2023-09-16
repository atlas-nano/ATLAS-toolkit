#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use MolData;
use Getopt::Std qw(getopt);
use General qw(GetSelections ShowSelectionInfo);

sub numerically { $a<=> $b; }

my ($fileName, $saveName, $fileType, $mutName);
my ($molStructure, $RES, $mRes, $LIBS);

$|++;
&init;

print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
print "Done\nSelecting atoms/residues...";
&parseAtomSelection;
print "Done\nMutating $mRes->{count} residue(s) to " . uc $mutName . "...";
&mutateRes($molStructure, $mRes, $mutName, $LIBS);
print "Done\nCreating $fileType file $saveName...";
$molStructure->write($saveName, $fileType);
system("rm -fr _tmpres.pdb leaprc mol.prmtop mol.rst7 leap.log _tmpres.bgf _out");
print "Done\n";

sub mutateRes {
# the idea here is to take the selected resid and mutate it to the new residue
# by: 
# 1. removing all but the backbone atoms (C O N CA CB) 
# 2. renaming the residue to the new name
# 3. writing out the new residue (plus 1 residue before and after) as a pdb file
# 4. using xleap to add back the hydrogens/side chain for the fragment
# 5. inserting the new residue back into the structure
    my ($struct, $resIDs, $resName, $amberLibs) = @_;
    my ($resAtoms, $count, $resStruct, $bondOpts); 
    my ($atom, $resList, $res, $newRes, $resNum, $atomC);

    for $resNum (sort numerically keys %{ $resIDs->{"resid"} }) {
	$count = 0;
	$resList = buildResList($struct, $amberLibs, $resNum);
	$resList->{res}{atoms} = $resIDs->{"resid"}{$resNum};
	$resList->{res}{resid} = $resNum;
	$resStruct =  Packages::MolData->new();
	for $res ("prev", "res", "next") {
	    next if (! $resList->{$res});
	    $resAtoms = $resList->{$res}{atoms};
	    for $atomC (keys %{ $resAtoms }) {
		$atom = $struct->cloneAtom($resAtoms->{$atomC});
		$atom->("resname", uc $resName) if ($res eq "res"); # update the resname if targed
		if ($res =~ /prev|next/) {
		    $atom->("resname", "HIS") if ($atom->resname =~ /HI/);
		    $atom->("resname", "CYS") if ($atom->resname =~ /CYX/);
		}
		if ($atom->atmname =~ /^\s*(C|N|CA|O)\s*$/) { #only write backbone atoms to pdb
		    $resStruct->insertAtom($atom);
		    $count++;
		}
		$atom = undef;
	    }
	    undef($resAtoms);
	}
	next if (! $count);
	$resStruct->count->{"atoms"} = $count;
	$resStruct->write("_tmpres.pdb", "pdb");
	&runLeap("_tmpres.pdb");
	$newRes = extractRes("_tmpres.bgf", $resList);
	$struct->deleteRes($resNum, \%{ $bondOpts });
	$struct->insertRes($resNum, $newRes, $bondOpts);
	#$struct->updateAtomIndex;

        #now do some garbage cleanup. Probably not necessary
	$resList = undef;
	undef($resStruct);
	undef($atom);
	undef($newRes);
	undef($bondOpts);
    }
}

sub extractRes {
    my ($bgfFile, $resList) = @_;
    my ($resStruct, $res, $resIDAtoms, $resID);

    $resStruct =  Packages::MolData->new();
    $resStruct->testFile($bgfFile);
    $resStruct->read($bgfFile, "bgf");
    $resID = 2;
    $resID = 1 if (! exists($resList->{prev}));
    $resIDAtoms = $resStruct->shash->{resid}{ $resID };
    $res = $resStruct->extractAtoms($resIDAtoms);
    return $res;
}

sub buildResList {
   my ($struct, $amberLibs, $resNum) = @_;
   my ($resList, $i, @atomList, $sHash, $resname, $atomId);

   $sHash = { ($resNum - 1) => "prev", ($resNum + 1) => "next" };
   for $i (($resNum - 1), ($resNum + 1)) {
	next if (! exists($struct->shash->{"resid"}{$i}));
	@atomList = keys %{ $struct->shash->{"resid"}{$i} };
	next if (! @atomList);
	$atomId = shift @atomList;
	next if (! $atomId or ! $struct->shash->{"resid"}{$i}{$atomId});
	$resname = $struct->shash->{"resid"}{$i}{$atomId}->resname;
	$resname = "HIE" if ($resname =~ /his|hse|hdd/i);
	$resname = "CYS" if ($resname =~ /cyx/i);
	next if (! exists($amberLibs->{uc $resname}));
	$resList->{ $sHash->{$i} }{atoms} = $struct->shash->{"resid"}{$i};
	$resList->{ $sHash->{$i} }{resid} = $i;
    }

    return $resList;
}

sub getSelectedResID {
    my ($struct, $selList, $libList, $resList, $mutateName) = @_;
    my ($atomList, $i, $j, $operator, $val, $foundAtoms, $k);
    my ($allResNames, $resname, $foundRes, @tmp, $resid);

    $atomList = $struct;
    # first find all the atoms corresponding to the selection
    for $i (keys %{ $selList }) {
	for $j (keys %{ $selList->{$i} }) {
	    if ($j =~ /^(>|<|^|=)(.+)/) {
		($operator, $val) = ($1, $2);
		$operator = "!=" if ($operator eq "^");
		if ($val =~ /[A-Za-z]/) {
		    $operator = "lt" if ($operator eq "<");
		    $operator = "gt" if ($operator eq ">");
		    $operator = "ne" if ($operator eq "^");
		    $operator = "eq" if ($operator eq "=");
		}
	    } else {
		$val = $j;
		$operator = "==";
		$operator = "eq" if ($val !~ /\-?\d+\.?\d*e?\-?\d*/);
	    }
	    $i = "RESID" if ($i eq "RESNUM");
	    $foundAtoms = $atomList->find($i, $val, $operator);
	    undef($atomList);
	    $atomList = Packages::MolData->new();
	    for $k (keys %{ $foundAtoms }) {
		$atomList->insertAtom($foundAtoms->{$k});
	    }
	}
    }
    undef($atomList);

    #now find all the residues associated with these atoms
    %{ $allResNames } = %{ $struct->shash->{"resname"} };
    for $i (keys %{ $foundAtoms }) {
	$resname = $foundAtoms->{$i}->resname;
	$resid = $foundAtoms->{$i}->resid;
	$foundRes->{$resid} = $resname;
    }
    for $resid (keys %{ $foundRes }) {
	$resname =  uc $foundRes->{$resid};
	delete $foundRes->{$resid} if (! exists($libList->{$resname}) or ($resname eq $mutateName));
	next if (exists($resList->{resid}{$resid}));
	%{ $resList->{resid}{$resid} } = %{ $struct->shash->{"resid"}{$resid} };
	$resList->{count}++;
    }

    die "ERROR: No atoms matched the search parameter!\n" if (! $resList->{count});
}

sub init {
    my (%OPTS, $resSel);

    getopt('fsmrt', \%OPTS);
    for ("f", "m") {
	&showusage() if (! exists($OPTS{$_})); 
    }

    print "Initialzing...";
    $molStructure =  Packages::MolData->new();
   ($fileName, $fileType, $saveName, $resSel, $mutName) = ($OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{r}, $OPTS{m});

    $resSel = "*" if (! defined($resSel));
    $molStructure->testFile($fileName);
    if (! defined($fileType)) {
	$fileType = "bgf";
	if ($fileName =~ /\.(\w+)$/) {
	    $fileType = lc $1;
	}
    }
    $saveName = $molStructure->getFileName($fileName) if (! defined($saveName));
    $RES = getSelOpts($resSel);
    $LIBS = parseAmberLib();
    die "RESIDUE \"$mutName\" is not a valid AMBER residue!\n" if (! exists($LIBS->{uc $mutName }));
    print "Done\n";
}

sub getSelOpts {
    my ($selection) = $_[0];
    my ($SELECT, $group, $count, $tmpStr);

    $count = 1;
    $selection =~ s/^\s*//;
    $selection =~ s/\s*$//;

    if ($selection !~ /.*\(.+\)/) {
	while ($selection =~ /(\S+)/g) {
	    $tmpStr .= "($1) ";
	}
	$selection = $tmpStr;
    }
    while ($selection =~/\((\S+)\)/g) { #everything that is group in brackets represents a and group
	@{ $group } = split /\s/, $1;
	$SELECT->{$count} = GetSelections($group,0);
	$count++;
    }

    return $SELECT;
}

sub parseAmberLib {
    my (%RESNAMES);
    my ($libfile) = "/exec/amber9/dat/leap/lib/all_amino02.lib";
    
    open LIB, $libfile or die "ERROR: Cannot open $libfile: $!\n";
    while (<LIB>) {
	chomp;
	if ($_ =~ /^\s+\"(\w+)\"$/) {
	    $RESNAMES{$1} = 1;
	}
    }
    close LIB;
    die "No valid data found in AMBER lib file $libfile!\n" if (! %RESNAMES);
    return \%RESNAMES;
}

sub runLeap {
    my ($pdbfile) = $_[0];
    my ($leapcmd, $leaprc, $cmdStr, $createBGFcmd);

    $leapcmd = "/ul/tpascal/programs/ambertools/exe/tleap";
    $createBGFcmd = "/ul/tpascal/scripts/amber2bgf.pl mol.prmtop mol.rst7 _tmpres.bgf >> _out";

    die "ERROR: Cannot copy /ul/tpascal/amber/leaprc to current directory!\n"
        if(system("cp /ul/tpascal/amber/leaprc ."));

    $cmdStr = <<DATA;
mol = loadpdb $pdbfile
saveamberparm mol mol.prmtop mol.rst7
quit
DATA

    die "ERROR: Cannot write to ./leaprc!\n"
        if (system("echo \"$cmdStr\" >> ./leaprc"));
    die "ERROR: Cannot execute tleap!\n"
        if (system("$leapcmd > _out"));
    die "ERROR: Cannot create bgf file _tmp.bgf!\n"
        if(system("$createBGFcmd"));
}

sub parseAtomSelection {
    $mRes = { resid => {}, count => 0 };
    if (exists($RES->{ALL})) {
	$mRes = {resid => $molStructure->shash->{"resid"}, count => "ALL" };
    } else {
	&getSelectedResID($molStructure, $RES->{$_}, $LIBS, $mRes, $mutName) for (keys %{ $RES });
    }
}

sub showusage {
    my ($usage) = "usage: $0 -f filename -m mutation name -r (atom selection) -t (filetype) -s (save name)\n" .
	"options:\n\t-f filename: location of file\n\t-m mutation name: the name of the amber residue to mutate too\n" .
	"\t-r atom selection: the list of atoms (residues) to mutate (see below for syntax). Optional. Default: all\n" .
	"\t-r filetype: bgf|mol2|pdb|msi. Optional. Will be calculated if omitted\n\t-s savename: name of output file\n" .
	&ShowSelectionInfo;
    die "$usage";
}
