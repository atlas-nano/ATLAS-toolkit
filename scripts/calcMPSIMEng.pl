#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use General qw(ReadFFs STDev FileTester CoM GetSelections TrjSelections ShowSelectionInfo execCmd);
use AMBER qw(ParseAmberTrj GetAmberByteOffset ConvertAmberBox);
use LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
use ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList GetSolvent SplitAtomsByMol GetAtmData);
use BOX qw(GetBox);

sub init;
sub showUsage;
sub boxConvert;
sub calcMPSimEng;
sub numerically { ($a<=>$b); }
sub writeEng;
sub getGroupData;
sub execMPSim;
sub parseMPSimFile;
sub readGroupFile;

my ($SOLUTE, $trjFile, $BGF, $BONDS, $printStr, $reImage, $ff);
my ($getSnapshot, $getByteOffset, $trjType, $LAMMPSOPTS, $GROUPS);
my ($SELECT, $saveFile, $field, $DATA, $molList, $prefix);
my ($HEADERS, $BBOX, $GRPENG);

$|++;
&init;

if (! $trjType) {
    print "Computing MPSIM singlepoint...";
    &calcMPSimEng($BGF, $BBOX, 1, undef);
    print "Done\n";
} else {
    $field = scalar keys %{ $BGF };
    $getByteOffset->($SELECT, $trjFile, $field);
    if ($trjType == 2) {
        &GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS });
        $field = "atom";
    }
    $printStr = "Computing MPSIM singlepoint from $trjFile...";
    $getSnapshot->($BGF, $trjFile, $SELECT, $field, \&calcMPSimEng, $printStr, undef);
}
print "Writing stats to $saveFile...";
&writeEng($DATA, $saveFile, 0);
print "Done\n";
if (defined($GROUPS)) {
    $saveFile =~ s/\.\w+$//;
    $saveFile .= "_grps.dat";
    print "Writing grp stats to $saveFile...";
    &writeEng($GRPENG, $saveFile, 1);
    print "Done\n";
}

sub writeEng {
   my ($data, $savename, $writeStdev) = @_;
   my ($i, $field, $STATS);

   open OUTDATA, "> $savename" or die "ERROR: Cannot create $savename: $!\n";
   printf OUTDATA "%-8s%12s%12s%12s%12s%12s%12s\n",qw(Atom VDW Coul HB HB_dre Bond Total) if (! $writeStdev);
   printf OUTDATA "%-8s%20s%20s%20s%20s%20s%20s\n",qw(Atom VDW Coul HB HB_dre Bond Total) if ($writeStdev);
   for $i (sort numerically keys %{ $data->{ATOMS} }) {
	for $field (qw(vdw coul hb hb_dre bond total)) {
	    chop $data->{ATOMS}{$i}{$field};
	    ($STATS->{$field}{avg}, $STATS->{$field}{stdev}, undef) = STDev($data->{ATOMS}{$i}{$field});
	}
        printf OUTDATA "%-8d", $i;
	for $field (qw(vdw coul hb hb_dre bond total)) {
	    printf OUTDATA "%12.3f", $STATS->{$field}{avg};
	    printf OUTDATA "%8.2f", $STATS->{$field}{stdev} if ($writeStdev);
	}
	printf OUTDATA "\n";
    }
    close OUTDATA;
}

sub calcMPSimEng {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($tot, $i, $j, @tmp, $MOLECULE, $CENTER, $engFile, $HEADERS);

    if ($trjType == 2) { #LAMMPS
        $frameNum = $ATOMS->{TIMESTEP}[0];
	$tot = $ATOMS->{"NUMBER OF ATOMS"}[0];
        $BOX = ConvertLammpsBox(\%{ $ATOMS->{"BOX BOUNDS"} });
        $ATOMS = $ATOMS->{ATOMS};
	if ($LAMMPSOPTS->{scaled} or $LAMMPSOPTS->{imaged}) {
	    UnwrapAtoms(\%{ $ATOMS }, \%{ $BBOX }, $LAMMPSOPTS->{scaled});
	}
	@tmp = grep {!/COORD|VISITED/i} keys %{ $BGF->{1} };
	for $i (keys %{ $BGF }) {
	    for $j (@tmp) {
		$ATOMS->{$i}{$j} = $BGF->{$i}{$j};
	    }
	    $ATOMS->{$i}{IS_SOLVENT} = 1 if (exists($BGF->{$i}{IS_SOLVENT}));
	}
    } elsif ($trjType == 1) {
        $BOX = ConvertAmberBox(\%{ $BOX });
	$tot = scalar(keys %{ $ATOMS });
    } 
    if ($reImage) {
        $MOLECULE = GetAtmData($ATOMS, $SOLUTE);
        $CENTER = CoM($MOLECULE);
        @tmp = ("XCOORD", "YCOORD", "ZCOORD");
        for $j (@tmp) {
            $BOX->{$j}{CENTER} = $BOX->{$j}{len}/2;
            $BOX->{$j}{hi} = $BOX->{$j}{len};
            $BOX->{$j}{lo} = 0;
            for $i (keys %{ $ATOMS }) {
                $ATOMS->{$i}{$j} += ($BOX->{$j}{CENTER} - $CENTER->{$j});
            }
        }
        for $i (keys %{ $molList }) {
            $MOLECULE = GetAtmData($ATOMS, $molList->{$i});
            $CENTER = CoM($MOLECULE);
            ImageAtoms($MOLECULE, $CENTER, $BOX);
        }
    }

    # now create a bgf file for this snapshot
    $HEADERS = createHeaders($BOX, "${prefix}_${frameNum}.bgf");
    &addHeader($ATOMS, $HEADERS);
    &createBGF($ATOMS, $BONDS, "${prefix}_${frameNum}.bgf");

    $engFile = execMPSim($prefix, $ff, "${prefix}_${frameNum}.bgf", "${prefix}_${frameNum}"); 
    &parseMPSimFile($engFile, \%{ $DATA }, $frameNum);
}

sub parseMPSimFile {
    my ($prefix, $data, $frame) = @_;
    my ($engFile, @vals, $field, $isvalid, $atomId); 
    my ($i, $snap, $grpID, $j, $grpEng, $tot);

    $engFile = "${prefix}.fin.ener";
    $isvalid = $tot = 0;
    if (open ENGFILE, $engFile) {
    while (<ENGFILE>) {
	chomp;
	if ($_ =~ /^\d+\s*\-?\d+\.\d+/) {
	    @vals = split /\s+/, $_;
	    $atomId = shift @vals;
	    for $field (qw(vdw coul hb hb_dre bond total)) {
		$snap->{$atomId}{$field} = shift(@vals);
	    }
	    $tot += $snap->{$atomId}{total};
	    $isvalid = 1;
	}
    }
    close ENGFILE;
    }
    if (! $isvalid) {
	print "ERROR: Invalid MPSIM energy file $engFile!\n";
	&execCmd("rm -fr ${prefix}.trj ${prefix}.traj1 ${prefix}.fin.ener ${prefix}.fin.bgf");
        return;
    }
    $data->{frames} += 1;
    &execCmd("rm -fr ${prefix}.bgf ${prefix}.trj ${prefix}.traj1 ${prefix}.fin.ener ${prefix}.fin.bgf");
    for $i (keys %{ $snap }) {
	for $field (qw(vdw coul hb hb_dre bond total)) {
	    $data->{ATOMS}{$i}{$field} .= $snap->{$i}{$field} . " ";
	}
	next if (! defined($GROUPS));
	undef($grpID);
	for $j (0 .. $#{ $GROUPS }) {
	    $grpID = $j if (exists($GROUPS->[$j]{$i}));
	}
	next if (! defined($grpID));
	for $field (qw(vdw coul hb hb_dre bond total)) {
	    $grpEng->{$grpID}{$field} += $snap->{$i}{$field};
	}
    }
    return if (! defined($GROUPS));
    for $i (keys %{ $grpEng }) {
	for $field (qw(vdw coul hb hb_dre bond total)) {
	    $GRPENG->{ATOMS}{$i}{$field} .= $grpEng->{$i}{$field} . " ";
	}
    }
    $GRPENG->{frames} += 1;
}
		    
sub execMPSim {
    my ($name, $ff, $bgfFile, $projName) = @_;
    my ($outstring, $outCmd);
    my ($mpsim_cmd) = "/ul/maiti/src/non_pme/build/linux/mpsim.20030918 ";

    $outstring = "";
    &execCmd("cp /ul/tpascal/mpsim/ewald.par .");
    open CTL, "/ul/tpascal/mpsim/1-pme.ctl" || die "Cannot open 1-pme.ctl: $!\n";
    while (<CTL>) {
        chomp;
        if ($_ =~ /^PROJECT/) {
            $outstring .= "PROJECT               $projName";
        } elsif ($_ =~ /^STRUCTURE/) {
            $outstring .= "STRUCTURE             $bgfFile";
	} elsif ($_ =~ /^FF/) {
	    $outstring .= "FF                    $ff";
        } else {
            $outstring .= $_;
        }
        $outstring .= "\n";
    }
    close CTL;
    open CTL, "> ${name}_pme.ctl" || die "Cannot modify control file: $!\n";
    print CTL $outstring;
    close CTL;
    $outCmd = "$mpsim_cmd ${name}_pme.ctl >& ${projName}_mpsim.out";
    if (system($outCmd)) {
        print "Mpsim Error!\n";
        return 0;
    } else {
        #&execCmd("rm -f ${name}_pme.ctl ewald.par ${projName}_mpsim.out");
    }
    return $projName;
}

sub init {
    my (%OPTS, $usage);
    my ($ffStr, $groupStr, $trjStr, $trjSelStr, $grpFile, $soluteStr); 
    my ($BOX, $trjTypeStr, $bgfFile, $list, $j, @tmp, $atmBondList);

    getopt('bfgtsnlrai',\%OPTS);
    $usage = &showUsage;
    die "$usage" if (! exists($OPTS{b}));

    print "Initializing...";
    ($bgfFile, $ff, $groupStr, $trjFile, $trjSelStr, $trjTypeStr, $reImage, $saveFile, $grpFile, $soluteStr) = 
	($OPTS{b}, $OPTS{f}, $OPTS{g}, $OPTS{t}, $OPTS{l}, $OPTS{n}, $OPTS{r}, $OPTS{s}, $OPTS{a}, $OPTS{i});
 
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_eng.dat";
    }
    $reImage = 0 if (! defined($reImage) or $reImage !~ /^\s*(1|y|yes)\s*$/);
    $reImage = 1 if ($reImage =~ /^\s*(1|y|yes)\s*$/);

    if (defined($ff)) {
	if (! -e $ff or ! -r $ff or ! -T $ff) {
	    print "invalid force field \"$ff\": $!. Ignoring...";
	    $ff = "/ul/tpascal/mpsim/AMBER03.ff";
	}
    } else { $ff = "/ul/tpascal/mpsim/AMBER03.ff" };

    if (defined($trjFile) and -e $trjFile and -r $trjFile and -T $trjFile) {
        if (! defined($trjTypeStr)) {
            if ($trjFile =~ /\.lammpstrj/) {
                $trjTypeStr = "lammps";
            } else {
                $trjTypeStr = "amber";
            }
        }

        if (lc($trjTypeStr) ne "lammps") {
            $trjType = 1;
            $getSnapshot = \&ParseAmberTrj;
            $getByteOffset = \&GetAmberByteOffset;
        } else {
            $trjType = 2;
            $getSnapshot = \&ParseLAMMPSTrj;
            $getByteOffset = \&GetLammpsByteOffset;
        }
        if (! defined($trjSelStr)) {
            $trjSelStr = "*";
        }
        $list = TrjSelections($trjSelStr);
        for $j (keys %{ $list }) {
            $SELECT->{$j} = $list->{$j};
        }

        die "ERROR: No valid frames selected with selection $trjSelStr!\n"
            if (! keys %{ $SELECT } and $trjSelStr ne "*");
    } else {
        $trjType = 0;
    }
    print "Done\nParsing BGF file $bgfFile...";
    ($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
    $GROUPS = getGroupData($groupStr, $BGF) if (defined($groupStr));
    $prefix = basename($saveFile);
    $prefix =~ s/\.\w+$//;

    if (defined($soluteStr)) {
        @tmp = split /\s+/, $soluteStr;
        $atmBondList = GetSelections(\@tmp, 0);
	$SOLUTE = GetAtmList($atmBondList, $BGF);
    } else {
	$SOLUTE = \%{ $BGF->{1}{MOLECULE}{MEMBERS} };
    }
    die "ERROR: Solute atoms does not correspond to any in BGF file!\n"
        if (! keys %{ $SOLUTE });
    $BOX = GetBox($BGF, undef, $HEADERS);
    %{ $BBOX } = %{ $BOX };
    &boxConvert($BBOX);
    for (keys %{ $BGF }) {
	if (! exists($SOLUTE->{$_})) {
	    $molList->{$_} = 1;
	}
    }
    $molList =  SplitAtomsByMol($BGF, $molList);

    print "Done\n";
    $GROUPS = readGroupFile($grpFile) if (! defined($GROUPS) and defined($grpFile));
}

sub getGroupData {
    my ($groupStr, $atoms) = @_;
    my (@list, $GROUPS, $i, $selection, $rec);
    
    if ($groupStr =~ /\(.+\)/) {
	@list = split /\)/, $groupStr;
    } else {
	$list[0] = $groupStr;
    }


    for $i (@list) {
	$i =~ s/^\s*\(//;
	if ($i =~ /\s+/) {
            @{ $selection } = split /\s+/, $i;
	} else {
            $selection->[0] = $i;
	}
	$rec = GetSelections($selection, 0);
	$rec = GetAtmList($rec, $atoms);
	push @{ $GROUPS }, $rec;
    }

    return $GROUPS;
}

sub readGroupFile {
    my ($grpFile) = $_[0];
    my ($DATA, @dat, $i, $j, $grpNum);

    print "Parsing Grp file $grpFile...";
    open GRPFILE, $grpFile or die "ERROR: Cannot read from $grpFile: $!\n";
    while (<GRPFILE>) {
	chomp;
	last if ($_ =~ /Constraints/);
	if ($_ =~ /^Group (\d+)/) {
	    $grpNum = $1;
	} elsif ($grpNum and $_ =~ /^\d+/) {
	    @dat = split /,/, $_;
	    for $i (@dat) {
		if ($i =~ /(\d+)\-(\d+)/) {
		    for $j ($1 .. $2) {
			$DATA->[$grpNum]{$j} = $grpNum;
		    }
		} elsif ($j =~ /(\d+)/) {
		    $DATA->[$grpNum]{$j} = $grpNum;
                }
            }
        }
    }
    close GRPFILE;
    if (! $DATA) {
	print "ERROR!\n";
	return undef;
    }

    print "Done\n";

    shift @{ $DATA };

    return $DATA;
}

sub boxConvert {
    my ($box) = $_[0];

    for ("X", "Y", "Z") {
        %{ $box->{$_ . "COORD"} } = %{ $box->{$_} };
    }
}

sub showUsage {
    my ($usage) = "usage: $0 -b bgf file -f (MPSIM forcefield) -g (groups) -t (trajectory)" . 
		   " -l (trj frames) -n (trj type) -s (savename) -r (reimage)\n" .
	"options:\n\t-b bgf file: location of bgf structure file. This is the only required options." .
	"\n\t-f MPSIM forcefield: location of MPSIM forcefield to use. Defaults to /ul/tpascal/mpsim/AMBER03.ff" .
	"\n\t-g groups: group selection (see below for format). Enclose different group in brackets. eg \"(Ir 2) (Ir 3)\"" .
	"\n\t-t trajectory: LAMMPS/AMBER trajectory file. If specified, will calc single points for specified snapshots" .
	"\n\t-l trj frames: Trajectory selection. Can be single number, range (:Ita-b:c) or any combination of the two" .
	"\n\t-n trj type: Trajectory type. Either AMBER or LAMMPS. If not specified, will attempt to guess" .
	"\n\t-s savename: Name of file to save. Will default to {bgf}_eng.dat" .
	"\n\t-r reimage: Specifies whether to wrap atoms back into the box. Default no\n" .
	"\n\t-a group file: 2pt formatted group file. The -g option overwrites this\n" .
	&ShowSelectionInfo;
    return $usage;
}
