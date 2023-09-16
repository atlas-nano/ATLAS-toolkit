#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Packages::General qw(FileTester Trim execCmd PrintProgress STDev TrjSelections GetSoluteAtoms CoP CenterOnMol GetSelections);
use Packages::GetParms;
use Packages::HelixLayout;
use Packages::FileFormats qw(GetBGFFileInfo createBGF addHeader addBoxToHeader createHeaders GetBondList);
use Packages::AMBER qw(ParseAmberTrj GetAmberByteOffset);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset GetLammpsTrjType);
use Packages::ManipAtoms qw(UnwrapAtoms ImageAtoms GetAtmList);
use Cwd 'abs_path';
use Cwd;
use Getopt::Std qw(getopt);

sub init;
sub determineHelixLayout;
sub setLayout;
sub numerically;
sub showHelp;
sub createMPSimFile;
sub execCmd;
sub getEnergies;
sub writeOutput;
sub calculateTotals;
sub createGraphs;
sub makeBox;
sub getAtoms;
sub convertBox;

my ($bgfFile, $parmFile, $filebase, $atomSelection, $getSnapshot, $reImageCoords);
my ($SELECT, $BGF, $BONDS, $Parm, $LAYOUT, $currDir, $HEADERS, @time_array, $RES); 
my ($printStr, %ENERGIES, $HC, $trjFile, @eng_header, @unit_array, $ATOMLIST);
my ($molList, $SOLUTEATMS, $getByteOffset, $LAMMPSOPTS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($BGF, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$SOLUTEATMS = GetAtmList($RES, $BGF);
$molList = GetBondList($BGF, $BONDS);
print "Done\n";
$getByteOffset->($SELECT, $trjFile, scalar keys %{ $BGF });
&GetLammpsTrjType($SELECT, $trjFile, "atom", \%{ $LAMMPSOPTS }) 
    if ($Parm->{Files}{trjType} == 2);
chdir $filebase;
open DUMPFILE, "> dumpfile.txt" or die "Cannot create dumpfile.txt: $!\n";
doMPSimAnal($Parm, $BGF, $SOLUTEATMS, $SELECT, $trjFile);
close DUMPFILE;
print "Calculating Total Energy...";
&writeOutput;
chdir $currDir;

sub doMPSimAnal {
    my ($PARMS, $bgfinfo, $atmSelection, $trjSelection, $tFile) = @_;
    my ($i, $engFile, $currBGF, $tmpSelect, $field, $tot);
    my ($strLen, $printStr, $isValid, $start, $index);

    $field = scalar keys %{ $bgfinfo };
    $field = "atom" if ($PARMS->{Files}{trjType} == 2);
    $tot = scalar keys %{ $trjSelection };
    $start = time();
    $index = 1;
    $isValid = 0;
    $printStr = "Analysing files...";
    print "${printStr}...calculating time remaining\r";
    $strLen = 30;

    for $i (sort numerically keys %{ $trjSelection }) {
	next if (! defined($trjSelection->{$i}));
	$engFile = "${currDir}/outfiles/" . $PARMS->{Molecule}{name} . "_${i}.out";
	if (! -e $engFile || ! -r $engFile || ! -T $engFile) { # then create the mpsim energy file
	    $currBGF = "${currDir}/bgffiles/" . $PARMS->{Molecule}{name} . "_${i}.bgf";
	    if (! -e $currBGF || ! -r $currBGF || ! -T $currBGF) { # then create the bgf file
		$tmpSelect = ();
		$tmpSelect->{$i} =  $trjSelection->{$i};
		$getSnapshot->($bgfinfo, $tFile, $tmpSelect, $field, \&saveBGF, undef, $currBGF);
	    }
	    if ( ! createMPSimFile($PARMS, $currBGF, $PARMS->{Molecule}{name} . "_$i", $i)) {
		$index++;
		next;
	    }
	}
	if (getEnergies($engFile, $atmSelection, $i)) {
	    push @time_array, $i;
	    $isValid = 1;
	}
	$strLen = PrintProgress($index, $tot, $start, $printStr);
	$index++;
    }
    printf "$printStr%-${strLen}s\n", "Done";
    die "ERROR: No valid files found\n" if (! $isValid);
    @unit_array = sort numerically keys %ENERGIES;
}

sub getAtoms {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);

    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = $allAtoms->{$i};
    }

    return \%ATOMS;
}

sub makeBox {
    my ($box) = $_[0];
    my (%BOX, $i, @dim);

    @dim = ("X", "Y", "Z");
    for $i (0 .. 2) {
	$BOX{$dim[$i]}{lo} = 0;
	$BOX{$dim[$i]}{hi} = $box->{$i + 2}{DATA};
	$BOX{$dim[$i] . "COORD"}{lo} = 0;
	$BOX{$dim[$i] . "COORD"}{hi} = $box->{$i + 2}{DATA};
	$BOX{$dim[$i] . "COORD"}{len} = $box->{$i + 2}{DATA};
    }

    return \%BOX;
}

sub convertBox {
    my ($lammpsBox) = $_[0];
    my (%amberBox, $i);

    $amberBox{1}{DATA} = 90;
    for $i (0 .. 2) {
	$amberBox{$i + 2}{DATA} = $lammpsBox->[$i]{hi} - $lammpsBox->[$i]{lo};
    }

    return \%amberBox;
}


sub saveBGF {
    my ($ATOMS, $BOX, $junk, $currBGF) = @_;
    my ($MOL, $CENTER, $i, $BBOX, $trjType, $headers, $j, @tmp, $tmp1);

    $trjType = 1;
    if (exists($ATOMS->{ATOMS})) { # then it's lammps
	$BOX = convertBox($ATOMS->{"BOX BOUNDS"});
	$BBOX = makeBox($BOX);
	$ATOMS = $ATOMS->{ATOMS};
	UnwrapAtoms($ATOMS, $BBOX, $LAMMPSOPTS->{scaled});
	$trjType = 2;
	$currBGF = $junk;
	@tmp = grep {!/COORD/i} keys %{ $BGF->{1} };
	for $i (keys %{ $ATOMS }) {
	    for $j (@tmp) {
		next if ($j =~ /COORD/i);
		$ATOMS->{$i}{$j} = $BGF->{$i}{$j};
	    }
	}
    } else {
	$BBOX = makeBox($BOX);
    }

    if ($reImageCoords) {
	$MOL = getAtoms($ATOMS, $SOLUTEATMS);
	$CENTER = CoP($MOL);
	CenterOnMol($ATOMS, $CENTER);

	for $i (keys %{ $molList }) {
	    $MOL = getAtoms($ATOMS, $molList->{$i});
	    $CENTER = CoP($MOL);
	    ImageAtoms($ATOMS, $molList->{$i}, $CENTER, $BBOX);
	}
    }
    $headers = createHeaders($BOX, $bgfFile);
    addHeader($ATOMS, $headers);
    createBGF($ATOMS, $BONDS, $currBGF);
}

sub doMPSimAnalOld {
    my ($PARMS, $bgfinfo, $headers, $trjSelection, $trjFileHandle, $bonds) = @_;
    my (@tmp, $frame, $currBGF, $DATA, $start, $tot, $index, $engFile);
    my ($lastFrame, $totAtms, $strLen, $printStr, $isValid);

    @tmp = sort numerically keys %{ $bonds };
    $totAtms = $tmp[$#tmp];

    @tmp = sort numerically keys %{ $trjSelection };
    $tot = $#tmp + 1;
    $start = time();
    $index = 1;
    $lastFrame = $isValid = 0;
    $printStr = "Analysing files...";
    print "${printStr}...calculating time remaining\r";
    for $frame (@tmp) {
	$engFile = "${currDir}/outfiles/" . $PARMS->{Molecule}{name} . "_${frame}.out";
	if (! -e $engFile || ! -r $engFile || ! -T $engFile) { # then create the mpsim energy file
	    $currBGF = "bgffiles/" . $PARMS->{Molecule}{name} . "_${frame}.bgf";
	    if (! -e $currBGF || ! -r $currBGF || ! -T $currBGF) { # then create the bgf file
		$DATA = getTrjFrame($trjFileHandle, $bgfinfo, $totAtms, ($frame - $lastFrame));
		if (! $DATA) {
		    $index++;
		    next;
		}
		updateBOX($DATA->{BOX}, $headers);
		addHeader($DATA->{ATOMS}, $headers);
		createBGF($DATA->{ATOMS}, $bonds, $currBGF);
		$lastFrame = $frame;
		$DATA = ();
	    }
	    if ( ! createMPSimFile($PARMS, $currBGF, $PARMS->{Molecule}{name} . "_$frame", $frame)) {
		$index++;
		next;
	    }
	}
	if (getEnergies($engFile, $PARMS->{Molecule}{total_bases}, $frame)) {
	    push @time_array, $frame;
	    $isValid = 1;
	}
	$strLen = PrintProgress($index, $tot, $start, $printStr);
	$index++;
    }
    printf "$printStr%-${strLen}s\n", "Done";
    die "ERROR: No valid files found\n" if (! $isValid);
    @unit_array = sort numerically keys %ENERGIES;

}

sub getEnergies {
    my ($curr_fle, $atomSelection, $hash_key) = @_;
    my ($data_start, $in_data, $unit_nm, $counter, $eng);
    my ($curr_val, $curr_total, $unit_label, $isValid);

    $isValid = 0;
    $data_start = 0;
    print DUMPFILE "Calculating energy from $curr_fle...";
    if (open INPFILE, $curr_fle) {
        while (<INPFILE>) {
            chomp;
            $in_data = $_;
            if ($in_data =~ /^\s*[Atom|GROUP]/ && $#eng_header <= 0) {
                while ($in_data =~ /\s+(\w+)/g) {
                    if ( ($1 !~ /GROUP/i) && ($1 !~ /TOTAL/i) ) {
                        push @eng_header, $1;
#                           print "GOT TOO GROUP HEADINGS: $1\n";
                    }
                }
            } elsif ($in_data =~ /^\s*(\d+)\s+/ && $#eng_header > 0) {
                $unit_label = $BGF->{$1}{ATMNAME};

#               print DUMPFILE "unit: $unit_nm, time: $hash_key\n";
                if (exists($atomSelection->{$1})) {
		    $unit_nm = $atomSelection->{$1};
                    $counter = 0;
                    $curr_val = $curr_total = 0;
                    while ($in_data =~ /(\-?\d+\.\d+)\s*/g && $counter <= $#eng_header) {
                        $eng = $eng_header[$counter];
                        $curr_val = $1;
                        $curr_total += $curr_val;
#                               print "counter: $counter eng: $eng\n";
                        if ($ENERGIES{$unit_nm}{$hash_key}{$eng}) {
                            $ENERGIES{$unit_nm}{$hash_key}{$eng} += $curr_val;
                        } else {
                            $ENERGIES{$unit_nm}{$hash_key}{$eng} = $curr_val;
                        }
                        $counter++;
                    }
                    if ($ENERGIES{$unit_nm}{$hash_key}{"TOTAL"}) {
                        $ENERGIES{$unit_nm}{$hash_key}{"TOTAL"} += $curr_total;
                    }else {
                        $ENERGIES{$unit_nm}{$hash_key}{"TOTAL"} = $curr_total;
                    }
                    $curr_val = $curr_total = 0;
                }
            }
        }

        if ($#eng_header > 0) {
            $isValid = 1;
            print DUMPFILE "Done\n";
        }else {
            print DUMPFILE "Failure. $curr_fle is invalid\n";
        }
        close INPFILE;
    }else {
        print DUMPFILE "\nWARNING: Invalid file: $curr_fle\n";
    }
    return $isValid;
}

sub createMPSimFile {
    my ($PARMS, $bgfFile, $engFile, $frame) = @_;
    my ($outstring, $curr_project, $outCmd);
    my ($mpsim_cmd) = "/ul/maiti/src/non_pme/build/linux/mpsim.20030918 ";

    $outstring = "";
    $curr_project = $PARMS->{Molecule}{name} . "_" . $frame;
    execCmd("cp /ul/tpascal/mpsim/ewald.par .");
    open CTL, "/ul/tpascal/mpsim/1-pme.ctl" || die "Cannot open 1-pme.ctl: $!\n";
    while (<CTL>) {
	chomp;
	if ($_ =~ /PROJECT/) {
	    $outstring .= "PROJECT               $engFile";
	} elsif ($_ =~ /STRUCTURE/) {
	    $outstring .= "STRUCTURE             $bgfFile";
	} else {
	    $outstring .= $_;
	}
	$outstring .= "\n";
    }
    close CTL;
    open CTL, "> 1-pme.ctl" || die "Cannot modify control file: $!\n";
    print CTL $outstring;
    close CTL;
    $outCmd = "$mpsim_cmd 1-pme.ctl >& mpsim.out";
    if (system($outCmd)) {
	print "Mpsim Error!\n";
	return 0;
    } else {
	execCmd("rm -f 1-pme.ctl ewald.par mpsim.out");
	execCmd("mv ${curr_project}.fin.ener ${currDir}/outfiles/${curr_project}.out", "${curr_project}.*");
    }

    return 1;
}

sub getTrjFrame {
    my ($TRJFILE, $ATOMS, $totAtms, $frame) = @_;
    my ($currFrame, $i, %DATA, $atomC, $offSet, $j, @dim, @tmp);

    @tmp = keys %{ $ATOMS->{1} };
    $atomC = $currFrame = 1;
    @dim = ("XCOORD","YCOORD","ZCOORD");
    $i = 0;
  MAINLOOP: while (<$TRJFILE>) {
      chomp;
      while ($_ =~ /(\-?\d+\.\d{3})/g) {
	  if ($atomC <= $totAtms) {
	      $DATA{ATOMS}{$atomC}{ $dim[$i] } = $1;
	      $i++;
	      if ($i == 3) {
		  $i = 0;
		  for $j (@tmp) {
		      next if ($j =~ /coord/i);
		      $DATA{ATOMS}{$atomC}{$j} = $ATOMS->{$atomC}{$j};
		  }
		  $atomC++;
	      }
	  } else {
	      $offSet = $atomC - $totAtms - 1;
	      if ($offSet < 2) { # store box information
		  $DATA{BOX}{ $dim[$offSet] } = $1;
		  $atomC++;
	      } else { # finished reading frame
		  $DATA{BOX}{ $dim[2] } = $1;
		  if ($currFrame == $frame) {
		      last MAINLOOP;
		  } else {
		      %DATA = ();
		      $currFrame++;
		      $atomC = 1;
		  }
	      }
	  }
      }
  }
    
    return \%DATA;
}

sub init {
    my ($i, $cmd, %OPTS, $selection);

    getopt('bpstaci',\%OPTS);
    ($parmFile, $bgfFile, $filebase, $selection, $atomSelection, $trjFile, $reImageCoords) = 
	($OPTS{p},$OPTS{b},$OPTS{s},$OPTS{t},$OPTS{a},$OPTS{c}, $OPTS{i});
    if (! defined($parmFile) || ! defined($bgfFile)) {
	&showHelp;
    }

    print "Initializing...";
    FileTester($parmFile);
    FileTester($bgfFile);

    print "parm: $parmFile...";
    $Parm = Packages::GetParms->new();
    die "Error in Paramater file\n" if (! $Parm->IsValidParams($parmFile));

    $trjFile = $Parm->{Files}{trajectory} if (! defined($trjFile));
    FileTester($trjFile);
    $trjFile = abs_path($trjFile);

    if (! $Parm->{isLayout} ) {
        $LAYOUT = determineHelixLayout($Parm);
    } else {
        $LAYOUT = setLayout($Parm);
    }

    if (! defined($filebase)) {
	$filebase = "./";
    } else {
	execCmd("mkdir -p $filebase");
    }

    $SELECT = TrjSelections($selection);
    $currDir = cwd();

    if ($Parm->{Files}{trjType} == 1) {
	$getSnapshot = \&ParseAmberTrj;
	$getByteOffset = \&GetAmberByteOffset; 
    } else {
	$getSnapshot = \&ParseLAMMPSTrj if ($Parm->{Files}{trjType} == 2);
	$getByteOffset = \&GetLammpsByteOffset;
    }

    for $i (1 .. $Parm->{Molecule}{total_bases}) {
	$RES->{RESNUM}{$i} = 1;
    }
  
    $reImageCoords = 0 if (! defined($reImageCoords) or $reImageCoords !~ /^1$/);
    execCmd("mkdir -p bgffiles outfiles");
    print "Done\n";
}

sub getAtomSelection {
    my ($selection, $resInfo, $parm) = @_;
    my (%ATMLIST, $i, $j, @tmp, $group, $rec, @tmp1, $inStr);

    if (! defined($selection)) {
	$selection = "all";
	for $i (1 .. $parm->{Molecule}{total_bases}) {
	    next if (! exists($resInfo->{$i}));
	    for $j (keys %{ $resInfo->{$i}{ATOMS} }) {
		$ATMLIST{$j} = $i;
	    }
	}
    } else {
	if ($selection =~ /\s+/) {
	    @tmp = split /\s+/, $selection;
	} else {
	    $tmp[0] = $selection;
	}
	for $i (@tmp) {
	    if (! -e $i) {
		$rec = GetSelections($i);
		for $j (keys %{ $rec }) {
		    $ATMLIST{$j} = $rec->{$j};
		}
	    } else {
		open SELECTFILE, $i || die "ERROR: Cannot open $i: $!\n";
		$group = 1;
		@tmp1 = ();
		while (<SELECTFILE>) {
		    chomp;
		    $inStr = $_;
		    if ($inStr =~ /^Group (\d+) Atoms \d+/) {
			$group = $1;
			@tmp1 = ();
		    }
		    next if ($inStr !~ /^\d+/);
		    while ($inStr =~ /(\d+\-\d+)/g) {
			push @tmp1, $1;
		    }
		    $inStr = $_;
		    while ($inStr =~ /\s+(\d+)\,/) {
			push @tmp1, $1;
		    }
		    for $i (@tmp1) {
			if ($i =~ /(\d+)\-(\d+)/) {
			    for $j ($1 .. $2) {
				$ATMLIST{$j} = $group;
			    }
			} else {
			    $ATMLIST{$i} = $group;
			}
		    }
		}
		close SELECTFILE;
	    }
	}
    }

    die "ERROR: \"$selection\" is not a valid atom list\n" if (! keys %ATMLIST);
    return \%ATMLIST;
}

sub determineHelixLayout {
    my ($P_File) = $_[0];
    my (@helix);

    my $hl = Packages::HelixLayout->spawn();
    $hl->DetermineHelixLayout(
                              $P_File->{"Molecule"}->{"major_groove"},
                              $P_File->{"Molecule"}->{"minor_groove"},
                              $P_File->{"Molecule"}->{"is3PrimeIn"},
                              $P_File->{"Molecule"}->{"bases_at_end"},
                              $P_File->{"Molecule"}->{"total_bases"},
                              $P_File->{"Molecule"}->{"crossovers"}
                              );

    @helix = $hl->GetHelixInfo();
    return \@helix;
}

sub setLayout {
    my ($P_File) = $_[0];
    my (@helix, $i, $j, $helixC, $chainC, $regionC, @tmp, $k, $rec, $len);

    $helixC = $chainC = $regionC = 0;

    for $i (sort numerically keys %{ $P_File->{Layout} }) {
        $chainC = 0;
        for $j (sort numerically keys %{ $P_File->{Layout}{$i} }) {
            $regionC = 0;
            @tmp = split /\,/, $P_File->{Layout}{$i}{$j};
            for $k (@tmp) {
                if ($k =~ /(\d+)\-(\d+)/) {
                    $rec = (
                            {
                                "StartUnit" => $1,
                                "EndUnit"   => $2,
                                "Strand"    => $chainC + 1,
                            }
                            );
                } else {
                    $k = Trim($k);
                    $rec = (
                            {
                                "StartUnit" => $k,
                                "EndUnit"   => $k,
                                "Strand"    => $chainC + 1,
                            }
                            );
                }
                $helix[$helixC][$chainC][$regionC] = $rec;
                $len = abs($rec->{EndUnit} - $rec->{StartUnit}) + 1;
                $HC->{$helixC}{LEN} += $len;
                $regionC++;
                $rec = ();
            }
            $chainC++;
        }
        $helixC++;
    }

    return \@helix;
}

sub numerically {
    ($a<=>$b);
}

sub showHelp {
    my ($error_str) = "usage: $0 -p parmFile -b bgf file -t \"trjSelect\" -s [save prefix] -a [atom selection] -c [trajectory] -i [reimage]\n";

    $error_str .= "\nOptions:\n";
    $error_str .= sprintf("%-30s%-50s\n", "-p parmfile:", "The parameter file");
    $error_str .= sprintf("%-30s%-50s\n", "-b bgf file:", "The BGF file");
    $error_str .= sprintf("%-30s%-50s\n", "-t trajSelection:", "trajectory selection. ex \":It1-10:2\" frames 1 - 10 every 2");
    $error_str .= sprintf("%-30s%-50s\n", " ", "can have multiple selections within the quotes. for all use \"*\"");
    $error_str .= sprintf("%-30s%-50s\n", "-s [save prefix]:", "(Optional) The directory for the saved files. Default is current");
    $error_str .= sprintf("%-30s%-50s\n", "-a [atom selection]:", "(Optional) The atoms to select as the output. Default is all");
    $error_str .= sprintf("%-30s%-50s\n", "-c [trajectory]:", "(Optional) Location of the trajectory. Overwrites the value in the parm file");
    $error_str .= sprintf("%-30s%-50s\n", "-i [reimage]:","(Optional) Reimage coordinates into unit cell. Default no");
    die "$error_str";
}

sub updateBOX {
    my ($BOX, $headers, $boxType) = @_;
    my ($i, $hasBox, $boxData);
    
    $hasBox = 0;
    if ($boxType == 1) { #amber
	$boxData = sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f", $BOX->{2}{DATA},
			   $BOX->{3}{DATA}, $BOX->{4}{DATA}, 90, 90, 90);
    } else { #lammps
	$boxData = sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f", ($BOX->[0]{hi} - $BOX->[0]{lo}),
			   ($BOX->[1]{hi} - $BOX->[1]{lo}), ($BOX->[2]{hi} - $BOX->[2]{lo}), 90, 90, 90);
    }
    for $i (0 .. $#{ $headers }) {
	if ($headers->[$i] =~ /^CRYSTX/) {
	    $hasBox = 1;
	    $headers->[$i] = $boxData;
	    last;
	}
    }
    
    if (! $hasBox) {
	push @{$headers}, "PERIOD 111";
	push @{$headers}, "AXES   ZYX",
	push @{$headers}, "SGNAME P 1",
	push @{$headers}, $boxData;
	push @{$headers}, "CELLS    -1    1   -1    1   -1    1";
    }
}
	
sub writeOutput() {

    my ($unit_key, $eng_key, $in_energy, $time_avg, $time_stdev);
    my ($avg_total_energy, $component_total_energy, $component_counter);
    my ($out_line, $avg_structure_out, $outfile, $time_total);

    &calculateTotals;

    $avg_total_energy = 0.0;

    mkdir "avg_structures", 0777 if (! -d "avg_structures");

    $out_line = sprintf("%8s", "GROUP");

    for (@eng_header) {
	$out_line .= sprintf("%10s", $_);
    }

    $out_line .= sprintf("%10s\n", "Total");
    $component_counter = 0;
    $avg_structure_out = "";
    for $unit_key ( @unit_array ) {
	$component_total_energy = 0.0;
	$out_line .= sprintf("%8s", $unit_key);
	for $eng_key (@eng_header) {
	    $in_energy = $ENERGIES{$unit_key}{Statistics}{$eng_key}{Average};
	    $component_total_energy += $in_energy;
	    $out_line .= sprintf("%10.2f", $in_energy);
	    $component_counter++;
	}
	$out_line .= sprintf("%12.3f", $component_total_energy) . "\n";
	$time_avg = $ENERGIES{$unit_key}{Statistics}{TIMEINFO}{Average};
	$time_stdev = $ENERGIES{$unit_key}{Statistics}{TIMEINFO}{StDev};
	$time_total = $ENERGIES{$unit_key}{Statistics}{TIMEINFO}{Total};
	$avg_structure_out .= sprintf("%5d%12.3f%10.3f%5d\n", 
				      $unit_key, $time_avg, $time_stdev, $time_total);
	$avg_total_energy += $time_avg;
	$component_total_energy = 0.0;
    }

    $outfile = "avg_structures/Total_all_Energies.dat";
    open AVGSTRU, "> $outfile" or die "Cannot create $outfile: $!\n"; 
    print AVGSTRU $out_line;
    close AVGSTRU;

    $outfile = "avg_structures/TotalEng.dat";
    open OUTFILE, "> $outfile" or die "Cannot create $outfile: $!\n";
    print OUTFILE "$avg_structure_out";
    close OUTFILE;

    $avg_total_energy = sprintf("%12.3f", $avg_total_energy);
    print "Done\n\tAverage Total Energy: $avg_total_energy kcal/mol\n";

    print "\tCreating graphs...";
    &createGraphs;
    print "Done\n";

}

sub createGraphs {

    my (@eng_comp_array, $this_unit, $total_energies, $avg);
    my ($i, $j, $outfile, $eng_total, $eng_data, $running_total);
    my ($graph_2d_out, $graph_3d_out, $my_total, $my_stdev, $my_num);

#   populate the arrays with the keys for the time and energy contributions


    mkdir "2d_graphs", 0777 if (! -d "2d_graphs");
    mkdir "3d_graphs", 0777 if (! -d "3d_graphs");

#    Create Per_Energy Graphs
    for $i (@eng_header) {
	$graph_3d_out = "";
	$graph_2d_out = "";
	for $j (@time_array) {
	    $eng_total = $eng_data = "";
	    for $this_unit ( @unit_array ) {
		$eng_data .= sprintf("%5d%5d%10.2f\n", $this_unit, $j, 
				     $ENERGIES{$this_unit}{$j}{$i});
		$eng_total .= "$ENERGIES{$this_unit}{$j}{$i} ";
	    }
	    ($avg, $my_stdev, $my_total) = STDev($eng_total);
	    $my_num = $my_total/$avg if ($avg > 0);
	    $graph_2d_out .= sprintf("%5d%10.2f %6.2f\n", $j,$my_total, $my_stdev);
	    $graph_3d_out .= $eng_data . "\n";
	}

	$outfile = "2d_graphs/${i}.dat";
	open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
	print OUTFILE "$graph_2d_out";
	close OUTFILE;
	
	
	$outfile = "3d_graphs/${i}.dat";
	open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
	print OUTFILE "$graph_3d_out";
	close OUTFILE;
    }

    $graph_2d_out = $graph_3d_out = "";

#    Create Total Energy per unit graphs
    for $j (@time_array) {
	$running_total = $eng_data = "";
	for $this_unit ( @unit_array ) {
	    $eng_total = 0.0;
	    for $i (@eng_header) {
		$eng_total += $ENERGIES{$this_unit}{$j}{$i};
	    }
	    $eng_data .= sprintf("%5d%5d%10.2f\n", $this_unit, $j, $eng_total);
	    $running_total .= "$eng_total ";
	}
	chop $running_total;
	($avg, $my_stdev, $my_total) = STDev($running_total);
	$my_num = $my_total/$avg;
	$graph_2d_out .= sprintf("%5d%10.2f\n", $j,$my_total);
	$total_energies .= "$my_total ";
	$eng_data .= sprintf("%-5s%5d%10.2f%6.2f\n", "#Avg", $j, $my_total, $my_stdev);
	$eng_data .= sprintf("%-5s%5d%10.2f\n", "#Tot", $j, $my_num);
	$graph_3d_out .= ($eng_data . "\n");
    }

    chop $total_energies;
    ($avg, $my_stdev, $my_total) = STDev($total_energies);
    $my_num = $my_total/$avg;
    $graph_2d_out .= sprintf("#%5s%10.2f %6.2f\n", "Avg", $avg, $my_stdev);

    $outfile = "2d_graphs/TOTAL.dat";
    open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
    print OUTFILE "$graph_2d_out";
    close OUTFILE;
        
    $outfile = "3d_graphs/TOTAL.dat";
    open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
    print OUTFILE "$graph_3d_out";
    close OUTFILE;
    
#    create average structures
    for $j (@eng_header) {
	$eng_data = "";
	for $this_unit ( @unit_array ) {
	    $eng_data .= sprintf("%5d", $this_unit);
	    $eng_data .= sprintf("%10.2f", 
				 $ENERGIES{$this_unit}{Statistics}{$j}{Average});
	    $eng_data .= sprintf("%10.2f", 
				 $ENERGIES{$this_unit}{Statistics}{$j}{StDev});
	    $eng_data .="\n";
	}
	$outfile = "avg_structures/${j}.dat";
	open OUTFILE, "> $outfile" || die "Cannot create $outfile: $!\n";
	print OUTFILE "$eng_data";
	close OUTFILE;
    }
}

sub calculateTotals {

    my ($unit_key, $time_key, $eng_key, $unit_val, $time_val, $eng_val, $counter);
    my ($avg, $StDev, $total);
    $avg = $StDev = $total = "";

# Reset the Statistics information
    for $unit_key (@unit_array) {
	if ($ENERGIES{$unit_key}{Statistics}) {
	    delete $ENERGIES{$unit_key}{Statistics};
	}
    }

# Calc the Energy per unit per Contributor
    for $time_key (@time_array) {
	for $unit_key (@unit_array) {
	    $total = 0;
	    for $eng_key (@eng_header) {
		$eng_val = $ENERGIES{$unit_key}{$time_key}{$eng_key};
		$ENERGIES{$unit_key}{Statistics}{$eng_key}{"DataVals"} .= "$eng_val ";
		$total += $eng_val;
	    }
	    $ENERGIES{$unit_key}{Statistics}{TIMEINFO}{DataVals} .= "$total ";
	    $total = 0;
	}
    }
		
    $total = 0;
# Store the Energy per unit per Contributor
    for $unit_key (@unit_array) {
	while ( ($eng_key, $eng_val) = each %{ $ENERGIES{$unit_key}{Statistics} }) {
	    chop $eng_val->{"DataVals"};
	    ($avg, $StDev, $total) = STDev($eng_val->{"DataVals"});
	    $eng_val->{"Average"} = $avg;
	    $eng_val->{"StDev"} = $StDev;
	    if ($avg < 0 || $avg > 0) {
		$eng_val->{"Total"} = ($total/$avg);
	    } else {
		$eng_val->{"Total"} = 0;
	    }
	    
	}
    }

}
