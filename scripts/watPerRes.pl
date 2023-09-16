#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo sortByRes GetBGFAtoms addHeader createBGF);
use General qw(GetSelections FileTester CoM GetBondLength PrintProgress);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetAtmList GetSolvent SplitAtomsByMol);

sub usage;
sub getCons;
sub numerically { ($a<=>$b); }
sub findWat;
sub writeGrpFile;
sub getList;

my ($bgfFile, $saveName, $selection);
my ($ATOMS, $BONDS, $SOLUTE, $SOLVENT, $SHELL, $RES, $GRPS, $HEADERS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Done\nParsing atom/residue selection...";
$SOLUTE = sortByRes($ATOMS, GetAtmList(GetSelections($selection, 0), $ATOMS));
die "ERROR: No valid atoms selected!\n" if (! keys %{ $SOLUTE });
$selection = ["NrWAT"];
$SOLVENT = SplitAtomsByMol($ATOMS, GetAtmList(GetSelections($selection, 0), $ATOMS));
print "Done\n";
$GRPS = findWat($ATOMS, $SOLUTE, $SOLVENT, $SHELL);
print "Writing vac grp file $saveName...";
&writeGrpFile($ATOMS, $SOLVENT, $GRPS, $saveName);
$saveName =~ s/\.\w+$//;
$saveName .= "_resgrps.bgf";
print "Done\nWriting compatible bgf file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub writeGrpFile {
    my ($atoms, $solv, $grpData, $outfile) = @_;
    my ($soluIons, $i, $tot, $count, $aList, $j, $tmp, $shellN);
    
    for $i (keys %{ $solv }) {
	for $j (keys %{ $solv->{$i} }) {
	    $tmp->{$j} = 1;
	}
    }

    for $i (keys %{ $atoms } ) {
	next if (exists($tmp->{$i}));
	$soluIons->{$i} = 1;
    }
    $soluIons = SplitAtomsByMol($atoms, $soluIons);

    $tot = scalar(keys %{ $soluIons }) + scalar(keys %{ $grpData->{SHELL} }) + scalar(keys %{ $grpData->{RES}} );
    open GRPFILE, "> $outfile" or die "ERROR: Cannot write to $outfile: $!\n";
    print GRPFILE "Total groups $tot\n";
    $count = 0;
    for $i (sort numerically keys %{ $soluIons }) {
	$count++;
	$tot = scalar(keys %{ $soluIons->{$i} });
	print GRPFILE "Group $count Atoms $tot\n";
	$aList = getList($atoms, $soluIons->{$i}, $count);
	print GRPFILE $aList;
	print GRPFILE "\n";
    }

    for $i ("RES", "SHELL") {
	next if (! keys %{ $grpData->{$i} });
	for $j (sort numerically keys %{ $grpData->{$i} }) {
	    $count++;
	    $shellN->{$count} = "${i} ${j}";
	    $tot = scalar(keys %{ $grpData->{$i}{$j} });
	    print GRPFILE "Group $count Atoms $tot\n";
	    $aList = getList($atoms, $grpData->{$i}{$j}, $count);
	    print GRPFILE $aList;
	    print GRPFILE "\n";
	}
    }
    close GRPFILE;
    $outfile =~ s/\.w+$//;
    $outfile .= "_resgrplist.dat";
    open OUTFILE, "> $outfile" or die "ERROR: Cannot create $outfile: $!\n";
    for $i (sort numerically keys %{ $shellN }) {
	print OUTFILE "$i $shellN->{$i}\n";
    }
    close OUTFILE;
}

sub getList {
    my ($atoms, $atomList, $resid) = @_;
    my ($i, $start, $prev, $groups, $counter);

    $start = $prev = -1;
    $groups = "";
    $counter = 0;
    for $i (sort numerically keys %{ $atomList }) {
	$atoms->{$i}{RESNUM} = $resid;
	if ($start == -1) {
	    $start = $i;
	} elsif (($i - $prev) > 1) {
	    $groups .= "${start}-${prev}, ";
	    $start = $i;
	    $counter++;
	}
	$prev = $i;
	if ($counter == 10) {
	    $counter = 0;
	    $groups = substr($groups, 0, -2);
	    $groups .= "\n";
	}
    }
    $groups .= "${start}-${prev}, ";
    $groups = substr($groups, 0, -2);

    return $groups;
}

sub findWat {
    my ($atoms, $resList, $solvent, $shellData) = @_;
    my ($i, $dist, $solvCenter, $mol, $j, $k, $min, $tot, $count); 
    my (%WATLOC, $soluCenter, @tmp, @shellList, $found, $start, $strLen);

    @shellList = sort numerically keys %{ $shellData };
    $tot = scalar(keys %{ $solvent });

    $count = 0;
    $start = time();
    for $i (keys %{ $solvent }) {
	$count++;
	$strLen = PrintProgress($count, $tot, $start, "Calculating solvent distance to solute...");
	undef $dist;
	$mol = ();
	for $j (keys %{ $solvent->{$i} }) {
	    $mol->{$j} = $atoms->{$j};
	}
	$solvCenter = CoM($mol);
	for $j (keys %{ $resList }) {
	    $mol = ();
	    for $k (keys %{ $resList->{$j}{ATOMS} }) {
		$dist->{ GetBondLength($solvCenter, $atoms->{$k}) } = $j;
	    }
	}
	@tmp = sort numerically keys %{ $dist };
	$min = shift @tmp;
	$found = 0;
	for $j (1 .. ($#shellList + 1)) {
	    if ($shellData->{$j} > $min) {
		$found = $j+2; # found in shell + 1
		$found = 1 if ($j == 1); # found in first shell
		last;
	    }
	}
	if (! $found ) {
	    for $j (keys %{ $solvent->{$i} }) {
		$WATLOC{SHELL}{0}{$j} = 1;
	    }
	} elsif ($found == 1) {
	    for $j (keys %{ $solvent->{$i} }) {
		$WATLOC{RES}{ $dist->{$min} }{$j} = 1;
	    }
	} else {
	    for $j (keys %{ $solvent->{$i} }) {
		$WATLOC{SHELL}{ $found -1 }{$j} = 1;
	    }
	}
    }
    $start = time();
    printf "Calculating distance for solvent molecule... %-${strLen}s\n", "${start}s elapsed..Done";
    return \%WATLOC;
}

sub init {
    my (%OPTS, $resSelect, @tmp, $shellStr, $counter, $increment, $j);
    getopt('brso',\%OPTS);
    ($bgfFile, $shellStr, $resSelect, $saveName) = ($OPTS{b},$OPTS{s},$OPTS{r}, $OPTS{o});
    for ($bgfFile, $resSelect) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    $SHELL = ();

    FileTester($bgfFile);
    if ($resSelect =~ /\s+/) {
	@{ $selection } = split /\s+/, $resSelect;
    } else {
	$selection->[0] = $resSelect;
    }
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }

    $saveName =~ s/_nvt_/_vac_/;

    if ($shellStr =~ /\s+/) {
        @tmp = split /\s+/, $shellStr;
    } else {
        $tmp[0] = $shellStr;
    }

    $counter = 1;
    $increment = 0;
    for $j (@tmp) {
        if ($j =~ /(\d+\.*\d+)/) {
            $SHELL->{$counter} = $1;
            if ($counter == 1) {
                $increment = $1;
            } elsif (abs($SHELL->{$counter} - $SHELL->{($counter - 1)}) > $increment) {
                $increment = abs($SHELL->{$counter} - $SHELL->{($counter - 1)});
            }
            $counter++;
        }
    }
    $SHELL->{1} = 3.6 if (! defined($SHELL));

    print "Done\n";
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -r residue_selection -s (solvent_shell_data) -o (output_name)
Arguments:
 bgf_file: name of bgf_file
 residue_selection: see below
 water_shell_data: distance from surface to group solvent - optional
 output_name: name of file to save - optional
 residue select options:
    [^][:][I|T|N][a|r]
    a   - atom
    r   - residue
    IaX - atom number X
    IrX - residue index X
    TaX - atom type X
    NaX - atom name X
    NrX - residue name X
    Use ":" to specify a range, eg. :Tr1-8 :Ia3-66
    Use "^" to exclude a selection. You can use multiple combinations
    range and exclusion enclosed in quotes, eg, "^:TrIP-IM ^:Ia23-45"
    to exclude residues of type IM and IP and atoms 23-45
DATA

die "\n";

}
