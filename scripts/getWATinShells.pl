#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
no warnings "recursion";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester GetBondLength CoM LoadFFs);
use BOX qw(GetCellFromHeader GetBox);;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use ManipAtoms qw(GetMols BuildAtomSelectionString SelectAtoms ReimageAtoms GetAtmData);

sub numerically { ($a<=>$b); };

my ($bgfFile, $shellStr, $saveName, $solvAtomList, $soluAtomList, $reImage, $reCenter, $saveBGF);
my ($ATOMS, $BONDS, $HEADERS, $BOX, $SOLUTE, $SOLVENT, $WSHELLS, $watShells, $MOLS);

$|++;
$watShells = &init;
print "Getting atom data from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS, $BONDS);
$SOLUTE = SelectAtoms($soluAtomList,$ATOMS);
die "ERROR: No valid solute atoms found!\n" if(! keys%{ $SOLUTE });
$SOLVENT = SelectAtoms($solvAtomList,$ATOMS);
die "ERROR: No valid solvent atoms found!\n" if(! keys%{ $SOLVENT });
$SOLUTE = GetAtmData($ATOMS,$SOLUTE);
print "\nReimaging atoms in box around solute...";
$BOX = GetBox($ATOMS, undef, $HEADERS);
&centerOnSolute($ATOMS, $SOLUTE) if ($reCenter);
&ReimageAtoms($ATOMS,$BONDS,$MOLS,GetCellFromHeader($HEADERS),$SOLVENT) if ($reImage);
print "\nGetting solvent molecules in shells around solute...";
$SOLVENT = GetMols($ATOMS,$BONDS,$SOLVENT);
&getCoM($ATOMS,$SOLVENT);
$WSHELLS = getSolvShells($SOLUTE,$SOLVENT,$watShells,$BOX);
print "Done\nCreating $saveName...";
&createFile($WSHELLS, $ATOMS, $SOLUTE, $SOLVENT, $saveName);
&updateBGF($ATOMS, $SOLUTE, $SOLVENT, $WSHELLS);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveBGF);
print "Done\n";
#&createShellBGFs($ATOMS, $WSHELLS, $SOLUTE, $SOLVENT, $saveName);

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
	$shellID += 1;
}

sub getSolvShells {
	my ($solu, $solv, $shells,$box) = @_;
	my ($i, $j, $dists, $val, $prev, @tmp, $shellList, $shellID);

	$prev = 0;
	for $i (keys %{ $solv }) {
		$dists = ();
		for $j (keys %{ $solu }) {
			$dists->{GetBondLength($solv->{$i}{CoM},$solu->{$j},$box)} = 1;
		}
		@tmp = sort numerically keys %{ $dists };
		$shellID = getShellID($tmp[0],$shells);
		push @{ $shellList->[$shellID] }, $i;
		$solv->{$i}{SHELL} = $shellID;
	}

	return $shellList;
}

sub getCoM {
	my ($atoms, $mols) = @_;
	my ($com, $moldata, $i, $ret);

	for $i (keys %{ $mols }) {
		$moldata = GetAtmData($atoms, $mols->{$i}{MEMBERS});
		$mols->{$i}{CoM} = CoM($moldata);
	}
}

sub createShellBGFs {
	my ($atoms, $shells, $solute, $solvent, $savePrefix) = @_;
	my ($MOLECULE, $CONS, $atmList, $i, $bgfName, @tmp, $atom, $printStr, $mol);

	$savePrefix =~ s/\.\w+$//;
	$printStr = "Creating BGF file";

	$bgfName = "${savePrefix}.solute.bgf";
	($MOLECULE, $CONS, undef) = GetBGFAtoms($solute, $ATOMS, $BONDS);
	 &addHeader($MOLECULE, $HEADERS);
	&createBGF($MOLECULE, $CONS, $bgfName);
	for $i (0 .. $#{ $shells }) {
		$MOLECULE = ();
		$CONS = ();
		$atmList = ();
		if ($i == $#{ $shells }) {
			$bgfName = "${savePrefix}.bulk.bgf";
		} else {
			$bgfName = "${savePrefix}.shell" . ($i+1). ".bgf";
		}
		print "${printStr} ${bgfName}...\r";
		for $mol (@{ $shells->[$i] }) {
			for $atom (keys %{ $solvent->{$mol}{MEMBERS} }) {
			$atmList->{$atom} = 1;
			}
		}
		($MOLECULE, $CONS, undef) = GetBGFAtoms($atmList, $ATOMS, $BONDS);
		&addHeader($MOLECULE, $HEADERS);
		&createBGF($MOLECULE, $CONS, $bgfName);
	}
	printf "${printStr}s...Done%50s\n", "";

}

sub createFile {
	my ($wshells, $atoms, $solute, $solvent, $save) = @_;
	my (@tmp, @tmp2, $i, $j, $dofStr, $moldof, $sol);
	
	open OUTDATA, "> $save" || die "ERROR: Cannot create $save: $!\n";
	print OUTDATA "Total Groups: " . ($#{ $wshells } + 1) . "\n";

	$dofStr = getShakeDOF($solute, $atoms);
	print OUTDATA "Group 1 Atoms " . scalar(keys %{ $solute }) . "\n";
	@{ $sol } = keys %{ $solute };
	print OUTDATA getRange($sol);
	
	$i = 0;
	while ($i <= $#{ $wshells }) {
		if (!defined($wshells->[$i])) {
			splice @{ $wshells }, $i, 1;
		} else {
			$i++;
		}
	}
	for $i (0 .. $#{ $wshells }) {
		@tmp = sort numerically @{ $wshells->[$i] };
		@tmp2 = ();
		$moldof = 0;
		for $j (@tmp) {
			$moldof += getShakeDOF($solvent->{$j}{MEMBERS}, $atoms);
			for (keys %{ $solvent->{$j}{MEMBERS} }) {
				push @tmp2, $_;
			}
		}
		$dofStr .= " $moldof";
		print OUTDATA "Group " . ($i + 2) . " Atoms " . ($#tmp2 + 1) . "\n";
		print OUTDATA getRange(\@tmp2);
	}
	print OUTDATA "Constraints\n$dofStr\n";
	close OUTDATA;
}

sub getRange {
	my ($atmList) = @_;
	my ($j, $start, $prev, $counter, $groups);
   
	$start = $prev = -1;
	$groups = "";
	$counter = 0;
	for $j (sort numerically @{ $atmList }) {
		if ($start == -1) {
			$start = $j;
		} elsif (($j - $prev) > 1) {
			$groups .= "${start} - ${prev} ";
			$start = $j;
			$counter++;
		}
		$prev = $j;
		if ($counter == 10) {
			$counter = 0;
			$groups = substr($groups, 0, -2);
			$groups .= "\n";
		}
	}
	$groups .= "${start} - ${prev} ";
	$groups = substr($groups, 0, -1);

	return "$groups\n";
}

sub updateBGF {
	my ($atoms, $solute, $solvent, $shells) = @_;
	my ($i, $j, $tot);

	$tot = scalar(@{ $shells }) + 1;
	for $i (keys %{ $atoms }) {
		$atoms->{$i}{CHAIN} = chr(64 + $tot);
	}
	for $i (keys %{ $solute }) {
		$atoms->{$i}{CHAIN} = "A";
	}
	for $i (keys %{ $solvent }) {
		for $j (keys %{ $solvent->{$i}{MEMBERS} }) {
		$atoms->{$j}{CHAIN} = chr(65 + $solvent->{$i}{SHELL});
		}
	}
}

sub init {
	my (@tmp, $i, $j, $tmp, $SHELL, %OPTS, $usage, $solvStr, $soluStr);

	getopt('bswuvci',\%OPTS);

	($bgfFile, $shellStr, $saveName, $solvStr, $soluStr, $reImage, $reCenter) = 
	($OPTS{b},$OPTS{w},$OPTS{s},$OPTS{v},$OPTS{u},$OPTS{i},$OPTS{c});
	$usage = &showUsage;
	for ("b", "w") {
		die "$usage\n" if (! $OPTS{$_});
	}
	
	print "Initializing...";
	FileTester($bgfFile);
	
	if ($shellStr =~ /\s+/) {
		@tmp = split /\s+/, $shellStr;
	} else {
		$tmp[0] = $shellStr;
	}

	for $i (@tmp) {
		$tmp->{$i} = 1 if ($i =~ /(\d+\.*\d*)/); 
	}
	die "Expected at least one decimal(integer) value for the shell limts. Got \"$shellStr\"\n" 
		if (! keys %{ $tmp });
	@{ $SHELL } = sort numerically keys %{ $tmp };

	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/.solvshell\.grps/;
	}

	print "total watershells: " . ($#tmp + 1) . "...Done\n";
	$solvStr = "resname eq 'WAT'" if(!defined($solvStr));
	$soluStr = "moleculeid==1" if (!defined($soluStr));
	$solvAtomList = BuildAtomSelectionString($solvStr);
	$soluAtomList = BuildAtomSelectionString($soluStr);
	$reImage = 1 if (! defined($reImage) or $reImage !~ /0|no/i);
	$reImage = 0 if ($reImage =~ /0|no/i);
	$reCenter = 0 if (! defined($reCenter) or $reCenter !~ /1|yes/i);
	$reCenter = 1 if ($reCenter =~ /1|yes/i);
	$saveBGF = $saveName;
	$saveBGF =~ s/\.\w+$//;
	$saveBGF .= ".bgf";
	return $SHELL;
}

sub getShakeDOF {
	my ($select, $atoms) = @_;
	my ($i, $dof);
	$dof = 0;
	for $i (keys %{ $select }) {
		if ($atoms->{$i}{FFTYPE} =~ /^H/) { #if this is a hydrogen
			$dof++;
			if (${ $atoms->{$i}{MOLSIZE} } == 3 ) {
				$dof += 0.5;
			}
		}
	}
	return $dof;
}

sub showUsage {

	my ($usage) = <<DATA;
This script will break up the bgf file into:  solute, ions, water shell 1...n, and bulk
usage: $0 -b bgf file le -w shell distance(s) -s [save name] -c [recenter] -i [reimage] -u [solute] -v [solvent]
	options:
	-b bgf file: the location of the bgf file
	-w water shell distances(s): the distance from the surface of the 1st water shell.
	   you can specify multiple water shell distances by enclosing them in "" qoutes
	-s [save name]: (optional) the name to save the group information generated by the script
	   if not specified, will be bgf_file.grps. Each water shell will be save_name_shellx.bgf
	-c [recenter]: (optional) place the center of mass of the solute at the box center. Default no
	-i [reimage]: (optional) wrap solvent atoms back into simulation cell. Default yes
DATA

	return $usage;
}
