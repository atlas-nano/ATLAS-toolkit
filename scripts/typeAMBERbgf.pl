#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(PrintProgress);

sub init;
sub showUsage;
sub typeMols;
sub execCmd;

my ($sfile, $rfile, $struct, $atomSel, $charge, $fftype, $sFF, $multi);
my ($aList, $ffStr);

$|++;
&init;
print "Gettting data from struct $rfile->{type} file $rfile->{name}...";
$struct->read($rfile->{name}, $rfile->{type});
print "Done\nSelecting atoms...";
$aList = $struct->findAtoms($atomSel,-1,1);
die "ERROR: No valid atoms found\n" if (! keys %{ $aList });
print "Done\n";
$ffStr = typeMols($struct,$aList,$charge,$fftype,$sFF);
print "Saving file...";
$struct->write($sfile->{name},$sfile->{type});
print "Done\n";
&pruneFF($ffStr) if($sFF);

sub typeMols {
	my ($s,$list,$chrg,$ff, $getFFparms) = @_;
	my ($i, $j, $k, $molFile, $molData, $bondCounter, $molN, $tmpMolData);
	my ($c, $tot, $start, $pStr, $index, $curr, $tmp, $ffStr, $ffref);

	$start = time();
	$c = 0;
	$tot = scalar keys %{ $list };
	$ffref = "$ENV{AMBERHOME}/dat/leap/parm/gaff.dat";
	$ffref = "$ENV{AMBERHOME}/dat/leap/parm/parm10.dat" if ($fftype eq "amber");
	$tmpMolData = MolData->new();
	print "Assigning $fftype atom types to system... Getting timing...\r";
	for $i (keys %{ $list }) {
		$j = $list->{$i};
		$molData = getMolAtoms($j);
		$molFile = MolData->new();
		$bondCounter = 0;
		$tmp = ();
		for $k (sort {$a <=> $b} keys %{ $molData }) {
			$molData->{$k}->("charge",0);
			$molFile->insertAtom($molData->{$k});
			push @{ $tmp }, $k;
			$molN = uc $molData->{$k}->resname;
			$bondCounter += $s->atoms($k)->bondlist->count;
		}
		$molFile->count->{atoms} = $j->count;
		$molFile->count->{molecule} = 1;
		$molFile->header->{descrp} = "mol_${i}\n";
		$molFile->header->{remark} = "SMALL\n";
		$molFile->header->{footer} = "END\n";
		$molFile->write("__${i}.mol2","mol2");
		&execCmd("$ENV{AMBERHOME}/bin/antechamber -i __${i}.mol2 -fi mol2 -fo mol2 -o __${i}.amber.mol2 -m 1 -rf $molN -rn $molN -pf y -nc $chrg -at $ff -c bcc -m $multi");
		if ($getFFparms) {
			&execCmd("$ENV{AMBERHOME}/bin/parmchk -i __${i}.amber.mol2 -f mol2 -p $ffref -o __${i}.frcmod");
			$ffStr .= "__${i}.frcmod ";
		}
		$tmpMolData->read("__${i}.amber.mol2", "mol2");
		for $index (0 .. $#{ $tmp }) {
			$k = $index + 1;
			$curr = $s->atoms($tmp->[$index]);
			$curr->("atmname", $tmpMolData->atoms($k)->atmname);
			$curr->("charge", $tmpMolData->atoms($k)->charge);
			$curr->("fftype", $tmpMolData->atoms($k)->fftype);
		}
		&execCmd("rm -fr __${i}.mol2 __${i}.amber.mol2 sqm.in sqm.out sqm.pdb");
		$molFile = ();
		$molData = ();
		undef $molFile;
		undef $molData;
		$c++;
		$pStr = "Assigning $fftype atom types to system... molecule $i...";
		PrintProgress($c, $tot, $start, $pStr);
	}
	print "Assigning $fftype atom types to system... molecules $tot of $tot...Done                           \n";
	return $ffStr;
}

sub pruneFF {
	my ($ffStr) = $_[0];
	my ($sname, $ffref, $refFF);

	$sname = $sfile->{name};
	$sname =~ s/\.\w+$//;

	$ffref = "$ENV{AMBERHOME}/dat/leap/parm/gaff.dat";
    	$ffref = "$ENV{AMBERHOME}/dat/leap/parm/parm10.dat" if ($fftype eq "amber");

    	$refFF = basename($ffref);
    	$refFF =~ s/\.\w+$//;
    	$refFF .= ".ff";

    	&execCmd("$Bin/amberff2cerius.pl $ffref $ffStr");
    	`cat $Bin/../ff/AMBERFF_header.ff $refFF > ${sname}.ff`;
    	&execCmd("$Bin/pruneCerius2FF.pl -b ${sname}.bgf  -f ${sname}.ff -s ${sname}.ff");
    	&execCmd("rm -fr $ffStr $refFF")
}

sub getMolAtoms {
	my ($i) = @_;
	my ($aList, $j, $c);

	for $j (keys %{ $i->atomlist }) {
		$c = $i->atomlist->{$j};
		$aList->{ $i->atom($c)->index } = $i->atom($c);
	}
	return $aList;
}

sub init {
	my (%OPTS);

	getopt('isocfpm', \%OPTS);
	die &showUsage . "\n" if (! exists($OPTS{i}));

	print "Initialzing...";
	$ENV{AMBERHOME} = "$Bin/../codes/amber12/" if (!exists($ENV{AMBERHOME}));
	($rfile->{name}, $sfile->{name}, $atomSel, $fftype, $charge, $sFF, $multi) = 
		($OPTS{i}, $OPTS{s}, $OPTS{o}, $OPTS{f}, $OPTS{c}, $OPTS{p}, $OPTS{m});
	$struct =  MolData->new();
	$struct->testFile($rfile->{name});
	if ($rfile->{name} =~ /\.(\w+)$/) {
		$rfile->{type} = lc $1;
	} else {
		die "ERROR: Cannot determine file type from $rfile->{name}!\n";
	}
	$atomSel = "index>0" if (! defined($atomSel));
	$charge = 0 if (! defined($charge) or $charge !~ /^\-?\d+\.?\d*$/);
	$fftype = "amber" if (! defined($fftype) or $fftype !~ /(amber|gaff)/i);
	$fftype =~ /(amber|gaff)/i;
	$fftype = lc $1;
	$sfile->{name} = $struct->getFileName($rfile->{name}) if (! defined($sfile->{name}));
	if ($sfile->{name} =~ /\.(\w+)$/) {
		$sfile->{type} = lc $1;
	} else {
		die "ERROR: Cannot determine file type from $sfile->{name}!\n";
	}
	system("rm -fr ambertype.log");
	$sFF = 0 if (! defined($sFF) or $sFF !~ /1|yes/i);
	$sFF = 1 if ($sFF =~ /1|yes/i);
	$multi = 1 if (! defined($multi) or $multi !~ /(\d+)/);
	$multi = $1 if ($multi =~ /(\d+)/);

	print "Done\n";
}

sub execCmd {
	my ($cmdStr, $output) = @_;
	my ($terminal);

	$terminal = 1;
	$terminal = 0 if (defined($output) && ! $output);
	if (! defined($output)) {
		$cmdStr .= ">> ambertype.log";
		$output = "ambertype.log";
	} elsif (defined($output) && ! $output) {
		$cmdStr .= ">& ambertype.log";
		$output = "ambertype.log";
	} else {
		$cmdStr .= ">& $output";
	}
	system("echo \"\ncommand:'$_[0]'\" >> $output");
	if (system($cmdStr)) {
		die "ERROR: Check ambertype.log\n" if ($terminal);
	}
}

sub showUsage {
	my $usage = <<DATA;
usage: $0 -i input structure -o (molecule selection = all) -c (molecule charge) -m (spin_multiplicity) -f (fftype) -s (save_name) -p (save ff parms)
OPTIONS:
		-i input structure:	REQUIRED. The following file formats are supported: BGF, PDB, MSI, MOL2
		-o atom selection	OPTIONAL. Any valid bgf field expression. E.g. resname eq 'WAT' will select
			all the "WAT" residues while index > 10 will select all indices > 10.
			combine multiple expressions to make complicated selections: e.g.
			(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
		-c molecule charge	OPTIONAL. charge on each molecule. Default 0
		-m spin multiplicity OPTIONAL. Multiplicity (2S+1). Default 1
		-f fftype			OPTIONAL. either amber (default) or gaff
		-p save ff parms	OPTIONAL. save the ff parms of the new system. Default is 0
		-s savename:	 OPTIONAL. Will assume \$prefix_mod.\$suffix if not specified.
DATA

	return $usage;
}
