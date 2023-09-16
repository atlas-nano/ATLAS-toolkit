#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use MolData;
use Getopt::Std qw(getopt);
use General qw(GetFileTypeStr GetFFTypeStr);

my ($sfile, $ffType, $solute, $solvent, $cellOpts, $randOpt, $randRotOpt, $ffields, $reversePlace);
my ($tot, $nMolu, $nMolv, $box);

$|++;
&init;
print "Gettting data from solute $sfile->{type} file $sfile->{name}...";
$solute->read($sfile->{name}, $sfile->{type},1);
print "Done\nCreating solvent box...";
&createSolventBox($solvent, $solute, $cellOpts, $randRotOpt);
$solute->write($sfile->{tmp}, "bgf");
print "Done\nEmbedding solute and removing bad contacts...";
($nMolu, $nMolv) = embedMols($solute, $sfile, "_solv_trim_centered.bgf", $ffields, $reversePlace);
print "Done\n";#Loading embeded system...";
$tot = removeSolv($sfile->{save}, $cellOpts, $nMolu, $nMolv, $solute, $reversePlace);
print "Added $tot->{atoms} solvent atoms ($tot->{molecule} molecules)...Done\n"; 
undef($solute);

sub getSoluData {
	my ($sol) = $_[0];
	my ($tot, @tmp);

	%{ $tot } = %{ $sol->count };
	@tmp = sort {$a<=>$b} keys %{ $sol->shash->{resid} };
	$tot->{resid} = pop @tmp;
	return $tot;
}

sub removeSolv {
	my ($bgffile, $cell, $nMsolu, $nMsolv, $solu, $rev) = @_;
	my ($tot, $count, $countStr, $getAtmStr, $i, $solvCount, $soluCount);

	print "Centering cell...";
	print "Done\n";
	print "Computing stats...";
	if ($rev) {
		$countStr = "$Bin/countAtoms.pl -f $bgffile -s \"moleculeid<=${nMsolv}\" -m 1";
	} else {
		$countStr = "$Bin/countAtoms.pl -f $bgffile -s \"moleculeid>${nMsolu}\" -m 1";
	}
	open COUNT, "$countStr |" or die "ERROR while executing $countStr\n";
	while (<COUNT>) {
		if ($_ =~ /^Found (\d+) (atoms|molecule)/) {
			$count->{$2} = $1;
		}
	}
	close COUNT;
	$tot = $count->{molecule}; # total number of solvent molecules
	if (! exists($cell->{total}) or $tot < $cell->{total}) {
		for $i ("atoms", "molecule") {
			$solvCount->{$i} = $count->{$i};
		}
		return $solvCount;
	}

	$tot -= $cell->{total};
	print "Done\nRemoving $tot solvent molecules to achieve $cell->{total}...";
	if ($rev) {
		$getAtmStr = "$Bin/removeMols.pl -b $bgffile -s $bgffile -a \"moleculeid<=${nMsolv}\" -m $tot -r $randOpt >> _out.dat";
	} else {
		$getAtmStr = "$Bin/removeMols.pl -b $bgffile -s $bgffile -a \"moleculeid>${nMsolu}\"  -m $tot -r $randOpt >> _out.dat";
	}
	die "\n" if (system($getAtmStr));
	$countStr = "$Bin/countAtoms.pl -f $bgffile -t bgf -m 1";
	$soluCount->{atoms} = $solu->count->{atoms};
	$soluCount->{molecule} = $solu->count->{molecule};
	open COUNT, "$countStr |" or die "ERROR while executing $countStr\n";
	while (<COUNT>) {
		if ($_ =~ /^Found (\d+) (atoms|molecule)/) {
			$count->{$2} = $1 - $soluCount->{$2};
		}
	}
	close COUNT;

	print "Done\n";
	return $count;
}
sub removeSolv_old {
	my ($structure, $cell, $count) = @_;
	my ($solvCount, $i, $tot, $mol, $molId, $start, @molIDs);
	
	$tot = $structure->count->{molecule} - $count->{molecule}; # total number of solvent molecules
	if (! exists($cell->{total}) or $tot < $cell->{total}) {
		for $i ("atoms", "molecule") {
			$solvCount->{$i} = $structure->count->{$i} - $count->{$i};
		}
		return $solvCount;
	}

	$start = $count->{molecule};
	$tot -= $cell->{total};
	for $i ($start .. $structure->count->{molecule}) {
	push @molIDs, $i;
	}
	$i = 1;
	while ($i <= $tot) {
		$molId = int(rand($#molIDs));
		print "Deleting molecule $molIDs[$molId] ($i of $tot)...\r";
		$mol = $structure->molecule($molIDs[$molId]);
		splice @molIDs, $molId, 1;
		$structure->deleteMol($mol);
		$i++;
	}
	print "Deleted $tot random solvent molecules to get to $cell->{total}..\n";
	for $i ("atoms", "molecule") {
		$solvCount->{$i} = $structure->count->{$i} - $count->{$i};
	}
	return $solvCount;
}

sub embedMols {
	my ($solute, $soluFile, $solvFile, $ffields, $rev) = @_;
	my ($embed, $countStr, $count);

	#center solute in box
	&recenter($soluFile->{tmp});
	$embed = "$Bin/embedMolecule.pl -m $solvFile -s $soluFile->{tmp} -f \"$ffields\" -c none -w $soluFile->{save} -r $rev";
	die "\n" if (system("$embed >> _out.dat"));
	system("rm -fr _out.dat _solv_1cell.bgf _solv_replicate.bgf _solv_trim.bgf _solv_trim_centered.bgf __solu.bgf");
	$countStr = "$Bin/countAtoms.pl -f $soluFile->{save} -m 1";
	open COUNT, "$countStr |" or die "ERROR while executing $countStr\n";
	while (<COUNT>) {
		if ($_ =~ /^Found (\d+) (atoms|molecule)/) {
			$count->{$2} = $1;
		}
	}
	close COUNT; 
	return ($solute->count->{molecule}, $count->{molecule}-$solute->count->{molecule});
}

sub createSolventBox {
	my ($solv, $solu, $cell, $rOpt) = @_;
	my ($i, $replicate, $blen, $trim, $bmin, $smin, $offset, $map, $solvBox, $cellScale);
	my ($tmp, $trans, $g, $rep_str);

	$map = ({ "a" => "x", "b" => "y", "c" => "z"});
	if (defined($cell->{density})) { #compress/expand the solvent cell to the new density. assume 1 g/cm3
		$cellScale = 1/($cell->{density}**(1/3)); #
		for $i ("x", "y", "z") {
			$solv->stressCell($i, $cellScale);
		}
	}
	$smin = $solv->getExtrema("min");
	if ($solv->cell->{valid}) {
		%{ $solvBox } = %{ $solv->cell };
	} else {
		for $i ("a","b","c") {
			$solvBox->{$i} = $solv->vdwbox->{$i}{max} - $solv->vdwbox->{$i}{min};
		}
	}
	if(defined($solu->cell->{a})) {
		$box->{a}{max}=$solu->cell->{a}; $box->{b}{max}=$solu->cell->{b}; $box->{c}{max}=$solu->cell->{c};
		$box->{a}{min}=$box->{b}{min}=$box->{c}{min}=0;
	} else {
		%{ $box } = %{ $solu->vdwbox };
	}
	for ("a", "b", "c") {
		$tmp = ();
		$tmp = $box->{$_}{max}-$box->{$_}{min};
		delete $box->{$_};
		$box->{$_}{max} = $tmp;
		$box->{$_}{min} = 0;
	}
	if (! defined($solu->cell->{a})) {
		($solu->cell->{a}, $solu->cell->{b}, $solu->cell->{c}) = ($box->{a}{max},$box->{b}{max},$box->{c}{max});
		$solu->cell->{alpha}=$solu->cell->{beta}=$solu->cell->{gamma}=90;
	}
	if (! defined($solu->cell->{image})) {
		$solu->cell->{image}{a}{p}=$solu->cell->{image}{b}{p}=$solu->cell->{image}{c}{p}=-1;
		$solu->cell->{image}{a}{n}=$solu->cell->{image}{b}{n}=$solu->cell->{image}{c}{n}=1;
	}

	$solv->write("_solv_1cell.bgf", "bgf");
	$replicate = "$Bin/replicate.pl -b ./_solv_1cell.bgf -d \"";
	$rep_str = "";
	for $i ("a", "b", "c") {
		$box->{$i}{max} += $cell->{cell}{$i}{max} 
			if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{max}));
		$box->{$i}{min} -= $cell->{cell}{$i}{min} 
			if (exists($cell->{cell}) and exists($cell->{cell}{$i}) and exists($cell->{cell}{$i}{min}));
		$box->{$i}{len} = $box->{$i}{max} - $box->{$i}{min};
		$bmin->{$i} = $box->{$i}{min};
		$blen .= "$box->{$i}{len} ";
		$g = int($box->{$i}{len}/$solvBox->{$i});
		$g = 1 if (! $g);
		$replicate .= sprintf("%d ",$g);
		$rep_str .= "$g"
	}
	$replicate .= "\" -s _solv_replicate.bgf";
	$replicate .= " -r 1 " if($rOpt == 1);
	#move the solvent box to the solute minima
	#for $i ("a", "b", "c") {
	#	$offset->{ $map->{$i }} = $bmin->{$i} - $smin->{$i};
	#}
	
	system("rm -fr _out.dat");
	if ($rep_str ne "111") {
		#replicate the cell by the replication vector calculated above
		die "ERROR while executing \"$replicate\"\n" if (system("${replicate} >> _out.dat"));
		# remove all molecules outside the solute (inflated) cell
		$trim = "$Bin/trimCell.pl -b _solv_replicate.bgf -c \"$blen $solv->cell->{alpha} $solv->cell->{beta} $solv->cell->{gamma}\" -s _solv_trim.bgf -m 1 -o 4";
		die "ERROR while executing \"$trim\"\n" if (system("${trim} >> _out.dat"));
	} else {
		$solv->write("_solv_trim.bgf", "bgf");
		$blen = "$solv->cell{a} $solv->cell->{b} $solv->cell->{c}";
	}
	#center solvent in box
	&recenter("_solv_trim.bgf", "_solv_trim_centered.bgf");

	#for $i (keys %{ $offset }) {
	#	$offset->{$i} = "+" . $offset->{$i} if ($offset->{$i}>=0);
	#}
	#die "ERROR: Cannot shift cell!\n" 
	#	if (system("$Bin/modifyAtomData.pl -s _solv_trim.bgf -w _solv_trim_centered.bgf -a \"index>0\" -f \"XCOORD:$offset->{x} YCOORD:$offset->{y} ZCOORD:$offset->{z}\" > _out.dat"));
	undef($solv);
}

sub recenter {
	my ($inFile, $outFile) = @_;
	my ($eStr, $center, $offset);

	$outFile = $inFile if (! defined($outFile));
	$box->{X} = $box->{a}; $box->{Y} = $box->{b}; $box->{Z} = $box->{c};
	open ECMD, "$Bin/getBounds.pl -b $inFile | " or die "ERROR: Cannot execute $Bin/getBounds.pl: $!\n";
	while (<ECMD>) {
		chomp;
		if($_ =~ /^(X|Y|Z)\s+(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)/) {
			$center->{$1} = ($2+$3)/2;
			$offset->{$1} = $box->{$1}{len}/2-$center->{$1};
			$offset->{$1} = "+$offset->{$1}" if ($offset->{$1}>=0);
		}
	}
	close ECMD;
	die "ERROR getting system bounds!"
		if (! exists($offset->{X}) or ! exists($offset->{Y}) or ! exists($offset->{Z}));
	system("$Bin/modifyAtomData.pl -s $inFile -w $outFile -a 'index>0' -f \"XCOORD:$offset->{X} YCOORD:$offset->{Y} ZCOORD:$offset->{Z}\" >> _out.dat");
	print "";
}	

sub init {
	my (%OPTS, $ffStr, $solvOptsStr, $solvTypeStr, $periodStr, $solvOpts);

	getopt('ifntwsroa', \%OPTS);
	for ("f", "i","n") {
		die &showUsage . "\n" if (! exists($OPTS{$_}));
	}

	print "Initialzing...";
	($sfile->{name}, $sfile->{type}, $ffields, $solvOptsStr, $solvTypeStr, $sfile->{save}, $randOpt,$reversePlace,$randRotOpt) = 
	($OPTS{i},       $OPTS{t},       $OPTS{f}, $OPTS{n},     $OPTS{w},     $OPTS{s},       $OPTS{r},$OPTS{o},     $OPTS{a});
	$solute =  MolData->new();
	if (! defined($sfile->{type})) {
		$sfile->{name} =~ /^\s*(\S+)/;
		$sfile->{type} = $1;
		$sfile->{type} =~ /\.(\w+)$/;
		$sfile->{type} = lc $1;
	}
	$sfile->{tmp} = "__solu.bgf";

	$cellOpts = getCellOpts($solvOptsStr);
	$solvOpts = getSolvOpts($solvTypeStr);
	if (! defined($sfile->{save})) {
		$sfile->{name} =~ /^\s*(\S+)/;
		$sfile->{save} = $solute->getFileName($1); 
	}
	$randOpt = 0 if (!defined($randOpt) or $randOpt !~ /yes|1/i);
	$randOpt = 1 if ($randOpt =~ /yes|1/i);
	print "...will randomize solvent deletions..." if ($randOpt);
	$reversePlace = 0 if (! defined($reversePlace) or $reversePlace !~ /yes|1/i);
	$reversePlace = 1 if ($reversePlace =~ /yes|1/i);
	print "...will save the structure as solvent:solute..." if ($reversePlace);
	print "Done\n";
	$solvent = MolData->new();
	print "Gettting data from solvent $solvOpts->{type} file $solvOpts->{file}...";
	$solvent->read($solvOpts->{file}, $solvOpts->{type});
	$randRotOpt = 0 if (! defined($randRotOpt) or $randRotOpt !~ /1|yes/i);
	$randRotOpt = 1 if ($randRotOpt =~ /1|yes/i);
	print "Done\n";
}

sub getCellOpts {
	my ($solventStr) = $_[0];
	my ($SOLVENT, $dim, $i, $operator, $val, $map);

	$map = ({ "x" => "a", "y" => "b", "z" => "c" });
	if ($solventStr =~ /total:\s*(\d+)/i) {
		$SOLVENT->{total} = $1;
	} elsif ($solventStr =~ /density:\s*(\d+\.?\d*)/i) {
		$SOLVENT->{density} = $1;
	}

	while ($solventStr =~ /(x|y|z|x.|y.|z.|x..|y..|z..):\s*(\+?|\-?|\+?\/?\-?|\=?)\s*(\d+\.?\d*)/gi) {
		($dim, $operator, $val) = (lc $1, "=", $2);
		($operator, $val) = ($2, $3) if ($3);
		if ($operator =~ /=/) {
			while ($dim =~ /(x|y|z)/g) {
				$SOLVENT->{cell}{ $map->{$1} }{max} = $val;
				$SOLVENT->{cell}{ $map->{$1} }{min} = 0;
			}
		} else {
			while ($operator =~ /(\+|\-)/g) {
				$i = "max" if ($1 eq "+");
				$i = "min" if ($1 eq "-");
				while ($dim =~ /(x|y|z)/g) {
					$SOLVENT->{cell}{ $map->{$1} }{$i} = $val;
				}
			}
		}
	}
	die "ERROR: Unindentified options obtained for solvent options: \"$solventStr\"\n"
	if (! $SOLVENT);
	return $SOLVENT;
}

sub getSolvOpts {
	my ($solventStr) = $_[0];
	my ($SOLVENT, $i);

	if (! defined($solventStr)) {
		$SOLVENT->{file} = "$Bin/dat/WAT/spc_box.bgf";
		$SOLVENT->{type} = "bgf";
	}elsif (-e $solventStr and -r $solventStr and -T $solventStr) {
		$SOLVENT->{file} = $solventStr;
		$SOLVENT->{type} = "bgf";
		if ($solventStr =~ /\.(\w+)$/) {
			$SOLVENT->{type} = lc $1;
		}
	} elsif ($solventStr =~ /(tip3|tip3_charmm|f3c|meso|spc|tip4)/i) {
		$i = lc $1;
		$SOLVENT->{file} = "$Bin/dat/WAT/${i}_box.bgf";
		$SOLVENT->{type} = "bgf";
	}

	if (! $SOLVENT) {
		$SOLVENT->{file} = "$Bin/dat/WAT/spc_box.bgf";
		$SOLVENT->{type} = "bgf";
	}

	return $SOLVENT;
}

sub showUsage {
	my ($fTypeStr) = GetFileTypeStr;
	my $ffTypeStr = GetFFTypeStr;
	my $usage = <<DATA;
usage: $0 -i input_structure -f "forcefield(s)" -n solvent_options -o (reverse_order_placement) -t (file type) -w (solvent type) -s (savename) -a (random_rotate_mol_opt)
Required arguments:
	-i input structure:	Location of input structure file.
$fTypeStr	-f forcefield(s): Forcefield file describing the atomic interactions
$ffTypeStr	-n solvent options:	Specifies the number of solvent molecules to add. Can be either:
		"density: x.x" - will add solvent atoms to the cell until density x.x is achieved. 
			Will assume that input solvent is equilibrated at denisty of 1 gm/cm3
		"total: x"	 - will add a total of x solvent molecules
		"[x|y|z]:[+|-|=] x.x" - will inflate|set the unit cell by x.x angstoms in specified direction 
			(=,+, -, +/- directions) for the specified dimension (any combo of x,y,z). 
			This option can be used with the previous 2. If none of the other 2 is 
			specified, will assume density 1 gm/cm3. Multiple entries should be
			enclosed in quotes. e.g. "x: +/- 10 y: =10 z: -12"
Optional arguments:
	-o reverse_order_placement: Flag for whether to reverse the order of the final structure, to be solvent:solute.
		Expected no|0 (default) or yes|1
	-t file type: Specifies the formatting of the input file. If not supplied, with default 
		to either the file suffix or to a BGF file if the file has no suffix
	-w solvent type: The following predetermined (equilibrated) solvents are available. Default F3C.
		Alternatively, you can provide the location of your solvent file.
		TIP3: the original Jorgenson TIP3 water model (rigid hydrogens)
		TIP4: TIP4P with massless pseudo-atom
		TIP3_CHARMM: TIP3 water model as implemented in CHARMM
		F3C: F3C water model (no rigid hydrogens)
		SPC: SPC water model (DEFAULT)
		MESO: Water model for Meso-scale simulations (Valeria Molinero)
	-r randomize: When removing molecules, randomize selection. Expected no|0 (default) or yes|1			
	-s savename: Will assume \$prefix_mod.\$suffix if not specified.
	-a rotate_mol_opt: Will randomely rotate the solvent molecule. Expected no|0 (default) or yes|1
DATA

	return $usage;
}
