#!/usr/bin/perl -w
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use General qw(ReadFFs GetSelections FileTester LoadFFs GetBondLength);
use Math::ematica qw(:PACKET :TYPE :FUNC);
use ManipAtoms qw(GetAtmList);

sub init;
sub MathLinkCmd;
sub parseEngFile;
sub createNBList;
sub pairMix;
sub numerically { ($a<=>$b); }
sub setVDWType;
sub getLine;
sub ShowError;

my ($selection, $BGFS, $engFile, $saveFile, $ffData);
my ($ATOMS, $SELECT, $DATA, $NLIST, $FF, $errFile, $MODEL);

$|++;
&init;
print "Parsing BGF File $BGFS->[0]...";
($ATOMS, undef, undef) = GetBGFFileInfo($BGFS->[0], 1);
print "Done\n";
$FF = LoadFFs($ffData, "", 0);
print "Parsing Energy file $engFile...";
$DATA = parseEngFile($engFile, $BGFS);
print "Done\nParsing atom/residue selection...";
$SELECT = GetSelections($selection, 0);
$SELECT = GetAtmList($SELECT, $ATOMS);
print "Done\nCreating linear system...";
open ERRFILE, "> $errFile" or die "ERROR: Cannot create error file $errFile: $!\n";
($NLIST, $MODEL) = createNBList($ATOMS, $SELECT, $FF);
&createLinearSystem($DATA, $NLIST, $MODEL);

close ERRFILE;
print "Done\n";

sub createLinearSystem {
    my ($data, $nlist, $model) = @_;
    my ($MathLinkResult, $holder, $link, $error);
    my ($i, $mathCmd, $count, $j, $k, $l);
    my ($atom1, $atom2, $fitCmd);

    $link = new Math::ematica '-linklaunch', '-linkname', 'math -mathlink';

    $count = 0;
    $fitCmd = "FindRoot[{";
    for $i (sort numerically keys %{ $data }) {
	$count++;
	($ATOMS, undef) = GetBGFFileInfo($data->{$i}{FILE}, 0);
	$mathCmd = "eqn${count} = -1*$data->{$i}{ENG} ";
	$fitCmd .= "eqn${count},";
	for $j (keys %{ $nlist }) {
	    for $k (keys %{ $nlist->{$j}} ) {
		for $l (@{ $nlist->{$j}{$k}{ATOMS} }) {
		    ($atom1, $atom2) = ($ATOMS->{ $l->{ATOM1} }, $ATOMS->{ $l->{ATOM2} });
		    $mathCmd .= getLine($atom1, $atom2, $model);
		}
	    }
	}
	$mathCmd .= ";";
	MathLinkCmd($mathCmd, 1, $link);
    }
    chop $fitCmd;
    $fitCmd .= "},{";
    for $i (keys %{ $model->{PARMS} }) {
	$fitCmd .= "{$i, $model->{PARMS}{$i}},";
    }
    chop $fitCmd;
    $fitCmd .= "}]";
    ($error, $holder) = MathLinkCmd($fitCmd, 0, $link);
    &ShowError($error, $link);
    print "";

}

sub getLine {
    my ($atom1, $atom2, $model) = @_;
    my ($dist, $vdwLine, $t1, $t2);

    ($t1, $t2) = ($atom1->{FFTYPE}, $atom2->{FFTYPE});
    ($t1, $t2) = ($t2, $t1) if ($t1 lt $t2);
    $dist = GetBondLength($atom1, $atom2);
    if ($model->{TYPE} =~ /EXPO_6/) {
	$vdwLine = "(" . $model->{"$t1 $t2 d"} . "/(" . $model->{"$t1 $t2 s"} . " - 6))*(6*Exp[" .
	    $model->{"$t1 $t2 s"} . "(1-(${dist}/" . $model->{"$t1 $t2 r"} . "))] - " .
	    $model->{"$t1 $t2 s"} . "*(${dist}/" . $model->{"$t1 $t2 r"} . ")^6)";
    } elsif ($model->{TYPE} =~ /LJ_12_6/) {
	$vdwLine = $model->{"${t1} ${t2} d"} . "*((${dist}/" . $model->{"${t1} ${t2} r"} . ")^(-12) " .
	    "-2*((${dist}/" . $model->{"${t1} ${t2} r"} . ")^(-6)))";
    }
    
    return " + $vdwLine";
}		    

sub createNBList {
    my ($atoms, $mol1, $ff) = @_;
    my ($mol2, $t1, $t2, $i, $j, %NBLIST, $curr, $rec, %TYPES);

    for $i (keys %{ $atoms }) {
	$mol2->{$i} = $atoms->{$i}{RESNUM} if (! exists($mol1->{$i}));
    }

    for $i (keys %{ $mol1 }) {
	$t1 = $atoms->{$i}{FFTYPE};
	die "ERROR: VDW for $t1 not found in forcefield!\n"
	    if (! exists($ff->{VDW}{$t1}));
	for $j (keys %{ $mol2 }) {
	    $t2 = $atoms->{$j}{FFTYPE};
	    $curr = ();
	    die "ERROR: VDW for $t2 not found in forcefield!\n"
		if (! exists($ff->{VDW}{$t2}));
	    die "ERROR: Incompatible VDW types between $t1 " . $ff->{VDW}{$t1}{$t1}{1}{TYPE} . 
		" and $t2 " . $ff->{VDW}{$t2}{$t2}{1}{TYPE} . "...Aborting\n"
		if ($ff->{VDW}{$t1}{$t1}{1}{TYPE} ne $ff->{VDW}{$t2}{$t2}{1}{TYPE});
	    $TYPES{ $ff->{VDW}{$t1}{$t1}{1}{TYPE} } = 1;
	    $TYPES{ $ff->{VDW}{$t2}{$t2}{1}{TYPE} } = 1;
	    die "ERROR: More than 1 VDW type recorded!\n" if (scalar(keys %TYPES) == 2);
	    if ($t1 gt $t2) {
		$curr = \%{ $NBLIST{$t1}{$t2} };
	    } else {
		$curr = \%{ $NBLIST{$t2}{$t1} };
	    }
	    if (exists($ff->{VDW}{$t1}{$t2})) {
		$curr->{FFDATA} = $ff->{VDW}{$t1}{$t2}{1};
	    } elsif (exists($ff->{VDW}{$t2}{$t1})) {
		$curr->{FFDATA} = $ff->{VDW}{$t2}{$t1}{1};
	    } else {
		$curr->{FFDATA} = pairMix($ff->{VDW}{$t1}{$t1}{1}, $ff->{VDW}{$t2}{$t2}{1});
	    }
	    $rec = (
		    {
			"ATOM1" => $i,
			"ATOM2" => $j,
		    }
		    );
	    push @{ $curr->{ATOMS} }, $rec;
	}
    }
    &setVDWType(\%TYPES,\%NBLIST);

    return (\%NBLIST, \%TYPES);
}

sub setVDWType {
    my ($types, $nlist) = @_;
    my ($i, $mType, $totParms, $offset, $j);

    $totParms = scalar(keys %{ $nlist });
    $totParms += $totParms/2;
    for $i (keys %{ $types }) {
	$mType = uc($i);
    }
    
    $offset = 1;
    $types->{TYPE} = $mType;
    for $i (keys %{ $nlist }) {
	for $j (keys %{ $nlist->{$i} }) {
	    if ($mType =~ /EXPO_6|MORSE/) {
		$types->{"${i} ${j} r"} = lc(chr($offset + 64));
		$types->{"${i} ${j} d"} = lc(chr($offset + 65));
		$types->{"${i} ${j} s"} = lc(chr($offset + 66));
		$types->{PARMS}{lc(chr($offset + 64))} = $nlist->{$i}{$j}{FFDATA}{VALS}[0];
		$types->{PARMS}{lc(chr($offset + 65))} = $nlist->{$i}{$j}{FFDATA}{VALS}[1];
		$types->{PARMS}{lc(chr($offset + 66))} = $nlist->{$i}{$j}{FFDATA}{VALS}[2];
		$offset += 3;
	    } elsif ($mType =~ /LJ/) {
		$types->{"${i} ${j} r"} = lc(chr($offset + 64));
		$types->{"${i} ${j} d"} = lc(chr($offset + 65));
		$types->{PARMS}{lc(chr($offset + 64))} = $nlist->{$i}{$j}{FFDATA}{VALS}[0];
		$types->{PARMS}{lc(chr($offset + 65))} = $nlist->{$i}{$j}{FFDATA}{VALS}[1];
		$offset += 2;
	    }
	}
    }
}

sub pairMix {
    my ($t1, $t2) = @_;
    my (%VDW);

    $VDW{TYPE} = $t1->{TYPE};
    $VDW{VALS}[0] = 0.5*($t1->{VALS}[0] + $t2->{VALS}[0]);
    $VDW{VALS}[1] = sqrt($t1->{VALS}[1] * $t2->{VALS}[1]);
    $VDW{VALS}[2] = sqrt($t1->{VALS}[2] * $t2->{VALS}[2]) if ($#{ $t1->{VALS} } > 1);

    return \%VDW;
}

sub parseEngFile {
    my ($inFile, $bgfs) = @_;
    my (%ENERGIES, $i, $j, $eng, $isValid);

    open ENGFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    while (<ENGFILE>) {
	chomp;
	if ($_ =~ /^\s*(\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	    $ENERGIES{$1}{ENG} = $2;
	}
    }
    close ENGFILE;
    die "ERROR: $inFile is not a valid energy file!\n"
	if (! %ENERGIES);

    for $i (keys %ENERGIES) {
	$eng = $ENERGIES{$i}{ENG};
	$isValid = 0;
	for $j (0 .. $#{ $bgfs }) {
            if ($bgfs->[$j] =~ /$i/) {
                $ENERGIES{$i}{FILE} = $bgfs->[$j];
                splice @{ $bgfs }, $j, 1;
                $isValid = 1;
                last;
            }
        }
	die "ERROR: No BGF file found for data point $i. Aborting\n"
	    if (! $isValid);
    }

    return \%ENERGIES;
}

sub init {
    my (%OPTS, $findCmd, $select);
    
    getopt('besmf', \%OPTS);
    for ("b", "e", "m", "f") {
	die "usage: $0 -b bgf files (location) -m mol1 atoms -e energy file -f force field(s) -s (save file)\n"
	    if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    $findCmd = "find $OPTS{b} -name '*.bgf' -print";
    if (-e $OPTS{b}) {
	push @{ $BGFS }, $OPTS{b};
    } elsif (open(FINDCMD, "$findCmd |")) {
	while (<FINDCMD>) {
	    chomp;
	    push @{ $BGFS }, $_;
	}
	close FINDCMD;
    } else {
	die "ERROR: No valid BGF files found while searching $OPTS{b}\n";
    }
    ($engFile, $saveFile, $select, $ffData) = ($OPTS{e}, $OPTS{s}, $OPTS{m}, $OPTS{f});
    if ($select =~ /\s+/) {
        @{ $selection } = split /\s+/, $select;
    } else {
        $selection->[0] = $select;
    }
    
    FileTester($engFile);
    if (! defined($saveFile)) {
	$saveFile = basename($engFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_opt_ff.dat";
    }
    $errFile = $saveFile;
    $errFile =~ s/\.\w+$//;
    $errFile .= "_errors.dat";

    ($ffData, undef) = ReadFFs($ffData);
    print "Done\n";
}

sub MathLinkCmd(@) {

    my ($inExpression, $isTerminal, $link) = @_;
    my ($tmp, $err, @holder);

    print ERRFILE "$inExpression\n";
    $link->PutFunction("EvaluatePacket",1);
    if (! $isTerminal) {
        $link->PutFunction("ToString", 1);
    }
    $link->PutFunction("ToExpression", 1);
    $link->PutString("$inExpression");
    $link->EndPacket;

    while ($link->NextPacket != RETURNPKT) {
        $link->NewPacket;
    }
    if ($isTerminal == 1) {
        $link->NewPacket;
        if (! $link) {
            print "ERROR: " . $link->ErrorMessage . "\n";
        }
    } else {
        $err = $link->ErrorMessage;
        if ($isTerminal == 2) {
            @holder = $link->GetRealList;
        } else  {
            $tmp = $link->GetString;
        }
        if ($err =~ /ok so far/) {
            $err = "";
        }
        if ($isTerminal == 2) {
            return ($err, @holder);
        } else {
            return ($err, $tmp);
        }
    }
}

sub ShowError {
    my ($err, $link) = @_;
    if ($err) {
        $link->PutFunction("Exit", 0);
        $link->EndPacket;
        $link = ();
        die "ERROR: $err\n";
    }
}
