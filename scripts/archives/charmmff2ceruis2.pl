#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::CERIUS2 qw(saveCeriusFF);
use Packages::General qw(FileTester LoadElements);
use File::Basename;
use strict;

sub init;
sub parseCharmmFF;
sub parseCharmmTop;
sub getFields;
sub getElementNum;
sub pruneAtomTypes;

die "usage: $0 charmmTop charmmPar [saveName]\n"
    if (! @ARGV || $#ARGV < 1);

my ($charmmTop, $charmmFF, $cerius2FF) = @ARGV;
my ($PARMS, $ELEMENTS);

$|++;
&init;
print "Parsing CHARMM force field $charmmFF...";
$PARMS = parseCharmmTop($charmmTop, $ELEMENTS);
$PARMS = parseCharmmFF($PARMS, $charmmFF);
pruneAtomTypes($PARMS);
print "Done\nCreating CERIUS2 force field $cerius2FF...";
saveCeriusFF($PARMS, $cerius2FF, $ELEMENTS);
print "Done\n";

sub pruneAtomTypes {
    my ($parms) = $_[0];
    my ($i);

    for $i (keys %{ $parms->{atoms} }) {
	delete $parms->{atoms}{$i} if (! exists($parms->{atoms}{$i}{VALS}));
    }
}

sub getElementNum {
    my ($elementName, $elements) = @_;
    my ($elementNum, $i);

    $elementName = uc($elementName);
    for $i (keys %{ $elements } ) {
	if (uc($elements->{$i}{SYMBOL}) eq $elementName) {
	    $elementNum = $i;
	    last;
	}
    }

    die "ERROR: Element not found in element table: $elementName!\n" if (! $elementNum);

    return $elementNum;

}

sub parseCharmmTop {
    my ($topFile, $elements) = @_;
    my (%DATA);
    
    open CHARMMTOP, $topFile || die "ERROR: Cannot open CHARMM topology file $topFile: $!\n";
    while (<CHARMMTOP>) {
	chomp;
	if ($_ =~ /^MASS\s+(\d+)\s+(\w+)\s+(\d+\.\d+)\s+(\w+)\s+\S\s+(.+)/) {
	    $DATA{atoms}{$2} = (
				{
				    "index"   => $1,
				    "mass"    => $3,
				    "element" => getElementNum($4, $elements),
				    "hybrid"  => 0,
				    "name"    => $5,
				}
				);
	} elsif ($_ =~ /^ATOM\s+(\w+)\s+(\w+)\s+(\-?\d+\.\d+)/) {
	    $DATA{atoms}{$1}{CHARGE} = $3 if (exists($DATA{atoms}{$1}));
	    $DATA{atoms}{$2}{CHARGE} = $3 if (exists($DATA{atoms}{$2}));
	}
    }
    close CHARMMTOP;

    die "ERROR: No valid data found while parsing $topFile!\n" if (! keys %DATA);

    return \%DATA;
}

sub parseCharmmFF {
    my ($PARMS, $ff) = @_;
    my (%DATA, %VAL, $counter, $header, $FIELDS, @atms); 
    my (@vals, $curr, $inStr, $i, $rec, @tmp, %TYPES, $j);
    my ($pi) = atan2(1,1) *4;

    %VAL = (
	    "BONDS"     => 2,
	    "ANGLES"    => 3,
	    "DIHEDRALS" => 4,
	    "IMPROPER"  => 4,
	    );
    %TYPES = (
	      "BONDS"      => "HARMONIC",
	      "ANGLES"     => "THETA_HARM",
	      "TORSIONS"   => "SHFT_DIHDR",
	      "INVERSIONS" => "IT_IJKL",
	      );

    %DATA = %{ $PARMS };
    open CHARMMPAR, $ff || die "ERROR: Cannot open CHARMM force field $ff: $!\n";
    while (<CHARMMPAR>) {
	chomp;
	next if ($_ !~ /^\w+/);
	if ($_ =~ /^(\w+)$/) {
	    $counter = 0;
	    if (exists($VAL{$1})) {
		$counter = $VAL{$1};
		$header = lc($1);
		$header = "torsions" if ($header eq "dihedrals");
		$header = "inversions" if ($header eq "improper");
	    }
	} elsif ($_ =~ /NONBONDED nbxmod/) {
	    $counter = 1;
	    $header = "atoms";
	} elsif ($counter > 0 && $_ =~ /^(\w+)(.+)\!*/) {
	    @atms = ();
	    @vals = ();
	    @tmp = split /\s+/, "$1 $2";
	    for $i (0 .. $#tmp) {
		if ($tmp[$i] !~ /([A-Z_a-z]+\w*)/) {
		    $j = $i;
		    last;
		}
		push @atms, $1 if (exists($DATA{atoms}{$1}) || $1 eq "X");
	    }

	    for $i ($j .. $#tmp) {
		last if ($tmp[$i] !~ /(\-?\d+\.*\d*)/);
		push @vals, $1;
	    }

	    if (($#atms + 1) == $counter) {
		@tmp = keys %{ $DATA{$header} };
		
		$FIELDS = getFileds(\@atms);
		$curr = \%{ $DATA{$header}{$atms[0]} };
		for $i (1 .. $#atms) {
		    $curr = \%{ $curr->{$atms[$i]} };
		}
		$curr->{counter} = ($#tmp + 2);
		$rec = ();
		if ($counter == 1) {
		    $vals[1] *= -1;
		    $rec->{mass} = $DATA{atoms}{$1}{mass};
		    $rec->{hybrid} = 0;
		    $rec->{element} = $DATA{atoms}{$1}{element};
		    $rec->{name} = $DATA{atoms}{$1}{name};
		    if ($#vals > 2) {
			$rec->{"1_4"}{e} = $vals[4];
			$rec->{"1_4"}{r} = $vals[5];
		    } else {
			$rec->{"1_4"}{e} = $vals[1];
			$rec->{"1_4"}{r} = $vals[2];
		    }
		}
		if ($counter == 4) {
		    if ($header eq "torsions") {
			$DATA{"torsionOrders"}{$curr->{counter}} = "0 1 2 3";
		    } else {
			$DATA{"inversionOrders"}{$curr->{counter}} = "3 1 2 0";
		    }
		    $rec->{type} = $TYPES{uc($header)};
		}
		for $i (0 .. $#{ $FIELDS }) {
		    if ($counter > 3 && $i == 1) {
			$vals[$#vals] *= ($pi/180);
		    } elsif ($counter == 3 && $i == 1) {
			$vals[1] *= ($pi/180);
		    }
		    $rec->{$FIELDS->[$i]} = $vals[$i];
		}
		push @{ $curr->{VALS} }, $rec;

		# UREY_BRADLEY
		if ($counter == 3 && $#vals == 3) {
		    @tmp = keys %{ $DATA{urey_bradley} };
		    $curr = \%{ $DATA{urey_bradley}{$atms[0]} };
		    for $i (1 .. $#atms) {
			$curr = \%{ $curr->{$atms[$i]} };
		    }
		    $curr->{counter} = ($#tmp + 2);
		    $rec = ();
		    $rec = (
			    {
				"ku" => $vals[2],
				"su" => $vals[3],
			    }
			    );
		    push @{ $curr->{VALS} }, $rec;
		}
	    }
	}
	
    }		
    close CHARMMPAR;

    return \%DATA;
	
}

sub getFileds {
    my ($atms) = $_[0];
    my ($i, @FIELDS);
    
    if ($#{ $atms } == 0) {
	@FIELDS = ("junk","e","r");
    } elsif ($#{ $atms } == 1) {
	@FIELDS = ("kb","r0");
    } elsif ($#{ $atms } == 2) {
	@FIELDS = ("kt","t0");
    } else {
	@FIELDS = ("kp","n","p0");
    }

    return \@FIELDS;
}

sub init {
    print "Initializing...";
    FileTester($charmmFF);
    FileTester($charmmTop);
    if (! defined($cerius2FF)) {
	$cerius2FF = basename($charmmFF);
	$cerius2FF =~ s/\.\w+$/\.ff/;
    }
    $ELEMENTS = &LoadElements;
    print "Done\n";

}
