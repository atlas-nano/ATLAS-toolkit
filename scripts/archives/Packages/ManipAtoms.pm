package Packages::ManipAtoms;
require Exporter;
use strict;
use Cwd;

our (@ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
$VERSION = "1.00";
@EXPORT = qw(ImageAtoms ScaleAtoms UnwrapAtoms FindElement CenterSystem GetAtmList);

sub CenterSystem {
    my ($atoms, $offset, $atomList) = @_;
    my ($i, $dim);

    for $i (keys %{ $atoms }) {
	next if (keys %{ $atomList } && ! exists($atomList->{$i}));
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $offset->{$dim} = $offset->{substr($dim,0,1)}{lo} if (! exists($offset->{$dim}));
	    next if (! $offset->{$dim});
	    $atoms->{$i}{$dim} -= $offset->{$dim};
	}
    }
}

sub FindElement {
    my ($fftype, $parms, $elements) = @_;
    my ($elementName, $elementNum, $i);

    die "ERROR: $fftype is not a valid force field type!\n" if (! exists($parms->{$fftype}));
    $elementName = $parms->{$fftype}{ATOM};
    $elementNum = 0;

    for $i (keys %{ $elements }) {
        if (uc($elements->{$i}{SYMBOL}) eq uc($elementName)) {
            $elementNum = $i;
            last;
        }
    }

    #die "ERROR: Element $elementName is not a valid element!\n" if (! $elementNum);

    return ($elementName, $elementNum);
}

sub ImageAtoms {
    my ($atoms, $moleculeAtoms, $CENTER, $box) = @_;
    my ($MOL, $i, $OFFSET, $dim, $index);

    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$OFFSET->{$dim} = sprintf("%.0f", ($CENTER->{$dim}/$box->{$dim}{len}));
	#$OFFSET->{$dim} = int($CENTER->{$dim}/$box->{$dim}{len});
    }
    
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$index = $dim;
	$index =~ s/COORD/INDEX/;
	for $i (keys %{ $moleculeAtoms }) {
	    $atoms->{$i}{$dim} -= ($OFFSET->{$dim} * $box->{$dim}{len});
	    $atoms->{$i}{$dim} -= $box->{$dim}{lo};
	    $atoms->{$i}{$index} = $OFFSET->{$dim};
	    $atoms->{$i}{IMAGED} = 1;
	}
    }
	    
}

sub ScaleAtoms {
    my ($atoms, $box) = @_;
    my ($i, $coord, $dim, $index, $pos);

    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$index = $dim;
	$index =~ s/COORD/INDEX/;
	for $i (keys %{ $atoms }) {
	    $pos = $atoms->{$i}{$dim};
	    $pos /= $box->{$dim}{len};
	    if ($pos =~ /(\d+)(\.\d+)/) {
		if ($1 > 0) {
		    $atoms->{$i}{$index} = $1 + 1;
		} else {
		    $atoms->{$i}{$index} = $1;
		}
		$atoms->{$i}{$dim} = sprintf("%.6f ", $pos);
	    }
	}
    }
}

sub UnwrapAtoms {
    my ($ATOMS, $BOX, $scaled) = @_;
    my ($atomC, $atom, $coord, $index, $pos);

    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	for $coord ("XCOORD", "YCOORD", "ZCOORD") {
	    $index = $coord;
	    $index =~ s/COORD/INDEX/;
	    $index = $atom->{$index};
	    if (! defined($scaled) or $scaled == 1) {
	        $atom->{$coord} *= $BOX->{$coord}{len};
	        $atom->{$coord} += $BOX->{$coord}{lo};
	    }
	    $atom->{$coord} += ($index * $BOX->{$coord}{len});
	}
    }
}

sub GetAtmList {
    my ($select, $ATOMS) = @_;
    my ($i, %LIST, $field, $val);

    for $i (keys %{ $ATOMS }) {
	for $field (keys %{ $select }) {
	    for $val (keys %{ $select->{$field} }) {
		$LIST{$i} = $val if ($val eq "*" or (exists($ATOMS->{$i}{$field}) and $ATOMS->{$i}{$field} eq $val) or 
				    ($field eq "INDEX" and $i == $val));
	    }
	}
    }
    
    return \%LIST;
}

1;
