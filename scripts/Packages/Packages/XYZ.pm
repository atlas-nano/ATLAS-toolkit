package XYZ;

require Exporter;
use strict;
use General qw(PrintProgress);
use constant q2e => 18.2223;
use Math::Trig qw(pi);

use strict;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(CreateXYZTrj);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub CreateXYZTrj {
    my ($ATOMS, $tstep, $OUTFILE) = @_;
    my ($atomC, $totC, @tmp, $dim);

    $totC = scalar(keys %{ $ATOMS  });
    @tmp = sort numerically keys %{ $ATOMS };
    print $OUTFILE "$totC\n";
    print $OUTFILE " i = $tstep\n";
    for $atomC (@tmp) {
	    printf $OUTFILE "%-2s ",$ATOMS->{$atomC}{ELEMENT};
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    printf $OUTFILE "%8.3f", $ATOMS->{$atomC}{$dim};
	}
	printf $OUTFILE "\n";
    }
}

1;	      
