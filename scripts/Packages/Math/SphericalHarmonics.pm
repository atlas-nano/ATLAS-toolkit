package Math::SphericalHarmonics;

require 5.005;

use Math::Complex;
use Carp;
use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use strict;
use constant PI => atan2(1,1) * 4;

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(plgndr ylm);
$VERSION = '0.1';

# Computes the spherical harmonic from the associated Legendre Polynomials
# use relation that
#	Ylm(theta,phi) = sqrt((2l+1*(l-m)!)/4pi/(l+m)! Plm(cos(theta))e^(imphi)

sub ylm {
	my ($l, $m, $theta, $phi) = @_;
	my ($pfactor, $imfactor, $plm, $val);

	$pfactor = $imfactor = 1;
	#first check if m < 0
	if ($m < 0) {
		# using relation Yl-m(theta,phi)=(-1)^m*complex_conjugate(Ylm(theta,phi)
		$m *= -1;
		$pfactor = (-1)**$m;
		$imfactor = -1;
	}
	$pfactor *= sqrt((2*$l+1)*factorial($l-$m)/4/PI/factorial($l+$m));
	$plm = plgndr($l,$m,cos($theta));
	$val = ();
	$val->{re} = $pfactor*$plm*cos($m*$phi);
	$val->{im} = $pfactor*$plm*$imfactor*sin($m*$phi);
	return $val;

}

#
# Computes the associated Legendre polynomial Pml(x). Here m and l are integers satisfying
# 0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1
# From numerical receipes in C
#
sub plgndr
{
	my ($l, $m, $x) = @_;
	my ($fact, $pll, $pmm, $pmmp1, $somx2);
	my ($i, $ll);

	die "Bad arguments passed to plgndr\n" if ($m < 0 or $m > $l or abs($x) > 1.0);

	$pmm = 1; # Compute Pmm
	if ($m > 0) {
		$somx2=sqrt((1-$x)*(1+$x));
		$fact=1;
		for $i (1 .. $m) {
			$pmm *= -$fact*$somx2;
			$fact += 2;
		}
	}
	return $pmm if ($l == $m);

	#Compute Pmm+1
	$pmmp1=$x*(2*$m+1)*$pmm;
	return $pmmp1 if($l==($m+1));

	for $ll ($m+2 .. $l) {
		$pll = ($x*(2*$ll-1)*$pmmp1-($ll+$m-1)*$pmm)/($ll-$m);
		$pmm=$pmmp1;
		$pmmp1=$pll;
	}
	return $pll;
}

# Recursive factorial function
sub factorial {
	my $arg = shift;

	return (($arg == 1) || ($arg == 0)) ? 1 : $arg * factorial($arg - 1);
}

1;
__END__

