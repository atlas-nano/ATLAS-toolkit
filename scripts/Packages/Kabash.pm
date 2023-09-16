package Kabash;
use strict;
use Math::MatrixReal;
use Math::GSL::Linalg::SVD;

use PDL;
use PDL::NiceSlice;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(SuperimposeAtomsPDL SuperimposeAtomsMR);

sub SuperimposeAtomsMR {
	my ($ref, $mol, $mol2) = @_;
	my ($parray, $qarray, $tarray, $tindx);
	my ($A, $B, $C, $R);
	my ($i, $j, $transform, $com_A);

   	$mol2 = $mol if (! defined($mol2));
	($parray, $qarray, $tarray, $tindx) = createKabashArrays($ref,$mol, $mol2);

	$A = Math::MatrixReal->new_from_rows([@{ $parray }]);
	$B = Math::MatrixReal->new_from_rows([@{ $qarray }]);
	$C = Math::MatrixReal->new_from_rows([@{ $tarray }]);

	#first remove center of mass displacement
	$com_A = col_sum($A);
	$A -= $com_A;
	$B -= col_sum($B);
	$C -= col_sum($C); 

	#now get rotation matrix
	$R = kabschMR($B, $A);

	#now apply best rotation
	$C *= $R;

	#finally place matched coordinates in original reference framework
	$C += $com_A;

	#and create the updated atom coordinate array
	for $i (0 .. $#{ $tarray }) {
		for $j (0 .. $#{ $tarray->[$i]}) {
			$transform->[$i][$j] = $C->element($i+1,$j+1);
		}
	}
	&updateCoords($mol2, $transform, $tindx);
}

sub SuperimposeAtomsPDL {
	my ($ref, $mol, $mol2) = @_;
	my ($A, $B, $C, $R);
	my ($parray, $qarray, $tarray, $tindx);
	my ($i, $j, $transform, $com_a);

    $mol2 = $mol if (! defined($mol2));
	($parray, $qarray, $tarray, $tindx) = createKabashArrays($ref,$mol, $mol2);
	
	$A = pdl(float, @{ $parray });
	$B = pdl(float, @{ $qarray });
	$C = pdl(float, @{ $tarray });

	#first remove center of mass displacement
	$com_a = centroid($A);
	$A -= $com_a;
	$B -= centroid($B);
	$C -= centroid($C);

	#now get rotation matrix
	$R = kabsch($B, $A);

	#now apply best rotation
	$C = $C x $R;

	#finally place matched coordinates in original reference framework
	$C += $com_a;

	#and create the updated atom coordinate array
	for $i (0 .. $#{ $tarray }) {
		for $j (0 .. $#{ $tarray->[$i]}) {
			$transform->[$i][$j] = $C->at($j,$i);
		}
	}
	&updateCoords($mol2, $transform, $tindx);
}

sub createKabashArrays {
	my ($ref,$mol, $mol2) = @_;
	my ($i, $parray, $qarray, $tarray, $tindx);

	for $i (keys %{ $ref}) {
		push @{ $parray }, ([$ref->{$i}{XCOORD}, $ref->{$i}{YCOORD}, $ref->{$i}{ZCOORD}]);
		push @{ $qarray }, ([$mol->{$i}{XCOORD}, $mol->{$i}{YCOORD}, $mol->{$i}{ZCOORD}]);
	}

	for $i (keys %{ $mol2 }) {
		push @{ $tarray }, ([$mol2->{$i}{XCOORD}, $mol2->{$i}{YCOORD}, $mol2->{$i}{ZCOORD}]);
		$tindx->{$i} = $#{$tarray};
	}

	return ($parray, $qarray, $tarray, $tindx);
}

sub updateCoords{
	my ($aList,$new_coords,$tindx) = @_;
	my ($i,$j);

	for $i (keys %{ $aList }) {
		$j = $tindx->{$i}; #get the array entry
		$aList->{$i}{XCOORD} = $new_coords->[$j][0];
		$aList->{$i}{YCOORD} = $new_coords->[$j][1];
		$aList->{$i}{ZCOORD} = $new_coords->[$j][2];
	}

}

sub rmsd {

	my ($u, $v) = @_;

    my ($dims, $n) = $u->dims;

    my $x = pdl(float, []);

    for (0..$dims-1) {
        $x = $x->append( ($u(($_),:) - $v(($_),:)) ** 2 );
    }
    return sqrt($x->sumover/$n);
}

sub centroid {
    my $pts = shift;

    my ($dims, $n) = $pts->dims;

    return $pts->xchg(0,1)->sumover/$n;
}

sub col_sum {
	my ($mat) = $_[0];
	my ($col_sum, $i, $j, $tmp, $ret);
	my ($rows,$cols) = $mat->dim();

	for $i (1 .. $cols) {
		$col_sum = 0;
		for $j (1 .. $rows) {
			$col_sum += $mat->element($j,$i);
		}
		$col_sum /= $rows;
		push @{ $tmp }, $col_sum;
	}

	$ret = Math::MatrixReal->new($rows, $cols);
	for $j (1 .. $rows) {
		for $i (1 .. $cols) {
			$ret->assign($j,$i,$tmp->[$i-1]);
		}
	}
	return $ret;
}

# https://en.wikipedia.org/wiki/Kabsch_algorithm
sub kabsch {
    my ($P, $Q) = @_;

    my $n = $P->dim(0);

    # Compute the covariance matrix
    my $H = ($Q->transpose) x $P;

    my ($U, $S, $V) = $H->svd;
    my $d = det($V x $U->transpose) < 0 ? -1 : 1;

    my $I = identity($n);
    $I($n-1,$n-1) .= $d;

    return $V x $I x $U->transpose;
}

sub kabschMR {
    my ($P, $Q) = @_;

    my ($n, $m) = $P->dim();

    # Compute the covariance matrix
    my $H = ~$Q * $P;

	#SVD 
	my $svd = Math::GSL::Linalg::SVD->new( { verbose=>0 });
	$svd->load_data( { data=> $H->[0] });
	$svd->decompose( { algorithm => q{gd} } );

    my ($S, $U, $V, $O) = $svd->results;
	my ($Uarray) = Math::MatrixReal->new_from_rows([@{ $U }]);
	my ($Varray) = Math::MatrixReal->new_from_rows([@{ $V }]);
	my ($Sarray) = Math::MatrixReal->new_from_rows([$S]);

	my $tmp = $Varray * ~$Uarray;
	my $d = $tmp->det < 0 ? -1 : 1;
	$tmp = ();
	for my $i (1 .. $m) {
		push @{ $tmp }, 1;
	}
    my $I = Math::MatrixReal->new_diag([@{ $tmp }]);
	$I->assign($m,$m,$d);

    return $Varray * $I * ~$Uarray;
}

;

