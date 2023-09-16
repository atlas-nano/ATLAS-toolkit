#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(addHeader createBGF insertHeaderRemark createHeaders);
use General qw(GetBondLength);
use Getopt::Std qw(getopt);

use constant PI => atan2(1,1) * 4;

sub usage;

my ($m, $n, $nwall, $saveName);
my ($ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Generating $n->[0] $m->[0] CNT...";
($ATOMS, $BONDS, $HEADERS) = createCNT($m, $n, $nwall);
print "Done\nCreating BGF File ${saveName}...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub createCNT {
	my ($m, $n, $nwall) = @_;
	my ($l, $il, $r, $ir, $Lz, $iLz, $id, $sid, $a0, $eps, $j, $i);
	my ($headers, $atoms, $latoms, $bonds, $lbonds, $box);

	$a0 = 2.46; $eps = 1e-08;
	%{ $atoms } = %{ $bonds } = ();
	$id = $r = $Lz = $l = 0;
	for $j (0 .. $nwall-1) {
		($latoms,$sid, $ir, $iLz, $il) = createAtomsXYZ($a0, $eps, $j, $id);
		$lbonds = createBonds($latoms, $a0, $iLz);
		&addMissingFields($latoms, $lbonds, $j);
		for $i (keys %{ $latoms }) {
			%{ $atoms->{$i} } = %{ $latoms->{$i} };
			$bonds->{$i} = ();
			@{ $bonds->{$i} } = @{ $lbonds->{$i} } if (defined($lbonds) and exists($lbonds->{$i}));
		}
		$id = $sid-1;
		$r = $ir if ($ir>$r);
		$Lz = $iLz if ($iLz>$Lz);
		$l = $il if ($il>$l);
	}

	#($r, $Lz, $l) = getLengthAndRadii($a0);
	$box = createBox($r, $Lz);
	$headers = newHeaders($a0, $l, $Lz, $box);

	return ($atoms, $bonds, $headers);
}

sub getLengthAndRadii {
	my ($a0) = $_[0];
	my ($l, $r, $v, $w, $s, $fac, $Lz);

	$l=sqrt($m->[0]**2+$n->[0]**2+$m->[0]*$n->[0]);
	$r=$l*$a0/PI+5;
	$v=$n->[0]+2*$m->[0];
	$w=2*$n->[0]+$m->[0];
	$s->[1]=$w/$l/2;                #don't write "w/2/l", w/2 is... 
	$s->[2]=$v/$l/2;                #...int,we want double
	$s->[3]=-sqrt(3)*$m->[0]/$l/2;
	$s->[4]=sqrt(3)*$n->[0]/$l/2;
	$fac=getFactor($v,$w);
	$v/=-$fac;
	$w/=$fac;
	$Lz=$a0*sqrt($v**2+$w**2+$v*$w);

	return ($r, $Lz, $l);
}

sub createBox {
	my ($r, $Lz) = @_;
	my ($box);

	$box->{X}{len}=$box->{XCOORD}{len}=$box->{X}{max}=$box->{XCOORD}{hi}=$r;
	$box->{Y}{len}=$box->{YCOORD}{len}=$box->{Y}{max}=$box->{YCOORD}{hi}=$r;
	$box->{Z}{len}=$box->{ZCOORD}{len}=$box->{Z}{max}=$box->{ZCOORD}{hi}=$Lz;
	$box->{X}{ANGLE}=$box->{Y}{ANGLE}=$box->{ALPHA}=$box->{BETA} = 90;
	$box->{Z}{ANGLE}=$box->{GAMMA}= 60;
	$box->{Z}{min}=$box->{Y}{min}=$box->{X}{min}=$box->{XCOORD}{lo}=$box->{YCOORD}{lo}=$box->{ZCOORD}{lo}=0;
	return $box

}

sub addMissingFields {
	my ($atoms, $bonds, $resnum) = @_;
	my ($i, $j, $nbond, $rec);

	$rec = (
			{
				LABEL     => "HETATM",
				RESNUM    => $resnum+1,
				RESNAME   => "CNT",
				FFTYPE    => "C_2G",
				LONEPAIRS => 0,
				CHARGE    => 0,
				CHAIN     => "A",
			}
		);

	for $i (keys %{ $atoms }) {
		for $j (keys %{ $rec }) {
			$atoms->{$i}{$j} = $rec->{$j};
		}
		$atoms->{$i}{ATMNAME} = "C${i}";
		$nbond = 0;
		$nbond = scalar(@{ $bonds->{$i} }) if (defined($bonds) and defined($bonds->{$i}));
		$atoms->{$i}{NUMBONDS} = $nbond;
	}
}

sub createBonds {
	my ($atoms, $a0, $Lz) = @_;
	my ($i, $j, $atom1,$atom2, $bonds, $dist, $dist_max, @tmp, $n, $m);

	$dist_max = $a0/sqrt(3);

	@tmp = sort { $a<=>$b } keys %{ $atoms };
	for $n (0 .. $#tmp) {
		$i = $tmp[$n];
		$atom1 = $atoms->{$i};
		for $m ($n+1 .. $#tmp) {
			$j = $tmp[$m];
			%{ $atom2 } = %{ $atoms->{$j} };
			$dist = GetBondLength($atom1, $atom2);
			if ($dist<=$dist_max) {
				push @{ $bonds->{$i} }, $j;
				push @{ $bonds->{$j} }, $i;
			}
			$atom2->{ZCOORD} += $Lz;
			$dist = GetBondLength($atom1, $atom2);
			if ($dist<=$dist_max) {
				push @{ $bonds->{$i} }, $j;
				push @{ $bonds->{$j} }, $i;
				$atoms->{$i}{DISPZ}[$#{ $bonds->{$i} } ] = 1;
				$atoms->{$j}{DISPZ}[$#{ $bonds->{$j} } ] = -1;
			} else {
				$atom2->{ZCOORD} -= 2*$Lz;
				$dist = GetBondLength($atom1, $atom2);
				if ($dist<=$dist_max) {
					push @{ $bonds->{$i} }, $j;
					push @{ $bonds->{$j} }, $i;
					$atoms->{$i}{DISPZ}[$#{ $bonds->{$i} } ] = -1;
					$atoms->{$j}{DISPZ}[$#{ $bonds->{$j} } ] = 1;
				}
			}
		}
	}
	return $bonds;
}

sub createAtomsXYZ {
	my ($a0, $eps, $j, $sid) = @_;
	my ($i, $l, $r, $v, $w, $s, $fac, $Lz, $id);
	my ($pmin, $pmax, $qmin, $qmax, $p, $q, $xx, $valid);
	my ($x, $y, $z, $atoms, $curr);

	$id = $sid + 1;
	%{ $atoms } = ();

	$l=sqrt($m->[$j]**2+$n->[$j]**2+$m->[$j]*$n->[$j]);
	$r=$l*$a0/PI+5;
	$v=$n->[$j]+2*$m->[$j];
	$w=2*$n->[$j]+$m->[$j];
	$s->[1]=$w/$l/2;
	$s->[2]=$v/$l/2; 
	$s->[3]=-sqrt(3)*$m->[$j]/$l/2;
	$s->[4]=sqrt(3)*$n->[$j]/$l/2;
	$fac=getFactor($v,$w);
	$v/=-$fac;
	$w/=$fac;
	$Lz=$a0*sqrt($v**2+$w**2+$v*$w);
	$pmin=minimum(minimum(0,$n->[$j]),minimum($v,$v+$n->[$j]));
	$pmax=maximum(maximum(0,$n->[$j]),maximum($v,$v+$n->[$j]));
	$qmin=minimum(minimum(0,$m->[$j]),minimum($w,$w+$m->[$j]));
	$qmax=maximum(maximum(0,$m->[$j]),maximum($w,$w+$m->[$j]));

	for $p ($pmin..$pmax) {
		for $q ($qmin..$qmax) {
			$xx = $a0*($s->[1]*$p+$s->[2]*$q);
			next if ($xx<-$eps or $xx>($l*$a0+$eps));
			$z->[$id] = ($s->[3]*$p+$s->[4]*$q)*$a0;
			next if ($z->[$id]<-$eps or $z->[$id]>($Lz+$eps));
			$y = -$l*$a0/PI/2*cos($xx*PI*2/$l/$a0);
			$x = $l*$a0/PI/2*sin($xx*PI*2/$l/$a0);
			$curr->{XCOORD} = $x; $curr->{YCOORD} = $y; $curr->{ZCOORD} = $z->[$id];
			next if (!checkAtomOverlap($atoms, $curr, $eps, $r, $Lz));
			$atoms->{$id}{ZCOORD} = $z->[$id]; $atoms->{$id}{YCOORD} = $y; $atoms->{$id}{XCOORD} = $x;
			$id++;

			$xx=$xx+$a0*($s->[1]+$s->[2])/3;     
			$z->[$id]=$z->[$id-1]+($s->[3]+$s->[4])*$a0/3;
			$y=-$l*$a0/PI/2*cos($xx*PI*2/$l/$a0);
			$x=$l*$a0/PI/2*sin($xx*PI*2/$l/$a0);
			$curr->{XCOORD} = $x; $curr->{YCOORD} = $y; $curr->{ZCOORD} = $z->[$id];
			next if (!checkAtomOverlap($atoms, $curr, $eps, $r, $Lz));
			$atoms->{$id}{ZCOORD} = $z->[$id]; $atoms->{$id}{YCOORD} = $y; $atoms->{$id}{XCOORD} = $x;
			$id++;
		}
	}
	return ($atoms, $id, $r, $Lz, $l);
}

sub checkAtomOverlap {
	my ($atoms, $curr, $eps, $r, $Lz) = @_;
	my ($i, $valid, $x, $y, $z);
	
	$valid = 1;
	$x = $curr->{XCOORD}; $y = $curr->{YCOORD}; $z = $curr->{ZCOORD};
	for $i (keys %{ $atoms }) {
		if ((abs($x-$atoms->{$i}{XCOORD})<$eps or abs(abs($x-$atoms->{$i}{XCOORD})-$r)<$eps) and
			(abs($y-$atoms->{$i}{YCOORD})<$eps or abs(abs($y-$atoms->{$i}{YCOORD})-$r)<$eps) and
			(abs($z-$atoms->{$i}{ZCOORD})<$eps or abs(abs($z-$atoms->{$i}{ZCOORD})-$Lz)<$eps)) {
			$valid = 0; 
			last;
		}
	}

	return $valid;
}

sub newHeaders {
	my ($a0, $l, $Lz, $box) = @_;
	my ($headers) = createHeaders($box, $saveName);

	&insertHeaderRemark($headers, "REMARK BGF file created by $0");
	&insertHeaderRemark($headers, sprintf("REMARK Bond Distance %10.5f",$a0/sqrt(3)));
	&insertHeaderRemark($headers, sprintf("REMARK C = %10.5f, with R = %10.5f", $l*$a0, $l*$a0/2/PI));
	&insertHeaderRemark($headers, sprintf("REMARK Lattice length in Z: %10.5f", $Lz));
	return $headers;
}

sub minimum {
	my ($i, $j) = @_;
	return $i if ($i<$j);
	return $j;
}

sub maximum {
	my ($i, $j) = @_;
	return $i if ($i>$j);
	return $j;
}

sub getFactor {
	my ($i, $j) = @_;
	my $fac = minimum($i, $j);
	while ($fac > 1) {
		return $fac if ($i%$fac==0 and $j%$fac==0);
		$fac--;
	}
	return $fac;
}

sub init {
	my (%OPTS, $ctheta, $n_str, $m_str, $list, $i);

	getopt('nms',\%OPTS);
	for ("n", "m") {
		die "usage: $0 -m 'm1 m2 m3...' -n 'n1 n2 n3...' -s [savename]\n" 
			if (! exists($OPTS{$_}));
	}
	($n_str, $m_str, $saveName) = ($OPTS{n}, $OPTS{m}, $OPTS{s});
	@{ $n } = split /\s+/, $n_str;
	@{ $m } = split /\s+/, $m_str;
	die "ERROR: Expect n ($n_str) and m ($m_str) to have the same number of entries!\n"
		if ($#{ $n } != $#{ $m });
	for $i (0 .. $#{ $m }) {
		die "ERROR: m ($m->[$i]) n ($n->[$i]) listed twice!\n"
			if (defined($list) and exists($list->{$m->[$i]}{$n->[$i]}));
		$list->{ $m->[$i] }{ $n->[$i] } = 1;
	}
	$nwall = scalar(@{ $n });
	print "Initializing...";
	$saveName = $n->[0] . "x" . $m->[0]. "cnt.bgf" if (! defined($saveName));
}
