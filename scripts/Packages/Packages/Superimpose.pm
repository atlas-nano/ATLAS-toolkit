use strict;
{
	package Kearsley;
	use Math::MatrixReal;
	use General qw(CoM);
	
	require Exporter;
	our @ISA = qw(Exporter);
	our @EXPORT = ();
	our @EXPORT_OK = qw(SuperimposeAtoms GetRotMatrix RotateMolecule TranslateAllAtoms);
	
	sub numerically { ($a<=>$b); }
	
	sub calcKearsleyMatrix {
	    my ($ref, $mov) = @_;
	    my (@km, $i, $j, @refMtchAtms, @movMtchAtms);
	    my ($Xm, $Xp, $Ym, $Yp, $Zm, $Zp);
	    my ($Xm2, $Xp2, $Ym2, $Yp2, $Zm2, $Zp2);
	    
	    for $i (0 .. 3) {
	        for $j (0 .. 3) {
	            $km[$i][$j] = 0;
	        }
	    }
	    
	    $j = 0;
	    for $i (sort numerically keys %{ $ref }) {
	        $refMtchAtms[$j][0] = $ref->{$i}{"XCOORD"};
	        $refMtchAtms[$j][1] = $ref->{$i}{"YCOORD"};
	        $refMtchAtms[$j][2] = $ref->{$i}{"ZCOORD"};
	        $j++;
	    }
	    
	    $j = 0;
	    for $i (sort numerically keys %{ $mov }) {
	        $movMtchAtms[$j][0] = $mov->{$i}{"XCOORD"};
	        $movMtchAtms[$j][1] = $mov->{$i}{"YCOORD"};
	        $movMtchAtms[$j][2] = $mov->{$i}{"ZCOORD"};
	        $j++;
	    }
	    
	    die "ERROR: Unequal number of atoms in ref and match\n" if ($#refMtchAtms != $#movMtchAtms);
	  
	    $j--;
	    for $i (0 .. $j) {
			$Xm = -1*($movMtchAtms[$i][0] - $refMtchAtms[$i][0]);
			$Ym = -1*($movMtchAtms[$i][1] - $refMtchAtms[$i][1]);
			$Zm = -1*($movMtchAtms[$i][2] - $refMtchAtms[$i][2]);
	
			$Xp = ($refMtchAtms[$i][0] + $movMtchAtms[$i][0]);
			$Yp = ($refMtchAtms[$i][1] + $movMtchAtms[$i][1]);
			$Zp = ($refMtchAtms[$i][2] + $movMtchAtms[$i][2]);
	
			$Xm2 = $Xm*$Xm;
			$Ym2 = $Ym*$Ym;
			$Zm2 = $Zm*$Zm;
	
			$Xp2 = $Xp*$Xp;
			$Yp2 = $Yp*$Yp;
			$Zp2 = $Zp*$Zp;
	
			$km[0][0] = $km[0][0] +  ($Xm2 + $Ym2 + $Zm2);
			$km[0][1] = $km[0][1] +  ($Yp*$Zm - $Ym*$Zp);
			$km[0][2] = $km[0][2] +  ($Xm*$Zp - $Xp*$Zm);
			$km[0][3] = $km[0][3] +  ($Xp*$Ym - $Xm*$Yp);
			$km[1][0] = $km[1][0] +  ($Yp*$Zm - $Ym*$Zp);
			$km[1][1] = $km[1][1] +  ($Yp2 + $Zp2 + $Xm2);
			$km[1][2] = $km[1][2] +  ($Xm*$Ym - $Xp*$Yp);
			$km[1][3] = $km[1][3] +  ($Xm*$Zm - $Xp*$Zp);
			$km[2][0] = $km[2][0] +  ($Xm*$Zp - $Xp*$Zm);
			$km[2][1] = $km[2][1] +  ($Xm*$Ym - $Xp*$Yp);
			$km[2][2] = $km[2][2] +  ($Xp2 + $Zp2 + $Ym2);
			$km[2][3] = $km[2][3] +  ($Ym*$Zm - $Yp*$Zp);
			$km[3][0] = $km[3][0] +  ($Xp*$Ym - $Xm*$Yp);
			$km[3][1] = $km[3][1] +  ($Xm*$Zm - $Xp*$Zp);
			$km[3][2] = $km[3][2] +  ($Ym*$Zm - $Yp*$Zp);
			$km[3][3] = $km[3][3] +  ($Xp2 + $Yp2 + $Zm2);        
	    }
	    
	    return \@km;
	}
	
	sub findDiagonalizingEuler {
	    my ($evec) = $_[0];
	    my ($e0, $e1, $e2, $e3) = (0, 0, 0, 0);
	    my ($ee0, $ee1, $ee2, $ee3);
	    my (@Rot_Mat);
	    
	    $e0 = $evec->[0][3] if($evec->[0][3]);
	    $e1 = $evec->[1][3] if($evec->[1][3]);
	    $e2 = $evec->[2][3] if($evec->[2][3]);
	    $e3 = $evec->[3][3] if($evec->[3][3]);
	
	    $ee0 = $e0*$e0;
	    $ee1 = $e1*$e1;
	    $ee2 = $e2*$e2;
	    $ee3 = $e3*$e3;
	
	
	    $Rot_Mat[0][0] = $ee0 + $ee1 - $ee2 - $ee3;
	    $Rot_Mat[0][1] = 2*($e1*$e2 + $e0*$e3);
	    $Rot_Mat[0][2] = 2*($e1*$e3 - $e0*$e2);
	    $Rot_Mat[1][0] = 2*($e1*$e2 - $e0*$e3);
	    $Rot_Mat[1][1] = $ee0 - $ee1 + $ee2 - $ee3;
	    $Rot_Mat[1][2] = 2*($e2*$e3 + $e0*$e1);
	    $Rot_Mat[2][0] = 2*($e1*$e3 + $e0*$e2);
	    $Rot_Mat[2][1] = 2*($e2*$e3 - $e0*$e1);
	    $Rot_Mat[2][2] = $ee0 - $ee1 - $ee2 + $ee3;
	
	    return \@Rot_Mat;   
	}
	
	
	sub SuperimposeAtoms {
	    my ($refMol, $movMol, $movMol2) = @_;
	    my ($ref_com, $mov_com, $kearsleyMatrix, $kearsleyEigenVec, $RotMat);
	    my ($start, $end);
	    
	    $ref_com = CoM($refMol);
	    $mov_com = CoM($movMol);
	    #translate both molecules to origin
	    &TranslateAllAtoms($refMol, $ref_com, -1);
	    &TranslateAllAtoms($movMol, $mov_com, -1);
	    
	    #get the rotation matrix
	    $RotMat = GetRotMatrix($refMol, $movMol);
		#printMatrix($RotMat);
		if (defined($movMol2)) {
			$mov_com = CoM($movMol2);
			&TranslateAllAtoms($movMol2, $mov_com, -1);
			&RotateMolecule($movMol2, $RotMat);
			&TranslateAllAtoms($movMol2, $ref_com, 1);
		} else {
			&RotateMolecule($movMol, $RotMat);
			&TranslateAllAtoms($movMol, $ref_com, 1);
		}
	    &TranslateAllAtoms($refMol, $ref_com, 1);
	}
	
	sub GetRotMatrix {
	    my ($refMol, $movMol) = @_;
	    my ($kearsleyMatrix, $kearsleyEigenVec, $RotMat);
	
	    $kearsleyMatrix = calcKearsleyMatrix($refMol, $movMol);
		#printMatrix($kearsleyMatrix);
	    $kearsleyEigenVec = diagMatrix($kearsleyMatrix, 4);
	    #printMatrix($kearsleyEigenVec);
	    $RotMat = findDiagonalizingEuler($kearsleyEigenVec);
	    return $RotMat;
	}
	
	sub printMatrix {
	    my ($matrix) = $_[0];
	    my ($i, $j, $outString);
	    
	    $outString = "\n{";
	    for $i (@{ $matrix }) {
	        $outString .= "{";
	        for $j (@{ $i }) {
	            $outString .= "$j,";
	        }
	        chop $outString;
	        $outString .= "},";
	    }
	    chop $outString;
	    $outString .= "}\n";
	    
	    print $outString;
	}
	
	sub RotateMolecule {
	    my ($atoms, $rotMatrix) = @_;
	    my (@angles, $PI, $i, $j, $atom, $xC, $yC, $zC, $dim);
	
	    for $i (keys %{ $atoms }) {
	        $atom = $atoms->{$i};
	        $xC = $atom->{"XCOORD"};
	        $yC = $atom->{"YCOORD"};
	        $zC = $atom->{"ZCOORD"};
	        $j = 0;
	        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	            $atom->{$dim} = $xC * $rotMatrix->[$j][0] + $yC * $rotMatrix->[$j][1] + $zC * $rotMatrix->[$j][2];
	            $j++;
	        }
	    }
	            
	}
	
	sub diagMatrix {
	    my ($matrixA, $size) = @_;
	    my ($new_matrix, $V, $E, $eigenVec);
	    
	    $new_matrix = Math::MatrixReal->new_from_rows([\@{ $matrixA->[0]}, \@{ $matrixA->[1]}, \@{ $matrixA->[2]}, \@{ $matrixA->[3]}]);
	    #print "Input Matrix:\n";
	    #print $new_matrix->as_yacas( ( format => "%.8f", align => "l",name => "A" ) );
	    
	    ($E, $V) = $new_matrix->sym_diagonalize();
	    #print "\nEigenvalues: \n";
	    #print $E->as_yacas( ( format => "%.5G", align => "l",name => "A" ) );    
	    #print "\nEigenvectors: \n";
	    #print $V->as_yacas( ( format => "%.5G", align => "l",name => "A" ) );
	    $eigenVec = eigenSort($E, $V, $size);
	    #print "\nBest Eigenvector: \n";
	    #printMatrix($eigenVec);
	    return $eigenVec;
	}
	
	sub eigenSort {
	    my ($eVal, $eVec, $size) = @_;
	    my ($i, $rec, %eigenVecs, $j, @tmp, @ret);
	    
	    for $i (1 .. $size) {
	        $rec = ();
	        for $j (1 .. $size) {
	            push @{ $rec }, $eVec->element($i,$j);
	        }
	        $eigenVecs{ $eVal->element($i,1) } = $rec;
	    }
	   
	    @tmp = sort numerically keys %eigenVecs;
	    
	    for $i (@tmp) {
	        push @ret, $eigenVecs{$i};
	    }
	    
	    return \@ret;
	}
	
	sub TranslateAllAtoms {
	    my ($atoms, $transVec, $direction) = @_;
	    my ($i, $dim);
	    
	    for $i (keys %{ $atoms }) {
	        for $dim (keys %{ $transVec }) {
	            $atoms->{$i}{$dim} += $direction * $transVec->{$dim};
	        }
	    }
	}
}

{
	#Superimpose atoms based on Quarternions
	#usage: QSuperimposeAtoms(reference,atoms,[moveAtoms])
	#note that the moveAtoms are optional and defaults to the atoms if not specified

	package Quarternion;
	use strict;
	use Math::MatrixReal;
	use General qw(CoM);
	
	require Exporter;
	our @ISA = qw(Exporter);
	our @EXPORT = ();
	our @EXPORT_OK = qw(SuperimposeAtoms QuarternionRotate QuarternionTransform);
	
	sub SuperimposeAtoms {
	    my ($refMol, $movMol, $movMol2) = @_;
	    my ($transform);
	
	    $movMol2 = $movMol if (! defined($movMol2));
	
	    #get quarternion transformation matrix
	    $transform = QuarternionTransform($refMol,$movMol);
	
	    #move molecule com to reference
	    &TranslateAllAtoms($movMol2, $transform->{Trans}, -1);
	
	    #rotate
	    &QuarternionRotate($movMol2, $transform->{Rot});

	}
	
	#do a superposition of an atom list to a reference (2nd argument)
	#selections must have the same number of atoms
	#Using method from: S. Kearsley, Acta Cryst. A45, 208-210 1989
	sub QuarternionTransform {
	    my ($a1, $a2) = @_;
	    my ($i, $j, $nrd1, $nrd2, $gc1, $gc2); 
	    my ($x1, $x2, $y1, $y2, $z1, $z2);
	    my ($iatom1, $iatom2);
	
	    my ($xm,$ym,$zm);
	    my ($xp,$yp,$zp);
	    my ($Sxmxm, $Sxpxp, $Symym, $Sypyp, $Szmzm, $Szpzp); 
	    my ($Sxmym, $Sxmyp, $Sxpym, $Sxpyp);
	    my ($Sxmzm, $Sxmzp, $Sxpzm, $Sxpzp);
	    my ($Symzm, $Symzp, $Sypzm, $Sypzp); 
	
	    $nrd1 = scalar(keys %{ $a1 }); 
	    $nrd2 = scalar(keys %{ $a2 });
	    die "superpose error: lists must have same number of atoms\n"  if ($nrd1 != $nrd2);
	
	    #get the geometric center of the molecules
	    $gc1 = CoM($a1);
	    $gc2 = CoM($a2);
	
	    #construct a 4X4 matrix in the quaternion representation
	    for $i (keys %{ $a1 }) {
	        $iatom1 = $a1->{$i};
	        $x1 =  $iatom1->{XCOORD} - $gc1->{XCOORD};
	        $y1 =  $iatom1->{YCOORD} - $gc1->{YCOORD};
	        $z1 =  $iatom1->{ZCOORD} - $gc1->{ZCOORD};
	    
	        $iatom2 = $a2->{$i};
	        $x2 =  $iatom2->{XCOORD} - $gc2->{XCOORD};
	        $y2 =  $iatom2->{YCOORD} - $gc2->{YCOORD};
	        $z2 =  $iatom2->{ZCOORD} - $gc2->{ZCOORD};
	
	        $xm = ($x1 - $x2);
	        $xp = ($x1 + $x2);           
	        $ym = ($y1 - $y2);
	        $yp = ($y1 + $y2);             
	        $zm = ($z1 - $z2);
	        $zp = ($z1 + $z2);
	
	        $Sxmxm  += $xm*$xm; 
	        $Sxpxp  += $xp*$xp;
	        $Symym  += $ym*$ym; 
	        $Sypyp  += $yp*$yp; 
	        $Szmzm  += $zm*$zm; 
	        $Szpzp  += $zp*$zp; 
	
	        $Sxmym  += $xm*$ym; 
	        $Sxmyp  += $xm*$yp; 
	        $Sxpym  += $xp*$ym; 
	        $Sxpyp  += $xp*$yp; 
	
	        $Sxmzm  += $xm*$zm; 
	        $Sxmzp  += $xm*$zp; 
	        $Sxpzm  += $xp*$zm; 
	        $Sxpzp  += $xp*$zp; 
	
	        $Symzm  += $ym*$zm; 
	        $Symzp  += $ym*$zp; 
	        $Sypzm  += $yp*$zm; 
	        $Sypzp  += $yp*$zp;
	    }
	    my @m;
	    $m[0]= $Sxmxm + $Symym + $Szmzm;
	    $m[1]= $Sypzm - $Symzp;
	    $m[2]= $Sxmzp - $Sxpzm;
	    $m[3]= $Sxpym - $Sxmyp;
	    $m[4]= $m[1];
	    $m[5]= $Sypyp + $Szpzp + $Sxmxm;
	    $m[6]= $Sxmym - $Sxpyp;
	    $m[7]= $Sxmzm - $Sxpzp;
	    $m[8]= $m[2];
	    $m[9]= $m[6];
	    $m[10]= $Sxpxp + $Szpzp + $Symym;
	    $m[11]= $Symzm - $Sypzp;
	    $m[12]= $m[3];
	    $m[13]= $m[7];
	    $m[14]= $m[11];
	    $m[15]=$Sxpxp + $Sypyp + $Szmzm;
	    #compute the egienvectors and eigenvalues of the matrix
	    my ( $revec, $reval ) = diagonalize(@m);
	    #the smallest eigenvalue is the rmsd for the optimal alignment
	    my $rmsd = sqrt(abs($reval->[0])/ $nrd1 );
	    #fetch the optimal quaternion
	    my @q;
	    $q[0]=$revec->[0][0];
	    $q[1]=$revec->[1][0];
	    $q[2]=$revec->[2][0];
	    $q[3]=$revec->[3][0];
	    #construct the rotation matrix given by the quaternion
	    my @mt;
	    $mt[0] = $q[0]*$q[0] + $q[1]*$q[1] - $q[2]*$q[2] - $q[3]*$q[3];
	    $mt[1] = 2.0 * ($q[1] * $q[2] - $q[0] * $q[3]);
	    $mt[2] = 2.0 * ($q[1] * $q[3] + $q[0] * $q[2]); 
	
	    $mt[3] = 2.0 * ($q[2] * $q[1] + $q[0] * $q[3]);
	    $mt[4] = $q[0]*$q[0] - $q[1]*$q[1] + $q[2]*$q[2] - $q[3]*$q[3];
	    $mt[5] =  2.0 * ($q[2] * $q[3] - $q[0] * $q[1]);
	
	    $mt[6] = 2.0 *($q[3] * $q[1] - $q[0] * $q[2]);
	    $mt[7] = 2.0 * ($q[3] * $q[2] + $q[0] * $q[1]);
	    $mt[8] = $q[0]*$q[0] - $q[1]*$q[1] - $q[2]*$q[2] + $q[3]*$q[3];
	
	    #compute the displacement vector
	    my @vt;
	    $vt[0] = $gc2->{XCOORD} - $mt[0]*$gc1->{XCOORD} - $mt[1]*$gc1->{YCOORD} - $mt[2]*$gc1->{ZCOORD};
	    $vt[1] = $gc2->{YCOORD} - $mt[3]*$gc1->{XCOORD} - $mt[4]*$gc1->{YCOORD} - $mt[5]*$gc1->{ZCOORD};  
	    $vt[2] = $gc2->{ZCOORD} - $mt[6]*$gc1->{XCOORD} - $mt[7]*$gc1->{YCOORD} - $mt[8]*$gc1->{ZCOORD};
	
	    my $ret = (
	                {
	                    "Trans" => {
	                                        "XCOORD" => $vt[0],
	                                        "YCOORD" => $vt[1],
	                                        "ZCOORD" => $vt[2],
	                    },
	                    "Rot" => \@mt,
	                }
	    );
	
	    #return the transformation as one list rotation first
	    return $ret;
	}
	
	sub diagonalize {
	   my ($onorm, $dnorm);
	   my ($b,$dma,$q,$t,$c,$s);
	   my ($atemp, $vtemp, $dtemp);
	   my ($i,$j,$k,$l);
	   my @a;
	   my @v;
	   my @d;
	   my $nrot = 30; #number of sweeps
	
	   for ($i = 0; $i < 4; $i++) 
	      {
	      for ($j = 0; $j < 4; $j++)
	         {
	         $a[$i][$j] =$_[4*$i + $j];
	         $v[$i][$j] = 0.0;
	         }
	      }
	
	   for ($j = 0; $j <= 3; $j++) 
	      {
	      $v[$j][$j] = 1.0;
	      $d[$j] = $a[$j][$j];
	      }
	
	   for ($l = 1; $l <= $nrot; $l++) 
	      {
	      $dnorm = 0.0;
	      $onorm = 0.0;
	      for ($j = 0; $j <= 3; $j++)
	         {
	         $dnorm +=  abs($d[$j]);
	         for ($i = 0; $i <= $j - 1; $i++)
	            {
	            $onorm += abs($a[$i][$j]);
	            }
	         }
	      last if(($onorm/$dnorm) <= 1.0e-12);
	      for ($j = 1; $j <= 3; $j++) 
	         {
	         for ($i = 0; $i <= $j - 1; $i++) 
	            {
	            $b = $a[$i][$j];
	            if(abs($b) > 0.0) 
	               {
	               $dma = $d[$j] - $d[$i];
	               if((abs($dma) + abs($b)) <=  abs($dma)) 
	                  {
	                  $t = $b / $dma;
	                  }
	               else 
	                  {
	                  $q = 0.5 * $dma / $b;
	                  $t = 1.0/(abs($q) + sqrt(1.0+$q*$q));
	                  $t *= -1.0 if($q < 0.0); 
	                  }
	               $c = 1.0/sqrt($t * $t + 1.0);
	               $s = $t * $c;
	               $a[$i][$j] = 0.0;
	               for ($k = 0; $k <= $i-1; $k++) 
	                  {
	                  $atemp = $c * $a[$k][$i] - $s * $a[$k][$j];
	                  $a[$k][$j] = $s * $a[$k][$i] + $c * $a[$k][$j];
	                  $a[$k][$i] = $atemp;
	                  }
	               for ($k = $i+1; $k <= $j-1; $k++)
	                  {
	                  $atemp = $c * $a[$i][$k] - $s * $a[$k][$j];
	                  $a[$k][$j] = $s * $a[$i][$k] + $c * $a[$k][$j];
	                  $a[$i][$k] = $atemp;
	                  }
	               for ($k = $j+1; $k <= 3; $k++) 
	                  {
	                  $atemp = $c * $a[$i][$k] - $s * $a[$j][$k];
	                  $a[$j][$k] = $s * $a[$i][$k] + $c * $a[$j][$k];
	                  $a[$i][$k] = $atemp;
	                  }
	               for ($k = 0; $k <= 3; $k++) 
	                  {
	                  $vtemp = $c * $v[$k][$i] - $s * $v[$k][$j];
	                  $v[$k][$j] = $s * $v[$k][$i] + $c * $v[$k][$j];
	                  $v[$k][$i] = $vtemp;
	                  }
	               $dtemp = $c*$c*$d[$i] + $s*$s*$d[$j] - 2.0*$c*$s*$b;
	               $d[$j] = $s*$s*$d[$i] + $c*$c*$d[$j] +  2.0*$c*$s*$b;
	               $d[$i] = $dtemp;
	               } 
	            }  
	         }
	      } 
	
	   $nrot = $l;
	   for ($j = 0; $j <= 2; $j++) 
	      {
	      $k = $j;
	      $dtemp = $d[$k];
	      for ($i = $j+1; $i <= 3; $i++) 
	         {
	         if($d[$i] < $dtemp) 
	            {
	            $k = $i;
	            $dtemp = $d[$k];
	            }
	         }
	
	      if($k > $j) 
	         {
	         $d[$k] = $d[$j];
	         $d[$j] = $dtemp;
	         for ($i = 0; $i <= 3; $i++) 
	            {
	            $dtemp = $v[$i][$k];
	            $v[$i][$k] = $v[$i][$j];
	            $v[$i][$j] = $dtemp;
	            }
	         }
	      }
	   
	   return (\@v,\@d)
	}
	
	#rotate an atom selection by the specified matrix, Ax = x'.
	#the matrix is specified by a flat list of the matrix elements.
	sub QuarternionRotate
	{
	    my ($atoms, $rotMat) = @_;
	    my ($x1,$y1,$z1);
	    my ($x2,$y2,$z2);
	    my ($nrd, $i, $iatom);
	
	    my ($u11, $u12, $u13, $u21, $u22, $u23, $u31, $u32, $u33) = @{ $rotMat };
	
	    $nrd = scalar(keys %{ $atoms });
	 
	    for $i (keys %{ $atoms })
	    {
	        $iatom = $atoms->{$i};
	        $x1 = $iatom->{XCOORD};
	        $y1 = $iatom->{YCOORD};
	        $z1 = $iatom->{ZCOORD};
	        $x2 = $u11*$x1 + $u12*$y1 + $u13*$z1;
	        $y2 = $u21*$x1 + $u22*$y1 + $u23*$z1;
	        $z2 = $u31*$x1 + $u32*$y1 + $u33*$z1;
	        $iatom->{XCOORD} = $x2;
	        $iatom->{YCOORD} = $y2;
	        $iatom->{ZCOORD} = $z2;
	    }
	    return
	}
	sub TranslateAllAtoms {
	    my ($atoms, $transVec, $direction) = @_;
	    my ($i, $dim);
	    
	    for $i (keys %{ $atoms }) {
	        for $dim (keys %{ $transVec }) {
	            $atoms->{$i}{$dim} += $direction * $transVec->{$dim};
	        }
	    }
	}
}
	
{
	package Kabash;
	#uses Math::MatrixReal and Linalg::SVD
	#see Kabash.pm alternative formulation using PDL
	use strict;
	use Math::MatrixReal;
	use Math::GSL::Linalg::SVD;

	require Exporter;
	our @ISA = qw(Exporter);
	our @EXPORT = ();
	our @EXPORT_OK = qw(SuperimposeAtoms);
	
	sub SuperimposeAtoms {
		my ($ref, $mol, $mol2) = @_;
		my ($parray, $qarray, $tarray, $tindx);
		my ($A, $B, $C, $R);
		my ($i, $j, $transform, $com_C);

		$mol2 = $mol if (! defined($mol2));
		($parray, $qarray, $tarray, $tindx) = createKabashArrays($ref,$mol, $mol2);

		$A = Math::MatrixReal->new_from_rows([@{ $parray }]);
		$B = Math::MatrixReal->new_from_rows([@{ $qarray }]);
		$C = Math::MatrixReal->new_from_rows([@{ $tarray }]);

		#first remove center of mass displacement
		$com_C = centroid($C);
		$A -= centroid($A);
		$B -= centroid($B);
		$C -= $com_C;

		#now get rotation matrix
		$R = kabsch($B, $A);

		#now apply best rotation
		$C *= $R;

		#finally place matched coordinates in original reference framework
		$C += $com_C;

		#and create the updated atom coordinate array
		for $i (0 .. $#{ $tarray }) {
			for $j (0 .. $#{ $tarray->[$i]}) {
				$transform->[$i][$j] = $C->element($i+1,$j+1);
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

	sub centroid {
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
}

1;
