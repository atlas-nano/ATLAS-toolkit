#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts");
}

use strict;
use Packages::FileFormats;
use Packages::General;
use File::Basename;

#   atomistic2Meso.pl: This script will open a pdb file and will obtain the
#   Mesoscale parameters (angles, torsions, bond_lengths)
#   for the appropriate "beads" defined in the parameter file
#   It will then output a .pqr file corresponding to the meso-scale file

sub Initialize();
sub GetParms(@);
sub CreateMesoModel(@);
sub CreateHBondBead(@);
sub WritePQRFile(@);
sub GetSaveName(@);
sub Numerically;
sub GET_COM(@);
sub GetMass(@);
sub IsEndOfChain(@);
sub isPresesnt(@);
sub IsMember(@);
sub GetName(@);
sub GetDiff(@);
sub GetIndex(@);
sub FixConnections(@);
sub IsIndex(@);
sub ExtractResidue(@);
sub GetChain(@);
sub ParsePDBFile(@);
sub GetConnections(@);
sub GetNameOrder(@);

die "usage: $0 parm_file pdb_file|directory [strand_length] [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($parm_file, $pdb, $strand_length, $save_name) = @ARGV;
my ($curr_file, $pqr_file, $PDB_FILES, $Connections);
my ($PARMS, $Atoms_Info, $Meso_Info, $counter, $BEADS);


$PDB_FILES = Initialize();
print "Obtaining Bead parameters ....";
($PARMS) = GetParms($parm_file);
print "Done\n\n";

$counter = 0;

for $curr_file (@{ $PDB_FILES }) {
    print "$curr_file:....READING";
    ($strand_length, $Atoms_Info) = ParsePDBFile($curr_file);
    next
	if (! $strand_length);
    print "....Creating MESO";
    $Meso_Info = CreateMesoModel($Atoms_Info, $PARMS);
    $pqr_file = GetSaveName($curr_file, $counter);
    print "....Writing Meso";
    WritePQRFile($Meso_Info, $pqr_file, $PARMS, $strand_length, $counter);
    print "Done\n";
}

sub Initialize() {
    my (@PDBS);

    if (! -e $parm_file or ! -r $parm_file or ! -T $parm_file) {
	die "ERROR: Invalid Parameter file: $!\n";
    }

    if (! -T $pdb) {
	if (-d $pdb) {

	    opendir PDBFILES, $pdb or die "ERROR: Cannot access directory $pdb: $!\n";
	    @PDBS = grep { /\.pdb$/ && -f} map { "$pdb/$_"} readdir PDBFILES;
	    closedir PDBFILES;
	}
    } else {
	$PDBS[0] = $pdb;
    }

    die "ERROR: No valid pdbfiles found!\n"
	if ($#PDBS == -1);

    if (! defined($save_name)) {
	$save_name = "./";
    }

    return \@PDBS;
}

sub GetParms(@) {
    my ($parm) = $_[0];
    my ($curr_id, $rec, $mykey, %PARMS_HOLDER);

    $curr_id = -1;

    open PARMFILE, $parm or die "Cannot open parameter file $parm: $!\n";
    while(<PARMFILE>) {
	chomp;
	if ($_ =~ /^BEAD_ID: (\d+)/) {
	    if ($curr_id > -1) {
		for $mykey (keys %{ $rec }) {
		    $PARMS_HOLDER{"BEADS"}{$curr_id}{$mykey} = $rec->{$mykey};
		}
	    }
	    $curr_id = $1;
	} elsif ($curr_id > -1 and $_ =~ /^BEAD_(\w+): (.+)/) {
	    $rec->{$1} = $2;
	} elsif ($_ =~ /^BOND_(\d+): (\w+\-\w+):(.+)/) {
	    $PARMS_HOLDER{"BONDS"}{$1} = (
					  {
					      "VAL"   => $2,
					      "RANGE" => $3,
					  }
					  );
	} elsif ($_ =~ /^ANGLE_(\d+): (\w+\-\w+\-\w+):(.+)/) {
	    $PARMS_HOLDER{"ANGLES"}{$1} = (
					   {
					       "VAL"   => $2,
					       "RANGE" => $3,
					   }
					   );
	} elsif ($_ =~ /^TORSION_(\d+): (\w+\-\w+\-\w+\-\w+):(.+)/) {
	    $PARMS_HOLDER{"TORSIONS"}{$1} = (
					     {
						 "VAL"   => $2,
						 "RANGE" => $3,
					     }
					     );
	} elsif ($_ =~ /^MASS_(\w+)\s+(\d+\.*\d*)/) {
	    $PARMS_HOLDER{"MASSES"}{$1} = $2;
	}
	    
    }
    
    for $mykey (keys %{ $rec }) {
	$PARMS_HOLDER{"BEADS"}{$curr_id}{$mykey} = $rec->{$mykey};
    }

    close PARMFILE;
    

    die "ERROR: Invalid parameter file $parm\n"
	if ($curr_id == -1);
    return (\%PARMS_HOLDER);
}

sub Numerically {
    ($a<=>$b);
}

sub CreateMesoModel(@) {
    my ($atoms, $parms) = @_;
    my ($atm_counter, $bead_counter, $atm_label, $bead_member, $curr_id, $res_id);
    my (%MESO_MODEL, $atm_res, $base_res, $hbond, $non_member, $counter);
    my ($result, $isAcceptor, @search_elements, @hblist, $index);
    my ($avg, $stdev, $total, $tmp, $hbond_nm, $companion);

    $res_id = $curr_id = -1;
    for $atm_counter (sort Numerically keys %{ $atoms }) {
	$res_id = $atoms->{$atm_counter}{"RES_ID"};
	$atm_label = $atoms->{$atm_counter}{"LABEL"};
	$atm_res = $atoms->{$atm_counter}{"RES_NAME"};
	$atm_res =~ s/\d+//g;

	for $bead_counter (sort Numerically keys %{ $parms->{"BEADS"} }) {
	    $bead_member = $parms->{"BEADS"}{$bead_counter}{"MEMBERS"};
	    $non_member = $parms->{"BEADS"}{$bead_counter}{"NONMEMBERS"};
	    $hbond = $parms->{"BEADS"}{$bead_counter}{"HBONDS"};
	    $base_res = $parms->{"BEADS"}{$bead_counter}{"RES"};

	    if ($bead_member and IsMember($bead_member,$atm_label) and IsMember($base_res,$atm_res)) {
		$MESO_MODEL{$res_id}{$bead_counter}{"ATOMLIST"} .= "$atm_label,";
		$MESO_MODEL{$res_id}{$bead_counter}{"ATOMS"} .= "$atm_counter ";
		$MESO_MODEL{$res_id}{$bead_counter}{"INDEX"} = $bead_counter;
		GET_COM(\%{ $MESO_MODEL{$res_id}{$bead_counter} },\%{ $atoms->{$atm_counter} });
		$MESO_MODEL{$res_id}{$bead_counter}{"TYPE"} = $parms->{"BEADS"}{$bead_counter}{"NAME"};
		
		for $counter ("CHARGE", "RADII", "NUMBONDS", "LONEPAIRS", "ELEMENT") {
		    $MESO_MODEL{$res_id}{$bead_counter}{$counter} = $parms->{"BEADS"}{$bead_counter}{$counter};
		}
		
	    }elsif ($non_member and IsMember($non_member,$atm_label) and IsMember($base_res,$atm_res)) {
		$MESO_MODEL{$res_id}{$bead_counter}{"ATOMLIST"} .= "$atm_label,";
	    }
	}
    }

    for $res_id (sort Numerically keys %MESO_MODEL) {
	for $curr_id (sort Numerically keys %{ $MESO_MODEL{$res_id} }) {
	    $bead_counter = (\%{ $MESO_MODEL{$res_id}{$curr_id} });
	    $bead_counter->{"XCOORD"} = $bead_counter->{"COM"}{"X"}/$bead_counter->{"COM"}{"TotalMass"};
	    $bead_counter->{"YCOORD"} = $bead_counter->{"COM"}{"Y"}/$bead_counter->{"COM"}{"TotalMass"};
	    $bead_counter->{"ZCOORD"} = $bead_counter->{"COM"}{"Z"}/$bead_counter->{"COM"}{"TotalMass"};
	    for $counter (keys %{ $bead_counter->{"HBONDS"} }) {
		$MESO_MODEL{$res_id}{$index} = CreateHBondBead($index, $counter, \%{ $bead_counter->{"HBONDS"}{$counter} }, $curr_id);
		$index++;
	    }
	}
    }
	
    return \%MESO_MODEL;	
}

sub CreateHBondBead(@) {
    my ($index, $nm, $hbond, $parent) = @_;
    my (%Curr_Bead, $counter);
    
    $Curr_Bead{"ATOMLIST"} = "$nm";
    $Curr_Bead{"INDEX"} = 100 + $index;
    
    for $counter ("XCOORD", "YCOORD", "ZCOORD", "ATOMS", "COMPANION") {
	$Curr_Bead{$counter} = $hbond->{$counter};
    }
    $Curr_Bead{"CHARGE"} = 0.0;
    $Curr_Bead{"RADII"} = 1.058;
    $Curr_Bead{"NUMBONDS"} = 1;
    $Curr_Bead{"LONEPAIRS"} = 0;
    if ($index eq "51") {
	$Curr_Bead{"ELEMENT"} = "O";
    } elsif ($index eq "52") {
	$Curr_Bead{"ELEMENT"} = "N";
    } elsif ($index eq "53") {
	$Curr_Bead{"ELEMENT"} = "S";
    }	
    $Curr_Bead{"TYPE"} = "H" . $index;
    $Curr_Bead{"PARENT"} = $parent;

    return \%Curr_Bead;
}

sub GET_COM(@) {
    my ($curr_bead, $atm_counter) = @_;
    my ($atm_mass);

    $atm_mass = GetMass($atm_counter->{"LABEL"});

    $curr_bead->{"COM"}{"X"} += $atm_counter->{"XCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"Y"} += $atm_counter->{"YCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"Z"} += $atm_counter->{"ZCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"TotalMass"} += $atm_mass;

#    return $curr_bead;
}

sub GetSaveName(@) {
    my ($curr_file, $file_counter) = @_;

    if (-d $save_name) {
	$curr_file = $save_name . "/" . basename($curr_file);
    } elsif (-e $save_name) {
	$curr_file = $save_name;
    }

    $curr_file =~ s/\.pdb$//i;
    $curr_file .= ".pqr";

    return $curr_file;

}

sub ExtractResidue(@) {
    my ($curr_res) = $_[0];
    my ($returnval, $tmp);

    for $tmp (keys %{ $curr_res }) {
	if ($curr_res->{$tmp}{"TYPE"} and $curr_res->{$tmp}{"TYPE"} ne "SUG" and $curr_res->{$tmp}{"TYPE"} ne "PHO"
	    and $curr_res->{$tmp}{"TYPE"} !~ /H\d+/) {
	    $returnval = $curr_res->{$tmp}{"TYPE"};
	    last;
	}
    }

    return $returnval;

}

sub GetChain(@) {
    my ($curr_chain, $strand_length) = @_;
    my ($counter, $returnval, $largest_val);

    $counter = 0;
    $counter = int($curr_chain/($strand_length + 1)) + 1;
    $returnval = uc(chr($counter + 64));
    
    return $returnval;
}

sub WritePQRFile(@) {
    my ($Meso, $file_nm, $parms, $strand_len, $file_no) = @_;
    my ($res_counter, $out_string, $bead_counter, $fmt1, $index, $res_type);
    my ($tmp, %CONNECTS, %MYLIST, $my_key, $j_index, $search_key);
    my ($BEADS, $tot_res, $tot_res1, $b_counter, @R_T, $j, $i_index, $parent);
    my ($r_counter, $diff, $is_Pres, $bgf_string, $fmt2, $old_chain);
    my ($chain, $element, $num_bonds, $lone_pairs, $bgf_file, $muser, $mdate, $radii);

    $index = 1;
    $bgf_file = $file_nm;
    $bgf_file =~ s/\.\w+$//g;
    $bgf_file .= ".bgf";

    $fmt1 = "%-6s%5d  %4s %3s  %4d    %8.3f%8.3f%8.3f %5.1f %5.1f\n";
    $fmt2 = "%-6s %5d %5s %3s %1s%5d %10.5f%10.5f%10.5f %-5s%3d%2d %8.5f%2d%4d%10.5f\n";

    $muser = $ENV{USER};
    $mdate = (localtime)[3] . "/" .  (localtime)[4] . "/" . ((localtime)[5] + 1900);
    $mdate .= " at " .  sprintf("%2d", (localtime)[2]) . ":" .   sprintf("%2d", (localtime)[1]) . 
	":" . sprintf("%2d", (localtime)[0]);
    $bgf_string = "BIOGRF  332\nDESCRP MesoDNA\nREMARK Bead Model for Atomistic DNA\n";
    $bgf_string .= "REMARK Created by $muser on $mdate\nFORCEFIELD DREIDING\n";
    $bgf_string .= "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)\n";

    for $res_counter (sort Numerically keys %{ $Meso } ) {
	for $tmp (sort Numerically keys %{ $Meso->{$res_counter} } ) {

	    $bead_counter  = (\%{ $Meso->{$res_counter}{$tmp} });
	    $bead_counter->{"ATOM_NUM"} = $index;
	    $res_type = ExtractResidue(\% {$Meso->{$res_counter} });
	    $out_string .= sprintf($fmt1, "ATOM",$index,$bead_counter->{"TYPE"},$res_type,$res_counter,
				   $bead_counter->{"XCOORD"},$bead_counter->{"YCOORD"},$bead_counter->{"ZCOORD"},
				   $bead_counter->{"CHARGE"},$bead_counter->{"RADII"});
	    $chain = GetChain($res_counter, $strand_len);
	    if (! $old_chain) {
		$old_chain = $chain;
	    } elsif ($old_chain ne $chain) {
		$old_chain = $chain;
		$bgf_string .= "TER\n";
		$out_string .= "TER\n";
	    }
	    $element = $bead_counter->{"ELEMENT"};
	    $num_bonds = $bead_counter->{"NUMBONDS"};
	    $lone_pairs = $bead_counter->{"LONEPAIRS"};
	    $radii = $bead_counter->{"RADII"};
	    $bgf_string .= sprintf($fmt2, "HETATM",$index,$bead_counter->{"TYPE"},$res_type,$chain,$res_counter,
				   $bead_counter->{"XCOORD"},$bead_counter->{"YCOORD"},$bead_counter->{"ZCOORD"},
				   $element,$num_bonds,$lone_pairs,$bead_counter->{"CHARGE"},0,0,$radii);
	    $my_key = $res_counter . "-" . $tmp;
	    $MYLIST{$index} = $my_key;
	    $index++;
	}
    }


    $bgf_string .= "TER\n";
    $out_string .= "TER\n";

    $index = 1;
    for $r_counter (sort Numerically keys %{ $Meso } ) {

	$BEADS = ();
	push @{ $BEADS }, (sort Numerically keys %{ $Meso->{$r_counter} });
	$tot_res = $tot_res1 = $#{ $BEADS };
	for $b_counter (0 .. $tot_res) {
	    $R_T[$b_counter] = $r_counter;
	}
	
	if (! IsEndOfChain($r_counter, $strand_len)) {
	    
	    push @{ $BEADS }, (sort Numerically keys %{ $Meso->{($r_counter + 1)} });
	    $tot_res1 = $#{ $BEADS };
	    for $b_counter (($tot_res + 1) .. $#{ $BEADS }) {
		$R_T[$b_counter] = $r_counter + 1;
	    }
	}
	
	for $tmp (0 .. $tot_res) {
	    $i_index = $Meso->{$R_T[$tmp]}{$BEADS->[$tmp]}{"INDEX"};
	    for $j (0 .. $tot_res1) {
		$j_index = $Meso->{$R_T[$j]}{$BEADS->[$j]}{"INDEX"};
		$search_key = $i_index . "-" . $j_index;
		$diff = GetDiff($R_T[$tmp], $R_T[$j]);
		($is_Pres) = isPresent($search_key, \$parms->{"BONDS"}, $diff);
		if ($is_Pres and $i_index < 100 and $j_index < 100) {
		    $my_key = $R_T[$j] . "-" . $BEADS->[$j];
		    $CONNECTS{$index} .= GetIndex(\%MYLIST, $my_key) . " ";
		} elsif ($R_T[$tmp] == $R_T[$j]) {
		    $parent = $Meso->{$R_T[$j]}{$BEADS->[$j]}{"PARENT"};
		    if (defined($parent) and $parent == $i_index) {
			$my_key = $R_T[$j] . "-" . $BEADS->[$j];
			$CONNECTS{$index} .= GetIndex(\%MYLIST, $my_key) . " ";
		    }			
		}
	    }
	    if (! exists($CONNECTS{$index})) {
		$CONNECTS{$index} = " ";
	    }
	    $index++;

	}
    }

    $Connections->{$file_no} = FixConnections(\%CONNECTS);
#    $Connections = \%CONNECTS;
    $bgf_string .= "FORMAT CONECT (a6,14i6)\n";
    for $index (sort Numerically keys %{ $Connections->{$file_no} }) {
	$bgf_string .= sprintf("CONECT%6d", $index);
	$bgf_string .= $Connections->{$file_no}{$index};
	$bgf_string .= "\n";
    }

    open OUTFILE, "> $file_nm" or die "Cannot write to $file_nm: $!\n";
    print OUTFILE $out_string;
    close OUTFILE;

    open OUTFILE, "> $bgf_file" or die "Cannot write to $bgf_file: $!\n";
    print OUTFILE $bgf_string;
    close OUTFILE;

}

sub FixConnections(@) {
    my ($connect) = $_[0];
    my ($counter, @c_list, $returnval, $index, $i, $Connects, @tmp);

    for $counter (keys %{ $connect }) {
	$connect->{$counter} = Trim($connect->{$counter});
	if ($connect->{$counter} ne "") {
	    @c_list = split /\s+/, $connect->{$counter};
	    for $index (0 .. $#c_list) {
		push @{ $returnval->{$c_list[$index]} }, $counter;
		push @{ $returnval->{$counter} }, $c_list[$index];
	    }
	} else {
	    push @{ $returnval->{$counter} }, "";
	}
    }

    for $counter (keys %{ $returnval }) {
	next
	    if ($returnval->{$counter}[0] eq "");
	@{ $returnval->{$counter} } = sort Numerically @{ $returnval->{$counter} };
    }

    for $counter (keys %{ $returnval }) {
	for $index (0 .. $#{ $returnval->{$counter} }) {
	    for $i (($index + 1) .. $#{ $returnval->{$counter} }) {
		if ($returnval->{$counter}[$index] eq $returnval->{$counter}[$i]) {
		    $returnval->{$counter}[$i] = "";
		}
	    }
	}
    }

    for $counter (sort Numerically keys %{ $returnval }) {
	for $index (0 .. $#{ $returnval->{$counter} }) {
	    if ($returnval->{$counter}[$index] ne "") {
		$Connects->{$counter} .= sprintf("%6d", $returnval->{$counter}[$index]);
	    } else {
		$Connects->{$counter} = " ";
	    }
	}
	chop $Connects->{$counter}; 
    }

    return $Connects;
}

sub GetIndex(@) {
    my ($my_list, $my_key) = @_;
    my ($counter, $returnval);

    $returnval = "";
    for $counter (keys %{ $my_list }) {
	if ($my_list->{$counter} eq $my_key) {
	    $returnval = $counter;
	    last;
	}
    }

    return $returnval;

}

sub GetMass(@) {
    my ($atm_label) = $_[0];
    my ($returnval, $in_label);

    $returnval = 0;
    
    if ($atm_label =~ /^\d*([a-z]+)/i) {
	$in_label = $1;
   
	if ($PARMS->{"MASSES"}{$in_label}) {
	    $returnval = $PARMS->{"MASSES"}{$in_label};
	}
    }
    return $returnval;
}


sub IsEndOfChain(@) {
    my ($curr_chain, $s_len) = @_;

    if ($curr_chain % $s_len == 0) {
	return 1;
    } else {
	return 0;
    }
}

sub isPresent(@) {
    my ($search_string, $parms, $difference) = @_;
    my ($counter, $index, $isvalid, $i_loop, @index_list, @compare_list, $my_num, $was_changed);

    $isvalid = 1;
    
    @index_list = split /\-/, $search_string;
    
    $my_num = $index_list[0];
    
    for $counter (0 .. ($#index_list - 1)) {
	if ($index_list[$counter] eq $index_list[$counter + 1]) {
	    $isvalid = 0;
	    last;
	}
    }
    
    if (! $isvalid) {
	return 0;
    }
    
    $isvalid = $was_changed = 0;
    
    for $i_loop (0 .. 0) { # do not use symmetry
	for $counter (keys %{ $$parms }) {
	    $search_string = "";
	    $was_changed = 0;
	    @compare_list = split /\-/, $$parms->{$counter}{"VAL"};
	    for $index (0 .. $#compare_list) {
		if ($compare_list[$index] ne "X") {
		    if ($compare_list[$index] ne $index_list[$index]) {
			$was_changed = 1;
		    }
		}
		if ($was_changed) {
		    last;
		} else {
		    $search_string .= $compare_list[$index] . "-";
		}
	    }
	    if (! $was_changed) {
		if (IsMember($$parms->{$counter}{"RANGE"}, $difference)) {
		    chop $search_string;
		    $isvalid = 1;
		    last;
		}
	    }
	}
	if ($isvalid) {
	    last;
	}
	@index_list = reverse @index_list;
    }

    return ($isvalid);
}


sub IsMember(@) {
    my ($search_host, $search_test) = @_;
    my (@search_items, $returnval, $counter);

    $returnval = 0;
    $search_host =~ s/\-\w+//g;
    @search_items = split /\,/, $search_host;

    for $counter (@search_items) {
	if ($counter eq $search_test) {
	    $returnval = 1;
	    last;
	}
    }

    return $returnval;
}
    
	
sub GetName(@) { 
    my ($id_string) = $_[0];
    my ($return_string, $counter);
    my (@ids) = split /\-/, $id_string;

    for $counter (@ids) {
	$return_string .= $PARMS->{"BEADS"}{$counter}{"NAME"} . "-";
    }

    chop $return_string;

    return $return_string;
}

sub GetDiff(@) {
    my (@in_vals) = @_;
    my (@pairs, $counter, $returnval, $curr_diff);

    for $counter (0 .. ($#in_vals - 1)) {
	$curr_diff = ($in_vals[$counter + 1] - $in_vals[$counter]);
	push @pairs, $curr_diff;
    }

    for $counter (@pairs) {
	$returnval += $counter;
    }

    return $returnval;

}


sub IsIndex(@) {
    my ($b_index, $a_index) = @_;
    my (@b_items, $counter, $result);

    @b_items = split /\s+/, $b_index;

    $result = 0;
    for $counter (@b_items) {
	if ($counter == $a_index) {
	    $result = 1;
	    last;
	}
    }

    return $result;
}


sub GetConnections(@) {
    my ($con_str, $bead_info, $file_nm) = @_;
    my ($res_counter, $bead_counter, @C_Index, %Conn, $con_counter, $index, $curr_bead);

    $con_counter = 0;
    $con_str = Trim($con_str);
    @C_Index = split /\s+/, $con_str;
    for $index (@C_Index) {
      RESLOOP: for $res_counter (keys %{ $bead_info->[$file_nm] } ) {
	  for $bead_counter (keys %{ $bead_info->[$file_nm]{$res_counter} }) {
	      $curr_bead = \%{ $bead_info->[$file_nm]{$res_counter}{$bead_counter} };
	      if ($curr_bead->{"ATOM_NUM"} == $index) {
		  $Conn{$con_counter} = (
					 {
					     "BEAD_ID" => $bead_counter,
					     "RES"     => $res_counter,
					 }
					 );
		  $con_counter++;
		  last RESLOOP;
	      }
	  }
      }
    }
    
    return \%Conn;    	
}

sub ParsePDBFile(@) {
    my ($pdb_file) = $_[0];
    my ($Atom_Info) = GetPDBFileInfo($pdb_file, 1, " ", " ");
    my ($counter, $curr_res_name, $curr_res_id, $strand_indicator);

    $curr_res_id = 0;
    $strand_indicator = "";
    for $counter (sort Numerically keys %{ $Atom_Info }) {
	if ($Atom_Info->{$counter}{"RES_ID"} > $curr_res_id) {
	    $curr_res_id = $Atom_Info->{$counter}{"RES_ID"};
	    $curr_res_name = $Atom_Info->{$counter}{"RES_NAME"};
	    if ($curr_res_name =~ /(\d+)/) {
		if ($strand_indicator eq "") {
		    $strand_indicator = $1;
		} elsif ($strand_indicator ne $1 and ! $strand_length) {
		    $strand_length = $curr_res_id;
		}
	    }
	}
    }

    if (! $strand_length) {
	print "ERROR processing $pdb_file: CANNOT find strand_length\n";
    }
    return ($strand_length, $Atom_Info);

}


sub strVal(@) {
    my ($inStr) = $_[0];
    my (@chars, $returnval);

    @chars = split //, $inStr;

    for (@chars) {
	$returnval .= ord($_);
    }

    return $returnval;
}
