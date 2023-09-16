#!/usr/bin/perl -w
#  This script will open a bgf file, with an accompanying cerius2 force field 
#  file and create a forcefield file that is compatible with lammps
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}
use strict;
use Packages::General;

sub ValidateFile(@);
sub ParseCerius2File(@);
sub ParseBGFFile(@);
sub AddUnitLabel(@);
sub StoreExtrema(@);
sub CreateOutputHeader();
sub Numerically;
sub PrintAtoms(@);
sub PrintValence(@);
sub GenerateValence(@);
sub GenerateAngles(@);
sub GenerateTorsions(@);
sub CreateInputFile(@);
sub FindAngleTypes(@);
sub FindBondAtomTypes(@);

die "usage: $0 bgf_file cerius_ff [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgf_file, $ff_file, $save_name) = @ARGV;
my (%ERRORS, $counter, $index, $PARMS, $ATOMS, $CONNECTIONS);
my ($BONDS, $ANGLES, $VDW, $TORSIONS, $dat_file, $in_file);

$dat_file = "data." . $save_name;
$in_file = "in." . $save_name;

for ($ARGV[0], $ARGV[1]) {
    ValidateFile($_);
}


if (! $save_name) {
    $save_name = $bgf_file;
    $save_name =~ s/\.bgf$/\.in/;
}

($PARMS) = ParseCerius2File($ff_file);
($ATOMS, $CONNECTIONS) = ParseBGFFile($bgf_file, $PARMS);

$ATOMS = AddUnitLabel($ATOMS, $PARMS);
$BONDS = GenerateValence($CONNECTIONS, "BONDS", 1);
$ANGLES = GenerateValence($CONNECTIONS, "ANGLES", 2);
$TORSIONS = GenerateValence($CONNECTIONS, "TORSIONS", 3);

open OUTFILE, "> $dat_file" or die "Cannot write to file $dat_file: $!\n";
print OUTFILE "LAMMPS Description\n\n";
CreateOutputHeader();
PrintAtoms($ATOMS);
PrintValence("Bonds", $BONDS);
PrintValence("Angles", $ANGLES);
PrintValence("Dihedrals", $TORSIONS);
close OUTFILE;

open OUTFILE, "> $in_file" or die "Cannot write to file $in_file: $!\n";
CreateInputFile($PARMS, $dat_file);
close OUTFILE;
for $counter (keys %ERRORS) {
    print "\n=====$counter ERRORS ====\n";
    for $index(keys %{ $ERRORS{$counter} }) {
	print "$index: " . $ERRORS{$counter}{$index} . "\n";
    }
}


sub ValidateFile(@) {
    my ($in_file) = $_[0];

    die "Error accessing file $in_file: $!\n"
	if (! -e $in_file or ! -r $in_file or ! -T $in_file);

}

sub ParseBGFFile(@) {
    my ($in_file, $PAR) = @_;
    my ($ATOMS, $BONDS, @CON, $counter, $in_data, $atm_patern, $total_atoms);

    $total_atoms = 0;
    $atm_patern = '^ATOM\s+(\d+)\s+(\w+)\s+(\w+)\s\w?\s+(\d+)\s+' . 
	'(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\w+)' .
	'\s+(\d+)\s+(\d+)\s+(\d+\.\d+)';  
    open BGFFILE, $in_file or die "Cannot open BGFfile $in_file: $!\n";
    while (<BGFFILE>) {
	chomp;
	$in_data = $_;
	$in_data =~ s/HETATM/ATOM/;
	if ($in_data =~ /$atm_patern/) {
	    $ATOMS->{$1} = (
			    {
				"ATMNAME"     => $8,
				"RESNAME"     => $3,
				"RESNUM"      => $4,
				"XCOORD"      => $5,
				"YCOORD"      => $6,
				"ZCOORD"      => $7,
				"REALATM"     => $8,
				"NUMBONDS"    => $9,
				"LONEPAIRS"   => $10,
				"CHARGE"      => $11,
				"OCCUPANCY"   => 0,
				"RESONANCE"   => 0,
				"RADII"       => 0,
			    }
			    );
	    StoreExtrema($5, $6, $7, $PAR);
	    $total_atoms++;
	} else {
	    if ($_ =~ /^CONECT\s+(.+)$/) {
		@CON = split /\s+/, $1;
		for $counter (1 .. $#CON) {
		    push @{ $BONDS->{$CON[0]} }, $CON[$counter];
		}
	    }
	}
    }
    close BGFFILE;

    die "ERROR: $in_file does not contain any ATOM/CONNECTION information!\n"
	if (! defined($BONDS));
    die "ERROR: IT's The Connections stupid!\n" 
	if (! defined ($ATOMS));

    print "TOTAL ATOMS FOUNDS IN FILE: $total_atoms\n";
    return ($ATOMS, $BONDS);
}

sub  ParseCerius2File(@) {
    my ($ff_file) = $_[0];
    my (%PARMS, $which_var, $in_data, $type_counter, @vdws);
    my ($type_id, $type_id2, @bonds, $type_id3, @angles);
    my ($type_id4, @torsions, $key_code);
    my ($bond_counter, $angle_counter, $torsion_counter);

    $which_var = $bond_counter = $angle_counter = $torsion_counter = $counter = 0;
    open FORCEFIELD, $ff_file or die "Cannot open force field file $ff_file: $!\n";
    while (<FORCEFIELD>) {
	chomp;
	$in_data = $_;
	if ($in_data =~ /^ATOMTYPES/) {
	    $which_var = 1;
	} elsif ($in_data =~ /HYDROGEN_BONDS/) {
		$which_var = 2;	    
	} elsif ($in_data =~ /^DIAGONAL_VDW/) {
		$which_var = 3;
	} elsif ($in_data =~ /^BOND_STRETCH/) {
	    $which_var = 4;
	} elsif ($in_data =~ /^ANGLE_BEND/) {
	    $which_var = 5;
	} elsif ($in_data =~ /^TORSIONS/) {
	    $which_var = 6;
	} elsif ($in_data =~ /^\s*(\w+)\s+(\w+)\s+(\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ and ($which_var == 1)) {
	    $type_counter += 1;
	    $PARMS{"ATOMTYPES"}{$1} = (
				       {
					   "TYPEID"    => $type_counter,
					   "ATOM"      => $2,
					   "MASS"      => $3,
					   "CHARGE"    => $4,
					   "NUMBONDS"  => $5,
					   "LONEPAIRS" => $6,
					   "OTHER"     => $7,
				       }
				       );
	    print "ATOMTYPES: $1: $type_counter\n";
	} elsif ($in_data =~ /^\s*(\w+)\s+(\w+)\s+(\w+)\s+(.+)/ and ($which_var == 2)) {
	    if (uc($2) ne "IGNORE") {
		$type_id = $PARMS{"ATOMTYPES"}{$1}{"TYPEID"};
		next
		    if (! defined($type_id));
		$type_id2 = $PARMS{"ATOMTYPES"}{"$2"}{"TYPEID"};
		next
		    if (! defined($type_id2));
		@vdws = split /\s+/,$4;
		$type_id3 = $3;
		$PARMS{"HBOND"}{$type_id . " " . $type_id2} = (
				     {
					 "TYPE"   => $type_id3,
					 "COMP"   => $type_id2,
					 "VALS"   => [@vdws],
				     }
				     );
		print "HBOND: $1 - $2\n";
	    }
	} elsif ($in_data =~ /^\s*(\w+)\s+(\w+)\s+(\d+\.\d+)\s*(.+)/ and ($which_var == 3)) {
	    if (uc($2) ne "IGNORE") {
		$type_id = $PARMS{"ATOMTYPES"}{$1}{"TYPEID"};
		next
		    if (! defined($type_id));
		@vdws = split /\s+/,$4;
		$type_id3 = $2;
		$type_id2 = ""; 
		$PARMS{"VDW"}{$type_id} = (
				     {
					 "TYPE"   => $type_id3,
					 "COMP"   => "",
					 "VALS"   => [$3, @vdws],
				     }
				     );
		print "VDW: $1: $type_id\n"
	    }
	} elsif ($in_data =~ /^\s*(\w+)\s+(\w+)\s+(\w+)\s+(.+)/ and ($which_var == 4)) {
	    $bond_counter++;
	    $type_id = $PARMS{"ATOMTYPES"}{$1}{"TYPEID"};
	    next
	        if (! defined($type_id));
		$type_id2 = $PARMS{"ATOMTYPES"}{$2}{"TYPEID"};
		next
		    if (! defined($type_id2));
		
		@bonds = split /\s+/, $4;
		$key_code = "$type_id - $type_id2";
	    if (uc($3) ne "IGNORE") {
		$PARMS{"BONDS"}{$bond_counter} = (
					      {
						  "KEY"      => $key_code,
						  "TYPE"     => $3,
						  "VALS"     => [@bonds],
					      }
					      );
	    } else {
		$PARMS{"BONDS"}{$bond_counter} = (
                                              {                                                                     
                                                  "KEY"      => $key_code,                                          
                                                  "TYPE"     => $3,                                                 
                                                  "VALS"     => ([0]),                                           
                                              }                                                                     
                                              );                                                                    
	    }
		print "BOND $bond_counter: $key_code\n";
	    
	} elsif ($in_data =~ /^\s*(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(.+)/ and ($which_var == 5)) {
		$angle_counter++;
		$type_id = $PARMS{"ATOMTYPES"}{$1}{"TYPEID"};
		next
		    if (! defined($type_id));
		$type_id2 = $PARMS{"ATOMTYPES"}{$2}{"TYPEID"};
		next
		    if (! defined($type_id2));
		$type_id3 = $PARMS{"ATOMTYPES"}{$3}{"TYPEID"};
		next
		    if (! defined($type_id3));
		
		@angles = split /\s+/, $5;
		$key_code = "$type_id - $type_id2 - $type_id3";
		if (uc($4) ne "IGNORE") {
		$PARMS{"ANGLES"}{$angle_counter} = (
					       {
						   "KEY"      => $key_code,
						   "TYPE"     => $4,
						   "VALS"     => [@angles],
					       }
					       );
		} else {
                $PARMS{"ANGLES"}{$angle_counter} = (
                                               {
                                                   "KEY"      => $key_code,
                                                   "TYPE"     => $4,
                                                   "VALS"     => ([0]),
                                               }
                                               );
		}
		print "ANGLE $angle_counter: $key_code\n";
	} elsif ($in_data =~ /^\s*(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(.+)/ and ($which_var == 6)) {
	    if (uc($5) ne "IGNORE") {
		$torsion_counter++;
		$type_id = $PARMS{"ATOMTYPES"}{$1}{"TYPEID"};
		next
		    if (! defined($type_id));
		$type_id2 = $PARMS{"ATOMTYPES"}{$2}{"TYPEID"};
		next
		    if (! defined($type_id2));
		$type_id3 = $PARMS{"ATOMTYPES"}{$3}{"TYPEID"};
		next
		    if (! defined($type_id3));
		$type_id4 = $PARMS{"ATOMTYPES"}{$4}{"TYPEID"};
		next
		    if (! defined($type_id4));
		
		@torsions = split /\s+/, $6;
		$key_code = "$type_id - $type_id2 - $type_id3 - $type_id4";
		if (uc($5) ne "IGNORE") {
		$PARMS{"TORSIONS"}{$torsion_counter} = (
					       {
						   "TYPE"     => $5,
						   "KEY"      => $key_code,
						   "VALS"     => [@torsions],
					       }
					       );
		} else {
                $PARMS{"TORSIONS"}{$torsion_counter} = (                                                            
                                               {                                                                    
                                                   "TYPE"     => $5,                                                
                                                   "KEY"      => $key_code,                                         
                                                   "VALS"     => ([0]),                                       
                                               }                                                                    
                                               );                                                                   
		}
		print "TORSION $torsion_counter: $key_code\n";
	    } else {
		$torsion_counter++;
		$PARMS{"TORSIONS"}{$torsion_counter}{"TYPE"} = "IGNORE";
	    }
	}
		
	    
    }
    close FORCEFIELD;

    die "ERROR: Force Field file $ff_file is invalid!\n"
	if (! %PARMS);

    $PARMS{"ATOMTYPES"}{"Lammps"}{"LABEL"} = "atom";
    $PARMS{"BONDS"}{"Lammps"}{"LABEL"} = "bond";
    $PARMS{"ANGLES"}{"Lammps"}{"LABEL"} = "angle";
    $PARMS{"TORSIONS"}{"Lammps"}{"LABEL"} = "dihedral";
    $PARMS{"VDW"}{"Lammps"}{"LABEL"} = "pair";
    
    return (\%PARMS);
}

sub AddUnitLabel(@) {
    my ($atoms, $parms) = @_;

    my ($counter, $type_id, $atm_name);

    for $counter (keys %{ $atoms } ) {
	$atm_name = $atoms->{$counter}{"ATMNAME"};
	$type_id = $parms->{"ATOMTYPES"}{$atm_name}{"TYPEID"};
	if (! defined($type_id) or ! $type_id) {
	    if (! $ERRORS{"ATOMTYPES"}{$atm_name}) {
		$ERRORS{"ATOMTYPES"}{$atm_name} = 1;
	    } else {
		$ERRORS{"ATOMTYPES"}{$atm_name} += 1;
	    }
	} else {
	    $atoms->{$counter}{"TYPEID"} = $type_id;
	}
    }

    return ($atoms);
}

sub StoreExtrema(@) {
    my ($x, $y, $z, $par) = @_;
    my ($counter);

    if (! defined($par->{"EXTREMA"})) {
	for $counter ("X", "Y", "Z") {
	    $par->{"EXTREMA"}{$counter} = (
					   {
					       "lo" => 99999.999,
					       "hi" => -99999.999,
					   }
					   );
	}
    } else {
	$par->{"EXTREMA"}{"X"}{"hi"} = $x
	    if ($x > $par->{"EXTREMA"}{"X"}{"hi"});
	$par->{"EXTREMA"}{"X"}{"lo"} = $x
	    if ($x < $par->{"EXTREMA"}{"X"}{"lo"});
	$par->{"EXTREMA"}{"Y"}{"hi"} = $y
	    if ($y > $par->{"EXTREMA"}{"Y"}{"hi"});
	$par->{"EXTREMA"}{"Y"}{"lo"} = $y
	    if ($y < $par->{"EXTREMA"}{"Y"}{"lo"});
	$par->{"EXTREMA"}{"Z"}{"hi"} = $z
	    if ($z > $par->{"EXTREMA"}{"Z"}{"hi"});
	$par->{"EXTREMA"}{"Z"}{"lo"} = $z
	    if ($z < $par->{"EXTREMA"}{"Z"}{"lo"});
    }
}

sub CreateOutputHeader() {
    my ($tot, $index, $tmp, @holder, $counter);
   
    @holder = keys %{ $ATOMS };
    printf OUTFILE "%12s  atoms\n", ($#holder + 1);
    @holder = keys %{ $BONDS };
    printf OUTFILE "%12s  bonds\n", ($#holder + 1);
    @holder = keys %{ $ANGLES };
    printf OUTFILE "%12s  angles\n", ($#holder + 1);
    @holder = keys %{ $TORSIONS };
    printf OUTFILE "%12s  dihedrals\n", ($#holder + 1);
    printf OUTFILE "%12s  impropers\n\n", 0;

    for $index ("ATOMTYPES", "BONDS", "ANGLES", "TORSIONS") {
	@holder = ();
	@holder = keys %{ $PARMS->{$index} };
	$tmp = $PARMS->{$index}{"Lammps"}{"LABEL"};
	printf OUTFILE "%12s  $tmp types\n", $#holder;
    }
    printf OUTFILE "%12s  improper types\n\n", 0;

    for $index (keys %{ $PARMS->{"EXTREMA"} }) {
	printf OUTFILE "%10.6f %10.6f ", $PARMS->{"EXTREMA"}{$index}{"lo"}, 
	$PARMS->{"EXTREMA"}{$index}{"hi"};
	print OUTFILE lc($index) . "lo " . lc($index) . "hi\n";
    }

    print OUTFILE "\nMasses\n\n";

    $tmp = ();
    for $index (keys %{ $PARMS->{"ATOMTYPES"} }) {
	next
	    if ($index eq "Lammps");
	$tmp->{$PARMS->{"ATOMTYPES"}{$index}{"TYPEID"}} = $PARMS->{"ATOMTYPES"}{$index}{"MASS"};
    }

    for $index (sort Numerically keys %{ $tmp }) {
	printf OUTFILE "%3d %8.4f\n", $index, $tmp->{$index};
    }

    for $counter ("BONDS", "ANGLES", "TORSIONS", "VDW") {
	$tmp = "";
	print OUTFILE "\n" . ucfirst($PARMS->{$counter}{"Lammps"}{"LABEL"}) . " Coeffs\n\n";
	@holder = keys %{ $PARMS->{$counter} };
	$tot = 0;
	while ($tot <= $#holder) {
	    if ($holder[$tot] =~ /[a_zA_Z]/) {
		splice @holder, $tot, 1;
	    } else {
		$tot++;
	    }
	}
	for $index (sort Numerically @holder) {
	    $tmp = sprintf("%3d", $index);
	    for (@{ $PARMS->{$counter}{$index}{"VALS"} }) {
		$tmp .= sprintf("%12.6f", $_);
	    }
	    print OUTFILE "$tmp\n";
	}
    }

    print OUTFILE "\n";

}

sub PrintAtoms(@) {
    my ($atoms) = $_[0];

    my ($counter, $type_id, $atm_name, $fmt, $out_string);
    $fmt = "%8d%8d%8d%10.4f%15.6f%15.6f%15.6f\n";

    print OUTFILE "Atoms\n\n";
    for $counter (sort Numerically keys %{ $atoms } ) {
	$out_string = sprintf($fmt, $counter, 1, $atoms->{$counter}{"TYPEID"}, $atoms->{$counter}{"CHARGE"},
			      $atoms->{$counter}{"XCOORD"}, $atoms->{$counter}{"YCOORD"}, $atoms->{$counter}{"ZCOORD"});
	print OUTFILE $out_string;
    }

    print OUTFILE "\n";

}

sub Numerically {
    ($a<=>$b);
}

sub GenerateValence(@) {
    my ($connections, $which_valence, $depth) = @_;
    my ($counter, $key_code, $bond_index, $result, %OUTPUT);
    my ($valence_counter, $error_code, $atm_code, $iter);
    my ($angle_index, $torsion_index);

    $valence_counter = 0;
    $iter = 1;
    for $counter (sort Numerically keys %{ $connections }) {
	for $bond_index (@{ $connections->{$counter} }) {
	    next
		if ($bond_index < $counter and $depth == $iter);
	    $atm_code = "$counter - $bond_index";
	    $key_code = $ATOMS->{$counter}{"TYPEID"} . " - " . $ATOMS->{$bond_index}{"TYPEID"};
	    $error_code = $ATOMS->{$counter}{"ATMNAME"} . " - " . $ATOMS->{$bond_index}{"ATMNAME"};
	    if ($iter == $depth) {
		($result, $valence_counter) = FindKeyCode(\%{ $PARMS->{$which_valence} }, $key_code, $atm_code, 
							  $valence_counter, $which_valence, $error_code);
		if ($result ne "") {
		    $OUTPUT{$valence_counter} = $result;
		}
	    } else { # At least it's an angle, so go one level deeper
		for $angle_index (@{ $connections->{$bond_index} }) {
		    next
			if ($angle_index == $bond_index or $angle_index == $counter);
		    $atm_code = "$counter - $bond_index - $angle_index";
		    $key_code = $ATOMS->{$counter}{"TYPEID"} . " - " . $ATOMS->{$bond_index}{"TYPEID"} . " - " . $ATOMS->{$angle_index}{"TYPEID"};
		    $error_code = $ATOMS->{$counter}{"ATMNAME"} . " - " . $ATOMS->{$bond_index}{"ATMNAME"} . " - " . $ATOMS->{$angle_index}{"ATMNAME"};
		    if (($iter + 1) == $depth) {
			($result, $valence_counter) = FindKeyCode(\%{ $PARMS->{$which_valence} }, $key_code, $atm_code, 
								  $valence_counter, $which_valence, $error_code);
			if ($result ne "") {
			    $OUTPUT{$valence_counter} = $result;
			}
			
		    } else { # Finally, it's a torsion 
			for $torsion_index (@{ $connections->{$angle_index} }) {
			    next
				if ($torsion_index == $angle_index or $torsion_index == $bond_index or $torsion_index == $counter);
			    $atm_code = "$counter - $bond_index - $angle_index - $torsion_index";
			    $key_code = $ATOMS->{$counter}{"TYPEID"} . " - " . $ATOMS->{$bond_index}{"TYPEID"} . 
				" - " . $ATOMS->{$angle_index}{"TYPEID"} . " - " . $ATOMS->{$torsion_index}{"TYPEID"};
			    $error_code = $ATOMS->{$counter}{"ATMNAME"} . " - " . $ATOMS->{$bond_index}{"ATMNAME"} . 
				" - " . $ATOMS->{$angle_index}{"ATMNAME"} . " - " . $ATOMS->{$torsion_index}{"ATMNAME"};
			    ($result, $valence_counter) = FindKeyCode(\%{ $PARMS->{$which_valence} }, $key_code, $atm_code, 
								      $valence_counter, $which_valence, $error_code);
			    if ($result ne "") {
				$OUTPUT{$valence_counter} = $result;
			    }
			}
			
		    }
		}
		
	    }
	}
    }
    
    die "ERROR: $which_valence in structure does not correspond to any bonds in the force field!\n"
	if (! %OUTPUT);

    return (\%OUTPUT);
}

sub FindKeyCode(@) {
    my ($valence_object, $key_code, $atm_code, $index, $which_valence, $error_code) = @_;
    my ($counter, $result, $i);
    my (@IDS) = split /\s+\-\s+/, $key_code;
    my (@ATOMS) = split /\s+\-\s+/, $atm_code;
    $i = 0;

    $result = "";
    while ($i < 2) {
	for $counter (keys %{ $valence_object }) {
	    next
		if ($counter eq "Lammps");
	    if ($valence_object->{$counter}{"KEY"} eq $key_code) {
		$index++;
		$result = sprintf("%6d%4d", $index, $counter);
		for (@ATOMS) {
		    $result .= sprintf("%7d", $_);
		}
		last;
	    }
	}
	$key_code = "";
	for $counter (reverse @IDS) {
	    $key_code .= "$counter - ";
	}
	$key_code = substr($key_code, 0, (length($key_code)-3));
	$i++;
    }

    if ($result eq "") {
	if (! $ERRORS{$which_valence}{$error_code}) {
	    $ERRORS{$which_valence}{$error_code} = 1;
	} else {
	    $ERRORS{$which_valence}{$error_code} += 1;
	}
    }

    return ($result, $index);
}

sub PrintValence(@) {
    my ($title, $which_obj) = @_;
    my ($counter);

    print OUTFILE "$title\n\n";
    for $counter (sort Numerically keys %{ $which_obj } ) {
	print OUTFILE $which_obj->{$counter} . "\n";
    }
    print OUTFILE "\n";
}

sub CreateInputFile(@) {
    my ($parm, $datafile) = @_;
    my ($counter, @out_array, $pair_style, $bond_style, $dihedral_style, $HBONDS);
    my ($angle_style, $shake, @hb_data, $index, $tmp, $hash_key, $shake_angle);
    my (%VDWLIST, %ANGLELIST, %BONDLIST, %TORSIONLIST);

    $shake = " ";
    @out_array = keys %{ $parm->{"VDW"} };
    for $hash_key (@out_array) {
    next
	if (! defined($parm->{"VDW"}{$hash_key}{"TYPE"}));
    if ($parm->{"VDW"}{$hash_key}{"TYPE"} =~ /MORSE/i) {
	$pair_style = "morse";
    } elsif ($parm->{"VDW"}{$hash_key}{"TYPE"} =~ /LJ/i) {
	$pair_style = "lj";
    } else {
	$pair_style = $parm->{"VDW"}{$hash_key}{"TYPE"};
    }
    $VDWLIST{$pair_style} .= " $hash_key";
    }

#    @out_array = keys %{ $parm->{"HBOND"} };
#    for $hash_key (@out_array) {
#	if (defined($parm->{"HBOND"}{$hash_key}{"TYPE"})) {
#	    if ($hash_key =~ /(\d+)\s+(\d+)/) {
#		$HBONDS->{$1} = 1;
#		$HBONDS->{$2} = 1;
#	    }
#	}
#    }
    
#    $hash_key  = $out_array[0];
#    if (defined($HBONDS)) {
#	$pair_style .= "/hb_" . $parm->{"HBOND"}{$hash_key}{"TYPE"};
#	$shake = "b" . FindBondAtomTypes(\%{ $parm->{"BONDS"} }, $HBONDS);
#	$shake_angle = " a" . FindAngleTypes(\%{ $parm->{"ANGLES"}}, $HBONDS);
#    }
	
    @out_array = keys %{ $parm->{"BONDS"} };
    for $hash_key (@out_array) {
    next
        if (! defined($parm->{"BONDS"}{$hash_key}{"TYPE"}));

    if ($parm->{"BONDS"}{$hash_key}{"TYPE"} =~ /HARMONIC/i) {
	$bond_style = "harmonic";
    } else {
	$bond_style = $parm->{"BONDS"}{$hash_key}{"TYPE"};
    }
    $BONDLIST{$bond_style} .= " $hash_key";
    }

    @out_array = keys %{ $parm->{"ANGLES"} };
    for $hash_key (@out_array) {
    next
        if (! defined($parm->{"ANGLES"}{$hash_key}{"TYPE"}));

    if ($parm->{"ANGLES"}{$hash_key}{"TYPE"} =~ /HARM/i) {
	$angle_style = "harmonic";
    } else {
	$angle_style = $parm->{"ANGLES"}{$hash_key}{"TYPE"};
    }
    $ANGLELIST{$angle_style} .= " $hash_key";
    }

    @out_array = keys %{ $parm->{"TORSIONS"} };
    for $hash_key (@out_array) {
    next
        if (! defined($parm->{"TORSIONS"}{$hash_key}{"TYPE"}));

    if ($parm->{"TORSIONS"}{$hash_key}{"TYPE"} =~ /DIHEDRAL/i) {
	$dihedral_style = "harmonic";
    } else {
	$dihedral_style = $parm->{"TORSIONS"}{$hash_key}{"TYPE"};
    }
    $TORSIONLIST{$dihedral_style} .= " $hash_key";
    }

    print "";
    @out_array = (
		  ["units", "real"],
		  ["atom_style", "full"],
		  ["boundary", "p p p"],
		  ["", ""],
		);
    for $hash_key (keys %VDWLIST) {
	push @out_array, (["pair_style", $hash_key . $VDWLIST{$hash_key}]);
    }
    for $hash_key (keys %BONDLIST) {
        push @out_array, (["bond_style", $hash_key  . $BONDLIST{$hash_key}]);
    }
    for $hash_key (keys %ANGLELIST) {
        push @out_array, (["angle_style", $hash_key  . $ANGLELIST{$hash_key}]);
    }
    for $hash_key (keys %TORSIONLIST) {
        push @out_array, (["dihedral_style", $hash_key  . $TORSIONLIST{$hash_key}]);
    }
    push @out_array, (
		  ["improper_style", "none"],
		  ["kspace_style", "none"],
		  ["",""],
		  ["read_data", $datafile],
		  ["",""],
		  ["neighbor", "2.0 bin"],
		  ["neigh_modify", "every 2 delay 5 check yes"],
		  ["",""],
		  ["timestep", "2.0"],
		  ["",""],
		  ["thermo_style", "multi"],
		  ["thermo", "50"],
		  ["",""],
		  ["fix", "1 all nvt 275.0 275.0 100.0"],
		  ["",""],
		  ["dump", "1 all atom 10 out.trj"],
		  ["",""],
		  ["run", "300"],
		  );

# This is a hack of a hack!! Have to fix hydrogen bonds sometime!    
    for $counter (0 .. $#out_array) {
	printf OUTFILE "%-15s", $out_array[$counter][0];
	print OUTFILE $out_array[$counter][1] . "\n";
#	if ($counter == 11) {
#	    if ($shake ne " ") {
#		for $index (@hb_data) {
#		    printf OUTFILE "%-15s", "pair_coeff";
#		    print OUTFILE $index . "\n";
#		}
#	    }
#	} elsif ($counter == 21) {
#	    if ($shake ne " ") {
#		$shake = Trim($shake);
#		printf OUTFILE "%-15s", "fix";
#		print OUTFILE"2 all shake 0.0001 20 10 " . $shake;
#		if (defined($shake_angle)) {
#		    print OUTFILE $shake_angle;
#		}
#		    
#		print OUTFILE "\n";
#	    }
#	}
    }
}

sub FindBondAtomTypes(@) {
    my ($bond_list, $hbonds) = @_;
    my ($found_bonds, $counter, $returnval, $curr_angle);

    $returnval = " ";
    for $counter (keys %{ $bond_list }) {
	next
	    if (! defined($bond_list->{$counter}{"KEY"}));
	if ($bond_list->{$counter}{"KEY"} =~ /(\d+)\s+\-\s+(\d+)/) {
	    if (defined($hbonds->{$1}) or defined($hbonds->{$2})) {
		$found_bonds->{$counter} = 1;
	    }
	}
    }

    for $counter (keys %{ $found_bonds }) {
	$returnval .= "$counter ";
    }

    chop $returnval;
    return $returnval;
}

sub FindAngleTypes(@) {
    my ($angle_list, $hbonds) = @_;
    my ($found_angles, $counter, $returnval, $curr_angle);

    $returnval = " ";
    for $counter (keys %{ $angle_list }) {
	next
	    if (! defined($angle_list->{$counter}{"KEY"}));
	if ($angle_list->{$counter}{"KEY"} =~ /(\d+)\s+\-\s+(\d+)\s+\-\s+(\d+)/) {
	    if (defined($hbonds->{$1}) and defined($hbonds->{$3})) {
		$found_angles->{$counter} = 1;
	    }
	}
    }

    for $counter (keys %{ $found_angles }) {
	$returnval .= "$counter ";
    }

    chop $returnval;
    return $returnval;
}
