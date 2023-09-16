#!/usr/bin/perl -w  
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester);

my ($lmp_data_file, $pdb_file, $saveName);

$|++;
&init;
&CharmmCmap($lmp_data_file, $saveName);

sub init {
	my (%OPTS);

	getopt('dps',\%OPTS);
	($lmp_data_file, $pdb_file, $saveName) = ($OPTS{d}, $OPTS{p}, $OPTS{s});
	&usage if (! exists($OPTS{d}) or ! exists($OPTS{p}));

	print "Initializing...";
	FileTester($lmp_data_file);
	if (! defined($saveName)) {
		$saveName = basename($lmp_data_file);
		$saveName =~ s/^\w+\./cmap\./;
	}	
	print "Done\n";
}

sub usage {
	die "usage: $0 -d lmp_data_file -p pdb_file -s (savename)\n";
}

  # ----------------------- DESCRIPTION: sub CharmmCmap ------------------------ #
  # This subroutine add a new section "CMAP" to the LAMMPS data file, which      #
  # a part of the implementation of "CHARMM CMAP" (see references) in LAMMPS.    #
  # The section "CMAP" contains a list of dihedral ID pairs from adjecent        #
  # peptide backtone dihedrals whose dihedral angles are corrresponding to PHI   #
  # and PSI. (PHI: C--N--C_aphla_C and PSI: N--C_alpha--C--N)                    #
  #                                                                              #
  # Initiated by:   Xiaohu Hu (hux2@ornl.gov)                                    #
  # May 2009                                                                     #
  #                                                                              #
  # Finalized Oct 2016 by: Robert Latour (latourr@clemson.edu),                  #
  #   David Hyde-Volpe, and Tigran Abramyan, Clemson University,                 #
  #   and Chris Lorenz (chris.lorenz@kcl.ac.uk                                   #
  #                                                                              #
  # References:                                                                  #
  # - MacKerell, A.D., Jr., Feig, M., Brooks, C.L., III, Improved Treatment of   #
  #   Protein Backbone Conformational in Empirical Force Fields, J. Am. Chem.    #
  #   Soc. 126(2004): 698-699.                                                   #
  # - MacKerell, A.D., Jr., Feig, M., Brooks, C.L., III, Extending the Treatment #
  #   of Backbone Energetics in Protein Force Fields: Limitations of Gas-Phase   #
  #   Quantum Mechnacis in Reproducing Protein Conformational Distributions in   #
  #   Molecular Dynamics Simulations, J. Comput. Chem. 25(2004): 1400-1415.      #
  # ---------------------------------------------------------------------------- #

sub CharmmCmap
 {
	my ($data_file, $sname) = @_;

	print "\nINITIATING CHARMM CMAP SUBROUTINE...\n\n";

	# Reread and analyse $lmp_data_file
	my @raw_data;
	open(LAMMPS, "< $data_file") or
	die "\"sub CharmmCmap()\" cannot open \"$data_file!\n";
	print "Analyzing \"$data_file\"...\n\n";
	@raw_data = <LAMMPS>;
	close(LAMMPS);

	# Locate and extract the sections "Masses" and "Atoms"
	my $line_number = 0;
	# Header infos, 0 by default
	my $natom_types = 0;
	my $natom_number = 0;
	my $ndihedral_number = 0;
	my $temp_string;
	# splice points, 0 by default
	my $splice_onset_masses = 0;
	my $splice_onset_atoms = 0;
	my $splice_onset_dihedrals = 0;

	foreach my $line (@raw_data) {
		$line_number++;
		chomp($line);
	# Extract useful informations from the header
		if ($line =~ m/atom types/) {
			($natom_types,$temp_string) = split(" ",$line);
			if ($natom_types == 0) {
				die "\nError: Number of atom types is 0!\n";
			}
			print "Total atom types: $natom_types\n";
		}
		if ($line =~ m/atoms/) {
			($natom_number,$temp_string) = split(" ",$line);
			if ($natom_number == 0) {
				die "\nError: Number of atoms is 0!\n";
			}
			print "Total atoms: $natom_number\n";
		}
		if ($line =~ m/dihedrals/) {
			($ndihedral_number,$temp_string) = split(" ",$line);
			if ($ndihedral_number == 0) {
				die "\nError: Number of dihedrals is 0\n";
			}
			print "Total dihedrals: $ndihedral_number\n";
		}
	# Locate and data from sections "Masses", "Atoms" and "Dihedrals"
		if ($line =~ m/Masses/) {
			$splice_onset_masses = $line_number + 1;
			if ($splice_onset_masses-1 == 0) {
				die "\nError: Can not find the section \"Masses\"\n";
			}
			print "Section \"Masses\" found: line $splice_onset_masses\n";
		}
		if ($line =~ m/Atoms/) {
			$splice_onset_atoms = $line_number +1;
			if ($splice_onset_atoms-1 == 0) {
				die "\nError: Can not find the section \"Atoms\"\n";
			}
			print "Section \"Atoms\" found: line $splice_onset_atoms\n";
		}
		if ($line =~ m/Dihedrals/) {
			$splice_onset_dihedrals = $line_number + 1;
			if ($splice_onset_dihedrals-1 == 0) {
				die "\nError: Can not find the section \"Dihedrals\"\n";
			}
			print "Section \"Dihedrals\" found: line $splice_onset_dihedrals\n";
		}
	}

	print "\nGenerating PHI/PSI dihedral pair list...\n\n";

	my @temp1 = @raw_data;
	my @temp2 = @raw_data;
	my @temp3 = @raw_data;

	# Extract the section "Masses", "Atoms" and "Dihedrals"
	my @temp_masses_data = splice(@temp1,$splice_onset_masses,$natom_types);
	my @temp_atoms_data = splice(@temp2,$splice_onset_atoms,$natom_number);
	my @temp_dihedrals_data = splice(@temp3,$splice_onset_dihedrals,$ndihedral_number);
	
	# Store @temp_masses_dat into a matrix
	my @masses_matrix;
	my $atom_type;
	my $mass;
	for (@temp_masses_data) {
		($atom_type, $mass) = split(" ");
		push(@masses_matrix,[$atom_type,$mass]);
	}

	# Store @temp_atoms_data into a matrix
	my @atoms_matrix;
	my $atom_ID;
	my $molecule_ID;
	my $atype;
	my $charge;
	my $atom_x_coor;
	my $atom_y_coor;
	my $atom_z_coor;
	for (@temp_atoms_data) {
		($atom_ID,$molecule_ID,$atype,$charge,$atom_x_coor,$atom_y_coor,$atom_z_coor) = split(" ");
	 	push(@atoms_matrix,
			[$atom_ID,$molecule_ID,$atype,$charge,$atom_x_coor,$atom_y_coor,$atom_z_coor]);
	}

	# Store @temp_dihedrals_data into a matrix
	my @dihedrals_matrix;
	my $dihedral_ID;
	my $dihedral_type;
	my $dihe_atom1;
	my $dihe_atom2;
	my $dihe_atom3;
	my $dihe_atom4;
	for (@temp_dihedrals_data) {
		($dihedral_ID,$dihedral_type,$dihe_atom1,$dihe_atom2,$dihe_atom3,$dihe_atom4) = split(" ");
		push(@dihedrals_matrix,
			[$dihedral_ID,$dihedral_type,$dihe_atom1,$dihe_atom2,$dihe_atom3,$dihe_atom4]);
	}
	
	# Find out and extract the peptide backbone dihedrals
	#
	# Definitions of peptide backbone dihedrals
	#
	# For dihedral angle PHI: C--N--CA--C
	# For dihedral angle PSI: N--CA--C--N
	#
	# ---------------------------------------------------------
	#  atom  |  mass  |partial charge|     amino-acid
	# ---------------------------------------------------------
	#   C    | 12.011 |     0.51     | all except GLY and PRO
	#   N    | 14.007 |    -0.29     | PRO
	#   N    | 14.007 |    -0.47     | all except PRO
	#   CA   | 12.011 |     0.07     | all except GLY and PRO
	#   CA   | 12.011 |    -0.02     | GLY
	#   CA   | 12.011 |     0.02     | PRO
	# ---------------------------------------------------------
	#
	#  Peptide backbone
	#        ...
	#         /
	#      O=C
	#         \
	#          N-H
	#         / -----> PHI (C-N-CA-C)
	#      H-CA-R
	#         \ -----> PSI (N-CA-C-N)
	#          C=O
	#         /
	#      H-N
	#         \
	#         ...
	#
	# Criteria to be a PHI/PSI dihedral pair:
	# 1. Atoms have to match with the mass/charge constellations as
	#    defined above.
	# 2. The atoms N--CA--C needs to be covalently bonded with each
	#    other.

	# Find which types do C, N and CA correspond to and store them
	# in lists
	my $mass_carbon = 12.011;
	my $mass_nitrogen = 14.007;

	my @carbon_list;
	my @nitrogen_list;
	my $carbon_counter = 0;
	my $nitrogen_counter = 0;

	for (my $i = 0; $i < $natom_types; $i++) {
		if (${$masses_matrix[$i]}[1] == $mass_carbon) {
			push(@carbon_list,${$masses_matrix[$i]}[0]);
			$carbon_counter++;
		}
		if (${$masses_matrix[$i]}[1] == $mass_nitrogen) {
			push(@nitrogen_list,${$masses_matrix[$i]}[0]);
			$nitrogen_counter++;
		}
	}
	# Quit if no carbons or nitrogens
	if ($carbon_counter == 0 or $nitrogen_counter == 0) {
		if ($carbon_counter == 0) {
			print "No carbon atoms exist in the system\n";
		}
		if ($nitrogen_counter == 0) {
			print "No nitrogen atoms exist in the system\n";
		}
		print "CMAP usage impossible\n";
		return;
	}

	print "Carbon atom type/s: @carbon_list\n";
	print "Nitrogen atom type/s: @nitrogen_list\n";

	# Determine the atom types of C, CA and N

	# Charges of the backbone atoms
	my $charge_C = 0.51;
	my $charge_CA = 0.07;
	my $charge_N = -0.47;

	# Special setting for PRO
	my $charge_N_PRO = -0.29;
	my $charge_CA_PRO = 0.02;

	# Special setting for GLY
	my $charge_CA_GLY = -0.02;

	# Peptide backbone atom types
	my $C_type;
	my $CA_type;
	my $CA_GLY_type;
	my $CA_PRO_type;
	my $N_type;
	my $N_PRO_type;

	my $C_counter = 0;
	my $CA_counter = 0;
	my $CA_GLY_counter = 0;
	my $CA_PRO_counter = 0;
	my $N_counter = 0;
	my $N_PRO_counter = 0;

	my $C_flag = 0;

	for (my $i = 0; $i <= $natom_number; $i++) {
		my $cur_type = ${$atoms_matrix[$i]}[2];
		my $cur_charge = ${$atoms_matrix[$i]}[3];
		for (my $j = 0; $j <= $#carbon_list; $j++) {
			if ($cur_type == $carbon_list[$j]) {
				$C_flag = 1;
				if ($cur_charge == $charge_C) {
					$C_type = $cur_type;
					$C_counter++;
				}
				if ($cur_charge == $charge_CA) {
					$CA_type = $cur_type;
					$CA_counter++;
				}
				if ($cur_charge == $charge_CA_GLY) {
					$CA_GLY_type = $cur_type;
					$CA_GLY_counter++;
				}
				if ($cur_charge == $charge_CA_PRO) {
					$CA_PRO_type = $cur_type;
					$CA_PRO_counter++;
				}
			}
		}
		if ($C_flag == 0) {
			for (my $k = 0; $k <= $#nitrogen_list; $k++) {
				if ($cur_type == $nitrogen_list[$k]) {
					if ($cur_charge == $charge_N) {
						$N_type = $cur_type;
						$N_counter++;
					}
					if ($cur_charge == $charge_N_PRO) {
						$N_PRO_type = $cur_type;
						$N_PRO_counter++;
					}
				}
			}
		}
		$C_flag = 0;
	}

	# Quit if one of the atom types doesn't exist
	if ( $C_counter == 0 or
	    ($CA_counter == 0 and $CA_GLY_counter == 0 and $CA_PRO_counter == 0) or
	    ($N_counter == 0 and $N_PRO_counter == 0) ) {
		if ($C_counter == 0) {
			print "\nCannot find the peptide backbone C atom type\n";
		}
		if ($CA_counter == 0 and $CA_GLY_counter == 0 and $CA_PRO_counter == 0) {
			print "\nCannot find the peptide backbone C-alpha atom type\n";
		}
		if ($N_counter == 0 and $N_PRO_counter == 0) {
			print "\nCannot find the peptide backbone N atom type\n";
		}
		print "CMAP usage impossible\n";
		return;
	}

	print "Peptide backbone carbon type: $C_type\n";
	print "Alpha-carbon type: $CA_type\n" if ($CA_counter > 0);
	print "Alpha-carbon type (GLY): $CA_GLY_type\n" if ($CA_GLY_counter > 0);
	print "Alpha-carbon type (PRO): $CA_PRO_type\n" if ($CA_PRO_counter > 0);
	print "Peptide backbone nitrogen type: $N_type\n" if ($N_counter >0);
	print "Peptide backbone nitrogen type (PRO): $N_PRO_type\n" if ($N_PRO_counter > 0);

	# Loop through the dihedral list to find the PHI- and PSI-dihedrals
	my @PHI_dihedrals;
	my @PSI_dihedrals;
	my $PHI_counter = 0;
	my $PSI_counter = 0;

	for (my $i = 0; $i < $ndihedral_number; $i++) {
		my $cur_dihe_ID = ${dihedrals_matrix[$i]}[0];
		my $cur_atom1_type = ${atoms_matrix[${dihedrals_matrix[$i]}[2]-1]}[2];
		my $cur_atom2_type = ${atoms_matrix[${dihedrals_matrix[$i]}[3]-1]}[2];
		my $cur_atom3_type = ${atoms_matrix[${dihedrals_matrix[$i]}[4]-1]}[2];
		my $cur_atom4_type = ${atoms_matrix[${dihedrals_matrix[$i]}[5]-1]}[2];

		next if (${dihedrals_matrix[$i]}[2] == ${dihedrals_matrix[$i-1]}[2] and
		         ${dihedrals_matrix[$i]}[3] == ${dihedrals_matrix[$i-1]}[3] and
	                 ${dihedrals_matrix[$i]}[4] == ${dihedrals_matrix[$i-1]}[4] and
                         ${dihedrals_matrix[$i]}[5] == ${dihedrals_matrix[$i-1]}[5]);

		# Determine PHI-dihedrals; If C-CA-N-C or C-N-CA-C, then save it in a list
		if ($cur_atom1_type == $C_type and $cur_atom4_type == $C_type) {
			if ( ( ($cur_atom2_type == $CA_type or
			        $cur_atom2_type == $CA_GLY_type or
				$cur_atom2_type == $CA_PRO_type) and
			       ($cur_atom3_type == $N_type or
			        $cur_atom3_type == $N_PRO_type) ) or
			     ( ($cur_atom3_type == $CA_type or
				$cur_atom3_type == $CA_GLY_type or
			        $cur_atom3_type == $CA_PRO_type) and
			       ($cur_atom2_type == $N_type or
				$cur_atom2_type == $N_PRO_type) ) ) {
				push (@PHI_dihedrals,$cur_dihe_ID);
				$PHI_counter++;
			}
		}

		# Determin PSI-dihedrals; If N-CA-C-N or N-C-CA-N (N can be both normal N or N proline),
		# then save it in a list
		if ( ($cur_atom1_type == $N_type and $cur_atom4_type == $N_type) or
		     ($cur_atom4_type == $N_PRO_type and $cur_atom1_type == $N_PRO_type) or
	     	     ($cur_atom1_type == $N_type and $cur_atom4_type == $N_PRO_type) or
                     ($cur_atom4_type == $N_type and $cur_atom1_type == $N_PRO_type) ) {
			if ( ( ($cur_atom2_type == $CA_type or
			        $cur_atom2_type == $CA_GLY_type or
		                $cur_atom2_type == $CA_PRO_type) and
			        $cur_atom3_type == $C_type) or
			     ( ($cur_atom3_type == $CA_type or
				$cur_atom3_type == $CA_GLY_type or
			        $cur_atom3_type == $CA_PRO_type) and
				$cur_atom2_type == $C_type) ) {
				push (@PSI_dihedrals,$cur_dihe_ID);
				$PSI_counter++;
			}
		}
	}

	# Quit if no PHI or PSI dihedrals
	if ($PHI_counter == 0 or $PSI_counter ==0) {
		if ($PHI_counter == 0) {
			print "Can not find the PHI backbone dihedrals\n";
		}
		if ($PSI_counter ==0) {
			print "Can not find the PSI backbone dihedrals\n";
		}
		print "CMAP usage impossible\n";
		return;
	}

	# Construct the PHI/PSI diheral pair list
	#
	# The algorithm:
	#         _____
	#        |     |
	#     1--2--3--4      PHI-dihedral
	#     4--3--2--1
	#   --C--N-CA--C--N-- Peptide backbone
	#        1--2--3--4
	#        4--3--2--1   PSI-dihedral
	#        |_____|
	#
	# For a certain PHI dihedral, following conditions have to be met:
	#
	#        PHI         PSI
	#  If (2--3--4) = (1--2--3)
	#  or
	#  if (2--3--4) = (4--3--2)
	#  or
	#  if (3--2--1) = (1--2--3)
	#  or
	#  if (3--2--1) = (4--3--2),
	#
	# then these 2 dihedrals are a PHI/PSI pair. If a pair is found, the
	# dihedral IDs will be stored in "@PHI_PSI_matrix".
	
	my @PHI_PSI_matrix;
	my $crossterm_CA_charge;
	my $crossterm_type;
	my $crossterm_counter = 0;
	my $crossterm_type1_flag = 0;
        my $crossterm_type2_flag = 0;
	my $crossterm_type3_flag = 0;
	my $crossterm_type4_flag = 0;
	my $crossterm_type5_flag = 0;
	my $crossterm_type6_flag = 0;

	for (my $i = 0; $i <= $#PHI_dihedrals; $i++) {
		my $cur_PHI_dihe = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[0];
		my $phi1 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[2];
		my $phi2 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[3];
		my $phi3 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[4];
		my $phi4 = ${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[5];
		for (my $j = 0; $j <= $#PSI_dihedrals; $j++) {
			my $cur_PSI_dihe = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[0];
			my $psi1 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[2];
			my $psi2 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[3];
			my $psi3 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[4];
			my $psi4 = ${dihedrals_matrix[$PSI_dihedrals[$j]-1]}[5];
			if ( ($phi2 == $psi1 and $phi3 == $psi2 and $phi4 == $psi3) or
			     ($phi2 == $psi4 and $phi3 == $psi3 and $phi4 == $psi2) or
			     ($phi3 == $psi1 and $phi2 == $psi2 and $phi1 == $psi3) or
			     ($phi3 == $psi4 and $phi2 == $psi3 and $phi1 == $psi2) ) {

			     # Find out to which amino acid the cross-term belongs

			     if ($phi3 == $psi2 or $phi3 == $psi3) {
				     $crossterm_CA_charge = ${atoms_matrix[${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[4]-1]}[3];
			     }
			     if ($phi2 == $psi2 or $phi2 == $psi3) {
				     $crossterm_CA_charge = ${atoms_matrix[${dihedrals_matrix[$PHI_dihedrals[$i]-1]}[3]-1]}[3];
			     }

			     # Def. the type of the crossterm per cmap.data file; If C_alpha of the crossterm is
			     # - ALA type, then $crossterm_type = 1;
                             # - ALA-PRO (ALA is the current AA), then $crossterm_type = 2;
			     # - PRO type, then $crossterm_type = 3;
                             # - PRO-PRO (First PRO is the current AA), then $crossterm_type = 4;
			     # - GLY type, then $crossterm_type = 5;
                             # - GLY-PRO (GLY is the current AA), then $crossterm_type = 6;

			     if ($crossterm_CA_charge == $charge_CA) { $crossterm_type = 1; $crossterm_type1_flag = 1; }
			     if ($crossterm_CA_charge == $charge_CA_GLY) { $crossterm_type = 5; $crossterm_type5_flag = 1; }
			     if ($crossterm_CA_charge == $charge_CA_PRO) {
					$crossterm_type = 3; $crossterm_type3_flag = 1;
					# Checking the last crossterm, re-assign the last crossterm type if needed
					if ($crossterm_counter-1 >= 0 and $PHI_PSI_matrix[$crossterm_counter-1][0] == 1) {
						$PHI_PSI_matrix[$crossterm_counter-1][0] = 2;
						$crossterm_type2_flag = 1;
					}
					if ($crossterm_counter-1 >= 0 and $PHI_PSI_matrix[$crossterm_counter-1][0] == 3) {
                                                $PHI_PSI_matrix[$crossterm_counter-1][0] = 4;
                                                $crossterm_type4_flag = 1;
                                        }
					if ($crossterm_counter-1 >= 0 and $PHI_PSI_matrix[$crossterm_counter-1][0] == 5) {
                                                $PHI_PSI_matrix[$crossterm_counter-1][0] = 6;
                                                $crossterm_type6_flag = 1;
                                        }
			     }	
			     push(@PHI_PSI_matrix,[$crossterm_type,$phi1,$phi2,$phi3,$phi4,$psi4]);
			     $crossterm_counter++;

			     $crossterm_CA_charge = 0;
			     $crossterm_type = 0;
		     }
		}
	}

	# Check whether the amino acid at the C-terminus is a PRO or not. If yes, the type of the last crossterm
	# should be set to its X-PRO form  instead of X, where X is ALA, PRO, or GLY.  X-PRO form = X type + 1.

	my @pdb_data;
	my $line;
	my $prefix;

	$prefix = basename($data_file);
	$prefix =~ s/^\w+\.//;

	open(PDB,"< $pdb_file")
		or die "WARNING: Cannot open file \"$pdb_file\"! (required if the -cmap option is used)\n";
	@pdb_data = <PDB>;
	close(PDB);

        my @ter_line;
        my $ter_AA_type = 0;
        my $ter_flag = 0;
	foreach $line (@pdb_data) {
		if ($line =~ m/TER/) {
		      @ter_line = split(" ",$line);
                      $ter_AA_type = $ter_line[2]; 	
                      print "Terminal amino acid type is: $ter_AA_type\n";
                      $ter_flag = 1;
	        }
        }
        if ($ter_flag == 0) {
                print "\n*** ERROR IN THE PDB FILE: ***\n";
                print "In order for the CMAP section to be generated, the pdb file must \n";
                print "identify the C-terminus amino acid in the file with 'TER'. \n";
                print "This line is missing from the pdb file that was used.\n";
                print "To correct this problem, open the pdb file in an editor,\n";
                print "find the last atom of the last amino acid residue in the peptide\n";
                print "chain and insert the following line immediately after that atom:\n";
                print "                  'TER   <#1>   <RES>   <#2>' \n";
                print "where '<#1> is the next atom number, <RES> is the three letter amino\n";
                print "acid abbreviation for that amino acid, and <#2> is the molecule number\n";
                print "of the terminal amino acid residue.\n\n";
                print "For example, if the last atom of the last amino acid in the peptide\n";
                print "sequence is listed in the pdb file as:\n\n";
                print "  'ATOM   853  O  GLU  P  56  12.089 -1.695 -6.543 1.00 1.03  PROA'\n\n";
                print "you would insert the following line after it:\n\n";
                print "  'TER    854     GLU     56'\n\n";
                print "If any additional atoms are listed in the pdb file (e.g., water, ions)\n";
                print "after this terminal amino acid residue, their atom numbers and\n";
                print "molecule numbers must be incremented by 1 to account for the new line\n";
                print "that was inserted.\n\n";
                die "Error: No terminating atom designated in pdb file!  See above note to correct problem.\n\n";
        }

        if ($ter_AA_type eq "PRO") {
                $PHI_PSI_matrix[$crossterm_counter-1][0] = $PHI_PSI_matrix[$crossterm_counter-1][0]+1;
        }

	# Print out PHI/PSI diheral pair list
	my $pair_counter = 0;
        # Don't presently use $ncrosstermtypes but have this available if wish to print it out
	my $ncrosstermtypes = $crossterm_type1_flag + $crossterm_type2_flag + $crossterm_type3_flag +
                              $crossterm_type4_flag + $crossterm_type5_flag + $crossterm_type6_flag;
	print "\nWriting to \"$data_file\" with section \"CMAP crossterms\" added at the end.\n";
	
	# Writing the new lammps data file
	open(REWRITE,"> $data_file")
		or die "Cannot write file \"$data_file\"!\n";
	foreach $line (@raw_data) {
		printf(REWRITE "$line\n");
		if ($line =~ m/impropers/) {
			printf(REWRITE "%12d  crossterms\n", $crossterm_counter);
		}
	}
	printf(REWRITE "CMAP\n\n");

        my $ref_line;
        my $column;
	foreach $ref_line (@PHI_PSI_matrix) {
		$pair_counter++;
		printf(REWRITE "%8d",$pair_counter);
		foreach $column (@$ref_line) {
			printf(REWRITE " %7d",$column);
		}
		printf(REWRITE "\n");
	}
	close(REWRITE);
	print "\nDone!\n\n";
}
