#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use File::Basename;
use Packages::General;
use Packages::FileFormats;
use Packages::GetParms;
use Packages::HelixLayout;

sub DetermineHelixLayout();
sub GetCorrectUnit(@);
sub MapPDBFile();
sub CreateHBondList();
sub CreateCarnalInput(@);
sub GetBase(@);
sub Numerically;
sub CarnalInput(@);
sub ObtainHBondData();
sub ProcessData(@);
sub WriteOutput();
sub MapSingleHelix();

die "usage: $0 parmfile\n"
    if (!@ARGV);

my ($paramfile) = $ARGV[0];
my (@helix, $P_File, $PDBFile, %HBond, %BOND_DATA, $helix_length);

{
    
    print "\n-==START HBOND ANALYSIS ==-\n\n";
    open DUMPFILE, "> dumpfile.txt" or die "Cannot create dumpfile.txt: $!\n";

    $P_File = Packages::GetParms->new();
    if (! $P_File->IsValidParams($paramfile)) {
	die "Error in Paramater file\n";
    }

    print "Step 1: Determining Structural Layout...";
    if ($#{ $P_File->{"Molecule"}->{"crossovers"} } > 0 ||
        ! $P_File->{"Molecule"}->{"crossovers"}->[0] == -1) {
	DetermineHelixLayout();
    } else {
	MapSingleHelix();
    }
    print "Done\nStep 2: Obtaining Structural Map...";
    MapPDBFile();
    CreateHBondList();
    print "Done\nStep 3: Calculating H-Bonds...";
    ObtainHBondData();
    print "Done\nStep 4: Creating Tables...";
    WriteOutput();
    print "Done\n\tCreated file hbond.dat with hydrogen bond data.\n";
    print "\n-==END HBOND ANALYSIS ==-\n";
    close DUMPFILE;
}

sub DetermineHelixLayout() {
    my $hl = Packages::HelixLayout->spawn();
    $helix_length = $P_File->{"Molecule"}->{"total_bases"}/4;

    $hl->DetermineHelixLayout(
			      $P_File->{"Molecule"}->{"major_groove"}, 
			      $P_File->{"Molecule"}->{"minor_groove"}, 
			      $P_File->{"Molecule"}->{"is3PrimeIn"}, 
			      $P_File->{"Molecule"}->{"bases_at_end"}, 
			      $helix_length, 
			      @{ $P_File->{"Molecule"}->{"crossovers"} }
			      );

    @helix = $hl->GetHelixInfo();
}

sub MapSingleHelix() {
    my ($region_counter, $h1_counter, $h2_counter, $counter);
    my ($rec_h1, $rec_h2);
    $helix_length = $P_File->{"Molecule"}->{"total_bases"}/2;

    $h1_counter = 1;
    $h2_counter = $helix_length;

    $counter = 0;
    while ($h1_counter < $helix_length) {
	$rec_h1 = (
		    {
			"EndUnit" => $h1_counter + $P_File->{"Molecule"}->{"minor_groove"} - 1,
			"StartUnit" => $h1_counter,
			"Strand" => 1,
		    }
		  );
	$h1_counter += $P_File->{"Molecule"}->{"minor_groove"};

        $rec_h2 = (
                    {
                        "StartUnit" => $h2_counter - $P_File->{"Molecule"}->{"minor_groove"} + 1,
                        "EndUnit" => $h2_counter,
                        "Strand" => 2,
                    }
                  );
        $h2_counter -=  $P_File->{"Molecule"}->{"minor_groove"};

	push @{ $helix[0][0] }, $rec_h1;
	push @{ $helix[0][1] }, $rec_h2;

	if ($h1_counter >= $helix_length) {
	    last;
	}

        $rec_h1 = (
                    { 
                        "EndUnit" => $h1_counter + $P_File->{"Molecule"}->{"major_groove"} - 1,
                        "StartUnit" => $h1_counter, 
                        "Strand" => 1,
                    } 
                  );
        $h1_counter += $P_File->{"Molecule"}->{"major_groove"};
 
        $rec_h2 = (
                    { 
                        "StartUnit" => $h2_counter - $P_File->{"Molecule"}->{"major_groove"} + 1,
                        "EndUnit" => $h2_counter,
                        "Strand" => 2,
                    }
                  );
	$h2_counter -=  $P_File->{"Molecule"}->{"major_groove"};

        push @{ $helix[0][0] }, $rec_h1;
        push @{ $helix[0][1] }, $rec_h2;
    }

#    @{ $helix[0][1] } = reverse @{ $helix[0][1] };
}

sub MapPDBFile() {
    my ($pdbfile, $out_cmd);

# create a pdbfile and get the base name

    open OUTFILE, "> ptraj_in" || die "Cannot create tmp file ptraj_in: $!\n";
    print OUTFILE "trajin " . $P_File->{"Files"}->{"trajectory"} . " 1 1 1\n";
    print OUTFILE "trajout test pdb\n";
    close OUTFILE;

    $out_cmd = "/ul/maiti/ptraj-6.3/linux/ptraj ";
    $out_cmd .= $P_File->{"Files"}->{"topology"};
    $out_cmd .= " < ptraj_in >& junk"; 

    if (system($out_cmd) && ! -e "test.1") {
	die "Error executing ptraj\n";
    }

    $PDBFile = GetPDBFileInfo("test.1", 0);

    system "rm -f test.1";
}

sub CreateHBondList() {
    my ($curr_helix, $curr_region, $tot_regions, $str_len);
    my ($s1_base, $s2_base, $s1_strand, $s2_strand);
    my ($curr_s1_base, $curr_s2_base, $base_counter, $counter);
    my ($s1_base_res, $s2_base_res, $hbond_pair, $carnal_data, $num_bonds);
    my ($pair_counter);

    $tot_regions = $#{ $helix[0][0] };

    $str_len = $helix_length;
    $pair_counter = 0;
    for $curr_helix (@helix) {
	for $curr_region (0 ..$tot_regions) {
	    $s1_base = $curr_helix->[0]->[$curr_region]->{"StartUnit"};
	    $s2_base = $curr_helix->[1]->[$curr_region]->{"EndUnit"};
	    $s1_strand = $curr_helix->[0]->[$curr_region]->{"Strand"} - 1;
	    $s2_strand = $curr_helix->[1]->[$curr_region]->{"Strand"} - 1;
	    
	    $base_counter = $curr_helix->[0]->[$curr_region]->{"EndUnit"} - $s1_base;
	    for $counter (0 .. $base_counter) {
		$curr_s1_base = ($s1_base + $counter) + ($s1_strand * $str_len);
		$curr_s2_base = ($s2_base - $counter) + ($s2_strand * $str_len);
#		$curr_s1_base = $s1_base + $counter;
#		$curr_s2_base = $s2_base - $counter;
		$hbond_pair = $curr_s1_base . "_" . $curr_s2_base;
		$s1_base_res = GetBase($curr_s1_base);
		$s2_base_res = GetBase($curr_s2_base);
		($num_bonds, $carnal_data) = 
		    CreateCarnalInput($s1_base_res, $s2_base_res, 
				      $curr_s1_base, $curr_s2_base);
		if ($carnal_data) {
		    $HBond{$pair_counter}->{"PAIR"} = $hbond_pair;
		    $HBond{$pair_counter}->{"CARNAL_IN"} = $carnal_data;
		    $HBond{$pair_counter}->{"NUM_BONDS"} = $num_bonds;
		    $HBond{$pair_counter}->{"HAS_DATA"} = 0;
		}
		$pair_counter++;
	    }
	}
    }
}

sub CreateCarnalInput(@) {
    my ($s1, $s2, $s1_num, $s2_num) = @_;
    my ($is_valid, $num_bonds, $result);


    $is_valid = 0;
    if ( ($s1 =~ "DG" && $s2 =~ "DC")) {
	$result = CarnalInput($s2_num, $s1_num, 1);
	$num_bonds = 3;
	$is_valid = 1;
    }elsif ($s1 =~ "DC" && $s2 =~ "DG") {
	$result = CarnalInput($s1_num, $s2_num, 1);
	$num_bonds = 3;
	$is_valid = 1;
    }elsif ($s1 =~ "DA" && $s2 =~ "DT") {
	$result = CarnalInput($s2_num, $s1_num, 0);
	$num_bonds = 2;
	$is_valid = 1;
    }elsif ($s1 =~ "DT" && $s2 =~ "DA") {
	$result = CarnalInput($s1_num, $s2_num, 0);
	$num_bonds = 2;
	$is_valid = 1;
    }

    if (! $is_valid) {
	print DUMPFILE "WARNING: Invalid base pair: $s1_num $s1 - $s2_num $s2\n";
	return;
    } else {
	return ($num_bonds, $result);
    }
}

sub GetBase(@) {
    my ($counter, $result);
    my ($in_base) = $_[0];

    for $counter (sort Numerically keys %{ $PDBFile }) {
	if ($in_base == $PDBFile->{$counter}->{"RES_ID"}) {
	    $result = $PDBFile->{$counter}->{"RES_NAME"};
	    last;
	}
    }

    if (! $result) {
	print DUMPFILE "WARNING: Invalid base id $in_base\n";
	return;
    }

    return $result;
}

sub Numerically {
    ($a <=> $b);
}

sub CarnalInput(@) {
    my ($s1_num, $s2_num, $is_GC) = @_;
    my ($result);

    $result = "FILES_IN\n";
    $result .= "\tPARM p1 " . $P_File->{"Files"}->{"topology"} . ";\n";
    $result .= "\tSTREAM s1 " . $P_File->{"Files"}->{"trajectory"} . ";\n";
    $result .= "FILES_OUT\n";
    $result .= "\tHBOND h1 hb1 TABLE;\n\tHBOND h2 hb2 TABLE;\n";
    if ($is_GC) {
	$result .= "\tHBOND h3 hb3 TABLE;\n";
    }
    $result .= "DECLARE\n";
    $result .= "\tGROUP g1 ( (RES $s1_num) & (ATOM TYPE O) );\n";
    $result .= "\tGROUP g2 ( (RES $s2_num) & (ATOM TYPE N2) );\n";
    if ($is_GC) {
	$result .= "\tGROUP g3 ( (RES $s1_num) & (ATOM TYPE NC) );\n";
	$result .= "\tGROUP g4 ( (RES $s2_num) & (ATOM TYPE NA) );\n";
	$result .= "\tGROUP g5 ( (RES $s2_num) & (ATOM TYPE O) );\n";
	$result .= "\tGROUP g6 ( (RES $s1_num) & (ATOM TYPE N2) );\n";
    } else {
	$result .= "\tGROUP g3 ( (RES $s2_num) & (ATOM TYPE NC) );\n";
	$result .= "\tGROUP g4 ( (RES $s1_num) & (ATOM TYPE NA) );\n";
    }

    $result .= "OUTPUT\n";
    $result .= "\tHBOND h1 DONOR g2 ACCEPTOR g1;\n";
    $result .= "\tHBOND h2 DONOR g4 ACCEPTOR g3;\n";
    if ($is_GC) {
	$result .= "\tHBOND h3 DONOR g6 ACCEPTOR g5;\n";
    }
    $result .= "END\n";

    return $result;
}

sub ObtainHBondData() {

    my ($curr_pair, $out_cmd);

    system "mkdir -p data";
    for $curr_pair (sort Numerically keys %HBond) {
	open OUTFILE, "> carnal_in" or die "Cannot create temporary files: $!\n";
	print OUTFILE $HBond{$curr_pair}->{"CARNAL_IN"};
	close OUTFILE;
	$out_cmd = "carnal -O < carnal_in >& junk";
	if (system($out_cmd)) {
	    print DUMPFILE "ERROR while executing CARNAL on ";
	    print DUMPFILE $HBond{$curr_pair}->{"PAIR"} . "\n";
	}else {
	    $HBond{$curr_pair}->{"HAS_DATA"} = 1;
	    ProcessData($HBond{$curr_pair}->{"NUM_BONDS"}, 
			$HBond{$curr_pair}->{"PAIR"}, $curr_pair);
	}
	system "rm -f hb*.tab";

    }
}

sub ProcessData(@) {
    my ($tot_files, $pair_nm, $pair) = @_;
    my ($file_counter, $curr_file, $inText, $bond_counter, $time_step);

    for $file_counter (1 .. $tot_files) {
	$curr_file = "hb" . $file_counter . ".tab";
	if (open INFILE, $curr_file) {
	    $time_step = 1;
	    while (<INFILE>) {
		chomp;
		$inText = $_;
		if ($inText =~ /^(\d+)/) {
		    $bond_counter = $1;
		    while ($inText =~ /\s+(\d+)/g) {
			$bond_counter += $1;
		    }
		    $BOND_DATA{$time_step}{$pair}{"BOND_LIST"}[$file_counter-1] = $bond_counter;
		    $BOND_DATA{$time_step}{$pair}{"PAIR"} = $pair_nm;	
		    $BOND_DATA{$time_step}{"AVAILABLE_BONDS"} += $bond_counter;
		    $BOND_DATA{$time_step}{"TOTAL_BONDS"} += 1;
		    $time_step++;
		}
	    }
	    close INFILE;
	} else {
	    print DUMPFILE "ERROR: No H-Bond Data found for file $curr_file in $pair_nm\n";
	}
	system "cat $curr_file >> data/RES_" . $pair_nm . ".dat";

    }
}

sub WriteOutput() {
    my ($time_step, $curr_percent, $outData);
    $outData = "";
    for $time_step (sort Numerically keys %BOND_DATA) {
	$curr_percent = 100 * $BOND_DATA{$time_step}{"AVAILABLE_BONDS"} / 
	    $BOND_DATA{$time_step}{"TOTAL_BONDS"};
	$outData .= sprintf("%6d %7.3f\n", $time_step, $curr_percent);
    }

    open OUTFILE, "> hbond.dat" || die "Cannot create file hbond.dat: $!\n";
    print OUTFILE $outData;
    close OUTFILE;
}
