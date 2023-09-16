#!/usr/bin/perl -w
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use Packages::HelixLayout;
use Packages::General;

sub GetFileInfo();
sub StoreInfo($);
sub ObtainBaseInfo(@);
sub DetermineUnit(@);
sub CreateHelices();
sub GetUnit(@);
sub GetOutName($);
sub ValidateInput();
sub ReverseRegions(@);

if (!@ARGV or $#ARGV <4) {
    die "usage: helix_from_strands.pl pdbfile majorgroove:minorgroove ", 
    "is3primeIn numBasesAtEnd crossovers\n";
}

my ($in_file, $strand_unit_holder, $is3primeIn, $numBasesAtEnd, @crossovers) = @ARGV;
my ($majorgroove, $minorgroove) = split /:/, $strand_unit_holder;

my ($i, @pdb_info, $strand_data, $total_base_pairs);
my (@helix , $unit_tracker, $unit_counter, $base_counter);

my $hl = Packages::HelixLayout->spawn();

#  --== START ==--

ValidateInput();
GetFileInfo();

print "Processing $in_file....";
open INFILE, $in_file or die "Cannot open $in_file: $!\n";
while (<INFILE>) {
    chomp;
    StoreInfo($_);
}
print "Done\n";
close INFILE;

print "Sending $#crossovers crossovers\n";

$hl->DetermineHelixLayout($majorgroove, $minorgroove, $is3primeIn, 
			  $numBasesAtEnd, $total_base_pairs, @crossovers);
@helix = $hl->GetHelixInfo();

CreateHelices;
 
sub CreateHelices() {
    my ($helices, $helix_strand, $helix_regions);
    my ($curr_strand, $start_unit, $end_unit);
    my ($lineout, $outfile, $no_cross_file, $file_unit_counter, $file_base_counter);

    $file_unit_counter = 1;
    $file_base_counter = 1;

    system "mkdir -p pdbfiles/nocrossovers";

    $no_cross_file = "pdbfiles/nocrossovers/$in_file";
    open CROSSFILE, "> $no_cross_file" or die "Cannot create $no_cross_file: $!\n";

    for $helices (0 .. $#helix) {
	print "\nHelix $helices\n";
	$outfile = GetOutName($helices);

	open OUTFILE, "> $outfile" or die "Cannot create $outfile: $!\n";
	print "Creating $outfile...";
	$unit_counter = 1;
	$base_counter = 1;

	for $helix_strand (@{$helix[$helices]}) {
	    ReverseRegions($helices, 1) if ($unit_counter > 1);
	    for $helix_regions (@{$helix_strand}) {
		$curr_strand = $helix_regions->{"Strand"} - 1;
		$start_unit = $helix_regions->{"StartUnit"};
		$end_unit = $helix_regions->{"EndUnit"};

		if ($start_unit > $end_unit) {
		    ($start_unit, $end_unit) = ($end_unit, $start_unit);
		}
		for $i ($start_unit .. $end_unit) {
		    ($unit_counter, $lineout) = ObtainBaseInfo($curr_strand, $i, $unit_counter,
							       $base_counter);
		    print OUTFILE "$lineout";

		    ($file_unit_counter, $lineout) = ObtainBaseInfo($curr_strand, $i, 
								    $file_unit_counter, 
								    $file_base_counter);
		    print CROSSFILE "$lineout";
		    $base_counter++;
		    $file_base_counter++;
		}
	    }
	    print OUTFILE "TER\n";
	    print CROSSFILE "TER\n";
	}
	close OUTFILE;
	print "Done\n";
    }
    close CROSSFILE;
# now type the no crossover file for xleap
    system "/ul/tpascal/scripts/fixxleappdb.pl $no_cross_file";
}

sub ObtainBaseInfo(@) {
    my ($strand_id, $unit_id, $u_counter, $b_counter) = @_;
    my ($i, $result, $atom_info, $update_unit);

    for $i (@{$pdb_info[$strand_id]->[$unit_id]}) {
	$atom_info = sprintf("ATOM%7d", $u_counter);
	$result .= $atom_info . "$i\n";
	$update_unit = sprintf("%7d", $b_counter);
	$result =~ s/_.+_/$update_unit/g;
	$u_counter++;
    }

    return ($u_counter, $result);
}



sub GetFileInfo() {
    
    $total_base_pairs = 0;

    open TMPFILE, "tail -1 $in_file |" or die "Cannot open $in_file\n";
    while (<TMPFILE>) {
	chomp;
	if ($_ =~ /^ATOM\s+\d+\s*\S*\s*\w+\s+(\d+)/) {
	    $total_base_pairs = $1/4;
	}
    }
    close TMPFILE;
    
    $total_base_pairs ?
	print "Found $total_base_pairs units per strand\n" :
	die "File: $in_file is invalid\n";
}

sub StoreInfo($) {
    my ($line_in) = $_[0];
    my ($current_strand, $current_unit, $unitline) = (0, 0, "");

    if ($line_in =~ /^ATOM\s+\d+(\s*\S*\s*)(\w+)\s+(\d+)\s+(.+)$/) {
	$current_strand = int(($3 - 1) / $total_base_pairs );
	$current_unit = GetUnit($current_strand, $3);

	$unitline = sprintf("%6s%-2s_%7d_%40s", $1, $2, $current_unit, $4);
	push @{$pdb_info[$current_strand][$current_unit]}, $unitline;
     }
}

sub GetUnit(@) {

    my ($curr_strand, $curr_unit) = @_;
    my ($returnval, $isNewStrand) = ("", 1);
    my ($total_units_in_strand, $total_units);

    my ($old_strand, $old_unit) = split /:/, $strand_unit_holder;

    if ($curr_strand == $old_strand) {
	if ($curr_unit > $old_unit or $curr_unit < $old_unit) {
	    $unit_tracker += 1;
	}
	$returnval = $unit_tracker; 
    } else {
	$returnval = 1; 
	$unit_tracker = 1;
    }

    $strand_unit_holder = $curr_strand . ":" . $curr_unit;

    return $returnval;

}

sub GetOutName($) {
    my ($helix_no) = $_[0];
    my ($tempstr);

    $tempstr = $in_file;
    print "INFILE: $in_file\n";
    if ($tempstr =~ /^(.+)_(\d+)\.pdb$/) {
	$tempstr = $1 . "_Helix" .($helix_no + 1) . "_" . $2 .".pdb";
    }else {
	$tempstr .= "_Helix" . $helix_no . ".pdb";
    }

    return $tempstr;
}

sub ValidateInput() {
    my (@tmp_info);

    if (! IsInteger($majorgroove) || ! IsInteger($minorgroove)) {
	die "Invalid major_groove:minor_groove specification. Expected integer\n";
    }

    die "Invalid is3primeIn value: Expected 1 or 0\n" 
	if (! IsInteger($is3primeIn) );
    die "Invalid number of bases at the ends. Expected Integer\n" 
	if (! IsInteger($numBasesAtEnd));

    -e $in_file or die "Cannot locate $in_file: $!\n";

    for $i (@crossovers) {
	push @tmp_info, $i 
	    if ( IsInteger($i) || $i <= $total_base_pairs);
    }

    $#tmp_info == -1 ?
	die "Invalid crossover specification. ", 
	"Please supply base numbers delimited by spaces\n" :
	print "" . ($#crossovers + 1) . " total crossovers\n";
    @crossovers = ();
    @crossovers = @tmp_info;
}

sub ReverseRegions(@) {
    my ($whichHelix, $whichStrand) = @_;
    my ($totalRegions) = $#{ $helix[$whichHelix]->[$whichStrand] };

    @{ $helix[$whichHelix]->[$whichStrand] } = reverse @{ $helix[$whichHelix]->[$whichStrand] };

}
