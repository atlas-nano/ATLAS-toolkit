#!/usr/bin/perl -w
# BaseComposition.pl - This script will create PX Molecules of the same length
# and same number of crossovers, with the same number of base pairs in the 
# major and minor grooves. These molecules will differ in the base pair
# composition, measured as the %GC. It will start at the starting value
# and continue untill the max %GC is achieved.
# Rules:
# In order too force the crossovers at a certain point, the opposite chains
# cannot complement

# -=== DECLARATION SECTION
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

use strict;
use File::Basename;
use Packages::General;
use Packages::GetParms;
use Packages::HelixLayout;

sub ValidateInput();
sub GetParameters();
sub CreateSequences();
sub SubmitJob();

my $VERSION = "1.00";

# -=== START MAIN PROGRAM
if (! @ARGV or $#ARGV < 2) {
    die "Usage: $0 parameter_file start_percent increment [end_percent]\n";
}
my ($parmfile, $start_percent, $increment, $end_percent) = @ARGV;
my ($P_File, @helices);
ValidateInput();
GetParameters();
CreateSequences();
SubmitJob();
# -=== END MAIN PROGRAM

# -=== SUBROUTINES
sub ValidateInput() {

    -e $parmfile or die "Cannot locate $parmfile: $!\n";
    die "Invalid starting percentage, expected decimal got $start_percent"
	if (! IsDecimal($start_percent) );
    die "Starting Percent Cannot be greater than 99%"
	if ($start_percent > 99);

    die "Invalid increment. Expected decimal got $increment"
	if (! IsDecimal($increment) );
    die "Increment cannot be greater than 99%"
	if ($increment > 99);

    if ($end_percent) {
	die "Invalid end percentage, expected decimal got $start_percent"
	    if (! IsDecimal($end_percent) );
	die "Ending Percent Cannot be less than 1%"
	    if ($end_percent < 1);
    }
}

sub GetParameters() {
    $P_File = Packages::GetParms->new();
    if (! $P_File->IsValidParams($parmfile)) {
	die "Error in Paramater file\n";
    }
}

sub CreateSequences() {
   my ($max_GC, $counter, $free_base_pairs, $constrained_base_pairs, $i);
   
   my ($base_pairs, @crossovers) = 
       ($P_File->{"Molecule"}->{"total_bases"}, 
	$P_File->{"Molecule"}->{"crossovers"} );
   
   $base_pairs = $base_pairs/4;
   
   $free_base_pairs = ($base_pairs - ($#crossovers * 2) );
   $max_GC = $free_base_pairs/$base_pairs;
   
   $end_percent = $max_GC if (! $end_percent);
   
   ($end_percent, $start_percent) = ($start_percent, $end_percent)
       if ($start_percent > $end_percent);
   
   $counter = $start_percent;

   while ($counter < $end_percent) {
       $counter += $increment;
       for $i (1 .. 2) {
	   @helix = GetRndBasePairs($free_base_pairs, $counter);
	   $sequence = CreateSequence(\@helix, $i);
	   WriteSequence($sequence);
       }
       CreateStructure($counter);
   }
}

sub GetRndBasePairs(@) {
    my ($total_bp, $percent_GC) = @_;
    my ($amt_GC) = int ($total_bp * $percent_GC/100);
    my (@bps, @result, $curr_val);

    for (0 .. ($amt_GC - 1)) {
	$bps[$_] = "GC";
    }

    for ($amt_GC .. ($total_bp - $amt_GC)) {
	$bps[$_] = "AT";
    }
	
    while ($#bps) {
	$curr_val = int (rand $total_bp); 
	push @result, $bps[$curr_val];
	delete $bps[$curr_val];
    }

    push @result, $bps[0];
    return @result;
}

sub CreateSequence(@) {
    my ($helix_info, $which_helix) = @_;
    my ($result, $i, $counter, $tot_base_counter);

#create bases at the end
    $counter = 0;
    for $i (1 .. $P_File->{"Molecule"}->{"bases_at_end"}) {
	$result .= $helix_info->[$counter];
	$counter++;
    }

    $tot_base_counter = ($P_File->{"Molecule"}->{"total_bases"}/4) - 
	$P_File->{"Molecule"}->{"bases_at_end"};

    while ($tot_base_counter > 0) {
}
sub SubmitJob() {

}
