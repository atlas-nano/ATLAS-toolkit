#!/usr/bin/perl -w
use strict;
BEGIN {
    push (@INC, "/ul/tpascal/scripts/");
}

# This script will open a PDB file, determine the sequence
# and determine the randomness index:
# index = sum(#occurances of pair - expectation of pair)/(total pairs)
use strict;
use Packages::General;
use Packages::FileFormats;

sub Initialize();
sub ParsePDBFILE(@);
sub DetermineRandomness(@);
sub Numerically;
sub GetValidFiles(@);
sub WriteOutput(@);

die "usage: $0 pdbfile|directory (strand_start) (strand_end)\n"
    if (! @ARGV);

my ($loc, $base_start, $base_end) = @ARGV;
my ($pdbfile, $SEQUENCE, $RI, $index, $MYFILES, @STATS, $rec, $mySequence, $longestFileName);

$base_start = 1
    if (! defined($base_start));
$base_end = 1000
    if (! defined($base_end));

{
    Initialize();
    ($MYFILES, $longestFileName) = GetValidFiles($loc);
    for $pdbfile (@{ $MYFILES }) {
    
	$SEQUENCE = ParsePDBFILE($pdbfile, $base_start, $base_end);
	$RI = DetermineRandomness($SEQUENCE);

	$mySequence = "";
	for $index (@{ $SEQUENCE }) {
	    $mySequence .= $index;
	}
	
	$rec = (
		{
		    "FILE"       => $pdbfile,
		    "SEQUENCE"   => $mySequence,
		    "RANDOMNESS" => $RI,
		}
		);
	push @STATS, $rec;
	
    }

    WriteOutput(\@STATS, $longestFileName);
}

sub Initialize() {
    if (! -d $loc) {
	die "ERROR accessing PDB file $loc: $!\n"
	    if (! -e $loc or ! -r $loc or ! -T $loc);
    } else {
	die "ERROR accessing directory $loc: $!\n"
	    if (! -e $loc);
    }

    die "ERROR: Expected integer for startbase. Got $base_start\n"
	if (! IsInteger($base_start));

    die "ERROR: Expected integer for endbase. Got $base_end\n"
	if (defined($base_end) and ! IsInteger($base_end));

}

sub GetValidFiles(@) {
    my ($loc) = $_[0];
    my (@holder, $name_len, $counter);
                                                                                                                     
    if (-e $loc) {
        if (! -d $loc) {
            $holder[0] = $loc;
        } else {
            opendir (WORKDIR, $loc) or die "Cannot access specified directory: $!\n";
            @holder = grep { /\w+\.pdb$/ } map {"$loc/$_"} readdir WORKDIR;
            closedir WORKDIR;
        }
    }
                                                                                                                     
    die "ERROR: Cannot find any valid pdbfile here: $loc\n"
        if (! @holder);
        
    $name_len = 0;
    for $counter (@holder) {
	if (length($counter) > $name_len) {
	    $name_len = length($counter);
	}
    }

    return (\@holder, $name_len);
}

sub ParsePDBFILE(@) {
    my ($pdb_file, $start, $end) = @_;
    my ($counter, $Atom_Info, $curr_res_id, $curr_res_name, $old_res_id, @sequence, @holder);

    

    $Atom_Info = GetPDBFileInfo($pdb_file, 0, " ", " ");
    if (! defined($ARGV[2])) {
	@holder = sort Numerically keys %{ $Atom_Info };
	$end = $old_res_id = $curr_res_id = 0;
	for $counter (@holder) {
	    $curr_res_id = $Atom_Info->{$counter}{"RES_ID"};
	    if ($old_res_id < $curr_res_id) {
		$old_res_id = $curr_res_id;
		$end++;
	    }
		
	}
	$end = ($end/2);
    }
    
    $old_res_id = $curr_res_id = 0;
    for $counter (sort Numerically keys %{ $Atom_Info }) {
	$curr_res_id = $Atom_Info->{$counter}{"RES_ID"};
	$curr_res_name = $Atom_Info->{$counter}{"RES_NAME"};
	next
	    if ($old_res_id == $curr_res_id);
	if ($curr_res_id >= $start and $curr_res_id <= $end) {
	    if ($curr_res_name =~ /D?(\w)/) {
		push @sequence, $1;
		$old_res_id = $curr_res_id;
	    } else {
		print "INVALID RESIDUE NAME: $curr_res_name\n";
	    }
	} else {
	    last;
	}
    }

    die "ERROR: Error obtaining sequence from $pdb_file\n"
	if (! @sequence);

    return \@sequence;
}

sub DetermineRandomness(@) {
    my ($helix_sequence) = $_[0];

    my ($sequence_expectation, $strand_length, $counter, $curr_base, $next_base);
    my ($hash_key, %SEQ_TRACKER, $random_index, $tmp);

    #There are 16 unique pairs, so round up to the nearest integer
    $sequence_expectation = int(($#{ $helix_sequence } + 1)/16); 
    if ((($#{ $helix_sequence } + 1) % 16) > 0) {
	$sequence_expectation++;
    }
    
    #Now march along the helix, record the occurance of each sequence pair
    $strand_length = $#{ $helix_sequence}; # Only march along one strand

    for $counter (0 .. ($strand_length - 1)) {
	$curr_base = $helix_sequence->[$counter];
	if (defined($helix_sequence->[$counter + 1])) {
	    $next_base = $helix_sequence->[$counter + 1];
	    $hash_key = $curr_base . $next_base;
	    if (! defined($SEQ_TRACKER{$hash_key})) {
		$SEQ_TRACKER{$hash_key} = 1;
	    } else {
		$SEQ_TRACKER{$hash_key} += 1;
	    }
	}
    }

    $random_index = 0;
    for $counter (keys %SEQ_TRACKER) {
	if ($SEQ_TRACKER{$counter} > $sequence_expectation) {
	    $random_index += ($SEQ_TRACKER{$counter} - $sequence_expectation);
	}
    }

    $random_index = 1 - ($random_index/$strand_length);
    
    return $random_index;
}

sub Numerically {
    ($a <=> $b);
}

sub WriteOutput(@) {
    my ($STATS, $col_2_len) = @_;
    my ($counter, $fmt, $fmt1);

    $col_2_len += 3;
    $fmt = "%-" . $col_2_len . "s%13s %20s\n";
    $fmt1 = "%-" . $col_2_len . "s%13.2f";
    
    printf $fmt, "File", "Randomness", "SEQUENCE";
    for $counter (1 .. ($col_2_len + 34)) {
	print "-";
    }

    print "\n";

    for $counter (@{ $STATS}) {
	printf $fmt1, $counter->{"FILE"}, $counter->{"RANDOMNESS"};
	print " " .  $counter->{"SEQUENCE"} . "\n";
    }


}
