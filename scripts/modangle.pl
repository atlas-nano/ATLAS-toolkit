#!/usr/bin/perl -w

for ($i=4; $i<10; $i++) {
    for ($j=4; $j<10; $j++) {
	open SEQUENCFLE, "sequence.txt" or die "ERROR: Cannot open seqyence.txt: $!\n";
	while (<SEQUENCFLE>) {
	    chomp($_);
	    push @outarray, $_;
	}
	close SEQUENCFLE;

	open NEWSEQUENCE, "> sequence.txt";
	print NEWSEQUENCE "$i:$j\n";
	for ($k = 1; $k <=$#outarray; $k++) {
	    print NEWSEQUENCE $outarray[$k] . "\n";
	}
	close NEWSEQUENCE;

	@outarray = ();
	$pxmol = $i . $j;
	$curr_fle = $pxmol . ".pdb";
	$my_cmd =  "~/scripts/createhelix.pl sequence.txt $curr_fle 0";
	system $my_cmd;
	$my_cmd = "~/scripts/prep_sander.pl $curr_fle $pxmol /ul/tpascal/scripts";
	system $my_cmd;
    }
}
 
