#!/usr/bin/perl -w

# This program will determine the load on the node and wolf machines
# and submit a job.

# usage submit_job.pl topfile coord_file no_of_bases scripts_dir [machinename]

my ($topfile, $coord_file, $no_bases, $scripts_dir, $node_name) = @ARGV;

if (! $topfile or ! $coord_file or ! $scripts_dir or !@ARGV) {
    die "usage: submit_job.pl topfile coord_file no_of_bases scripts_dir [machinename]\n";
}

-e $topfile or die "Cannot find $topfile: $!\n";
-e $coord_file or die "Cannot find $coord_file: $!\n";
-e $scripts_dir . "/run_md.pl" or die "Cannot find $scripts_dir: $!\n";

if (! $no_bases =~ /^\d+$/) {
    die "Invalid no_of_bases. Got $no_bases, expected integer\n";
} 

#first figure out which is the best computer to use

if (! $node_name) { 
    $my_cmd = "mscload -m\"node wolf\"";
    open FINDPC, "$my_cmd |" or die "Cannot execute mscload command\n";
    
    @loadhash = (
		 {
		     "Name" =>"",
		     "Load" =>0
		     }
		 
		 );
    
    my $my_counter = 0;
    while (<FINDPC>) {
	$myin = $_;
	if ($myin =~ /\s+(\w+)\s+(\d+)/) {
	    if ($2 <3) {
		$loadhash[$my_counter]{"Name"} = $1;
		$loadhash[$my_counter]{"Load"} = $2;
		$my_counter += 1;
	    }
	}
    }
    close FINDPC;

    $low_load = 99;
    $node_name = "";

    if ($my_counter > 0) {
	for ($i=0; $i<$my_counter; $i++) {
	    if ($loadhash[$i]{"Load"} < $low_load and $loadhash[$i]{"Load"} < 3) {
		$node_name = $loadhash[$i]{"Name"};
		$low_load = $loadhash[$i]{"Load"};
	    }
	}
    } else {
	die "Cannot find any machine with a load less than 3\n";
    }
}

if ($node_name ne "") {
    print "\nRunning simulation on $node_name\n";
    #$my_cmd = "| rlogin $node_name | $scripts_dir/run_md.pl $topfile $coord_file $no_bases";
    #open TEST, $my_cmd or die "Cannot fork: $!\n";
    #while (<TEST>) {
    #}
    #close TEST;

    $scr_nm = "> " . $ENV{"HOME"} . "/". $topfile . ".md.script";
    open WRITEFLE, $scr_nm or die "Cannot create $scr_nm: $!\n";
    #print WRITEFLE "#!/bin/csh -f\n";
    print WRITEFLE "cd " . $ENV{"PWD"} . "\n";
    print WRITEFLE $scripts_dir . "/run_md.pl $topfile $coord_file $no_bases";
    close WRITEFLE;

    system "chmod +x " . $ENV{"HOME"} . "/". $topfile . ".md.script";
    $my_cmd = "rsh $node_name " . $topfile . ".md.script";
    #print $my_cmd . "\n";
    system $my_cmd; 
} else {
    die "Cannot find any node with a load less than 2\n";
}

$new_dir = $topfile . "-results";
system "rm -fr $new_dir";
mkdir "$new_dir", 0777;
chdir "$new_dir";

$my_cmd = "rcp " . $node_name . ":/temp1/tpascal/" . $topfile . "-dir/products/finalstruct* .";
if (system $my_cmd) {
    print "\n\nSimulations Completed\n";
    print "The ouputfile can be found in the $new_dir directory\n";
    chdir "../";
} else {
    die "There was an error\n";
}
