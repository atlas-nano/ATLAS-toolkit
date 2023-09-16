#!/usr/bin/perl -w
# Finds the number of jobs that the specified user has on the specified
# machines
#
# usage: myjobs.pl [username] [mscload string]

if ($#ARGV <0 ){
    $selected_nm = $ENV{"LOGNAME"};
} else {
    $selected_nm = $ARGV[0];
}

if ($ARGV[1]) {
    $mscload_string = "mscload " . $ARGV[1];
} else {
    $mscload_string = "mscload";
}

print "Current jobs for $selected_nm ";
if ($ARGV[1]) {
    print "using: $ARGV[1]";
}
print "\n";

print " Host     Load  Rating  Mb Free     %CPU    Mb    User      Command\n";
print "-------------------------------    ---------------------------------------\n";

$no_jobs = 0;

open LOADLIST, "$mscload_string |" or die "Cannot execute mscload command\n";
while(<LOADLIST>) {
    $my_in = $_;
    chomp($my_in);
    if ($my_in =~ /(\s+\w+\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+)\s\d+\.\d+\s+\d+\.\d+\s+(\w+)/) {
	$machine_nm = $1;
	$user_nm = $2;
	if ($user_nm eq $selected_nm) {
	    print "$my_in\n";
	    $no_jobs += 1;
	}
    } else {
	if ($my_in =~ /^\s+(\d+\.\d+\s+\d+\.\d+\s+)(\w+)(.+)/) {
	    $user_nm = $2;
	    if ($user_nm eq $selected_nm) {
		print "$machine_nm $1" . $2 . $3 . "\n";
		$no_jobs += 1;
	    }
	}
    }
}

close LOADLIST;

if (! $no_jobs) {
    print "No jobs found\n";
} else {
    print "$no_jobs total jobs\n";
}
