#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Net::SSH qw(ssh);

sub GetNodeList();
sub GetDirs(@);

# This script will access the borg cluster and look for jobs under the
# /temp1/$username

my (@nodes, $i, %dirHash, $outdata, $usrname);
my ($largest_f_len, $largest_list_len) = (0, 0);
if (@ARGV) {
    $usrname = $ARGV[0];
}

if (! $usrname) {
    $usrname = $ENV{USER};
}

print "Getting Node List..";
GetNodeList();
print "Done\n";

print "Obtaining Directory Listings...";
for $i (@nodes) {
    GetDirs($i);
}
print "Done\n";

@nodes = sort (keys %dirHash);

for $i (@nodes) {
    if (length($i) > $largest_list_len) {
	$largest_list_len = length($i);
    }
}

if ($#nodes < 0) {
    die "No files found for $usrname on the borg machines\n";
}

print "\nListing:\n";
printf "%-" . ($largest_f_len + 10) . "s Located on\n", "File";
for $i (0 .. ($largest_f_len + $largest_list_len + 20)) {
    print "-";
}
print "\n";


for $i (@nodes) {
    printf "%-" . ($largest_f_len + 10) . "s " . $dirHash{$i} . "\n", $i;
}

sub GetNodeList() {

    my ($i, $node_list_file);
    if (@ARGV and $#ARGV > 1) {
	$node_list_file = $ARGV[0];
    } else {
	if (system("pbsnodes -a >& _node_file")) {
	    die "Error: No node file was specified and cannot execute pbsnodes file: $!\n";
	}
	$node_list_file = "./_node_file";
    }
    open NODELIST, $node_list_file or die "Cannot open $node_list_file: $!\n";
    while (<NODELIST>) {
	chomp;
	if ($_ =~ /^(node\d+\-\d+)\.wag\.caltech\.edu/) {
	     push @nodes, $1;
	}
    }
    close NODELIST;

    if ($#nodes < 0) {
	die "ERROR: No valid node specification found in $node_list_file\n";
    }

}

sub GetDirs(@) {

    my ($curr_node) = $_[0];
    my ($tmp_file, $cmd, $returnval);

    $tmp_file = "/ul/" . $ENV{USER} . "/tmp_file";
    $cmd = "ls /temp1/";
    $cmd .= $usrname . "/ > $tmp_file";

   my $username = $ENV{USER};
#    print "$cmd\n";
    ssh("$username\@$curr_node", "$cmd");
 

    if (open INFILE, $tmp_file) {
	while (<INFILE>) {
	    chomp;
	    if ($_ =~ /\w+/) {
		$dirHash{$_} .= "$curr_node ";
		if (length($_) > $largest_f_len) {
		    $largest_f_len = length($_);
		}
	    }
	}
	close INFILE;
	system "rm -f $tmp_file";
    }
}
