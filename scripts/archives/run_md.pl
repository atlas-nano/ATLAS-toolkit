#!/usr/bin/perl -w

# This program will run molecular dynamics on a structure
# It does the following:
#   1) runs an initial minimization
#   2) run an equilibration
#   3) runs a production

# usage: run_md.pl topfile coord_file no_of_bases
use Sys::Hostname;

$top_file = "";
$coord_file = "";
$no_bases = 0;

if (!@ARGV || $#ARGV < 2) {
    die "usage: do_min.pl topfile coord-file no_of_bases\n";
} else {
    $top_file = $ARGV[0];
    $coord_file = $ARGV[1];
    $no_bases = $ARGV[2];
}

#check for files
if (! -e $top_file) {
die "Cannot open $top_file: File does not exist\n";
}

if(! -e $coord_file) {
die "Cannot open $coord_file: File does not exist\n";
}

$no_bases = $no_bases *2;
print "\n\n--====Start Program====\n\n";

$curr_dir = "/temp1/tpascal";
if (! -e $curr_dir) {
    $my_cmd = "mkdir /temp1/tpascal";
    system $my_cmd;
}

$curr_dir = "/temp1/tpascal/" . $coord_file . "-min";
my $curr_host = hostname;

if (open LOGFILE, ">> logfile.txt") {
    my $date = localtime;
    print LOGFILE "$date\n";
    print LOGFILE "Running Dynamics using the following:\n";
    print LOGFILE "Topology file: $top_file\n";
    print LOGFILE "Coordinate file: $coord_file\n";
    print LOGFILE "Results stored in $curr_dir on $curr_host\n\n";
    print LOGFILE "\n";
    close LOGFILE;
}

$my_cmd = "rm -fr $curr_dir";
system $my_cmd;

$my_cmd = "cp -R ~/md_simulations/procedure $curr_dir/";
system $my_cmd;

$curr_top_file = $top_file;
$curr_coord_file = $coord_file;

$my_cmd = "cp $coord_file $curr_dir";
system $my_cmd;

$my_cmd = "cp $top_file $curr_dir";
system $my_cmd;

chdir $curr_dir;

$curr_top_file = $curr_dir . "/" . $top_file;
$curr_coord_file = $curr_dir . "/" . $coord_file;

chdir "minimize";
$my_cmd = "do_min.pl . $curr_top_file $curr_coord_file $no_bases";
system $my_cmd;

$my_cmd = "cp step3.restrt ../dynamics/md_ntr.restrt";
system $my_cmd;

chdir "../dynamics";

$my_cmd = "do_md.pl $curr_top_file";
system $my_cmd;

$my_cmd = "cp md.restrt ../prduction/md_ntr.restrt";
system $my_cmd;

chdir "../production";
$my_cmd = "do_prod.pl $curr_top_file";
system $my_cmd;

print "\n\nSimulation Completed\n";

