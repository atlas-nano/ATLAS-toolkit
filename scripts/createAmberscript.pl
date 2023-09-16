#!/usr/bin/perl -w
# This script will take a topology and coordinate file and create a
# Parm and Script file for use on the Matrix or Borg clusters

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester IsInteger);
use File::Basename qw(basename dirname);
use Cwd qw(realpath);

sub Initilize;
sub CreateParmFile;
sub CreateScriptFile;

die "usage: $0 topology_file coordinate_file solute_bases [save_name]\n"
    if (! @ARGV or $#ARGV < 2);

my ($top_file, $coord_file, $base_num, $save_name) = @ARGV;

Initilize;

print "Creating Parameter file $save_name" . ".parm...";
CreateParmFile($save_name . ".parm");
print "Done\nCreating cluster script file $save_name" . ".script...";
CreateScriptFile($save_name . ".script", $save_name . ".parm");
print "Done\n";

sub CreateParmFile {
    my ($parmFile) = $_[0];

    open PARMFILE, "> $parmFile" or die "Cannot create parameter file $parmFile: $!\n";
    print PARMFILE <<DATA;
Mol: 6:5
Total bases: $base_num
Bases at end: 5
Crossovers: -1
Topology file: $top_file
Trajectory file: $coord_file
3 prime in: 1
Directory Name: $save_name
Production Steps: 250000
Production Traj Dump: 2500

DATA

close PARMFILE;
    



}


sub CreateScriptFile {
    my ($fileName, $parmFile) = @_;
    my ($err_file, $out_file);

    $err_file = dirname($top_file) . "/" . $save_name . ".err";
    $out_file = dirname($top_file) . "/" . $save_name . ".out";
    $parmFile = realpath("./" . $parmFile);

    open SCRIPTFILE, "> $fileName" or die "Cannot create script file $fileName: $!\n";
    
    print SCRIPTFILE <<DATA;
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1000:00:00
#PBS -q workq
#PBS -e $err_file
#PBS -o $out_file
#PBS -N $save_name
#!/bin/tcsh
#!/bin/csh
                                                                                                                  
echo "Nodes:"
cat \$PBS_NODEFILE
echo "JOB: $save_name"
echo

set curr_dir = $ENV{PWD}
set parm_file = $parmFile
set nodes = 1
set nprocs = `wc -l < \$PBS_NODEFILE`

setenv DO_PARALLEL ""
setenv executable "/exec/amber9/exe/sander"

cat \$PBS_NODEFILE > \${curr_dir}/nodelist.\$PBS_JOBID

/ul/tpascal/scripts/dna_sim.pl  \$parm_file \$executable || goto error

exit:
echo " "
echo No errors detected
exit(0)
                                                                             
exit:
error:
echo "  FAILED.  Sorry bud."
cho " "
exit(1)

DATA

close SCRIPTFILE;
    


}

sub Initilize {
    FileTester($top_file);
    FileTester($coord_file);
    
    $save_name = $top_file
	if (! defined($save_name));
    
    if ($save_name =~ /(.+)\.\w{3}$/) {
	$save_name = $1;
    }

    die "Error: Expected integer for variable \"solute_bases\". Got $base_num\n"
	if (! IsInteger($base_num));

    $top_file = realpath($top_file);
    $coord_file = realpath($coord_file);

}
	
