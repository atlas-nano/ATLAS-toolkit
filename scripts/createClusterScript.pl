#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester);
use Cwd;
use Cwd 'abs_path';

sub init;
sub createScriptFile;
sub showUsage;

my ($OPTS, $scriptFile, $reaxFF);

$reaxFF = "";
$|++;
&init;
print "Create Cluster script file $scriptFile...";
&createScriptFile($OPTS, $scriptFile);
print "Done\n";

sub createScriptFile {
    my ($opts, $fileName) = @_;
    my ($currDir) = &getcwd();
    my ($prefix, $bgffile, $USER, $forcefield, $inFile, $dataFile);
    
    ($prefix, $bgffile, $forcefield, $inFile, $dataFile) = 
	($opts->{p}, abs_path($opts->{b}), $opts->{f}, $opts->{i}, $opts->{d});
    $USER = $ENV{USER};

    open SCRIPTFILE, "> $fileName" or die "ERROR: Cannot write to $fileName: $!\n";
    print SCRIPTFILE <<DATA;
#!/bin/bash
#Job name:
#SBATCH --job-name=${prefix}
#
#Partition:
#SBATCH --partition=shared
#
#Account
#SBATCH --account=csd626
#
#Nodes
#SBATCH --nodes=1
#
#Processors:
#SBATCH --ntasks-per-node=16
#
#Wall clock limit:
#SBATCH --time=10:0:0
#
#SBATCH --mail-type=ALL
#

prefix=${prefix}
rtemp=298
press=1
curr_dir=\$SLURM_SUBMIT_DIR
temp_dir=/expanse/lustre/scratch/\$USER/temp_project/md/lammps/\${prefix}/\${rtemp}K

lmp_equil_file=in.\${prefix}
lmp_data_file=data.\${prefix}

module purge
#module load cpu/0.15.4 slurm gcc/10.2.0 openmpi cmake gsl intel-mkl amdfftw
module load cpu/0.17.3b  gcc/10.2.0/npcyll4 cmake/3.21.4/teqow32 openmpi/4.1.3/oq3qvsv fftw/3.3.10/bmvqsbr python/3.8.12/7zdjza7 intel-mkl/2020.4.304/ghfk3mu gsl/2.7/wtlsmyy gnuplot/5.4.2/mfinpvw sdsc slurm

nprocs=\$(( \$SLURM_NTASKS_PER_NODE * \$SLURM_NNODES / 2 ))
PARALLEL="mpirun -n \$nprocs -mca btl vader,self"

LMP="/home/tpascal/codes/bin/lmp_expanse -screen none -var rtemp \$rtemp -var press \$press"

mkdir -p \$temp_dir/analysis \${curr_dir}/results
cd \$temp_dir
for i in \$lmp_data_file \$lmp_equil_file 
do
	cp \${curr_dir}/\${i} ./
done

echo "LAMMPS dynamics of \${prefix} at \${rtemp}K"
echo "running in \$temp_dir"
\$PARALLEL \$LMP -in \${lmp_equil_file} -log \${prefix}.\${rtemp}K.equil.lammps.log
cp *.log *.lammpstrj \${curr_dir}/results

DATA

close SCRIPTFILE;
}

sub init {
    my ($i);
    getopt('pbfidsr', \%{ $OPTS });
    for $i ("p", "b", "f", "i", "d") {
	die &showUsage . "\n" if (! exists($OPTS->{$i}));
    }
    print "Initializing...";
    for $i ("b", "f", "i", "d") {
	if ($i eq "f") {
	    if ($OPTS->{$i} =~ /\s/) {
		while ($OPTS->{$i} =~ /(\S+)/g) {
		    FileTester($1);
		}
	    } else {
		FileTester($OPTS->{$i});
	    }
	} else {
	    FileTester($OPTS->{$i});
	}
    }
    print "Done\n";
    $scriptFile = $OPTS->{s};
    $reaxFF = "cp $OPTS->{r} ./ffield.reax" if (exists($OPTS->{r}));
    if (! defined($scriptFile)) {
	$scriptFile = basename($OPTS->{b});
	$scriptFile =~ s/\.\w+$//;
	$scriptFile .= "_cluster.script";
    }
    print "Done\n";
}

sub showUsage {
    return "usage: $0 -p name of job -b bgf file -f forcefield file -i lammps input file -d lammps data file\n" .
	"Options\n\tprefix: This will be the name of the job. The files will be stored in /temp1/{USER}/{name of job}\n";
}
