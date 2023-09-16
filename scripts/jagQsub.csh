#!/bin/tcsh
#!/bin/csh

if ($#argv < 1) then
echo "usage: $0 [jaguar prefix| input file] [jagtype=run] [num nodes = 1] [ppn=1]"
exit(1)
endif

set curr_dir = `echo $PWD`

set mol = `/ul/tpascal/scripts/getFileAbsPath.pl $1`
set molname = `basename $mol .in`
set inputfile = ${molname}.in
set maefile = ${molname}.mae
if (((-e $inputfile) == 0)) then
    echo "ERROR: Cannot find $inputfile or $maefile"
    exit(1)
endif

set jagType = "run"
if ($#argv > 1) then
    set jagType = $2
endif

set nodes = 1
set ppn = 1
set host = ""
set hostfile = ""
set tot_procs = 1
if ($#argv > 2) then
    set nodes = $3
    set ppn = 1
    set tot_procs = $3
    if ($#argv > 3) then
        set ppn = $4
        @ tot_procs = $nodes * $ppn
    endif
    set  jagType = "$jagType" 
endif

set procs_str = ""
if ($tot_procs > 1) then
    set procs_str = "-PROCS $tot_procs"
endif

cat > ${molname}.jaguar.pbs <<EOF
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -l walltime=96:00:00
#PBS -q workq
#PBS -j oe
#PBS -N ${molname}
#!/bin/csh
#!/bin/tcsh

# For jaguar on cluster uses a special license server
setenv LM_LICENSE_FILE @10.0.0.1

# borg's nodes aren't in jaguar.hosts file so define temp dir
setenv JAGUAR_TEMP /temp1/`echo $USER`

echo "Nodes:"
cat \$PBS_NODEFILE
echo "Jaguar run of ${molname}"
echo

cd \$PBS_O_WORKDIR

set procs = \$PBS_NUM_PPN
set procs_str = ""
set hosts = \`cat \$PBS_NODEFILE | tr '\n' ' '\`
set nodes = \`echo \$PBS_NODEFILE | wc -l | awk '{print \$1}'\`
if(\$nodes>1) then
  set procs = \`echo \$procs \$nodes | awk '{print \$1*\$2}'\`
endif
if(\$procs>1) then
  set procs_str = "-PARALLEL \$procs -HOST '\$hosts'"
  set procs_str = "-PARALLEL \$procs"
endif
jaguar $jagType \$procs_str ${molname}.in -WAIT || goto error

echo " No errors detected - All Tasks Completed "
exit(0)

error:
echo "Error occurred. See ${molname}.out"
exit(1)

EOF
