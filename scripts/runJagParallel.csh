#!/bin/csh
#!/bin/tcsh

if ($#argv < 1) then
echo "usage: $0 [jaguar prefix| input file] [nprocs]"
exit(1)
endif

if ($#argv > 1) then
    set nprocs = $2
else
    set nprocs = 8
endif

set molname = `basename $1 .in`

cat >>${molname}.script <<DATA;
#PBS -l nodes=1:ppn=${nprocs}
#PBS -l walltime=960:00:00
#PBS -q workq
#PBS -e ${PWD}/${molname}_jag.err
#PBS -o ${PWD}/${molname}_jag.out
#PBS -N ${molname}_jag

echo "Nodes:"
cat \$PBS_NODEFILE
echo "Jaguar run of ${molname}"
echo

# For jaguar on cluster uses a special license server
setenv LM_LICENSE_FILE \@10.0.0.1

# borg's nodes aren't in jaguar.hosts file so define temp dir
setenv JAGUAR_TEMP /temp1/${USER}

#set up variables for parallel run
set mpi_used = "/opt/mpich-1.2.5.10-ch_p4-gcc"
set port = 2345
set procs = ${nprocs}
set path=(\${mpi_used}/bin \$path)

setenv SCHRODINGER_NODEFILE \$PBS_NODEFILE
setenv MPI_USEP4SSPORT yes
setenv MPI_P4SSPORT \$port
setenv SCHRODINGER_MPI_START yes

\$SCHRODINGER/utilities/mpich start -m \$PBS_NODEFILE -p \$port

cd ${PWD}

jaguar run ${molname} -PROCS \$procs -WAIT || goto error

echo " No errors detected - All Tasks Completed "
exit(0)

error:
echo "Error occurred."
exit(1)
DATA
