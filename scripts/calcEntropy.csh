#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
    echo "usage: $0 [vel file list] [datafile] [2pttype] [numprocs=8]"
    exit(1);
endif

set flist = "$1"
set datafile = $2
set nprocs = 8

if !((-e $datafile) && (-r $datafile)) then
    echo ERROR: Cannot find LAMMPS datafile $datafile
    exit(1)
endif
set type2pt = 1
if ($#argv > 2) then
    set type2pt = $3
endif
set nprocs = 8
if ($#argv > 3) then
    set nprocs = $4
endif

set list = `ls $flist |grep .lammps`
if ($#list == 0) then
    echo ERROR: No valid files found while searching $1
    exit(1)
endif

set done = 0
set counter = 0
set offset = 0;
set iter = 1;

setenv LD_LIBRARY_PATH /ul/tpascal/vac:$LD_LIBRARY_PATH

echo Iteration $iter
while (! $done)
    if ($counter < $nprocs) then
        @ counter = $counter + 1;
        @ index = $counter + $offset
        if ($index <= $#list) then
	    set mol = `basename $list[$index] .lammps`
	    set dir = `dirname $list[$index]`
            set trj = $list[$index]
            set log = ${dir}/${mol}.eng
            if !((-e $trj) && (-r $trj) && (-e $log) && (-r $log)) then
		echo ERROR: Cannot find either logfile $log or traj file $trj
		exit(1)
	    endif
	    echo vac for $mol
	    csh -f /ul/tpascal/scripts/createinputfile.csh ${dir} ${mol} ${datafile} ${type2pt}
            /ul/tpascal/programs/md_analysis/driver ${mol}.in >& ${mol}.screen.out &
	else
	    wait;
	    set done = 1           
        endif
    else
        @ offset = $index
        set counter = 0
        if ($offset >= $#list) then
	    wait;
            set done = 1
        else
            @ iter = $iter + 1
	    wait;
            echo Iteration $iter
        endif
    endif
end

