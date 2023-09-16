#!/bin/tcsh
#!/bin/csh

if ($#argv < 1) then
    echo "usage: $0 [log file list] [numprocs=8]"
    exit(1);
endif

set flist = $1
set nprocs = 8
if ($#argv > 1) then
    set nprocs = $2
endif

set list = `ls $flist |grep .lammps`
if ($#list == 0) then
    echo ERROR: No valid files found while searching $2
    exit(1)
endif

set done = 0
set counter = 0
set offset = 0;
set iter = 1;

echo Iteration $iter
while (! $done)
    if ($counter < $nprocs) then
        @ counter = $counter + 1;
        @ index = $counter + $offset
        if ($index <= $#list) then
	    set mol = `basename $list[$index] _log.lammps`
	    set log = $list[$index]
	    set save = ${mol}_thermo.dat
	    echo calculating thermo for $mol
	    /ul/tpascal/scripts/getLammpsLogInfo.pl -f "temp poteng volume press" -t "*" -a 0 -s $save -l $log > /dev/null &
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

