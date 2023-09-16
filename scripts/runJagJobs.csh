#!/bin/tcsh
#!/bin/csh

if ($#argv < 1) then
    echo "usage: $0 filelist [procs_per_node=1] [max_processes=8]"
    exit(1)
endif

set fileList = `ls $1`
if ($#{fileList} == 0) then
    echo "Error: No valid files found while searching $1"
    exit(1)
endif

set ppn = 1
set nprocs = 8
if ($2) set ppn = $2
if ($3) set nprocs = $3
set myhost = `basename $HOST .local`
set running = `jaguar jobs running | grep -c ${myhost}`

setenv LM_LICENSE_FILE @10.0.0.1

foreach i ($fileList)
    set mol = `basename $i .01.in`
    set mol = `basename $mol .02.in`
    set mol = `basename ${mol} .in`
    set outfile = ${mol}.0.out
    if (-e $outfile) continue;
    set running = `jaguar jobs running | grep -c ${myhost}`
    #echo "${running} of ${nprocs} jobs running on ${myhost}"
    while ($running >= $nprocs)
	sleep 30s
	set running = `jaguar jobs running | grep -c $myhost`
    end

    echo "running jaguar singlepoint on ${mol}"
    if ($ppn == 1) then
	jaguar run ${mol} > /dev/null
    else
	jaguar run ${mol} -PROCS $ppn > /dev/null
    endif
    sleep 5s
end
