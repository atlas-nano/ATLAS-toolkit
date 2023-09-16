#!/bin/tcsh
#!/bin/csh

if ($#argv < 4) then
	echo "usage: $0 'pwr_files' 'thermo_files' pwr_column thermo_column"
	exit(1)
endif

set pwr_files = (`ls $1`)
set thermo_files = (`ls $2`)
set pcol = $3
set tcol = $4

if !($#pwr_files == $#thermo_files) then
	echo "ERROR: #pwr_files ($#pwr_files) != #thermo_file ($#thermo_files)"
	exit(1)
endif

if ($#pwr_files < 2) then
	echo "ERROR: Need at least 2 pwr_files to average"
	exit(1)
endif

if !(`echo $pcol | egrep -c '^[0-9]+$'`) then
	echo "ERROR: Expected positive integer for column. Got '$pcol'"
	exit(1)
endif
if !(`echo $tcol | egrep -c '^[0-9]+$'`) then
	echo "ERROR: Expected positive integer for column. Got '$tcol'"
	exit(1)
endif

set nf = `tail -1 $pwr_files[1] | awk '{print NF}'`
set nmols = (`grep nmolecules $2 | awk -v col="$tcol" '{print $col}'`)
echo "#nmols: $nmols nfields $nf nfiles $#pwr_files"
awk -v nmols="$nmols" -v col="$pcol" -f ~/scripts/calc_col_avg.awk $1
