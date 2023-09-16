#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
	echo "usage: $0 pmf_lmp_file window_size [force_constant=5] [save_prefix]"
	exit(1)
endif

if !(-e $1 && -r $1) then
	echo "ERROR: Cannot access $1"
	exit(1)
endif

set pmf_lmp_file = $1
set window_size = $2
if !(`echo $2 | awk '{if($1 ~ /^[0-9]*$/) print 1; else print 0}'`) then
	echo "ERROR: Expected positive number for window size. Got '$2'"
	exit(1)
endif

set force_constant = 5
if ($#argv > 2) set force_constant = $3
set savePrefix = `basename $1`
set savePrefix = $savePrefix:r
if ($#argv > 3) set savePrefix = $4

echo "FILENAME: $pmf_lmp_file"
set tot = `wc -l $pmf_lmp_file | awk '{print $1}'`
set comments = `grep -c '^#' $pmf_lmp_file`
@ lines = $tot - $comments
set step = (`head -10 $pmf_lmp_file | egrep -v '^#' | sed 's/:.*$//' | awk '{print $1}'`)
@ incr = $step[2] - $step[1]
set num = `echo $lines $incr $window_size | awk '{print int($1*$2/$3)}'`
cat <<DATA
STATS:
	#lines in ${1}: $lines
	TSTEP increment: $incr
	Window Size: $2
	Num windows: $num
DATA

if ($num < 2) then
	echo "ERROR: Need at least two windows!" 
	exit(1)
endif

echo
echo "1. Getting timeseries data from $pmf_lmp_file"
set bounds = (`egrep -v '^#' $pmf_lmp_file | awk 'BEGIN{min=999999999999999; max=-99999999999999999}{if($2>max)max=$2;if($2<min)min=$2}END{printf "%d %d",int(min)-1,int(max)+1}'`)
@ incr = $2 / $incr
egrep -v '^#' $pmf_lmp_file | awk -v incr=$incr -v prefix="$savePrefix" 'BEGIN{c=0; i=1; outfile=sprintf("%s.1.timeseries.dat",prefix)}{c++; if(c>=incr/2) printf "%f %f %f\n",$1,$2,$5 >> outfile; if(c==incr) { c = 0; i++; outfile=sprintf("%s.%d.timeseries.dat",prefix,i); }}'

echo "2. Creating histograms"
rm -fr ${savePrefix}.wham.in
rm -fr ${savePrefix}.hist.dat
foreach i (`seq 1 $num`)
	set fle = ${savePrefix}.${i}.timeseries.dat
	set min = `sort -n -k 2 ${savePrefix}.${i}.timeseries.dat | head -1 | awk '{printf "%f\n",$2}'`
	set min = `tail -1 ${savePrefix}.${i}.timeseries.dat | awk '{print $3}'`
	awk -v bins=200 -v min=$bounds[1] -v max=$bounds[2] -v col=2 -f ~tpascal/scripts/histogram.awk ${savePrefix}.${i}.timeseries.dat > __tmp.dat
	#set min = `cat __tmp.dat | awk 'BEGIN{y_max=0}{if($1 ~ /^[0-9]/ && $2>y_max){y_max=$2;x_max=$1}}END{print x_max}'`
	echo "#$i" >> ${savePrefix}.hist.dat 
	cat __tmp.dat >> ${savePrefix}.hist.dat 
	echo >> ${savePrefix}.hist.dat
	echo >> ${savePrefix}.hist.dat
	echo "${savePrefix}.${i}.timeseries.dat $min $force_constant" >> ${savePrefix}.wham.in
end

echo "3. Running wham analysis on data bounds $bounds"
set bounds = (`egrep -v '^#' ${savePrefix}.wham.in | awk 'BEGIN{min=999999999999999; max=-99999999999999999}{if($2>max)max=$2;if($2<min)min=$2}END{printf "%d %d",int(min),int(max)+1}'`)
~tpascal/codes/bin/wham $bounds 100 0.01 298 0 ${savePrefix}.wham.in ${savePrefix}.pmf.wham.dat  > /dev/null

echo "4. Creating gnuplot file ${savePrefix}.hist.plt"
cat > ${savePrefix}.hist.plt <<DATA
datafile = "${savePrefix}.hist.dat"
stats datafile
pl for [ind=0:(STATS_blocks-1)] datafile i ind u 1:2:-2 w l t columnheader(1) lc variable
DATA

sed -i '$ d' ${savePrefix}.hist.dat
sed -i '$ d' ${savePrefix}.hist.dat
rm -fr __tmp.dat
