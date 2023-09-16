#!/bin/bash

if [ $# -lt 2 ]; then
	echo "usage: $0 pmf_lmp_file window_size [force_constant=5] [save_prefix]"
	exit 1
fi

if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot access $1"
	exit 1
fi
pmf_lmp_file=$1

if [ $(echo $2 | awk '{if($1 ~ /^[0-9]*$/) print 1; else print 0}') -eq 0 ]; then
	echo "ERROR: Expected positive number for window size. Got '$2'"
	exit 1
fi
window_size=$2

force_constant=5
if [ $# -gt 2 ]; then
	force_constant=$3
fi

savePrefix=$(basename $1)
savePrefix=${savePrefix%.*}
if [ $# -gt 3 ]; then
	savePrefix=$4
fi

echo "FILENAME: $pmf_lmp_file"
tot=$(wc -l $pmf_lmp_file | awk '{print $1}')
comments=$(grep -c '^#' $pmf_lmp_file)
lines=$(( $tot - $comments ))
step=($(head -10 $pmf_lmp_file | egrep -v '^#' | sed 's/:.*$//' | awk '{print $1}'))

incr=$(( ${step[2]} - ${step[1]}))
num=$(echo $lines $incr $window_size | awk '{print int($1*$2/$3)}')
if [ $num -lt 2 ]; then
	echo "ERROR: Need at least two windows!" 
	exit 1
fi

cat <<DATA
STATS:
    #lines in ${1}: $lines
    TSTEP increment: $incr
    Window Size: $2
    Num windows: $num
DATA

echo
echo "1. Getting timeseries data from $pmf_lmp_file"
bounds=($(egrep -v '^#' $pmf_lmp_file | awk 'BEGIN{min=999999999999999; max=-min}{if($2>max)max=$2;if($2<min)min=$2}END{printf "%d %d",int(min)-1,int(max)+1}'))
incr=$(( $2 / $incr ))
echo "    bounds: ${bounds[*]}"
echo "    incr: $incr"
awk -v incr=$incr -v prefix="$savePrefix" '
BEGIN{
	c=0; 
	i=1; 
	outfile=sprintf("%s.1.timeseries.dat",prefix)
}
{
	if($1 !~ /^#/) { 
		c++; 
		if(c>=incr/2) 
			printf "%f %f %f\n",$1,$2,$5 >> outfile; 
		if(c==incr) { 
			c = 0; 
			i++; 
			outfile=sprintf("%s.%d.timeseries.dat",prefix,i); 
		}
	}
}' $pmf_lmp_file

echo "2. Creating histograms"
rm -fr ${savePrefix}.wham.in
rm -fr ${savePrefix}.hist.dat
for i in $(seq 1 $num)
do
	fle=${savePrefix}.${i}.timeseries.dat
	#min=$(sort -n -k 2 ${savePrefix}.${i}.timeseries.dat | head -1 | awk '{printf "%f\n",$2}')
	#min=$(tail -1 ${savePrefix}.${i}.timeseries.dat | awk '{print $3}')
	awk -v bins=200 -v min=${bounds[0]} -v max=${bounds[1]} -v col=2 -f ~tpascal/scripts/histogram.awk $fle > __tmp.dat
	min=$(cat __tmp.dat | awk 'BEGIN{y_max=0}{if($1 ~ /^[0-9]/ && $2>y_max){y_max=$2;x_max=$1}}END{print x_max}')
	cat >> ${savePrefix}.hist.dat <<DATA


#${i}
DATA
	echo >> ${savePrefix}.hist.dat
	cat __tmp.dat >> ${savePrefix}.hist.dat 
	echo "${savePrefix}.${i}.timeseries.dat $min $force_constant" >> ${savePrefix}.wham.in
done

echo "3. Running wham analysis on data bounds ${bounds[*]}"
bounds=($(egrep -v '^#' ${savePrefix}.wham.in | awk 'BEGIN{min=999999999999999; max=-min}{if($2>max)max=$2;if($2<min)min=$2}END{printf "%d %d",int(min),int(max)+1}'))
n=$(echo ${bounds[*]} | awk '{print int(($2-$1)*10)}')
~tpascal/codes/bin/wham ${bounds[*]} $n 0.01 298 0 ${savePrefix}.wham.in ${savePrefix}.pmf.wham.dat  > /dev/null

echo "4. Creating gnuplot file ${savePrefix}.hist.plt"
cat > ${savePrefix}.hist.plt <<DATA
datafile = "${savePrefix}.hist.dat"
stats datafile
pl for [ind=0:(STATS_blocks-1)] datafile i ind u 1:2:-2 w l t columnheader(1) lc variable
DATA

sed -i '$ d' ${savePrefix}.hist.dat
sed -i '$ d' ${savePrefix}.hist.dat
rm -fr __tmp.dat
