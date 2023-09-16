#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 datafile (savename) (xrange)"
	exit(1)
endif

if !(-e $1 && -r $1) then
	echo "ERROR: Cannot access $1"
	exit(1)
endif

set sname = `basename $1`
set sname = $sname:r".gaussfit.dat"
if ($#argv > 1) set sname = $2

set xrange = ""
set vars = (`cat $1 | awk 'BEGIN{max=-999999999999999999999999; start = 0}{if(start==0&&$2>0) { start = $1; } if($2>max){max=$2;xval=$1;}if($2>0) end=$1;}END{print xval,max,(xval-start)/4,start-.1,end+.1}'`)
set xrange = "[$vars[4] : $vars[5]]"
if ($#argv>2) set xrange = "$3"
cat > fit.plt <<DATA;
set fit quiet
gauss(x)=a*exp(-(x-b)**2/2/c/c)
a = $vars[2]
b = $vars[1]
c = $vars[3]
set xrange $xrange
fit gauss(x) "$1" via a,b,c
set samples 1000
set table "$sname"
pl gauss(x)
unset table
DATA

gnuplot < fit.plt

set a = `egrep '^a\s*=' fit.log | tail -1 | awk '{print $3}'`
set b = `egrep '^b\s*=' fit.log | tail -1 | awk '{print $3}'`
set c = `egrep '^c\s*=' fit.log | tail -1 | awk '{print $3}'`
echo $b +/- $c $a

exit(0)

