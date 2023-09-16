#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
	echo "usage: $0 n m length[nm] [savename]"
	exit(1)
endif

set n = $1
set m = $2
set l = $3
set savename = ${1}.${2}.${l}nm.bgf
if ($#argv > 3) set savename = $4

if !(`echo $n | egrep '^[0-9]+$'`) then
	echo "Expected digit for n, got '$n'"
	exit(1)
endif

if !(`echo $m | egrep '^[0-9]+$'`) then
	echo "Expected digit for n, got '$m'"
#exit(1)
endif

set num_test = `echo $l | egrep '^[0-9.]+$'`
if ($num_test == "") then
	echo "Expected numeric for length, got '$l'"
	exit(1)
endif

cat > __cnt.in <<DATA
1
$1 $2
DATA

#/home/tpascal/codes/bin/mw_nanotube < __cnt.in > __cnt.out
/home/tpascal/scripts/createNtube.pl -m $1 -n $2 -s ${1}_${2}_mw.bgf > /dev/null || goto error

if !(-e ${1}_${2}_mw.bgf) then
	echo "Error while create ${1} ${2} ntube unit cell. See __cnt.out"
	exit(1)
endif

mv ${1}_${2}_mw.bgf ${savename}
set d = `grep CRYSTX ${savename} | awk -v l=$l '{printf "%.0f\n", l*10/$4}'`
echo "replicating cell $d"
if ($d > 1) /home/tpascal/scripts/replicate.pl -b ${savename} -s ${savename} -d "1 1 $d" > /dev/null
#/home/tpascal/scripts/removeBond.pl -b ${savename} -s ${savename} -i "index>0" -j "index>0" > /dev/null
#/home/tpascal/scripts/bondByDistance.pl -b ${savename} -s ${savename} -d 1.6 > /dev/null
rm -fr __cnt.in __cnt.out

exit(0)

error:
echo "ERROR occurred"
exit(1)

