#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
  echo "usage: $0 dens_file (bulk_dens) (col=2)"
  exit(1)
endif

if !(-e $1) then
  echo "ERROR: Cannot locate $1"
  exit(1)
endif

#source /etc/profile.d/modules.csh
#module load cpu/0.15.4  gcc/10.2.0 gnuplot
set col = 2
if ($#argv > 2) set col = $3
set flen = `wc -l $1 | awk '{print $1-5}'`
set extrema = (`tail -${flen} $1 | awk '{if($'$col'>0.01) { if(NR==1) { tot = 0; } else if (tot ==0) { min = $1; tot = 1;} else if (tot == 1) { max = $1; } } } END{ print min,max,min+(max-min)/2}'`)
set mid = $extrema[3]
set hi = `echo $extrema[2] | awk '{print $1+10}'`
set lo = `echo $extrema[1] | awk '{print $1-10}'`
set mol = `basename $1`
set mol = $mol:r

cat > ${mol}.gnuplot.script <<DATA
fitP(x) = a1*(1+tanh(-((x-z1)/s1)))
fitN(x) = a2*(1-tanh(-((x-z2)/s2)))
a1 = 0.5
s1 = 0.15
z1 = $extrema[2]
fit [${mid}:${hi}] fitP(x) "$1" u 1:$col via a1,z1,s1

a2 = a1
s2 = s1
z2 = $extrema[1]
fit [${lo}:${mid}] fitN(x) "$1" u 1:$col via a2,z2,s2
#set term post enhanced color solid "Helvetica,30" lw 2
set term png size 1400,1050 enhanced color font "Heveltica,40" rounded crop lw 3
set out "${mol}.fit.png"
pl "$1" u 1:$col w lp lt -1 title "data", fitN(x) w l title "fitN", fitP(x) w l title "fitP"
DATA

if ($#argv > 1) then
cat >> ${mol}.gnuplot.script <<DATA
delta = 0.1
int1a(x,d) = (x<=d*.1) ? 0 : (int1a(x-d,d)+(f(x-d)+4*f(x-d*.5)+f(x))*d/6.)
int1b(x,d) = (x>=-d*.1) ? 0 : (int1b(x+d,d)+(f(x+d)+4*f(x+d*.5)+f(x))*d/6.)
int2(x,y,d) = (x>y-d*.5) ? 0 : (int2(x+d,y,d) + (f(x)+4*f(x+d*.5)+f(x+d))*d/6.)
integral_f(x) = (x>0)?int1a(x,x/ceil(x/delta)):-int1b(x,-x/ceil(-x/delta))
integral2_f(x,y) = (x<y)?int2(x,y,(y-x)/ceil((y-x)/delta)): -int2(y,x,(x-y)/ceil((x-y)/delta))

f(x) = fitN(x) - $2
set table '__fitN.dat'
pl [${lo}:${mid}] integral_f(x)
unset table

f(x) = fitP(x) - $2
set table '__fitP.dat'
pl [${mid}:${hi}] integral_f(x)
unset table
DATA

gnuplot < ${mol}.gnuplot.script >& ${mol}.gnuplot.out
set a1 = `egrep "a1\s*=" ${mol}.gnuplot.out | tail -1 | awk '{print $3}'`
set a2 = `egrep "a2\s*=" ${mol}.gnuplot.out | tail -1 | awk '{print $3}'`
set s1 = `egrep "s1\s*=" ${mol}.gnuplot.out | tail -1 | awk '{print $3}'`
set s2 = `egrep "s2\s*=" ${mol}.gnuplot.out | tail -1 | awk '{print $3}'`
set z1 = `egrep "z1\s*=" ${mol}.gnuplot.out | tail -1 | awk '{print $3}'`
set z2 = `egrep "z2\s*=" ${mol}.gnuplot.out | tail -1 | awk '{print $3}'`

if !(`echo $a1 $a2 | awk '{val = 1; if ($1 <= 0) val = 0; if ($2 <=0) val = 0; if (($1/$2)<0.9 || ($1/$2)>1.1) val = 0; print val}'`) then
else if !(`echo $s1 $s2 | awk '{val = 1; if ($1 <= 0) val = 0; if ($2 <=0) val = 0; if (($1/$2)<0.9 || ($1/$2)>1.1) val = 0; print val}'`) then
else if !(`echo $z1 $z2 | awk '{val = 1; if (($1-$2)<0) val = 0; print val}'`) then
endif

echo $a1 $s1 $z1 $a2 $s2 $z2
#rm -fr __fitP.dat __fitN.dat
