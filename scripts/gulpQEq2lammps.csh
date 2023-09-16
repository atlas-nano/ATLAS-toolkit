#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 gulp_qeq_input"
	exit(1)
endif

if !(-e $1 && -r $1) then
	echo "ERROR: Cannot access $1"
	exit(1)
endif

if !(`grep -c '^qelec' $1`) then
	echo "ERROR: Cannot locate QEq section (searching for 'qelec') in $1"
	exit(1)
endif

set tot = `wc -l $1 | awk '{print $1}'`
set start = `grep -n '^qelec' $1 | tail -1 | sed 's/:.*$//'`
@ offset = $tot - $start
set c = 0
foreach i (`tail -${offset} $1 | awk '{print $1}'`)
	if !(`grep -ic " ${i} " /home/tpascal/scripts/Packages/elementList.txt`) then
		break
	endif
	@ c = $c + 1
end

set n_str = ()
foreach i (`tail -${offset} $1 | head -${c} | awk '{print $1}'`)
	set ele = `echo $i | sed 's/^.*\(\w+\).*/\1/' | sed 's/[0-9]//g'`
	set n_str = ($n_str `grep " $ele " /home/tpascal/scripts/Packages/elementList.txt | awk '{i=NF-1; print $i}' | sed 's/^\([0-9]*\).*$/\1/'`)
end

tail -${offset} $1 | head -${c} | awk -v n_str="$n_str" 'BEGIN{split(n_str,n," ")}{printf "%-3s %8.3f %8.3f %8.3f %8.3f\n",$1,$2*23.06,$3*23.06*2,(2*n[NR]+1)/$4/4,0.00}'
