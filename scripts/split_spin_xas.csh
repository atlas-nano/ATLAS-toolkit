#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 xas_file"
	exit(1)
endif

foreach i ($1)
	if !(-e $i && -r $i) continue
	cp $i _tmp.dat
	echo " #SPIN" >> _tmp.dat
	set spin_loc  = (`grep -n "^\s*#SPIN" _tmp.dat | sed 's/:.*//'`)
	if ($#spin_loc < 2) continue
	foreach j (`seq 2 $#spin_loc`)
		@ k = $j - 1
		set stop = $spin_loc[$j]
		set start = $spin_loc[$k]
		@ offset = $stop - $start
		head -${stop} _tmp.dat | tail -${offset} > ${i}_${k}
	end
end
