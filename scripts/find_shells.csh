#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 data_file (col=2) (num_shells=3) (symmetric=no)"
	exit(1)
endif

if !(-e $1 && -r $1) then
	echo "ERROR: Cannot access $1"
	exit(1)
endif

set col = 3
if ($#argv > 1) set col = $2


set num_shells = 3
if ($#argv > 2) set num_shells = $3

set is_sym = 0
if ($#argv > 3) set is_sym = $4

set nlines = `wc -l $1 | awk '{print $1}'`
set nskip = `cat $1 | tail -10 | head -2 | awk '{printf "%s ",$1}END{printf "\n"}' | awk '{print $2-$1}' | awk '{printf "%.0f\n",2/$1}'`

#left side
set start = `cat $1 | awk '$1 ~ /[0-9]/ {print NR}' | head -1`
set mid = $nlines
if ($is_sym) set mid = `echo $start $nlines | awk '{val = $1+($2-$1)/2; printf "%.0f\n",val}'`
@ start_block = $mid - $start

@ valley = 0
@ count = 0
set p_str = (`cat $1 | head -${mid} | tail -${start_block} | awk '{print $'$col',$1,NR+'$start'}' | sort -rn | awk '{print $2,$1,$3}' | head -1`)
while($count < $num_shells)
	#now find nearest peak after current peak
	set next_peak = (`cat $1 | head -${mid} | awk '{if((NR>'$p_str[3]'+'$nskip')&&(NR<'$p_str[3]'+'$nskip'*2)) print $'$col',$1,NR}' | sort -rn | head -1 | awk '{print $2,$1,$3}'`)
	#echo "peak at $p_str next peak at $next_peak"
	@ toffset = $next_peak[3] - $p_str[3]
	#now find valley
	set v_str = (`cat $1 | head -${next_peak[3]} | tail -${toffset} | awk '{ print $'$col',$1,NR}' | sort -n | awk '{print $2,$1,$3}' | head -1`)
	set valley = `echo $v_str[1]`
	set peak = $p_str[1]
	echo $peak $valley
	set p_str = ($next_peak)
	@ count = $count + 1
end

if !($is_sym) exit 0

#right side
set start = `tac $1 | awk '$1 ~ /[0-9]/ {print NR}' | head -1`
set mid = `echo $start $nlines | awk '{val = $1+($2-$1)/2; printf "%.0f\n",val}'`
@ start_block = $mid - $start

@ valley = 0
@ count = 0
set p_str = (`tac $1 | head -${mid} | tail -${start_block} | awk '{print $'$col',$1,NR+'$start'}' | sort -rn | awk '{print $2,$1,$3}' | head -1`)
while($count < $num_shells)
	#now find nearest peak after current peak
	set next_peak = (`tac $1 | head -${mid} | awk '{if((NR>'$p_str[3]'+'$nskip')&&(NR<'$p_str[3]'+'$nskip'*2)) print $'$col',$1,NR}' | sort -rn | head -1 | awk '{print $2,$1,$3}'`)
	#echo "peak at $p_str next peak at $next_peak"
	@ toffset = $next_peak[3] - $p_str[3]
	#now find valley
	set v_str = (`tac $1 | head -${next_peak[3]} | tail -${toffset} | awk '{ print $'$col',$1,NR}' | sort -n | awk '{print $2,$1,$3}' | head -1`)
	set valley = `echo $v_str[1]`
	set peak = $p_str[1]
	echo $peak $valley
	set p_str = ($next_peak)
	@ count = $count + 1
end
