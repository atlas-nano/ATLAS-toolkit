#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 'data_files' [column=2] [save_name]"
	exit(1)
endif

set list = ($1)
if ($#list == 1) then
	set list = (`cat $1`)
endif
set tot = `tail -1 $list[1] | awk '{print NF}'`

set col = 2
if ($#argv > 1) set col = $2

set save_name = `basename $list[1]`
set save_name = $save_name:r
set save_name = "${save_name}.avg.dat"
if ($#argv > 2) set save_name = $3

echo "col $col tot $tot"
paste $1 | egrep -v "#" | awk '{if ($1 ~ /^[0-9]+/) { tot = c = tot2 = 0; i = '$col'; while(i<NF) {tot += $i; tot2 += $i*$i; i+='$tot'; c++;} tot /= c; tot2 /= c; std=tot2-tot*tot; if(std>0) std=sqrt(std); else std = 0; print $1,tot,std}}' > $save_name
