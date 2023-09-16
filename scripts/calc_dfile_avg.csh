#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 data_files (save_name)"
	exit(1)
endif

set list = ($1)
if ($#list == 1) then
	set list = (`cat $1`)
endif
if ($#list == 0) then
	echo "ERROR: No valid files found while searching $1"
	exit(1)
endif
set fn_prefix = `set string1="$list[1]"; set string2="$list[$#list]"; printf "%s\n%s\n" "$string1" "$string2" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' | sed 's/\.\w*$//'`
set fn_prefix = `basename $fn_prefix`
set save_name = "${fn_prefix}.pwr"
if ($#argv > 1) set save_name = $2

set nfields = `tail -1 $list[1] | awk '{print NF}'`
echo $nfields
paste $list | awk '{if($2==($2+0)) { printf "%s ", $1; for(i=1;i<'$nfields';i++) {x1=0; x2=0; j = i + 1; for(k=0; k<'$#list'; k++) { l = j + k*'$nfields'; x1 += $l; x2 += $l*$l;} x1 /= '$#list'; x2 /= '$#list'; sigma = x2-x1*x1; if(sigma>0) sigma=sqrt(sigma); else sigma = 0; printf "%s %s ",x1,sigma} printf "\n"}}' > $save_name
