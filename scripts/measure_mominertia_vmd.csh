#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
	echo "usage: $0 bgf_file 'lmp_trj_file(s)' (selection=all) (Savename)"
	exit(1)
endif

if !(-e $1) then
	echo "ERROR: Cannot locate $1"
	exit(1)
endif
set trj_list = ($2)
foreach i ($trj_list)
	if !(-e $i) then
		echo "ERROR: Cannot locate $i"
		exit(1)
	endif
end

set sel = "all"
if ($#argv > 2) set sel = "$3"

set savename = `basename $1`
set savename = $savename:r
set savename = "${savename}.rgyr.dat"
if ($#argv > 3) set savename = $4

echo "mol new $1 type bgf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all" > __vmd.cmd
foreach i ($trj_list)
	set suf = $i:e
	set type = "lammpstrj"
	if ($suf == "mdcrd") set type = "crdbox"
	echo "mol addfile $i type $type first 1 last -1 step 1 filebonds 1 autobonds 1 waitfor all" >> __vmd.cmd
end

cat >> __vmd.cmd <<DATA;
package require pbctools
animate delete beg 0 end 0
pbc wrap -all -compound fragment -center com -centersel "resid 134"
set mol [molinfo top] 
set sel [atomselect top "$sel"]
set n [molinfo top get numframes]
set output [open "__pI.dat" w]
for {set i 1} {\$i <\$n} {incr i} {
	molinfo top set frame \$i
	set inertia [lindex [measure inertia \$sel] 0]
	puts \$output "\$i \$inertia"
}
close \$output
DATA

vmd -dispdev none -nt < __vmd.cmd  > /dev/null
if !(-e __pI.dat) then
	echo "ERROR occurred"
endif

mv __pI.dat $savename
egrep '^[0-9]' ${savename} | awk -f /home/tpascal/scripts/calcStats.csh
