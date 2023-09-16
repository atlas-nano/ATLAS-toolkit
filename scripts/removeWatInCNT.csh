#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 cnt_water_file cnt_atom_types [save_name]"
	exit(1)
endif

if !(-e $1 && -r $1) then
	echo "ERROR: Cannot access $1"
	exit(1)
endif
set savename = `basename $1`
set savename = $savename:r
set savename = "${savename}.new.bgf"
if ($#argv > 2) set savename = $3

set com = (`/home/tpascal/scripts/bgfCoM.pl -b $1 -a "fftype eq '$2'" | grep '^Center' | awk '{print $4,$5,$6}'`)
if !(${?com}) then
	echo "ERROR: Could not determine the center of mass of the CNT. Are you sure that the FFTYPE is '$2'?"
	exit(1)
endif
set bounds = (`/home/tpascal/scripts/getBounds.pl -b $1 -o "moleculeid==1" | egrep '^[X|Y|Z]' | awk '{print $2,$3}'`)
echo $com
echo "removing all water molecules enclosed by CNT in $1"
/home/tpascal/scripts/getBGFAtoms.pl -m 1 -b $1 -s $savename -o "fftype eq '$2' or zcoord<($bounds[5]) or zcoord>($bounds[6]) or sqrt((xcoord-($com[1]))**2+(ycoord-($com[2]))**2)>($bounds[2]-($bounds[1]))/2" > /dev/null || goto error
echo "Done"
exit(0)

error:
echo "ERROR Occurred"
exit(1)
