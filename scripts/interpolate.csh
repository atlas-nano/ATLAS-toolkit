#!/bin/tcsh
#!/bin/csh

if ($#argv < 4) then
	echo "usage: $0 file x_start x_stop n_pts (savename)"
	exit(1)
endif

if !(-e $1 && -r $1) then
	echo "ERROR: Cannot access $1"
	exit(1)
endif

set savename = `basename $1`
set savename = $savename:r
if ($#argv > 4) set savename = $5

cat > __xmgrace.script <<DATA;
read "$1"
interpolate(S0,mesh($2,$3,$4),SPLINE,ON)
write G0.S1 file "$savename"
exit
DATA

/usr/bin/grace -nosafe -batch __xmgrace.script
rm -fr __xmgrace.script
