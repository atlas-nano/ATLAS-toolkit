#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
    echo "usage: $0 prefix datafile [slab options] [orient shells] [orient increment]"
    exit(1)
endif

set prefix = $1
if (! -e ${prefix}.lammps) then
    if !( -e ${prefix}.lammpstrj) then
	echo "ERROR: Cannot find ${prefix}.lammps or ${prefix}.lammpstrj"
	exit(1)
    endif
    set mol = `basename $prefix`
    ln -s ${prefix}.lammpstrj ${mol}.lammps
    set prefix = $mol
endif
set mol = `basename $prefix`

set datafile = $2
if (! -e $2) then
    echo "ERROR: Cannot find data file $2"
    exit(1)
endif

set slab_opt = ""
set orient_shells = ""
set orient_ainc = ""
if ($#argv>2) then
    set slab_opt = "ANALYSIS_SLAB_SHAPE	$3"
endif
if ($#argv>3) then
    set orient_shells = "ANALYSIS_SLAB_ORIENT_SHELLS	$4"
    if($#argv>4) then
	set orient_ainc = "ANALYSIS_SLAB_ORIENT_AINCR	$5"
    endif
endif

cat > ${mol}.in <<DATA;
IN_LMPDATA                    ${datafile}
IN_LMPTRJ                     ${prefix}
ANALYSIS_FRAME_INITIAL        1
ANALYSIS_FRAME_FINAL          0
ANALYSIS_FRAME_STEP           1
ANALYSIS_OUT                  ${mol}
ANALYSIS_SLAB_THICKNESS       0.5
ANALYSIS_SLAB_EWALD_ACC	      0.001
ANALYSIS_SLAB_CENTER_COM      1
${slab_opt}
${orient_shells}
${orient_ainc}
DATA
