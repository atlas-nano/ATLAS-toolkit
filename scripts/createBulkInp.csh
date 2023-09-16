#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
    echo "usage: $0 prefix datafile [grpfile] [tstep]"
    exit(1)
endif

set prefix = $1
set mol = `basename $prefix`

#rm -fr ${mol}.lammps ${mol}.eng

if (! -e ${prefix}.lammps) then
    if(! -e ${prefix}.lammpstrj) then
	echo "ERROR: Cannot find ${prefix}.lammps"
	exit(1)
    endif
    ln -s ${prefix}.lammpstrj ${mol}.lammps
else
    ln -s ${prefix}.lammps ${mol}.lammps
endif
if (! -e ${prefix}.eng) then
    set nm = $prefix:r
    if(! -e ${nm}.equil.lammps.log) then
	echo "WARNING: Cannot find ${prefix}.eng or ${nm}.equil.lammps.log"
	#exit(1)
    else
        ln -s ${nm}.equil.lammps.log ${mol}.eng
    endif
else
    ln -s ${prefix}.eng ${mol}.eng 
endif

set datafile = $2
if (! -e $datafile) then
    echo "ERROR: Cannot find data file $datafile"
    exit(1)
endif

set grpString = ""
if ($#argv > 2) then
    if(-e $3) then
	set grpString = "IN_GROUPFILE                  $3"
    endif
endif

set tStep = 0.001
if ($#argv > 3) then
    set tStep = $4
endif

cat > ${mol}.stat.in <<DATA;
IN_LMPDATA                    ${datafile}
IN_LMPTRJ                     ${mol}
ANALYSIS_FRAME_INITIAL        1
ANALYSIS_FRAME_FINAL          0
ANALYSIS_FRAME_STEP           1
ANALYSIS_OUT                  ${mol}.stat
ANALYSIS_LMP_TSTEP            ${tStep}
ANALYSIS_STAT                 2
ANALYSIS_MSD_CORLENGTH        0.5
ANALYSIS_MSD_MEMORYMB         2500
ANALYSIS_MSD_CMCORRECTION     1
ANALYSIS_DIELECTRIC_EWALD     0.0001
${grpString}

DATA
