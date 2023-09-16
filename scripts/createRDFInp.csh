#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
    echo "usage: $0 prefix datafile [use mol com] [grpfile] [grp ids]"
    exit(1)
endif

set prefix = $1
if (! -e ${prefix}.lammps) then
    echo "ERROR: Cannot find ${prefix}.lammps"
    exit(1)
endif
set mol = `basename $prefix`

set datafile = $2
if (! -e $2) then
    echo "ERROR: Cannot find data file $2"
    exit(1)
endif

set molStr = ""
if ($#argv > 2 && $3 == 1) then
    set molStr = "ANALYSIS_RDF_MOL_COM"
endif

set grpStr = ""
if ($#argv > 3) then
	if (-e $4) then
	    set grpStr = "IN_GROUPFILE                  ${4}"
	endif
endif

set grpidStr = "ANALYSIS_RDF_GRPIDS          0 0"
if ($#argv > 4) then
    set grpidStr = "ANALYSIS_RDF_GRPIDS          $5"
endif
cat > ${mol}_rdf.in <<DATA;
IN_LMPDATA                    ${datafile}
IN_LMPTRJ                     ${prefix}
ANALYSIS_OUT                  ${mol}_rdf
ANALYSIS_FRAME_INITIAL        1
ANALYSIS_FRAME_FINAL          0
ANALYSIS_FRAME_STEP           1
ANALYSIS_RDF_RCUT            10
ANALYSIS_RDF_RINCRMNT        0.1
${grpStr}
${grpidStr}
${molStr}

DATA
