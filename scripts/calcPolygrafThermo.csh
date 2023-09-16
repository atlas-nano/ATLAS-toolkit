#!/bin/tcsh
#!/bin/csh

if ($#argv < 1) then
  echo "usage: $0 bgf_file (forcefield=dreiding) (rotational symmetry number = 1)"
  exit(1)
endif

set bgf_file = $1
if !(-e $bgf_file) then
  echo "ERROR: Cannot locate bgf file $bgf_file"
  exit(1)
endif

set forcefield = "dreiding"
if ($#argv > 1) then
  set forcefield = $2
  if (-e $forcefield) then
    cp $forcefield ./
    set ext = $forcefield:e
    set forcefield = `basename $forcefield .${ext}`
  endif
endif

set rot_sym = 1
if ($#argv > 2) then
  set rot_sym = $3
endif

set mol = `basename $bgf_file .bgf`
/ul/tpascal/scripts/createPolygrafThermoMacro.pl -b $bgf_file -r $rot_sym -s ${mol}.vacthermo.macro >& /dev/null || goto error
setenv BIOVER 400 
source ${membstruk}/Biograf/bg_login.csh
/biograf/bio_msc batpoly $forcefield ${mol}.vacthermo.macro >& ${mol}.vacthermo.log || goto error

exit:
echo "All tasks completed"
exit(0)

error:
echo "Error occurred. Check ${mol}.vacthermo.log"
exit(1)
