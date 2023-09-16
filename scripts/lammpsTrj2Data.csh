#!/bin/tcsh
#!/bin/csh

if ($#argv < 5) then
  echo "usage: $0 bgf_file force_field lammps_trj_file lammps_vel_file number_frames [prefix]"
  exit(1)
endif

set bgf_file = $1
if !(-e $bgf_file) then 
  echo "ERROR: Cannot locate BGF file $bgf_file"
  exit(1)
endif

set force_field = $2

set lmp_trj = $3
if !(-e $lmp_trj) then
  echo "ERROR: Cannot locate lammps trj file $lmp_trj"
  exit(1)
endif

set lmp_vel = $4
if !(-e $lmp_vel) then
  echo "ERROR: Cannot locate lammps vel file $lmp_vel"
  exit(1)
endif

set num = $5
set mol = `basename $bgf_file .bgf`
if ($#argv > 5) then
  set mol = $6
endif

set tot = `grep -c TIMESTEP $lmp_trj`
if ($num > $tot) then
  set num = $tot
endif
set tot1 = `grep -c TIMESTEP $lmp_vel`
if($tot != $tot1) then
  echo "ERROR: Number of frames in LAMMPS trajectory ($tot) not equal to LAMMPS velocity file ($tot1)"
  exit(1)
endif

set inc = `echo $tot $num | awk '{printf "%.0f",($1/$2)}'`
echo "selected frame 1 to $tot every $inc"
echo "velocities..."
/ul/tpascal/scripts/getLammpsVel.pl -v $lmp_vel -f ":It1-${tot}:${inc}" -b $bgf_file -t data -s ${mol} > /dev/null || goto error
echo "coordinates..."
/ul/tpascal/scripts/convertLammpsTrj.pl -t ":It1-${tot}:${inc}" -b $bgf_file -o bgf -s ${mol} -l $lmp_trj > /dev/null || goto error
foreach b (${mol}.*.bgf)
  set sname = $b:r
  echo "converting $sname to lammps data"
  /ul/tpascal/scripts/createLammpsInput.pl -b $b -f "$force_field" -s $sname -l 0 > /dev/null || goto error
  cat data.${sname}.vel >> data.${sname}
  rm -fr in.${sname}* ${sname}* data.${sname}.vel
end

exit:
echo "All tasks completed"
exit(0)

error:
echo "Error occurred"
exit(1)
