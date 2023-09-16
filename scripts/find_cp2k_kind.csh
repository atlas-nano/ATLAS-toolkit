#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 xyz_file [save_prefix]"
	exit(1)
endif

if !(-e $1) then
	echo "ERROR: Cannot locate $1"
	exit(1)
endif

set xyz_file = $1
set prefix = `basename $xyz_file`
set prefix = $prefix:r
set prefix = ${prefix}.kind
if ($#argv > 1) set prefix = $2

set basis_set_file = "~tpascal/codes/cp2k/trunk/cp2k/tests/QS/GTH_BASIS_SETS ~tpascal/codes/cp2k/trunk/cp2k/tests/QS/BASIS_MOLOPT"
set pseudo_file = "~tpascal/codes/cp2k/trunk/cp2k/tests/QS/GTH_POTENTIALS ~tpascal/codes/cp2k/trunk/cp2k/tests/QS/POTENTIAL"
set clist = ()
foreach i (`cat ${xyz_file} | awk '{val[$1]=1}END{for (ii in val) print ii}'`)
	set list = (`grep -c ${i} $basis_set_file | grep GTH | sed 's/.*://'`)
	if ($#list == 0) then
		echo "ERROR: Cannot find DZVP-GTH entry for $i in $basis_set_file"
		exit(1)
	endif
	set list = (`grep -c ${i} $pseudo_file | grep GTH | sed 's/.*://'`)
	if ($#list == 0) then
		echo "ERROR: Cannot finf GTH-PBE entry for $i in $pseudo_file"
		exit(1)
	endif
	set pseudo_str = `grep "${i} GTH-PBE-" $pseudo_file | head -1 | awk '{print $2}'`
	cat > ${prefix}.${i}.dat <<DATA;

  &KIND $i
   BASIS_SET DZVP-GTH
   POTENTIAL $pseudo_str
  &END KIND
DATA

	set clist = ($clist ${prefix}.${i}.dat)
end
cat $clist > ${prefix}.tot.dat
