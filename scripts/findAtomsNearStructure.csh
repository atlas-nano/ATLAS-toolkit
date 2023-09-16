#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
	echo "usage: $0 bgf_file 'reference_atom(s)' 'search_atom(s)' (distance=2.8)"
	exit(1)
endif

if !(-e $1) then
	echo "ERROR: Cannot locate $1\n"
	exit(1)
endif

set ref_atoms = (`~tpascal/scripts/queryAtom.pl -f "ATMNAME" -b $1 -a "$2" | grep "^ATOM" | awk '{print $2}'`)
if ($#ref_atoms == 0) then
	echo "ERROR: No valid reference atoms found while searching $2"
	exit(1)
endif

set dist = 2.8
if ($#argv > 3) set dist = $4

set searchStr = "($3) and ("
set atm_list = ()
foreach i ($ref_atoms)
    set coords = (`egrep "^(ATOM|HETATM)\s*$i " $1 | awk '{print $7,$8,$9}'`)
    set searchStr = "$searchStr (dist($i)<=$dist) or"
end
set searchStr = `echo $searchStr | sed 's/ or$/\)/'`
set atm_list = ($atm_list `~tpascal/scripts/queryAtom.pl -f "ATMNAME" -b $1 -a "$searchStr" | grep "^ATOM" | awk '{print $2}'`)
echo $atm_list
