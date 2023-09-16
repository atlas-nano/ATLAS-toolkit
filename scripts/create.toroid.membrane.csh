#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 cnt_file membrane_file '[ff_file(s)]' [savename]"
	exit(1)
endif

set cnt_file = $1
set mem_file = $2
if !(-e $cnt_file && -r $cnt_file) then
	echo "ERORR: Cannot access $cnt_file"
	exit(1)
endif
if !(-e $mem_file && -r $mem_file) then
	echo "ERORR: Cannot access $mem_file"
	exit(1)
endif

set savename = `basename $cnt_file`
set savename = $savename:r
set savename = "${savename}.mem.bgf"
if ($#argv > 3) set savename = $4

set bl = 1.4
echo "Step 1. Updating system dimensions..."
set new_cell = (`grep CRYSTX $mem_file | awk '{print $2,$3}'` `grep CRYSTX $cnt_file | awk '{print $4+20,$5,$6,$7}'`)
~tpascal/scripts/updateBGFBox.pl -b $cnt_file -s __cnt.bgf -c "$new_cell" > /dev/null || goto error
~tpascal/scripts/updateBGFBox.pl -b $mem_file -s __mem.bgf -c "$new_cell" > /dev/null || goto error
~tpascal/scripts/centerBGF.pl -b __cnt.bgf -s __cnt.bgf -c com_center -d xy > /dev/null || goto error
set cnt_zbounds = (`~tpascal/scripts/getBounds.pl -b __cnt.bgf -o "index>0" | grep '^Z ' | awk '{print $2,$3}'`)
set mem_zbounds = (`~tpascal/scripts/getBounds.pl -b __mem.bgf -o "index>0" | grep '^Z ' | awk '{print $2,$3}'`)
set zoffset = `echo $cnt_zbounds[1] $mem_zbounds[1] | awk -v bl=$bl '{val = $2-$1+bl; if(val>0) printf "+%s\n",val; else print val;}'`
echo "	cnt: $cnt_zbounds mem: $mem_zbounds zoffset: ${zoffset}"
~tpascal/scripts/modifyAtomData.pl -s __cnt.bgf -a "index>0" -f "ZCOORD:${zoffset}" -w __cnt.bgf > /dev/null || goto error
set radius = `~tpascal/scripts/getBounds.pl -b __cnt.bgf -o "index>0" | egrep '(X|Y)' | awk 'BEGIN{val=0; id=0}{val += $3-$2; id++}END{print val/id/2}'`
echo "	cnt radius: $radius"
set num_unsaturated_c = `~tpascal/scripts/queryAtom.pl -b __mem.bgf -a "numbonds==2" -f index | grep -c '^ATOM'` 
echo "	num_unsaturated_c: $num_unsaturated_c"
set mem_xycenter = (`echo $new_cell | awk '{print $1/2,$2/2}'`)
echo "	mem_xycenter: $mem_xycenter"

echo "Step 2. Preparing BGF files..."
set zcenter = `~tpascal/scripts/getBounds.pl -b __cnt.bgf -o "index>0" | grep '^Z ' | awk '{print ($3-$2)/2+$2}'`
set num_atms = `~tpascal/scripts/queryAtom.pl -b __cnt.bgf -a "zcoord<$zcenter and numbonds==2" -f index | grep -c '^ATOM'`
~tpascal/scripts/modifyAtomData.pl -s __cnt.bgf -w __cnt.bgf -a "index>0" -f "CHAIN:X" > /dev/null || goto error
~tpascal/scripts/modifyAtomData.pl -s __cnt.bgf -w __cnt.bgf -a "zcoord<$zcenter and numbonds==2" -f "CHAIN:A" > /dev/null || goto error
~tpascal/scripts/modifyAtomData.pl -s __cnt.bgf -w __cnt.bgf -a "zcoord>$zcenter and numbonds==2" -f "CHAIN:C" > /dev/null || goto error
echo "	zcenter: $zcenter"
echo "	num_atms: $num_atms"
@ valid = 0
set inc = `echo 0.355`

echo "Step 3: Aligning CNT and sheet"
foreach k (`seq 0 5`)
	set nrad = `echo $radius $k $inc | awk '{print $1+$2*$3}'`
	foreach i (0 -${inc} ${inc})
		set x = `echo $i $mem_xycenter[1] | awk '{print $1+$2}'`
		foreach j (0 -${inc} ${inc})
			set y = `echo $j $mem_xycenter[2] | awk '{print $1+$2}'`
			~tpascal/scripts/getBGFAtoms.pl -b __mem.bgf -o "sqrt((xcoord-$x)**2+(ycoord-$y)**2)>$nrad" -s __mem.hole.bgf > /dev/null || goto error
			~tpascal/scripts/modifyAtomData.pl -s __mem.hole.bgf -w __mem.hole.bgf -a "index>0" -f "CHAIN:X" > /dev/null || goto error
			set offset = `grep -c '^CONECT ' __mem.hole.bgf`
			set unsaturated_c_mem = `~tpascal/scripts/queryAtom.pl -b __mem.hole.bgf -a "numbonds<3" -f index | grep -c '^ATOM'`
			~tpascal/scripts/modifyAtomData.pl -s __mem.hole.bgf -w __mem.hole.bgf -a "numbonds<3" -f "CHAIN:B" > /dev/null || goto error
			@ diff = $unsaturated_c_mem - $num_atms - $num_unsaturated_c
			if ($diff == 0) then
				@ valid = 1
				break
			endif
		end
		if ($diff == 0) break
	end
	if ($diff == 0) break
end
if !($valid) then
	echo "ERROR: Cannot align sheet and cnt. Exiting"
	exit(1)
endif
set x = `echo $x $mem_xycenter[1] | awk '{val=$1-$2; if(val>=0) printf "+%s\n",val; else print val}'`
set y = `echo $y $mem_xycenter[2] | awk '{val=$1-$2; if(val>=0) printf "+%s\n",val; else print val}'`
~tpascal/scripts/modifyAtomData.pl -s __cnt.bgf -a "index>0" -f "XCOORD:${x} YCOORD:${y}" -w __cnt.bgf > /dev/null || goto error

#bottom
echo "Step 4. Bottom Sheet"
~tpascal/scripts/combineBGF.pl __mem.hole.bgf __cnt.bgf __mem.hole.cnt.bgf > /dev/null || goto error
~tpascal/scripts/createBestBonds.pl -b __mem.hole.cnt.bgf -i "chain eq 'A'" -j "chain eq 'B'" -s __mem.hole.cnt.bond.bgf > /dev/null || goto error
#top
echo "Step 5. Top Sheet"
set zoffset = `~tpascal/scripts/getBounds.pl -b __mem.hole.bgf -o "index>0" | grep '^Z ' | awk '{print $3}'`
set zoffset = `~tpascal/scripts/getBounds.pl -b __mem.hole.cnt.bgf -o "index>0" | grep '^Z ' | awk -v offset=$zoffset -v bl=$bl '{print $3-offset+bl}'`
~tpascal/scripts/modifyAtomData.pl -s __mem.hole.bgf -a "index>0" -f "ZCOORD:+${zoffset}" -w __mem.top.bgf > /dev/null || goto error
~tpascal/scripts/modifyAtomData.pl -s __mem.top.bgf -a "chain eq 'B'" -f "CHAIN:D" -w __mem.top.bgf > /dev/null || goto error
~tpascal/scripts/combineBGF.pl __mem.hole.cnt.bond.bgf __mem.top.bgf __2mem.hole.cnt.bgf > /dev/null || goto error
~tpascal/scripts/createBestBonds.pl -b __2mem.hole.cnt.bgf -i "chain eq 'C'" -j "chain eq 'D'" -s __2mem.hole.cnt.bonds.bgf > /dev/null || goto error
#~tpascal/scripts/bondByDistance.pl -b __2mem.hole.cnt.bgf -c 1.5 -m 3 -s __2mem.hole.cnt.bgf > /dev/null || goto error

echo "Step 6. Creating $savename"
~tpascal/scripts/centerBGF.pl -b __2mem.hole.cnt.bonds.bgf -s $savename -c com_center > /dev/null || goto error
~tpascal/scripts/modifyAtomData.pl -s $savename -w $savename -a "index>0" -f "RESNUM:1" > /dev/null || goto error
~tpascal/scripts/modifyAtomData.pl -s $savename -w $savename -a "chain ne 'X'" -f "RESNUM:2" > /dev/null || goto error
rm -fr __mem.*bgf __cnt.bgf __2mem.*
if ($#argv == 2) exit(0)

echo "Step 7. Minimizing structure"
set lammps = ~tpascal/codes/bin/lmp_mpi
~tpascal/scripts/createLammpsInput.pl -b ${savename} -f "$3" -s test -t min -o "no cross" > /dev/null || goto error
$lammps -in in.test -screen none || goto error
~tpascal/scripts/convertLammpsTrj.pl -b $savename -l test.min.lammps -o bgf -t last -s $savename > /dev/null || goto error
~tpascal/scripts/modifyAtomData.pl -s $savename -w $savename -a "index>0" -f "RESNUM:1 RESNAME:GRA CHAIN:A" > /dev/null || goto error
rm -fr in.test in.test_singlepoint data.test test.lammps.pbs log.lammps log.cite test.min.lammps
exit(0)

error:
echo "ERROR occurred"
exit(1);
