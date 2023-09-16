#!/bin/tcsh
#!/bin/csh

if ($#argv < 5) then
  echo "usage: $0 cnt_m cnt_n cnt_length reservior1 reservior2 [savename]"
  exit(1)
endif

set scripts_dir = ~tpascal/scripts

foreach i ($4 $5)
	if !(-e $i) then
		echo "ERROR: Cannot access $i"
		exit(1)
	endif
end

echo "Options:"
echo "-----------------------------------------------------"
#set cnt_bgf = $1
echo "CNT:				$1,$2"

set delz = $3
echo "membrane seperation length: 	$delz Angstroms"

set reser1 = $4
set reser2 = $5
echo "Top reserviour :		$reser1"
echo "Bottom reservior:		$reser2"

set savename = `basename $4`
set savename = $savename:r
set savename = "${savename}.${delz}membrane."`basename $4 .bgf`"."`basename $5 .bgf`".bgf"
if ($#argv > 5) set savename = $6

set dims = (`grep "^CRYSTX " $reser1 $reser2 | awk 'BEGIN{x=y=0;}{if($2>x)x=$2; if($3>y)y=$3}END{print x,y}'`)
echo "membrane surace area:		$dims[1] x $dims[2] Angstroms^2"
set z = `grep "^CRYSTX " $reser1 $reser2 | awk -v zlen=$3 'BEGIN{z=0}{z+=$4}END{print z+zlen+20}'`

echo
echo "Step 1. Preparing membrane"
set gra_bgf = __gra.bgf
set graphene_bgf = ${scripts_dir}/dat/graphene_ortho.bgf
set gra_mol = `basename $graphene_bgf .bgf`
set gra_cell = (`grep CRYSTX $graphene_bgf | awk '{print $2,$3}'`)
set rep_str = (`echo $dims $gra_cell | awk '{print int($1/$3)+1,int($2/$4)+1,1}'`)
${scripts_dir}/splitAtomsByMol.pl -b $graphene_bgf -s $gra_bgf -f resnum > /dev/null
${scripts_dir}/groupAtoms.pl -b $gra_bgf -s $gra_bgf -f resnum > /dev/null
${scripts_dir}/replicate.pl -b $gra_bgf -d "$rep_str" -s $gra_bgf  > /dev/null || goto error
set dims = (`grep "^CRYSTX " $gra_bgf | awk '{print $2,$3}'`)

echo "Step 2. Creating CNT"
set cnt_bgf = __${1}x${2}cnt.bgf
set zlen = `echo $3 | awk '{print $1/10}'`
echo csh -f ${scripts_dir}/createCNT.csh $1 $2 $zlen $cnt_bgf
csh -f ${scripts_dir}/createCNT.csh $1 $2 $zlen $cnt_bgf > /dev/null || goto error
set radius = (`grep "^REMARK C =" $cnt_bgf | awk '{print $8}'`)
${scripts_dir}/decreaseDim.pl -b $cnt_bgf -s $cnt_bgf -d z > /dev/null || goto error
${scripts_dir}/updateBGFBox.pl -b $cnt_bgf -c "$dims $z 90 90 90" -s $cnt_bgf > /dev/null || goto error
${scripts_dir}/centerBGF.pl -b $cnt_bgf -c com_center -s $cnt_bgf > /dev/null || goto error
set delz = `${scripts_dir}/getBounds.pl -b $cnt_bgf | grep '^Z ' | awk '{print $3-$2-3}'`

echo "Step 3. Creating toroid membrane"
csh -f ${scripts_dir}/create.toroid.membrane.csh $cnt_bgf $gra_bgf ${scripts_dir}/../ff/graphite.ff _test.bgf >& /dev/null|| goto error
set z_max = `${scripts_dir}/getBounds.pl -b _test.bgf -o "index>0" | grep "^Z" | awk '{print $3}'`
set z_min = `${scripts_dir}/getBounds.pl -b _test.bgf -o "index>0" | grep "^Z" | awk '{print $2}'`

echo "Step 4a. Preparing $reser1"
set solv_z = `${scripts_dir}/getBounds.pl -b $reser1 -o "index>0" | grep '^Z' | awk '{print $2}'`
set solv_mol = `basename $reser1 .bgf`
set dz = `echo $z_max $solv_z | awk '{dz=$1-$2+3; if (dz<0) print dz; else printf "+%s\n",dz;}'`
${scripts_dir}/modifyAtomData.pl -s $reser1 -a "index>0" -f "ZCOORD:${dz}" -w _reser1.ready.bgf > /dev/null || goto error
${scripts_dir}/updateBGFBox.pl -b _reser1.ready.bgf -c "$dims $z 90 90 90" -s _reser1.ready.bgf > /dev/null || goto error
${scripts_dir}/centerBGF.pl -b _reser1.ready.bgf -s _reser1.ready.bgf -d xy -c com_center > /dev/null || goto error

echo "Step 4b. Preparing $reser2"
set solv_z = `${scripts_dir}/getBounds.pl -b $reser2 -o "index>0" | grep '^Z' | awk '{print $3}'`
set solv_mol = `basename $reser1 .bgf`
set dz = `echo $z_min $solv_z | awk '{dz=$1-$2-3; if (dz<0) print dz; else printf "+%s\n",dz;}'`
${scripts_dir}/modifyAtomData.pl -s $reser2 -a "index>0" -f "ZCOORD:${dz}" -w _reser2.ready.bgf > /dev/null || goto error
${scripts_dir}/updateBGFBox.pl -b _reser2.ready.bgf -c "$dims $z 90 90 90" -s _reser2.ready.bgf > /dev/null || goto error
${scripts_dir}/centerBGF.pl -b _reser2.ready.bgf -s _reser2.ready.bgf -d xy -c com_center > /dev/null || goto error

#echo "Step 5: Creating capping sheets"
#set offset = `${scripts_dir}/getBounds.pl -b $reser1 -o "index>0" | grep '^Z' | awk '{print $3}'`
#set dz = `${scripts_dir}/getBounds.pl -b $gra_bgf -o "index>0" | grep '^Z' | awk -v offset=$offset '{dz=$3-offset-1; if(dz<0) print dz; else printf "+%s\n",dz}'`
#${scripts_dir}/modifyAtomData.pl -s $gra_bgf -w _${gra_mol}.cap1.bgf -a "index>0" -f "ZCOORD:$dz" > /dev/null || goto error
#set offset = `${scripts_dir}/getBounds.pl -b $reser2 -o "index>0" | grep '^Z' | awk '{print $3}'`
#set dz = `${scripts_dir}/getBounds.pl -b $gra_bgf -o "index>0" | grep '^Z' | awk -v offset=$offset '{dz=$3+offset+1; if(dz<0) print dz; else printf "+%s\n",dz}'`
#${scripts_dir}/modifyAtomData.pl -s $gra_bgf -w _${gra_mol}.cap2.bgf -a "index>0" -f "ZCOORD:$dz" > /dev/null || goto error
#${scripts_dir}/combineBGF.pl _${gra_mol}.cap1.bgf _${gra_mol}.cap2.bgf _${gra_mol}.caps.bgf > /dev/null || goto error

echo "Step 6. Combining systems"
#${scripts_dir}/combineBGF.pl _test.bgf _${gra_mol}.caps.bgf _test.bgf > /dev/null || goto error
${scripts_dir}/combineBGF.pl _test.bgf _reser1.ready.bgf _test.bgf > /dev/null || goto error
${scripts_dir}/combineBGF.pl _test.bgf _reser2.ready.bgf _test.bgf > /dev/null || goto error
${scripts_dir}/removeBGFBox.pl -b _test.bgf -s _test.bgf > /dev/null || goto error
${scripts_dir}/addBoxToBGF.pl _test.bgf  _test.bgf > /dev/null || goto error
${scripts_dir}/updateBGFBox.pl -b _test.bgf -c "$dims $z 90 90 90" -s ${savename} > /dev/null || goto error
${scripts_dir}/centerBGF.pl -b ${savename} -s ${savename} -c com_center -d z > /dev/null || goto error
${scripts_dir}/updateResNum.pl -b ${savename} -s ${savename} -a "index>0" > /dev/null || goto error

rm -fr $gra_bgf $cnt_bgf _${gra_mol}*.bgf _test.bgf _reser*.bgf

exit:
rm -fr core.*
echo "Created ${savename}"
echo "All Tasks Completed"
exit(0)

error:
echo "ERROR occurred"
exit(1)
