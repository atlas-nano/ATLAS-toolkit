#!bin/bash

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e

if [ $# -lt 5 ]; then
	echo "usage: $0 graphene_bgf x y z slit_size [savename]"
	exit 1
fi

if [ ! -e $1  ] ||  [ ! -r $1 ]; then
	echo "ERROR: Cannot access $1"
	exit 1
fi

scripts_dir=/home/tpascal/scripts/

bgf=$1
x=$2
y=$3
z=$4
slit_size=$5
savename=`basename $bgf`
savename=${savename%.*}
savename="${savename}.nanoslit.${x}x${y}Sheet.${z}Wall.${slit_size}Slit.bgf"
if [ $# -gt 5 ]; then
	savename=$6
fi

cat <<DATA;
==============
OPTIONS
==============
Graphene file: $bgf
Graphene Sheet X length: $x
Graphene Sheet Y length: $y
Graphene Wall Z length: $z
Slit Size: $slit_size
DATA

wat_bgf=${scripts_dir}/dat/WAT/tip3_box.bgf
wat_z=`grep CRYSTX $wat_bgf | awk '{print $4}'`

bgf_crystx=(`grep CRYSTX $bgf | awk '{print $2,$3,$4}'`)

#create sheet
sheet_rep=(`echo $x $y ${bgf_crystx[*]} | awk -v z=$4 '{print int(($1*2+z*2)/$3)+1,int($2/$4)+1}'`)
if [ ${sheet_rep[0]} -gt 1 ] || [ ${sheet_rep[1]} -gt 1 ]; then
	${scripts_dir}/replicate.pl -b $bgf -s __sheet.bgf -d "${sheet_rep[*]} 1" > /dev/null 
else
	cp $bgf __sheet.bgf
fi

${scripts_dir}/groupAtoms.pl -b __sheet.bgf -s __sheet.bgf -f "resnum" -g xcoord > /dev/null
#determine cuts
cuts=(`grep CRYSTX __sheet.bgf | awk -v z=$z -v x=$x '{tot=x*2+z*2; cx=$2; fac=cx/tot; printf "%.2f %.2f %.2f %.2f\n",fac*x/2,fac*(x/2+z),fac*(3*x/2+z),fac*(3*x/2+2*z)}'` 100000000)
echo "cuts: ${cuts[*]}"
#get pieces
cp __sheet.bgf __tmp.bgf
rot=180
for i in `seq 1 ${#cuts[*]}`
do
	v=${cuts[i-1]}
	echo cut $i $v
	${scripts_dir}/getBGFAtoms.pl -b __tmp.bgf -s __cut${i}.bgf -o "xcoord<$v" > /dev/null
	#rotate pieces
	if [ $rot -ne 0 ]; then
		echo "	rotating by $rot degrees"
		${scripts_dir}/rotatemol.pl -b __cut${i}.bgf -r "y:${rot}" -s __cut${i}.bgf > /dev/null
		#${scripts_dir}/centerBGF.pl -b __cut${i}.bgf -s __cut${i}.bgf -c box_origin > /dev/null
	fi
	rot=$(( $rot - 90 ))
	if [ $i -lt ${#cuts[*]} ]; then
		${scripts_dir}/getBGFAtoms.pl -b __tmp.bgf -s __tmp.bgf -o "xcoord>$v" > /dev/null
	fi
done
#position
##get the x vals for cut 3 (top center)
xvals=(`${scripts_dir}/getBounds.pl -b __cut3.bgf | grep '^X' | awk '{print $2,$3}'`)
zvals=(`${scripts_dir}/getBounds.pl -b __cut3.bgf | grep '^Z' | awk '{print $2,$3}'`)

##cut 1
dx=`${scripts_dir}/getBounds.pl -b __cut1.bgf | grep '^X' | awk -v lo=${xvals[0]} '{d=lo-$2; if(d>0) { printf "+%s",d} else {print d}}'`
dz=`${scripts_dir}/getBounds.pl -b __cut1.bgf | grep '^Z' | awk -v lo=${zvals[0]} -v off=$z '{d=lo-off-$2-2; if(d>0) { printf "+%s",d} else {print d}}'`
echo "shifting cut1 by xcoord $dx zcoord $dz"
${scripts_dir}/modifyAtomData.pl -s __cut1.bgf -w __cut1.bgf -a "index>0" -f "XCOORD:$dx ZCOORD:$dz" > /dev/null

##cut 2
dx=`${scripts_dir}/getBounds.pl -b __cut2.bgf | grep '^X' | awk -v lo=${xvals[0]} '{d=lo-$2; if(d>0) { printf "+%s",d} else {print d}}'`
dz=`${scripts_dir}/getBounds.pl -b __cut2.bgf | grep '^Z' | awk -v lo=${zvals[0]} -v off=$z '{d=$2-lo-1; if(d>0) { printf "+%s",d} else {print d}}'`
echo "shifting cut2 by xcoord $dx zcoord $dz"
${scripts_dir}/modifyAtomData.pl -s __cut2.bgf -w __cut2.bgf -a "index>0" -f "XCOORD:$dx ZCOORD:$dz" > /dev/null

##cut 4
dx=`${scripts_dir}/getBounds.pl -b __cut4.bgf | grep '^X' | awk -v hi=${xvals[1]} '{d=hi-$2; if(d>0) { printf "+%s",d} else {print d}}'`
dz=`${scripts_dir}/getBounds.pl -b __cut4.bgf | grep '^Z' | awk -v lo=${zvals[0]} -v off=$z '{d=$2-lo-1; if(d>0) { printf "+%s",d} else {print d}}'`
echo "shifting cut4 by xcoord $dx zcoord $dz"
${scripts_dir}/modifyAtomData.pl -s __cut4.bgf -w __cut4.bgf -a "index>0" -f "XCOORD:$dx ZCOORD:$dz" > /dev/null

##cut 5
dx=`${scripts_dir}/getBounds.pl -b __cut5.bgf | grep '^X' | awk -v hi=${xvals[1]} '{d=hi-$3; if(d>0) { printf "+%s",d} else {print d}}'`
dz=`${scripts_dir}/getBounds.pl -b __cut5.bgf | grep '^Z' | awk -v hi=${zvals[1]} -v off=$z '{d=hi-off-$3-2; if(d>0) { printf "+%s",d} else {print d}}'`
echo "shifting cut5 by xcoord $dx zcoord $dz"
${scripts_dir}/modifyAtomData.pl -s __cut5.bgf -w __cut5.bgf -a "index>0" -f "XCOORD:$dx ZCOORD:$dz" > /dev/null

#now reassemble
${scripts_dir}/combineBGF.pl __cut1.bgf __cut2.bgf __temp.bgf > /dev/null
mv __temp.bgf __test.bgf
${scripts_dir}/combineBGF.pl __test.bgf __cut3.bgf __temp.bgf > /dev/null
mv __temp.bgf __test.bgf
${scripts_dir}/combineBGF.pl __test.bgf __cut4.bgf __temp.bgf > /dev/null
mv __temp.bgf __test.bgf
${scripts_dir}/combineBGF.pl __test.bgf __cut5.bgf __temp.bgf > /dev/null
mv __temp.bgf __test.bgf

#reconnect bonds
sed '1,/FORMAT CONECT/d' __sheet.bgf > __bonds.dat
sed -i '/FORMAT CONECT/,$d' __test.bgf
echo "FORMAT CONECT (a6,12i6)" >> __test.bgf
cat __bonds.dat >> __test.bgf
${scripts_dir}/decreaseDim.pl -b __test.bgf -s __test.bgf -d x > /dev/null

x=`echo $x $slit_size | awk '{print $1+$2}'`
y=`grep CRYSTX __sheet.bgf | awk '{print $3}'`
dz=`echo $z | awk '{print $1*3+10}'`
${scripts_dir}/updateBGFBox.pl -b __test.bgf -s __test.bgf -c "$x $y $dz" > /dev/null 
${scripts_dir}/centerBGF.pl -b __test.bgf -s __test.bgf -c com_center  > /dev/null 

#make water box
wat_rep=(`grep CRYSTX $wat_bgf | awk -v x=$x -v y=$y -v z=$z '{print int(x/$2)+1,int(y/$3)+1,int(z/$4)+1}'`)
if [ ${wat_rep[0]} -gt 1 ] || [ ${wat_rep[1]} -gt 1 ]; then
	${scripts_dir}/replicate.pl -b $wat_bgf -s __watbox.bgf -d "${wat_rep[*]}" > /dev/null 
else
	cp $wat_bgf __watbox.bgf
fi
wat_z=$z
${scripts_dir}/trimCell.pl -b __watbox.bgf -s __watbox.bgf -c "$x $y $wat_z" -m 1 -u 1 > /dev/null 
tz=`echo $z $wat_z 6.0 | awk '{print $1+$2+$3+4}'`
${scripts_dir}/transmol.pl __watbox.bgf $tz z __watbox.bgf > /dev/null 
${scripts_dir}/updateBGFBox.pl -b __watbox.bgf -s __watbox.bgf -c "$x $y $dz" > /dev/null 
${scripts_dir}/centerBGF.pl -b __watbox.bgf -s __watbox.bgf -c com_center  > /dev/null 

#create final system
${scripts_dir}/combineBGF.pl __test.bgf __watbox.bgf $savename > /dev/null 
${scripts_dir}/splitAtomsByMol.pl -b $savename -s $savename -f resnum > /dev/null 

trap : 0
rm -fr __*.bgf __bonds.dat
echo "Created $savename"
exit 0
