#!/bin/bash

if [ $# -lt 2 ]; then
	echo "usage: $0 bgf_file jag_output_file [save_bgf_name]"
	exit 1
fi

bgf=$1
jag=$2
scripts_dir=`dirname $0`

for i in $bgf $jag
do
	if [ ! -e $i ] || [ ! -r $i ]; then
		echo "ERROR: Cannot access $i"
		exit 1
	fi
done

savename=`basename $1 .bgf`".nboQ.bgf"
if [ $# -gt 2 ]; then
	savename=$3
fi

echo "savename: $savename"
#first, search jaguar output file for NBO charges
charges=(`awk '
BEGIN{
	valid=0;
	i=1;
}
{
	if($1 ~ /Atom/ && $2 ~ /No/ && $3 ~ /Charge/ && $4 ~ /Core/ && $5 ~ /Valence/ && $6 ~ /Rydberg/ && $7 ~ /Total/) {
		valid = 1;
	} else if (valid && $1 ~ /^==/) {
		valid = 0;
	} else if (valid && $1 !~ /^-/ && NF> 6) {
		val[i++]=$3
	}
}
END{
	for(j=1;j<i;j++) {
		print val[j]
	}
}' $jag`)
if [ ${#charges[*]} -eq 0 ]; then
	echo "ERROR: No valid NBO charges found in Jaguar NBO file '$nbo'. Aborting..."
	exit 1
fi
#now test whether we have the correct number of charges and atoms
natoms=`egrep -c '^(ATOM|HETATM) ' $bgf`
if [ $natoms -ne ${#charges[*]} ]; then
	echo "ERROR: Incorrect number of NBO charges (${#charges[*]}) in $jag compared to number of atoms ($natoms) in $bgf"
	exit 1
fi

gawk -v q_str="${charges[*]}" '
BEGIN{
	split(q_str,q," ");
	c=1;
}
{
	if($1 ~ /^HETATM/ || $1 ~ /^ATOM/) {
		n=split($0,a," ",b)
		a[NF]=sprintf("%8.5f",q[c++])
		line=b[0]
		for (i=1;i<=n; i++)
			line=(line a[i] b[i])
		print line
	} else {
		print
	}
}' $bgf > $savename
~/scripts/getBGFAtoms.pl -b $savename -s $savename -o "index>0" > /dev/null
