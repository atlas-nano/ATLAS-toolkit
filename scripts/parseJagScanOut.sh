#!/bin/bash

if [ $# -eq 0 ]; then
	echo "usage: $0 jag_outfile [savePrefix]"
	exit 1
fi

if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot access $1"
	exit 1
fi

sprefix=`basename $1`
sprefix=${sprefix%.*}
if [ $# -gt 1 ]; then
	sprefix=$2
fi

jag_out=$1
if [ `egrep -c ' scan: ' $jag_out` -eq 0 ] || [ `egrep -c ' geometry:$' $jag_out` -eq 0 ]; then
	echo "ERROR: Invalid Jaguar scan file $1"
	exit 1
fi

echo "will create ${sprefix}.scan.dat and xyz files ${sprefix}"
rm -fr ${sprefix}.scan.dat
gawk -v prefix=$sprefix '
BEGIN{
	save_coords=n=0;
}
{
	if($2 ~ /^geometry:/){
		save_coords=1;
		c=0;
	}else if(NF==0 && save_coords==1){
		save_coords=0;
	}else if(save_coords==1 && NF==4){
		c++;
		if(c>1) {
			gsub(/[^a-zA-Z]+/,"",$1);
			coords[c][1]=$1;
			coords[c][2]=$2;
			coords[c][3]=$3;
			coords[c][4]=$4;
		}
	}else if($1 ~ /^scan:/) {
		n++;
		save_coords=0;
		indx=$4;
		sval=$7*627.509; #kcal/mol
		scanfile=sprintf("%s.scan.dat",prefix);
		xyzfile=sprintf("%s.%d.xyz",prefix,n);
		printf "%d %f # %f\n",n,sval,indx >> scanfile
		printf "%-d\nenergy %f index %f\n",c-1,sval,indx > xyzfile
		for(i=2;i<=c;i++) {
			printf " %-10s %10.5f %10.5f %10.5f\n",coords[i][1],coords[i][2],coords[i][3],coords[i][4] >> xyzfile
		}
	}
}' $jag_out
