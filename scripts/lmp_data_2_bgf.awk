#!/bin/gawk
BEGIN{
	natom=-1
	has_masses=0
	atoms_list=0
	bonds_list=0
	atmprev=0
}
{
	if ($_ ~ /Masses/)
		has_masses=1
	else if (has_masses && $0 ~ /^\s*[0-9]/ && NF == 4 && $3 ~ /#/)
		fftype[$1]=$NF
	else if (has_masses && $0 ~ /(Coeffs|Atoms)/)
		has_masses=0
	else if ($1 ~ /^Atoms/) {
		atoms_list=1
		print "BIOGRF  332\nFORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)";
	} else if (atoms_list && $1 ~ /^Bonds/) {
		atmprev=1
		atoms_list=0
		print "FORMAT CONECT (a6,14i6)"
		bonds_list=1
	} else if (atmprev==0 && atoms_list && NF > 6 && $1 ~ /^[0-9]/ && $2 ~ /^[0-9]/ && $3 ~ /^[0-9]/ && $4 ~ /^\-?[0-9]*\.[0-9]*/) {
		printf "%-6s %5d %5s %3s %1s %5s %9.5f %9.5f %9.5f %5s%3d%2d %8.5f\n","ATOM",$1,fftype[$3],"RES","X",1,$5,$6,$7,fftype[$3],0,0,$4;
		natom=$1
	} else if (bonds_list && NF == 4) {
		nbond[$3]++
		bonds[$3][nbond[$3]]=$4
		nbond[$4]++
		bonds[$4][nbond[$4]]=$3
	}
}END{
	for(i=1;i<=natom;i++) {
		printf "%-6s%6d","CONECT",i;
		for(j=1;j<=nbond[i];j++)
			printf "%6d",bonds[i][j];
		printf "\n"
	}
	print "END"
}
