#!/usr/bin/gawk -f
BEGIN{
	if(ARGC < 3) {
		print "usage: $0 -v chrg=[chelpg(default)|esp|resp|mulliken] qchem_out bgf_file"
		exit(1)
	}
	if(chrg == "") {
		chrg = "chelpg"
	}
	chrg=tolower(chrg)

	if(chrg ~ /chelp/) { 
		chrgStr = "Ground-State ChElPG Net Atomic Charges"
	} else if (chrg ~ /resp/) {
		chrgStr = "Merz-Kollman RESP Net Atomic Charges"
	} else if (chrg ~ /esp/) {
		chrgStr = "Merz-Kollman ESP Net Atomic Charges"
	} else if (chrg ~ /mul/) {
		chrgStr = "Ground-State Mulliken Net Atomic Charges"
	} else {
		print "ERROR: Cannot figure out charge type from '" chrg "'"
		exit(1)
	}
	c=0
	chrg_start = chrg_valid = rem = 0
	sprintf("egrep -c '^(ATOM|HETATM)' %s",ARGV[2]) | getline nAtoms
}
{
	if(FNR==1)c++
	if(c==1 && chrg_start == 0 && $0 ~ chrgStr) {
		chrg_start = 1
		chrg_valid = 1
	} else if (c == 1 && chrg_start > 0 && $1 ~ /^[0-9]/) {
		aChrg[$1]=$3
		nChrg = $1
	} else if (c == 1 && chrg_start > 0 && $0 ~ /Sum of atomic charges/) {
		if (nChrg != nAtoms) {
			print "ERROR: Number of values from " chrg " file " ARGV[1] " (" nChrg ") is not equal to number of atoms (" nAtoms ") in " ARGV[2] " ... Cannot continue"
			exit(1)
		}
		chrg_start = 0
	} else if (c == 2) {
		if(! chrg_valid) {
			print "ERROR: No valid charges read from " ARGV[1]
			exit(1)
		}
		if($1 !~ /^[ATOM|HETATM]/) {
			if($1 ~ /FORCEFIELD/ && rem == 0) {
				print "REMARK updated charges to " chrgStr
				rem = 1
			}
			print
		} else {
			n=split($0,a," ",b)
			tmpStr=sprintf("%10.5f",aChrg[$2])
			split(tmpStr,tmp1," ",tmp2)
			b[12]=tmp2[0]
			a[13]=tmp1[1]
			line=b[0]
			for (i=1;i<=n; i++)
				line=(line a[i] b[i])
			print line
		}
	}
}
