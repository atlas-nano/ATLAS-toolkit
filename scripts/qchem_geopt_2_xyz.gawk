#!/bin/gawk
function basename(file) {
     sub(".*/", "", file)
     return file
}
function prefix(file) {
	n=split(file,a,".")
	ofile = a[1]
	for(i=2;i<n;i++)
		ofile=sprintf("%s.%s",ofile,a[i])
	return ofile
}
BEGIN{
	if(length(xyz)>0) {
		system("rm -fr " xyz)
		system("touch " xyz)
	}
}
{
	if(FNR==1) {
		c++
		read_xyz=0
		read_natom=0
		natom=-1
		if(length(xyz) == 0) {
			out_xyz = sprintf("%s.pes_scan.xyz",prefix(basename(FILENAME)))	
			print "writing to " out_xyz	
		}else {
			out_xyz = xyz
		}
	}
	if(read_natom==0 && $0 ~ /NAtoms,    NIC,     NZ,/) {
		read_natom=1
	} else if(read_natom == 1 && $1 ~ /^[0-9]/) {
		natom=$1
		read_natom=-1
	} else if($0 ~ /OPTIMIZATION CONVERGED/) {
		read_xyz=1
		print natom >> out_xyz
		print "step ",c >> out_xyz
		c++
	} else if(read_xyz == 1 && $0 ~ /Z-matrix Print/) {
		read_xyz = 0
	} else if(read_xyz == 1 && $1 ~ /^[0-9]/ && NF==5) {
		$1 = "";
		print >> out_xyz
	}
}
