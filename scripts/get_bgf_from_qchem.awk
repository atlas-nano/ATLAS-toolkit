#!/bin/gawk

BEGIN{
	if(prefix == "") {
		print "usage: awk -v prefix={prefix} $0 {qcout(s)}"
		exit 1
	}
	mode=-1;
	fc=0
	eng_file = sprintf("%s.eng.dat",prefix);
	print "#distance(A) energy(hartree)" > eng_file;
}
{
	if(FNR==1){
		if(fc>0) {
			xyz_data=sprintf("%d\ndistance %f energy %f\n%s",ac,indx,eng,xyz_data);
			print xyz_data>xyz_out;
			print "babel -ixyz " xyz_out " -obgf " bgf_out " > /dev/null 2>&1" | "/bin/sh"
			print "/home/tpascal/scripts/addBoxToBGF.pl " bgf_out " " bgf_bout " > /dev/null 2>&1" | "/bin/sh"
			print "/home/tpascal/scripts/updateBGFBox.pl -b " bgf_out " -s " bgf_out " -c '100 100 100 90 90 90' > /dev/null 2>&1" | "/bin/sh"
			close("/bin/sh")
			printf("%s %f\n",indx,eng) >> eng_file
		}
		xyz_data=""; 
		sub(/\/.*/,"",FILENAME);
		sub(/R_/,"",FILENAME); 
		indx=FILENAME; 
		xyz_out=sprintf("%s.dimer.%sA.xyz",prefix,indx); 
		bgf_out=sprintf("%s.dimer.%sA.bgf",prefix,indx); 
		fc++;
		ac=0;
	} else if($1 ~ /^\$molecule/)
		mode=1;
	else if (mode==2 && $1 ~ /^\$end/)
		mode = -1; 
	else if (mode==1) 
		mode++; 
	else if (mode==2)  {
		ac++;
		xyz_data = sprintf("%s%s\n",xyz_data,$0);
	} else if ($1 ~ /MP2/ && $2 ~ /total/ && $3 ~ /energy/) {
		eng = $5
	}
}
END{
	xyz_data=sprintf("%d\ndistance %f energy %f\n%s",ac,indx,eng,xyz_data);
	print xyz_data>xyz_out;
	print "babel -ixyz " xyz_out " -obgf " bgf_out " > /dev/null 2>&1" | "/bin/sh"
	print "/home/tpascal/scripts/addBoxToBGF.pl " bgf_out " " bgf_bout " > /dev/null 2>&1" | "/bin/sh"
	print "/home/tpascal/scripts/updateBGFBox.pl -b " bgf_out " -s " bgf_out " -c '100 100 100 90 90 90' > /dev/null 2>&1" | "/bin/sh"
	close("/bin/sh")
	printf("%s %f\n",indx,eng) >> eng_file
}
