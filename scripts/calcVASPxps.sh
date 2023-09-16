#!/bin/bash

gawk -v delE=0.5 'a
function basename(file) {
	sub(".*/", "", file)
	return file
}
BEGIN{
	i=0;
	fname="";
}
{
	cfname=gensub(/(.+)\.nvt(.+)/,"\\1","g",basename(FILENAME)); 
	if($1 == 1 && $2 ~ /^F=/) {
		if(FILENAME ~ /xps\.OUT/){
			icount=amap[cfname]
			dE=gs[icount]-$3; 
			idx=int(dE/delE); 
			val[icount][idx]++l
			#print "XPS file: ",FILENAME," energy ",$3," icount ",icount," GS file: ",nm[icount]," GS eng: ",gs[icount]," dE ",dE
		} else {
			if(cfname != fname){
				i++; 
				fname=cfname; 
				nm[i]=cfname;
				amap[cfname]=i
			};
			gs[i]=$3;
		}
	}
}
END{
	for (i in val){
		outfile=sprintf("%s.xps.hist.dat",nm[i]); 
		print "" > outfile; 
		for(j in val[i]){
			print j*delE,val[i][j] >> outfile;
		}
	}
}' *.GS.OUT *xps.OUT
