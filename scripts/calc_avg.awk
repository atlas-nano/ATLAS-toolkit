#!/bin/awk -f
BEGIN{
	if(nmols=="") norm=0
	else {
		if(norm=="") norm=1;
		split(nmols,mols," ");
	}
	if(norm=="") norm = 0;
	fi=0;
}
{ 
	if (FNR==1) fi++;
	if ($1 ~ /^\s*\-?[0-9]+/) {
		for(i=2;i<=NF;i++) {
			a[FNR, i] = (a[FNR, i] ? a[FNR, i] FS : "" ) ((norm>0) ? ((norm>1) ? $col/mols[fi] : $col*mols[fi]) : $col) 
			b[FNR] = $1+0; 
		}
		nframe=NF;
	}
} 
END { 
	for (k in a) {
		split(a[k],coldata,SUBSEP);
		for (i in coldata) {
			split(coldata[i],vals,FS); 
				avg=avg2=c=0; 
			for (j in vals) {
				v = vals[j];
				avg += v; 
				avg2 += v*v; 
				c++
			} 
			avg /= c; 
			avg2 /= c; 
			std=avg2-avg*avg; 
			if(std>0) 
				std=sqrt(std); 
			else 
				std = 0; 
			print b[i],avg,std | "sort -n"; 
		}
	} 
}
