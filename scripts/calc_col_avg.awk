#!/bin/awk -f
BEGIN{
	if(col=="") col=2
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
	if ($1 ~ /^\s*\-?[0-9]+/ && $col ~ /^[0-9]+/) {
		a[FNR] = (a[FNR] ? a[FNR] FS : "" ) ((norm>0) ? ((norm>1) ? $col/mols[fi] : $col*mols[fi]) : $col) 
		b[FNR] = $1+0; 
	}
} 
END { 
	for (i in a) {
		split(a[i],vals,FS); 
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
