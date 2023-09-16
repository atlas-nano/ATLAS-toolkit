BEGIN {
	if(norm=="")
		norm=1
	if(mode=="")
		mode=1
}
{ 
	for(i=2;i<=NF;i++) {
	       $i /= norm	
		x1[i] += $i; 
		x2[i] += $i*$i; 
	} 
	j = NF;
} 
END { 
	for(i=2;i<=j;i++) { 
		x1[i] /= NR;
		x2[i] /= NR;
		sigma = 0; 
		if((x2[i] - x1[i]*x1[i])>0) { 
			sigma = sqrt(x2[i]-x1[i]*x1[i]); 
		}
		if(mode==2) {
			#standard error
			sigma /= sqrt(NR)
		}
		printf "%s %s ",x1[i],sigma; 
	} 
	printf "\n"; 
}
