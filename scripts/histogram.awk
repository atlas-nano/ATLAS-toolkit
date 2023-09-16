#!/bin/awk -f
BEGIN {
# parameters
# can also be set by "-v bins=40 -v min=10 -v max=30 -v col=2 -v norm=1" command line options
	if(bins == "")
		bins = 40;
	if(min == "" || min == max)
		min = 0;
	if(max == "" || min == max)
		max = 1;
	if(factor == "")
		factor=1;
	if(col == "")
		col=2;
    if(norm == "")
		norm=0;

	less = 0;
	more = 0;
	tot = 1;
	for(i=0; i<bins; i++) {
		a[i] = 0;
	}
} ; 
{
	if($1 ~ /^[-0-9]/) {
		y = int((($col-min)*bins)/(max-min));
		if(y < 0)
			less += 1;
		else if(y >= bins)
			more += 1;
		else  {
			a[y] += 1;
			if(norm) tot+=1;
	}
	}
} ; END {
	bwidth = (max-min)/bins;
	for(i=0; i<bins; i++) {
		print (min + bwidth*i + bwidth/2)*factor, a[i]/tot;
	}
}
