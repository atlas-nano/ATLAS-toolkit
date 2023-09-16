#!/bin/gawk
BEGIN{
	if(length(c)==0)
		c=1
}
{
	if(length(p) == 0) {
		s = $c
		p = $c
	} else {
		d=$c-p
		if(d>1) {
			if((p-s)>0) 
				out_str = sprintf("%s%d-%d ",out_str,s,p)
			else
				out_str = sprintf("%s%d ",out_str,s)
			s=$c
		}
		p=l=$c
	}
}
END {
	if((l-s)>0)
		printf "%s%d-%d\n",out_str,s,l
	else if (length(s)>0)
		printf "%s%d\n",out_str,s
}
