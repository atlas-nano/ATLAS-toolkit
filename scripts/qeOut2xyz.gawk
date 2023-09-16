#!/bin/gawk
{
	if($0 ~ /number of atoms\/cell/)
		natom=$NF
	if($0 ~ /^!    total energy/) {
		ceng=$5
		nstep++
	}
	if($0 ~ /ATOMIC_POSITIONS/) {
		print natom
		print "step ",nstep," energy ",ceng," Ry"
		write_coords=1
	} else if (write_coords && NF != 4)
		write_coords = 0
	else if	(write_coords)
		print
}
