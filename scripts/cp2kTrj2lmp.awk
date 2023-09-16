BEGIN{
    if(cell==""){
        print "ERROR: Expected a b c [alpha=90] [beta=90] [gamma=90] for cell variable, but not defined. Cannot continue..."
        exit(1)
    }
    if(ff_str>"") {
        n=split(ff_str,tmp," ")
        for(i=1;i<=n;i++) {
            fftype[tmp[i]]=i
        }
    }
    c[4]=c[5]=c[6]=90
    n=split(cell,c," ")
    if(n<3){
        print "ERROR: Expected a b c [alpha=90] [beta=90] [gamma=90] for cell variable. Got '" cell "'. Cannot continue..."
        exit(1)
    }
    valid=0
    s=21.86681818181818181817; #Bohr/au -> Angstrom/ps
}
{
    if(NR==1) {
        num_atoms=$1; 
        snap_size=num_atoms+2
        atom_counter=0
        valid=0
    }
    if(NR%snap_size==1) {
        tstep++
        valid=0
        atom_counter=0
        print "ITEM: TIMESTEP"
        printf "%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n0 %f\n0 %f\n0 %f\n",tstep,num_atoms,c[1],c[2],c[3]
        print "ITEM: ATOMS id type xu yu zu vx vy vz"
    }
    if(NR%snap_size==3){
        valid=1
    }
    if(valid==1) {
        atom_counter++
        if(!fftype[$1]) {
            fftype[$1]=++type_counter
        }
        printf "%d %d %.10f %.10f %.10f %.10f %.10f %.10f\n",atom_counter,fftype[$1],$2,$3,$4,$6*s,$7*s,$8*s
    }   
}
