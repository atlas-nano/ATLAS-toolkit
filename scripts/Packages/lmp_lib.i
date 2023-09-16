%include "exception.i"
%include cpointer.i              // Grab the SWIG pointer library

%exception {
    try {
        $action
        } 
    catch (const std::exception &e) {
        SWIG_exception_fail(SWIG_RuntimeError, e.what());
    }
}
%module lmp_lib
%include "std_vector.i"
%include "std_string.i"
%include "typemaps.i"

// =====various function redefinitions=====
void lammps_compute_get_dim(void *ptr, char *id, int *OUTPUT, int *OUTPUT);
void lammps_fix_get_dim(void *ptr, char *id, int *OUTPUT, int *OUTPUT);

// =====TYPEMAPS=======
//This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
    AV *tempav;
    I32 len;
    int i;
    SV  **tv;
    if (!SvROK($input))
        croak("Argument $argnum is not a reference.");
    if (SvTYPE(SvRV($input)) != SVt_PVAV)
        croak("Argument $argnum is not an array.");
    tempav = (AV*)SvRV($input);
    len = av_len(tempav);
    Newx($1,len+2,char *);
    for (i = 0; i <= len; i++) {
        tv = av_fetch(tempav, i, 0);
        $1[i] = (char *) SvPV(*tv,PL_na);
    }
    $1[i] = NULL;
};
// This cleans up the char ** array after the function call
%typemap(freearg) char ** {
    Safefree($1);
}
//Creates a new Perl array and places a NULL-terminated char ** into it
%typemap(out) char ** {
    AV *myav;
    SV **svs;
    int i = 0,len = 0;
    /* Figure out how many elements we have */
    while ($1[len])
       len++;
    Newx(svs,len,SV *);
    for (i = 0; i < len ; i++) {
        svs[i] = sv_newmortal();
        sv_setpv((SV*)svs[i],$1[i]);
    };
    myav =  av_make(len,svs);
    Safefree(svs);
    $result = newRV_noinc((SV*)myav);
    sv_2mortal($result);
    argvi++;
}

// Convert Perl array reference to double *
%typemap(in) double * {
    AV *tempav;
    I32 len;
    int i;
    SV  **tv;
    if (!SvROK($input))
        croak("Argument $argnum is not a reference.");
    if (SvTYPE(SvRV($input)) != SVt_PVAV)
        croak("Argument $argnum is not an array.");
    tempav = (AV*)SvRV($input);
    len = av_len(tempav);
    Newx($1,len+2,double);
    for (i = 0; i <= len; i++) {
        tv = av_fetch(tempav, i, 0);
        $1[i] = (double) SvNV(*tv);
    }
    $1[i] = 0;
};
// This cleans up the double * array after the function call
%typemap(freearg) double * {
    Safefree($1);
}
//Convert double * array to Perl array reference
%typemap(out) double * {
    AV *myav;
    SV **svs;
    int i = 0,len = 0;
    /* Figure out how many elements we have */
    while ($1[len])
       len++;
    Newx(svs,len,SV *);
    for (i = 0; i < len ; i++) {
        svs[i] = sv_newmortal();
        sv_setnv((SV*)svs[i],$1[i]);
    };
    myav =  av_make(len,svs);
    Safefree(svs);
    $result = newNV_noinc((SV*)myav);
    sv_2mortal($result);
    argvi++;
}

// Convert Perl array reference to int *
%typemap(in) int * {
    AV *tempav;
    I32 len;
    int i;
    SV  **tv;
    if (!SvROK($input))
        croak("Argument $argnum is not a reference.");
    if (SvTYPE(SvRV($input)) != SVt_PVAV)
        croak("Argument $argnum is not an array.");
    tempav = (AV*)SvRV($input);
    len = av_len(tempav);
    Newx($1,len+2,int);
    for (i = 0; i <= len; i++) {
        tv = av_fetch(tempav, i, 0);
        $1[i] = (int) SvIV(*tv);
    }
    $1[i] = 0;
};
// This cleans up the int * array after the function call
%typemap(freearg) int * {
    Safefree($1);
}
//Convert int * array to Perl array reference
%typemap(out) int * {
    AV *myav;
    SV **svs;
    int i = 0,len = 0;
    /* Figure out how many elements we have */
    while ($1[len])
       len++;
    Newx(svs,len,SV *);
    for (i = 0; i < len ; i++) {
        svs[i] = sv_newmortal();
        sv_setnv((SV*)svs[i],$1[i]);
    };
    myav =  av_make(len,svs);
    Safefree(svs);
    $result = newIV_noinc((SV*)myav);
    sv_2mortal($result);
    argvi++;
}

//Convert a Perl reference to a void **
%typemap(in, numinputs=0) void ** (void * retval) {
  $1 = (void**)&retval;
}
%typemap(argout) void ** {
    %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum),$1_descriptor, 0));
}

//=====Perl inline helper functions======
%inline %{
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "library.h"

//Convert C void * to int
extern int c_void_p_int(void *ptr) 
{
    return *(static_cast<int *>(ptr));
}

//Convert C void * to double
extern double c_void_p_dbl(void *ptr) 
{
    return *(static_cast<double *>(ptr));
}

//Convert C void * to double array and return Perl reference
extern SV* c_void_p_array(void *ptr, size_t i, int type) 
{
    int j;
    AV *myav;
    SV *result;
    SV **svs;
    double *dptr = static_cast<double *>(ptr); //cast void as double pointer to allow arithmatic
    int *iptr = static_cast<int *>(ptr); //cast void as int pointer to allow arithmatic

    Newx(svs,i,SV *);
    for(j=0;j<i;j++) {
        svs[j] = sv_newmortal();
        if(type == 0) {
            sv_setiv((SV*)svs[j], (*iptr));
            iptr++;
        } else {
            sv_setnv((SV*)svs[j], (*dptr));
            dptr++;
        }
    }
    myav = av_make(i,svs);

    Safefree(svs);
    dptr = NULL;
    iptr = NULL;

    result = newRV_noinc((SV*)myav);
    return result;
}

//Convert C void ** to 2d array and return Perl reference
extern SV* c_void_p_2d_array(void *ptr, size_t i, size_t j, int type) 
{
    int k,l,t,idx;
    AV *myav;
    SV *result;
    SV **svs;
    double **dptr = (double **)((void *)ptr); //resolve void * to double ** pointer
    int **iptr = (int **)((void *)ptr); //resolve void * to double ** pointer
    
    t=i*j;
    Newx(svs,t,SV *);
    idx=0;
    
    for(k=0;k<i;k++) {
        for(l=0;l<j;l++) {
            svs[idx] = sv_newmortal();
            if(type == 0) 
                sv_setiv((SV*)svs[idx], iptr[k][l]);
            else
                sv_setnv((SV*)svs[idx], dptr[k][l]);
            idx++;
        }
    }
    
    myav =  av_make(t,svs);
    Safefree(svs);
    iptr=NULL;
    dptr=NULL;
    result = newRV_noinc((SV*)myav);
    return result;
}

extern SV* _perl_wrap_extract_box(void *ptr) 
{
    double *boxlo = new double [3];
    double *boxhi = new double [3];
    double xy,yz,xz;
    int *p = new int [3];
    int c;

    lammps_extract_box(ptr,boxlo,boxhi,&xy,&xz,&xz,p,&c);

    //convert to LAMMPS array reference
    AV *myav;
    SV **svs;
    SV *result;
    int len = 13;
    int i;

    Newx(svs,len,SV *);
    for(i=0;i<len;i++) 
        svs[i] = sv_newmortal();

    sv_setnv((SV*)svs[0], boxlo[0]);
    sv_setnv((SV*)svs[1], boxlo[1]);
    sv_setnv((SV*)svs[2], boxlo[2]);
    sv_setnv((SV*)svs[3], boxhi[0]);
    sv_setnv((SV*)svs[4], boxhi[1]);
    sv_setnv((SV*)svs[5], boxhi[2]);
    sv_setnv((SV*)svs[6], xy);
    sv_setnv((SV*)svs[7], yz);
    sv_setnv((SV*)svs[8], xz);
    sv_setnv((SV*)svs[9], (double)p[0]);
    sv_setnv((SV*)svs[10],(double)p[1]);
    sv_setnv((SV*)svs[11],(double)p[2]);
    sv_setnv((SV*)svs[12],(double)c);

    myav = av_make(len, svs);
    Safefree(svs);
    delete [] boxlo;
    delete [] boxhi;
    delete [] p;
    result = newRV_noinc((SV*)myav);
    sv_2mortal(result);
    return result;
}

extern SV* _perl_wrap_gather_atoms(void *ptr, char *name, int type, int len, int tot) 
{
    double *data_dbl = new double[tot];
    int *data_int = new int[tot];
    void *ptr_data;

    for(int j = 0; j<len; j++) {
        data_dbl[j] = 0.0;
        data_int[j] = 0;
    }
    if(type == 0) {
        ptr_data = static_cast<void *>(data_int);
    } else {
        ptr_data = static_cast<void *>(data_dbl);
    }
    lammps_gather_atoms(ptr,name,type,len,ptr_data);

    //convert to LAMMPS array reference
    AV *myav;
    SV **svs;
    SV *result;
    double *dval = static_cast<double *>(ptr_data);
    int *ival = static_cast<int *>(ptr_data);
    int i;

    Newx(svs,tot,SV *);
    for(i=0;i<tot;i++) {
        svs[i] = sv_newmortal();
        if(type == 0) {
            sv_setiv((SV*)svs[i],*ival);
            ival++;
        } else {
            sv_setnv((SV*)svs[i],*dval);
            dval++;
        }
    }

    myav = av_make(tot, svs);
    result = newRV_noinc((SV*)myav);
    sv_2mortal(result);

    Safefree(svs);
    delete [] data_dbl;
    delete [] data_int;
    dval = NULL;
    ival = NULL;
    ptr_data = NULL;

    return result;
                          
}

extern void _perl_wrap_lammps_scatter_atoms(void *ptr, char *name, int type, int count, int len, SV* array_ref)
{
    int i;
    double *data_dbl = new double [len];
    int *data_int = new int [len];
    void *vptr;

    AV *tempav;
    SV  **tv;

    //Error checking
    if (!SvROK(array_ref))
        croak("Argument $argnum is not a reference.");
    if (SvTYPE(SvRV(array_ref)) != SVt_PVAV)
        croak("Argument $argnum is not an array.");

    //Convert perl array ref -> C array and create void ptr to it
    tempav = (AV*)SvRV(array_ref);
    for (i = 0; i < len; i++) {
        tv = av_fetch(tempav, i, 0);
        if(type == 0) 
            data_int[i] = (int) SvIV(*tv);
        else
            data_dbl[i] = (double) SvNV(*tv);
    }
    if(type == 0)
        vptr = (void *) data_int;
    else
        vptr = (void *) data_dbl;

    //call scatter_atoms
    lammps_scatter_atoms(ptr,name,type,count,vptr);

    //cleanup
    delete [] data_dbl;
    delete [] data_int;
    vptr = NULL;
}

%}
/* Put headers and other declarations here */
%{
#include "library.h"
#include <iostream>
%}
%include "library.h"

