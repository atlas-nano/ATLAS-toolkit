#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

void dgesvd_(char *JOBU,char *JOBVT,int *M,int *N,double *A,int *LDA,
             double *S, double *U,int *LDU,double *VT,int *LDVT,
             double *WORK,int *LWORK,int *INFO);

MODULE = PDBsub		PACKAGE = PDBsub		

int
do_svd3(R_temp,e_temp,V_temp)
AV *R_temp
AV *e_temp
AV *V_temp
CODE:
  int i,three=3,lwork=256,info;
  double work[256],R[9],e[3],U[9],V[9];
  char A='A';
  
  for(i=0;i<9;i++) {
    R[i]=SvNV(*av_fetch(R_temp,i,0));
  }
  dgesvd_(&A,&A,&three,&three,R,&three,e,U,&three,V,&three,work,&lwork,&info);
  RETVAL=info;
  for(i=0;i<3;i++) {
    av_store(e_temp,i,newSVnv(e[i]));
  }
  for(i=0;i<9;i++) {
    av_store(V_temp,i,newSVnv(V[i]));
  }
OUTPUT:
  RETVAL
