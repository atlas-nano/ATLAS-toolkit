%module FSearch
%{
#include "FSearch.h"
%}
%include "typemaps.i"
%include "typemaps.i"
%typemap(perl5,in) double * (double dvalue) {
  SV* tempsv;
  if (!SvROK($input)) {
    croak("expected a reference\n");
  }
  tempsv = SvRV($input);
  if ((!SvNOK(tempsv)) && (!SvIOK(tempsv))) {
    croak("expected a double reference\n");
  }
  dvalue = SvNV(tempsv);
  $1 = &dvalue;
}

%typemap(perl5,argout) double * {
  SV *tempsv;
  tempsv = SvRV($arg);
  sv_setnv(tempsv, *$input);
}

%typemap(perl5,in) int * (int dvalue) {
  SV* tempsv;
  if (!SvROK($input)) {
    croak("expected a reference\n");
  }
  tempsv = SvRV($input);
  if ((!SvIOK(tempsv))) {
    croak("expected a int reference\n");
  }
  dvalue = SvIV(tempsv);
  $1 = &dvalue;
}
                                                                                                                       
void setResidual(int *, int, int, int, double *, int *);
int getResidual(int, double **, double *, int *, int, int);
