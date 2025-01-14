#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "FSearch.h"

static int
not_here(char *s)
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

static double
constant(char *name, int arg)
{
    errno = 0;
    switch (*name) {
    }
    errno = EINVAL;
    return 0;

not_there:
    errno = ENOENT;
    return 0;
}


MODULE = FSearch		PACKAGE = FSearch		


double
constant(name,arg)
	char *		name
	int		arg

void
setResidual(solvAtoms,rMax,j,tot,residuals,vRes)
	int *		solvAtoms
	double		rMax
	int		j
	int		tot
	double *	residuals
	int *		vRes

int
getResidual(arrayIndex,distArray,residuals,validRes,tot,boxLen)
	int		arrayIndex
	double **	distArray
	double *	residuals
	int *		valiRes
	int		tot
	double		boxlen

