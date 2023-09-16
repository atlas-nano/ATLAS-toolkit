%module p5namot
%include "typemaps.i"
%{
#include "cmds_h.h"
#include "scriptinter.h"
#include "brokenlinks.h"
#include "snapshot.h"
#include "flags.h"
#include "defs.h"
#include "general.h"

%}

%init %{
    Inits();
    FLAGS=(FLAGS|Europa)-Europa;
%}

%constant int  NAMOTERR=-1;
%constant int  NAMOTOK=0;

%rename(QueryAngle) SIQueryAngle;
%rename(Cmd) SICmd;
%rename(QueryCm) SIQueryCm;
%rename(QueryAtomDist) SIQueryAtomDist;
%rename(QueryAtomPos) SIQueryAtomPos;
%rename(QueryScale) SIQueryScale;
%rename(QueryTors) SIQueryTors;
%rename(QueryBaseParameter) SIQueryBaseParameter;
%rename(QueryPhopParameter) SIQueryPhopParameter;
%rename(QuerySugParameter) SIQuerySugParameter;
%rename(QueryUnitParameter) SIQueryUnitParameter;
%rename(QuerySSParameter) SIQuerySSParameter;
%rename(QueryMCG) SIQueryMCG;
%rename(QueryHUB) SIQueryHUB;
%rename(GroupTagFind) SIGroupTagFind;
%rename(FuranoseRingGenerate) SIFuranoseRingGenerate;
%rename(ChainNumberOfGroups) SIChainNumberOfGroups;
%rename(HelicesNumberOf) SIHelicesNumberOf;
%rename(ChainNumberOf) SIChainNumberOf;
%rename(CalcPenalty) FFCalcEnergy;

extern int Inits();
extern int SIQueryAngle(int,char **,double *OUTPUT);
extern int SIQueryCm(int,char **,double *OUTPUT,double *OUTPUT,double *OUTPUT);
extern int SIQueryAtomDist(int,char **,double *OUTPUT);
extern int SIQueryAtomPos(int,char **,double *OUTPUT,double *OUTPUT,
	    double *OUTPUT);
extern int SIQueryScale(double *OUTPUT);
extern int SIQueryTors(int,char **,double *OUTPUT);
extern int SICmd(char *);
extern int SIQueryBaseParameter(int,char **,double *OUTPUT);
extern int SIQueryPhopParameter(int,char **,double *OUTPUT);
extern int SIQuerySugParameter(int,char **,double *OUTPUT);
extern int SIQueryUnitParameter(int,char **,double *OUTPUT);
extern int SIQuerySSParameter(int,char **,double *OUTPUT);
extern char *SIQueryMCG();
extern char *SIQueryHUB();
extern int SIFuranoseRingGenerate(double,double,double[],double[],double[]);
extern int SIGroupTagFind(char *,int  *OUTPUT,int  *OUTPUT,int  *OUTPUT);

extern int BrokenLinkNum(int *OUTPUT);
extern int BrokenLinkGet(int,int *OUTPUT,int *OUTPUT,int *OUTPUT,
	    int *OUTPUT,int *OUTPUT,int *OUTPUT,
	    double *OUTPUT,double *OUTPUT,double *OUTPUT,double *OUTPUT);

extern int SnapShotAdd(char *);
extern int SnapShotNumber(int *OUTPUT);
extern int SnapShotClear(char *);
extern int SnapShotRevert(char *);
extern char *SnapShotName(int);
extern int SIChainNumberOfGroups(int, int, int *OUTPUT);
extern int SIHelicesNumberOf(int  *OUTPUT);
extern int SIChainNumberOf(int,int *OUTPUT);
extern int FFCalcEnergy(double *OUTPUT,double *OUTPUT,double *OUTPUT,
    double *OUTPUT,double *OUTPUT);









