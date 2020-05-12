#!/usr/bin/env python
from cffi import FFI
def env_header():
    header = """
    typedef struct{
        int iread, ivs, vpvs;
        int istacor, ivpvs, doshot, dotel, normal, havepo,
		icon, dmax, iaddcon, ipop, ittnum;
    }CommonEnv;

    typedef struct {
        int iread, ivs, nthres, kmin, ndiv, ndiv2;
        float vpvs, resthres, resthrep, stdmax;
        char leqsfil[MAXSTRLEN + 1], fsumfil[MAXSTRLEN + 1], outlfil[MAXSTRLEN + 1],
                    fhedfil[MAXSTRLEN + 1], fdatfil[MAXSTRLEN + 1];
    }LocEnv;

    typedef struct{
        int iray, iraystat, idatout, nomat, dmean, kmin, kmax, ido1d;
        float vpvsscale, resflag;
        //msut files
        char locdfil[MAXSTRLEN + 1], telrerr[MAXSTRLEN + 1], 
            dtdsfil[MAXSTRLEN + 1], resfile[MAXSTRLEN + 1], hitfile[MAXSTRLEN + 1], 
            dtdhfil[MAXSTRLEN + 1], bookfil[MAXSTRLEN + 1], sclefil[MAXSTRLEN + 1];
        //optinal files
        char telefil[MAXSTRLEN],  pbasfil[MAXSTRLEN], sbasfil[MAXSTRLEN], shotfil[MAXSTRLEN],
            elipfil[MAXSTRLEN], raystat[MAXSTRLEN],  dotfile[MAXSTRLEN], 
            headfil[MAXSTRLEN], entfile[MAXSTRLEN], stcfile[MAXSTRLEN], specfile[MAXSTRLEN];
    }SphraydervEnv;

    typedef struct{
        int intlim;
        float damper;
        char nmodfil[MAXSTRLEN + 1], fresfil[MAXSTRLEN + 1];
    }RunlsqrEnv;

    typedef struct{
        int mavx, mavy, mavz, nsmooth, limitu, ipscflg, ido1d;
        float dvperc, pertscl;
    }MakenewmodEnv;
    """

    func = """
    CommonEnv setCommonEnv(char *);
    void setCommonVariables(CommonEnv *, char *);
    void setLocFiles(LocEnv *, char *);
    void setLocVariables(LocEnv *, char *);
    LocEnv setLocEnv(char *);
    SphraydervEnv setSphraydervEnv(char *);
    void setSphrayderVariables(SphraydervEnv *, char *);
    void setSphrayderFiles(SphraydervEnv *, char *);
    RunlsqrEnv setRunlsqrEnv(char *);
    void setRunlsqrVariables(RunlsqrEnv *,char *);
    void setRunlsqrFiles(RunlsqrEnv *, char *);
    MakenewmodEnv setMakeNewmodEnv(char *);
    void setMakeNewmodVariables(MakenewmodEnv *, char *);
    """

    return header + func