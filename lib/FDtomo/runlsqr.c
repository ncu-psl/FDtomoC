// c 
// c---------------------------------------------------------------------------
// c
// c  Current version: 2004.0924
// c
// c  Original code: Steve Roecker, RPI, January, 1993
// c
// c  This is the driver routine for running lsqr, modified from the program
// c  main.f I got from Gary Pavlis, which he modified from a program by
// c  Sally Szpakowski.
// c
// c  The three arrays a,
// c  ja, and na are created and stored in common. 
 
// c  Note that the subroutine aprod must be declared external if it is used
// c  only in the call to lsqr.
// c
// c
// c
// c  The following is a list of input parameters that must be dealt
// c  with before running the program:
// c
// c	damp - the damping factor which controls lsqr (see lsqr comments)
// c	m - the number of rows in the  matrix a
// c	n - the number of columns in the matrix a
// c	SIZEOFA - the size of a
// c	atol,btol,conlim,intlim,nout - all used in lsqr see  comments there
// c
// c	MMAX is the maximum number of rows, which is the same as the number of phases
// c	NMAX is the maximum number of columns, which is the same as the number of variables
// c
// c       Input File attachments:
// c
// c       fn      Name   Uname   U#  Comments
// c       11     dtdsfil lundts  19  Existing A Matrix
// c       12     resfile lunres  20  Existing data vector
// c
// c       Output File attachments:
// c
// c       fn      Name   Uname   U#  Comments
// c                      lout     2  runlsqr.log Log File
// c       36     nmodfil lunfmd  45  Perturbation file (output)
// c       37     fresfil lunfrs  46  Residual file (output)
// c
// c
// c--------------------------------------------------------------------------
// c
// c  functions and subroutines
// c     main    lsqr,remakea,aprod
// c     blas    scopy,snrm2,sscal
// c     fortran   sqrt,max
// c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <memory.h>
#include "common/string_process.h"
#include "common/geographic_method.h"
#include "common/parameter.h"
#include "common/gridspec.h"
#include "common/parseprogs.h"
#include "common/shared_variables.h"
#include "common/read_spec.h"
#include "FDtomo/runlsqr.h"

#include "runlsqr/aprod.h"
#include "runlsqr/makea.h"
#include "runlsqr/normlz.h"
#include "runlsqr/scopy.h"
#include "runlsqr/snrm2.h"
#include "runlsqr/sscal.h"

void lsqr(int , int , float , int , int ,int *, float *, float *, float *, float *, float *, float *, float , float , float , int , FILE *, int , float , float , float *, float *, float *, FILE *);

// c---gfortran objects to nout being declared here, so it is initialized below
// c	parameter(nout=6)
// c--mustv is the number of variables that must be assigned in the parameter
// c       file in order for this program to run correctly.  The names of the
// c       variables are specified in the mvals string array below.
// c--mustf is the number of files that must be attached; see the files array below.
#define MUSTV 4
#define MUSTF 4
#define MAXSTRLEN 132

int m, n, istop, intlim, j;

FILE *fp_nout = NULL;

// c  functions and local variables
int *ja, *na, *jndx;
float *b, u[MMAX];
float damp, v[NMAX], w[NMAX], x[NMAX], se[NMAX];
float atoL, btol, conlim, anorm;
float acond, rnorm, arnorm, dampsq, xnorm;
float *a;


char logfile[80 + 1];
float one = 1.0f;

int runlsqr(SPEC spec, SPHRAYDERV_DATA *SPHRAYDERV) {
	nxc = spec.nxc; nyc = spec.nyc; nzc = spec.nzc; 
	nx = spec.nx;   ny = spec.ny;   nz = spec.nz;
	
	h = spec.h; x0 = spec.x0; y = spec.y; 
	z0 = spec.z0; dq = spec.dq; df = spec.df; x00 = spec.x00; y00 = spec.y00;
	igridx = spec.igridx; igridy = spec.igridy; igridz = spec.igridz;

	float damper = spec.damper;
	int intlims = spec.intlims, ittnum = spec.ittnum;
	char dtdsfil[MAXSTRLEN + 1], resfile[MAXSTRLEN + 1], nmodfil[MAXSTRLEN + 1],
		fresfil[MAXSTRLEN + 1];
	char VERSION[10] = "2018.0114\0";
	
	strcpy(dtdsfil, spec.dtdsfil);
	strcpy(resfile, spec.resfile);
	strcpy(nmodfil, spec.nmodfil);
	strcpy(fresfil, spec.fresfil);

	int i, len, ierr;

	sprintf(logfile, "runlsqr.log%d\n", ittnum);
	FILE *fp_log = fopen(logfile, "w");
	if (!fp_log) {
		printf("(Error in runlsqr.c)create fp_log file error.\n");
		assert(0);
	}
	fprintf(fp_log, "  \n");
	fprintf(fp_log,
			" *************************************************************** \n");
	fprintf(fp_log, "          Parameters Set For This Run of runlsqr.c\n");
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " VERSION: %s\n", VERSION);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Current parameter specification file: %-40s\n",
			spec.spec_file);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " iteration counter: %d \n", ittnum);
	fprintf(fp_log, "  \n");

	fprintf(fp_log, " Damper: %f\n", damper);
	fprintf(fp_log, " intlim: %d\n", intlims);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Input file attachments: \n");
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " A matrix: %-60s \n", dtdsfil);
	fprintf(fp_log, " Data Vector: %-60s \n", resfile);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Output file attachments: \n");
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Perturbations: %-60s \n", nmodfil);
	fprintf(fp_log, "  \n");
	fprintf(fp_log,
			" *************************************************************** \n");

	damp = damper;
	FILE *fp_fmd = fopen(nmodfil, "w");
	if (!fp_fmd) {
		printf("file create error: %s\n", nmodfil);
		assert(0);
	}
	FILE *fp_frs = fopen(fresfil, "w");
	if (!fp_frs) {
		printf("file create error: %s\n", fresfil);
		assert(0);
	}
// c----------------------------------------------------------------------
// c
// c  generate a and b.  the vector r, ja, and na will define the matrix a.
// c
// c------------------------------------------------------------------------
	
	printf("  Reading in a and b ... \n");
	fprintf(fp_log, "  Reading in a and b ... \n");
	int nbl;

	m = SPHRAYDERV->mat->number_rows;
	n = nbl = SPHRAYDERV->mat->number_columns;
	jndx = SPHRAYDERV->mat->jndx;
	a = SPHRAYDERV->mat->elements;
	na = SPHRAYDERV->mat->elements_row;
	ja = SPHRAYDERV->mat->column_elements;
	b = SPHRAYDERV->b;

	
// c---------------------------------------------------------------------
// c
// c  solve the problem defined by aprod, damp, and b
// c
// c----------------------------------------------------------------------

// c  copy the rhs vector b into u (lsqr will overwrite u) and set the
// c  other input parameters for lsqr.  change the value of damp as needed.
	
	scopy(m, b, 1, u, 1);

// c	write(*,*) 'enter the damping factor'
// c	read(*,*) damp
// c	damp = 0.
// c--Documentation for lsqr says this sets atoL to eps
// c	atoL = 0.0
// c	write(*,*) 'Enter data relative precision (BTOL)'
// c	read(*,*) btol
// c
// c	write(*,*) 'enter the conlim'
// c	read(*,*) conlim
// c
// c
// c	relpr = 1.0e-6
	float relpr = 1.0e-12;
	atoL = relpr;
	btol = relpr;
	conlim = 1. / (10. * sqrt(relpr));
// c  This value is of intlim is appropriate for ill conditioned systems
	intlim = 4 * n;
	printf("  Suggested intlim = %12d\n", intlim);
	fprintf(fp_log, "  Suggested intlim = %12d\n", intlim);
	if (intlims != 0) {
		intlim = intlims;
	}
	if(fp_nout) {
		fprintf(fp_nout, " least-squares test problem      p(%8d%8d  %10.6f%10.2f%8d\n\n\n", m, n, damp, conlim, intlim);
	} else {
		printf(" least-squares test problem      p(%8d%8d  %10.6f%10.2f%8d\n\n\n", m, n, damp, conlim, intlim);
	}
	fprintf(fp_log, " %d %d %f %f %d \n", m, n, damp, conlim, intlim);

	lsqr(m, n, damp, 1, 1, ja, a, u, v, w, x, se, atoL, btol, conlim, intlim, fp_nout, istop, anorm, acond, &rnorm, &arnorm, &xnorm, fp_log);
// c-----------------------------------------------------------------------
// c
// c  examine the results.
// c  we print the residual norms rnorm and arnorm given by lsqr,
// c  and then compute their true values in terms of the solution x
// c  obtained by lsqr.  at least one of them should be small.
// c
// c------------------------------------------------------------------------

	one = 1.0f;
	dampsq = damp * damp;

	printf("                  residual norm (abar*x - bbar)    residual norm (normal eqns)    solution norm (x)\n");
	printf("                  -----------------------------    ---------------------------    -----------------\n\n");
	printf(" estimated by lsqr          %.6E                   %.6E              %.6E\n", rnorm, arnorm, xnorm);

// c  compute u = a*x
// c  Do this using aprod by setting u to zero then forming the product
// c  u <- A*x - u
// c  Thus u will contain predicted travel times

	for (int i = 0; i < m; i++) {
		u[i] = 0.;
	}
	int *iw = NULL;
	float *rw = NULL;
	aprod(1, m, n, x, u, 1, 1, iw, rw, a, na, ja);
	for (int i = 0; i < m; i++) {
		fprintf(fp_frs, "%lf %lf %lf", b[i], u[i], u[i] - b[i]);
	}
// c  compute v = a(transpose)*u + damp**2 * x
// c  this will be close to zero in all cases if x is close to
// c  a solution.
// c
// c--Sally's old code calculated several residual norms.
// c--since I removed the residual calculation, these are no longer
// c--correct.
// c      call scopy(n,x,1,v,1)
// c      call sscal(n,dampsq,v,1)
// c      call aprod(2,m,n,v,u,SIZEOFA,SIZEOFA,iw,rw)
// c
// c  compute the norms associated with x,u,v.
// c
// c      xnorm = snrm2(n,x,1)
// c      rnorm = sqrt(snrm2(m,u,1)**2 + dampsq*xnorm**2)
// c      arnorm = snrm2(n,v,1)
// c      write(nout,2200) rnorm,arnorm,xnorm
// c	write(9,3500) (x(j), j=1,n)
// c	do 200 j = 1, n
// c	  sc = scale(j)
// c	  write(19,3500) x(j)/sc, se(j)/sc, x(j), sc , jndx(j)
// c	  x(j) = x(j)/sc
// c	  se(j) = se(j)/sc
// c200	continue
	fwrite(&n, sizeof(n), 1, fp_fmd);
	fwrite(x, sizeof(x[0]), n, fp_fmd);
	fwrite(jndx, sizeof(jndx[0]), n, fp_fmd);
	fwrite(se, sizeof(se[0]), n, fp_fmd);
	return 0;
}

// c     ------------------------------------------------------------------
// c
// c     lsqr  finds a solution  x  to the following problems...
// c
// c     1. unsymmetric equations --    solve  a*x = b
// c
// c     2. linear least squares  --    solve  a*x = b
// c                                    in the least-squares sense
// c
// c     3. damped least squares  --    solve  (   a    )*x = ( b )
// c                                           ( damp*i )     ( 0 )
// c                                    in the least-squares sense
// c
// c     where  a  is a matrix with  m  rows and  n  columns, b  is an
// c     m-vector, and  damp  is a scalar (all quantities real).
// c     the matrix  a  is intended to be large and sparse.  it is accessed
// c     by means of subroutine calls of the form
// c
// c                call aprod( mode,m,n,x,y,leniw,lenrw,iw,rw )
// c
// c     which must perform the following functions...
// c
// c                if mode = 1, compute  y = y + a*x.
// c                if mode = 2, compute  x = x + a(transpose)*y.
// c
// c     the vectors x and y are input parameters in both cases.
// c     if mode = 1, y should be altered without changing x.
// c     if mode = 2, x should be altered without changing y.
// c     the parameters leniw, lenrw, iw, rw may be used for workspace
// c     as described below.
// c
// c     the rhs vector  b  is input via  u,  and subsequently overwritten.
// c
// c
// c     note.  lsqr uses an iterative method to approximate the solution.
// c     the number of iterations required to reach a certain accuracy
// c     depends strongly on the scaling of the problem.  poor scaling of
// c     the rows or columns of  a  should therefore be avoided where
// c     possible.
// c
// c     for example, in problem 1 the solution is unaltered by
// c     row-scaling.  if a row of  a  is very small or large compared to
// c     the other rows of  a,  the corresponding row of  ( a  b )  should
// c     be scaled up or down.
// c
// c     in problems 1 and 2, the solution  x  is easily recovered
// c     following column-scaling.  in the absence of better information,
// c     the nonzero columns of  a  should be scaled so that they all have
// c     the same euclidean norm (e.g.  1.0).
// c
// c     in problem 3, there is no freedom to re-scale if  damp  is
// c     nonzero.  however, the value of  damp  should be assigned only
// c     after attention has been paid to the scaling of  a.
// c
// c     the parameter  damp  is intended to help regularize
// c     ill-conditioned systems, by preventing the true solution from
// c     being very large.  another aid to regularization is provided by
// c     the parameter  acond,  which may be used to terminate iterations
// c     before the computed solution becomes very large.
// c
// c
// c     notation
// c     --------
// c
// c     the following quantities are used in discussing the subroutine
// c     parameters...
// c
// c     abar   =  (   a    ),          bbar  =  ( b )
// c               ( damp*i )                    ( 0 )
// c
// c     r      =  b  -  a*x,           rbar  =  bbar  -  abar*x
// c
// c     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
// c            =  norm( rbar )
// c
// c     relpr  =  the relative precision of floating-point arithmetic
// c               on the machine being used.  for example, on the ibm 370,
// c               relpr  is about 1.0e-6 and 1.0d-16 in single and double
// c               precision respectively.
// c
// c     lsqr  minimizes the function  rnorm  with respect to  x.
// c
// c
// c     parameters
// c     ----------
// c
// c     m       input      the number of rows in  a.
// c
// c     n       input      the number of columns in  a.
// c
// c     aprod   external   see above.
// c
// c     damp    input      the damping parameter for problem 3 above.
// c                        (damp  should be 0.0 for problems 1 and 2.)
// c                        if the system  a*x = b  is incompatible, values
// c                        of  damp  in the range 0 to sqrt(relpr)*norm(a)
// c                        will probably have a negligible effect.
// c                        larger values of  damp  will tend to decrease
// c                        the norm of  x  and to reduce the number of
// c                        iterations required by lsqr.
// c
// c                        the work per iteration and the storage needed
// c                        by lsqr are the same for all values of  damp.
// c
// c     leniw   input      the length of the workspace array  iw.
// c     lenrw   input      the length of the workspace array  rw.
// c     iw      workspace  an integer array of length  leniw.
// c     rw      workspace  a real array of length  lenrw.
// c
// c             note.  lsqr does not explicitly use the previous four
// c             parameters, but passes them to subroutine aprod for
// c             possible use as workspace.  if aprod does not need
// c             iw  or  rw,  the values  leniw = 1  or  lenrw = 1  should
// c             be used, and the actual parameters corresponding to
// c             iw  or  rw  may be any convenient array of suitable type.
// c
// c     u(m)    input      the rhs vector  b.  beware that  u  is
// c                        over-written by lsqr.
// c
// c     v(n)    workspace
// c     w(n)    workspace
// c
// c     x(n)    output     returns the computed solution  x.
// c
// c     se(n)   output     returns standard error estimates for the
// c                        components of  x.  for each  i,  se(i)  is set
// c                        to the value  rnorm * sqrt( sigma(i,i) / t ),
// c                        where  sigma(i,i)  is an estimate of the i-th
// c                        diagonal of the inverse of abar(transpose)*abar
// c                        and  t = 1      if  m .le. n,
// c                             t = m - n  if  m .gt. n  and  damp = 0,
// c                             t = m      if  damp .ne. 0.
// c
// c     atol    input      an estimate of the relative error in the data
// c                        defining the matrix  a.  for example,
// c                        if  a  is accurate to about 6 digits, set
// c                        atol = 1.0e-6 .
// c
// c     btol    input      an estimate of the relative error in the data
// c                        defining the rhs vector  b.  for example,
// c                        if  b  is accurate to about 6 digits, set
// c                        btol = 1.0e-6 .
// c
// c     conlim  input      an upper limit on  cond(abar),  the apparent
// c                        condition number of the matrix  abar.
// c                        iterations will be terminated if a computed
// c                        estimate of  cond(abar)  exceeds  conlim.
// c                        this is intended to prevent certain small or
// c                        zero singular values of  a  or  abar  from
// c                        coming into effect and causing unwanted growth
// c                        in the computed solution.
// c
// c                        conlim  and  damp  may be used separately or
// c                        together to regularize ill-conditioned systems.
// c
// c                        normally,  conlim  should be in the range
// c                        1000  to  1/relpr.
// c                        suggested value --
// c                        conlim = 1/(100*relpr)  for compatible systems,
// c                        conlim = 1/(10*sqrt(relpr))  for least squares.
// c
// c             note.  if the user is not concerned about the parameters
// c             atol, btol  and  conlim,  any or all of them may be set
// c             to zero.  the effect will be the same as the values
// c             relpr, relpr  and  1/relpr  respectively.
// c
// c     itnlim  input      an upper limit on the number of iterations.
// c                        suggested value --
// c                        itnlim = n/2     for well conditioned systems,
// c                        itnlim = 4*n     otherwise.
// c
// c     nout    input      file number for printer.  if positive,
// c                        a summary will be printed on file  nout.
// c
// c     istop   output     an integer giving the reason for termination...
// c
// c                0       x = 0  is the exact solution.
// c                        no iterations were performed.
// c
// c                1       the equations  a*x = b  are probably
// c                        compatible.  norm(a*x - b)  is sufficiently
// c                        small, given the values of  atol  and  btol.
// c
// c                2       the system  a*x = b  is probably not
// c                        compatible.  a least-squares solution has
// c                        been obtained which is sufficiently accurate,
// c                        given the value of  atol.
// c
// c                3       an estimate of  cond(abar)  has exceeded
// c                        conlim.   the system  a*x = b  appears to be
// c                        ill-conditioned.  otherwise, there could be an
// c                        an error in subroutine aprod.
// c
// c                4       the equations  a*x = b  are probably
// c                        compatible.  norm(a*x - b)  is as small as
// c                        seems reasonable on this machine.
// c
// c                5       the system  a*x = b  is probably not
// c                        compatible.  a least-squares solution has
// c                        been obtained which is as accurate as seems
// c                        reasonable on this machine.
// c
// c                6       cond(abar)  seems to be so large that there is
// c                        not much point in doing further iterations,
// c                        given the precision of this machine.
// c                        there could be an error in subroutine aprod.
// c
// c                7       the iteration limit  itnlim  was reached.
// c
// c     anorm   output     an estimate of the frobenius norm of  abar.
// c                        this is the square-root of the sum of squares
// c                        of the elements of  abar.
// c                        if  damp  is small and if the columns of  a
// c                        have all been scaled to have length  1.0,
// c                        anorm  should increase to roughly  sqrt(n).
// c                        a radically different value for  anorm  may
// c                        indicate an error in subroutine aprod (there
// c                        may be an inconsistency between modes 1 and 2).
// c
// c     acond   output     an estimate of  cond(abar),  the condition
// c                        number of  abar.  a very high value of  acond
// c                        may again indicate an error in aprod.
// c
// c     rnorm   output     an estimate of the final value of norm(rbar),
// c                        the function being minimized (see notation
// c                        above).  this will be small if  a*x = b  has
// c                        a solution.
// c
// c     arnorm  output     an estimate of the final value of
// c                        norm( abar(transpose)*rbar ), the norm of
// c                        the residual for the usual normal equations.
// c                        this should be small in all cases.  (arnorm
// c                        will often be smaller than the true value
// c                        computed from the output vector  x.)
// c
// c     xnorm   output     an estimate of the norm of the final
// c                        solution vector  x.
// c
// c
// c     subroutines and functions used
// c     ------------------------------
// c
// c     user       aprod
// c     lsqr       normlz
// c     blas       scopy,snrm2,sscal  (see lawson et al. below)
// c                (snrm2 is used only in normlz)
// c     fortran    abs,mod,sqrt
// c
// c
// c     precision
// c     ---------
// c
// c     the number of iterations required by lsqr will usually decrease
// c     if the computation is performed in higher precision.  to convert
// c     lsqr  and  normlz  between single- and double-precision, change
// c     the words
// c                scopy, snrm2, sscal
// c                abs, real, sqrt
// c     to the appropriate blas and fortran equivalents.
// c
// c
// c     references
// c     ----------
// c
// c     paige, c.c. and saunders, m.a.  lsqr: an algorithm for sparse
// c        linear equations and sparse least squares.
// c        acm transactions on mathematical software 8, 1 (march 1982).
// c
// c     lawson, c.l., hanson, r.j., kincaid, d.r. and krogh, f.t.
// c        basic linear algebra subprograms for fortran usage.
// c        acm transactions on mathematical software 5, 3 (sept 1979),
// c        308-323 and 324-325.
// c
// c
// c     lsqr.      this version dated 22 february 1982.
// c     ------------------------------------------------------------------
// c
// c     functions and local variables
void lsqr(int m, int n, float damp, int leniw, int lenrw, int *iw, float *rw,
		float *u, float *v, float *w, float *x, float *se, float atoL,
		float btol, float conlim, int itnlim, FILE *fp_nout, int istop,
		float anorm, float acond, float *rnorm, float *arnorm, float *xnorm,
		FILE *fp_out) {

	if (fp_nout != NULL) {
		fprintf(fp_nout, "%25slsqr   --   least-squares solution of  a*x = b\n\n", "");
		fprintf(fp_nout, "%25sthe matrix  a  has%8d rows   and%8d cols\n", "", m, n);
		fprintf(fp_nout, "%25sthe damping parameter is    damp   =%10.2E\n\n", "", damp);
		fprintf(fp_nout, "%25satol   =%10.2E          conlim =%10.2E\n", "", atoL, conlim);
		fprintf(fp_nout, "%25sbtol   =%10.2E          itnlim =%10d\n", "", btol, itnlim);
	} else {
		printf("%25slsqr   --   least-squares solution of  a*x = b\n\n", "");
		printf("%25sthe matrix  a  has%8d rows   and%8d cols\n", "", m, n);
		printf("%25sthe damping parameter is    damp   =%10.2E\n\n", "", damp);
		printf("%25satol   =%10.2E          conlim =%10.2E\n", "", atoL, conlim);
		printf("%25sbtol   =%10.2E          itnlim =%10d\n\n\n", "", btol, itnlim);
	}
	fprintf(fp_out, "%25slsqr   --   least-squares solution of  a*x = b\n\n", "");
	fprintf(fp_out, "%25sthe matrix  a  has %8d rows   and %8d cols\n", "", m, n);
	fprintf(fp_out, "%25sthe damping parameter is    damp   =%10.2E\n\n", "", damp);
	fprintf(fp_out, "%25satol   =%10.2E          conlim =%10.2E\n", "", atoL, conlim);
	fprintf(fp_out, "%25sbtol   =%10.2E          itnlim =%10d\n\n\n", "", btol, itnlim);


	float ctol = 0;
	float one = 1.0f;
	if (conlim > 0) {
		ctol = 1 / conlim;
	}
	float dampsq = damp * damp;
	anorm = 0;
	acond = 0;
	float bbnorm = 0;
	float ddnorm = 0;
	float res2 = 0;
	*xnorm = 0;
	float xxnorm = 0;
	float cs2 = -1;
	float sn2 = 0;
	float z = 0;
	int itn = 0;
	istop = 0;
	int nstop = 0;

	for (int i = 0; i < n; i++) {
		v[i] = 0;
		x[i] = 0;
		se[i] = 0;
	}

// c     set up the first vectors for the bidiagonalization.
// c     these satisfy   beta*u = b,   alfa*v = a(transpose)*u.
	float beta = 0;
	normlz(m, u, &beta);
	aprod(2, m, n, v, u, leniw, lenrw, iw, rw, a, na, ja);
	float alfa = 0;
	normlz(n, v, &alfa);
	scopy(n, v, 1, w, 1);

	float rhobar = alfa;
	float phibar = beta;
	float bnorm = beta;
	*rnorm = beta;
	*arnorm = alfa * beta;
	if (*arnorm <= 0) {
		goto a800;
	}
	if(fp_nout != NULL) {
		goto a100;
	}
	if (dampsq <= 0) {
		if (fp_nout) {
			fprintf(fp_nout, "   itn         x(1)           function       compatible incompatible norm(a) cond(a)\n\n");
		} else {
			printf("   itn         x(1)           function       compatible incompatible norm(a) cond(a)\n\n");

		}
	} else { // if (dampsq > 0) {
		if (fp_nout) {
			fprintf(fp_nout, "   itn         x(1)              function       compatible incompatible norm(abar) cond(abar)\n\n");
		} else {
			printf("   itn         x(1)              function       compatible incompatible norm(abar) cond(abar)\n\n");
		}
	}

	float test1 = one;
	float test2 = alfa / beta;
	if (fp_nout) {
		fprintf(fp_nout, "%6d    %.10E   %.10E    %.3E    %.3E\n\n", itn,x[0], *rnorm, test1, test2);
	} else {
		printf("%6d    %.10E   %.10E    %.3E    %.3E\n\n", itn,x[0], *rnorm, test1, test2);
	}
	
// c     ------------------------------------------------------------------
// c     main iteration loop.
// c     ------------------------------------------------------------------
	a100: itn++;
//      perform the next step of the bidiagonalization to obtain the
//      next  beta, u, alfa, v.  these satisfy the relations
//                 beta*u  =  a*v  -  alfa*u,
//                 alfa*v  =  a(transpose)*u  -  beta*v.
//
	sscal(m, (-alfa), u, 1);
	aprod(1, m, n, v, u, leniw, lenrw, iw, rw, a, na, ja);
	normlz(m, u, &beta);
	bbnorm = bbnorm + alfa * alfa + beta * beta + dampsq;
	sscal(n, (-beta), v, 1);
	aprod(2, m, n, v, u, leniw, lenrw, iw, rw, a, na, ja);
	normlz(n, v, &alfa);

//     use a plane rotation to eliminate the damping parameter.
//     this alters the diagonal (rhobar) of the lower-bidiagonal matrix.
	float rhbar2 = rhobar * rhobar + dampsq;
	float rhbar1 = sqrtf(rhbar2);
	float cs1 = rhobar / rhbar1;
	float sn1 = damp / rhbar1;
	float psi = sn1 * phibar;
	phibar = cs1 * phibar;

//      use a plane rotation to eliminate the subdiagonal element (beta)
//      of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

	float rho = sqrtf(rhbar2 + beta * beta);
	float cs = rhbar1 / rho;
	float sn = beta / rho;
	float theta = sn * alfa;
	rhobar = -cs * alfa;
	float phi = cs * phibar;
	phibar = sn * phibar;
	float tau = sn * phi;

//     update  x, w  and the standard error estimates.
	float t1 = phi / rho;
	float t2 = -theta / rho;
	float t3 = one / rho;
	for (int i = 0; i < n; i++) {
		float t = w[i];
		x[i] += t1 * t;
		w[i] = t2 * t + v[i];
		t = (t3 * t) * (t3 * t);
		se[i] += t;
		ddnorm += t;
	}

//      use a plane rotation on the right to eliminate the
//c     super-diagonal element (theta) of the upper-bidiagonal matrix.
//c     then use the result to estimate  norm(x).

	float delta = sn2 * rho;
	float gambar = -cs2 * rho;
	float rhs = phi - delta * z;
	float zbar = rhs / gambar;
	*xnorm = sqrtf(xxnorm + zbar * zbar);
	float gamma = sqrtf(gambar * gambar + theta * theta);
	cs2 = gambar / gamma;
	sn2 = theta / gamma;
	z = rhs / gamma;
	xxnorm = xxnorm + z * z;

//      test for convergence.
//      first, estimate the norm and condition of the matrix  abar,
//      and the norms of  rbar  and  abar(transpose)*rbar.

	anorm = sqrtf(bbnorm);
	acond = anorm * sqrtf(ddnorm);
	float res1 = phibar * phibar;
	res2 = res2 + psi * psi;
	*rnorm = sqrtf(res1 + res2);
	*arnorm = alfa * fabs(tau);
	
//      now use these norms to estimate certain other quantities,
//      some of which will be small near a solution.

	test1 = *rnorm / bnorm;
	test2 = *arnorm / (anorm * *rnorm);
	float test3 = one / acond;
	t1 = test1 / (one + anorm * *xnorm / bnorm);
	float rtol = btol + atoL * anorm * *xnorm / bnorm;

//      the following tests guard against extremely small values of
//      atoL, btol  or  ctol.  (the user may have set any or all of
//      the parameters  atoL, btol, conlim  to zero.)
//      the effect is equivalent to the normal tests using
//      atoL = relpr,  btol = relpr,  conlim = 1/relpr.

	t3 = one + test3;
	t2 = one + test2;
	t1 = one + t1;

	if (itn >= itnlim)
		istop = 7;
	if (t3 <= one)
		istop = 6;
	if (t2 <= one)
		istop = 5;
	if (t1 <= one)
		istop = 4;

//     allow for tolerances set by the user.
	
	if (test3 <= ctol)
		istop = 3;
	if (test2 <= atoL)
		istop = 2;
	if (test1 <= rtol)
		istop = 1;
//     ==================================================================
//     see if it is time to print something.
//     ==================================================================
	if ((m <= 40 || n <= 40) || (itn <= 10) || (itn >= itnlim - 10)
			|| (itn % 10 == 0) || (test3 <= 2.0 * ctol)
			|| (test2 <= 10.0 * atoL) || (test1 <= 10.0 * rtol)) {
		if(fp_nout){
			fprintf(fp_nout, "%6d    %.10E   %.10E    %.3E    %.3E   %.2E   %.2E\n", itn,x[0],*rnorm,test1,test2, anorm, acond);
		}else {
			printf("%6d    %.10E   %.10E    %.3E    %.3E   %.2E   %.2E\n", itn,x[0],*rnorm,test1,test2, anorm, acond);
		}
	}
//     print a line for this iteration.
//       ==================================================================
//
//      stop if appropriate.
//      the convergence criteria are required to be met on  nconv
//      consecutive iterations, where  nconv  is set below.
//      suggested value --   nconv = 1, 2  or  3.

	if (istop == 0)
		nstop = 0;
	if (istop == 0)
		goto a100;
	int nconv = 1;
	nstop++;
	if (nstop < nconv && itn < itnlim)
		istop = 0;
	if (istop == 0)
		goto a100;
// c     ------------------------------------------------------------------
// c     end of iteration loop.
// c     ------------------------------------------------------------------

// c     finish off the standard error estimates.

	float t = one;
	if (m > n)
		t = m - n;
	if (dampsq > 0)
		t = m;

	t = *rnorm / sqrtf(t);
	for (int i = 0; i < n; i++) {
		se[i] = t * sqrtf(se[i]);
	}
	a800:
	printf("\n no. of iterations =%6d         stopping condition =%3d\n\n", itn, istop);
	switch(istop) {
		case 0:
			printf(" the exact solution is  x = 0.                      \n\n\n");
			break;
		case 1:
			 printf(" a*x - b  is small enough, given  atol=%f, btol=%f        \n\n\n", atoL, btol);
			break;
		case 2:
			printf(" the least-sqrs soln is good enough, given  atol=%f    \n\n\n", atoL);
			break;
		case 3:
			printf(" the estimate of  cond(abar)  has exceeded  conlim  \n\n\n");
			break;
		case 4:
			printf(" a*x - b   is small enough for this machine\n\n\n");
			break;
		case 5:
			printf(" the least-sqrs soln is good enough for this machine\n\n\n");
			break;
		case 6:
			printf(" cond(abar)   seems to be too large for this machine\n\n\n");
			break;
		case 7:
			printf(" the iteration limit has been reached               \n\n\n");
			break;
	}
	if(fp_nout) {
		switch(istop) {
			case 0:
				printf(" the exact solution is  x = 0.                      \n\n\n");
				break;
			case 1:
				 printf(" a*x - b  is small enough, given  atol=%f, btol=%f        \n\n\n", atoL, btol);
				break;
			case 2:
				printf(" the least-sqrs soln is good enough, given  atol=%f    \n\n\n", atoL);
				break;
			case 3:
				printf(" the estimate of  cond(abar)  has exceeded  conlim  \n\n\n");
				break;
			case 4:
				printf(" a*x - b   is small enough for this machine\n\n\n");
				break;
			case 5:
				printf(" the least-sqrs soln is good enough for this machine\n\n\n");
				break;
			case 6:
				printf(" cond(abar)   seems to be too large for this machine\n\n\n");
				break;
			case 7:
				printf(" the iteration limit has been reached               \n\n\n");
				break;
		}
	}
	return;
}
