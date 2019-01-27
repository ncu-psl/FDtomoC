/*
 c   These parameters set the general array dimensions among
 c   common arrays by the various tomofd programs.
 c
 c	Fine grid specs:
 c
 c	nxm 	Maximum number of grid poconst ints in
 c	nym	the x,y,z directions
 c	nzm
 c
 c	Coarse grid specs:
 c
 c	nxcm 	Maximum number of coarse grid poconst ints in
 c	nycm	the x,y,z directions
 c	nzcm
 c
 c	Note that the programs calculate the number of find grid poconst ints
 c	based on the coarse grid specs and node locations.
 c
 c	Bookkeeping:
 c
 c       maxsta	This is the total number of travel time tables that
 c		routines like telrayderv can store while doing ray
 c		calculations.  It will be more efficient to make this
 c		number as large as your machines memory will allow, but
 c		need not be any larger than the total number of tables
 c		you have ( = number of stations in the station list or
 c		twice that if you are using both P and S waves)
 c
 c	maxobs	The maximum number of observations allowed per event
 c
 c	maxlst	The maximum number of stations in the station list
 c
 c       maxnbk  The maximum number of variables accumulated by a single ray
 c
 c	maxkbl  The maximun number of variables accumulated by a single event
 c
 c	maxmbl  The maximun number of variables accumulated by the entire
 c		dataset
 c
 c----Large matrices of normal equations
 c
 c       NMAX 	Maximum number of columns of A (total number of variables)
 c
 c       MMAX	Maximum number of rows of A (observations + constraconst ints)
 c
 c       SIZEOFA Maximum number of elements of A
 c
 c--------------------------------------------------------------------
 */
//c----fine grid specs
//c ****** crustal scale (1km)
#define nxm 301
#define nym 301
#define nzm 61
#define nxym nxm*nym
#define nxyzm nxym*nzm
#define nxyzm2 nxyzm*2

//c---coarse grid specs
#define nxcm 203
#define nycm 203
#define nzcm 105
#define nxycm nxcm*nycm
#define nxyzcm nxycm*nzcm
#define nxyzcm2 nxyzcm*2
#define nxcm1 nxcm-1
#define nycm1 nycm-1
#define nzcm1 nzcm-1
#define nblkcm nxcm1*nycm1*nzcm1

//c----bookkeeping
#define maxsta 2000
#define maxobs 800
#define maxlst 4000
#define maxlst2 2*maxlst
#define maxnbk 200000
#define maxkbl 100000
#define maxmbl 3000000

//c----large matrices of normal equations
#define NMAX 2500000
#define MMAX 50000000
#define SIZEOFA 200000000
//c----openmp
#define ncpu 20
#define maxtrd 1024
