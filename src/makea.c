// ---------------------------------------------------------------------------
//
//   a is the sparse matrix whose nonzero coefficients are stored by
//   rows.  let a have m rows, n columns, and nz nonzeroes.  we need
//   three arrays dimensioned as real a(nz), and integer ja(nz),
//   na(m) where:
//
//    a(l) is the lth nonzero of a, counting across row 1, then
//     row 2, and so on.
//
//    ja(l) is the column in which the lth nonzero of a lies.
//
//    na(i) is the number of nonzero coefficients  in the ith row
//     of a.
//
//   (see lsqr alg. paige and saunders p.19)
//   these re made available to aprod through common.  the actual array
//   dimensions are:
//
//    a(SIZEOFA),ja(SIZEOFA), and na(m)
//
//   SIZEOFA is the probable largest number of nonzero coefficients in a:
//   this is dependent upon details of how a is formed
//
//   other input dimensional parameters are:
//
//    m - the row dimension of calling program.  If the number
//        of rows read in exceed this the program is aborted.
//    n - the number of columns in a.  The program is aborted if
//        any elements of ja() exceed n
//
//   on output:
//
//    mm - counts the number of rows processed, used in error check
//    l - counts the total number of nonzero coefficients processed
//        also used as an error check
//    ja - as described above
//    a - as described above
//    na - as described above
//    b - the rhs vector
//
//    lundts = velocity element matrix
//    lunres = data vector
// ----------------------------------------------------------------------------

#include "../include/runlsqr/makea.h"
#include "../include/parameter.h"
#include "../include/environment_setting.h"

#define TYPE_INT 4
#define TYPE_LONG 8
int get_size(FILE *);
int get_namm(int *, FILE *, FILE *);
int get_b(float *, int, FILE *);

void makea(int *m, int *nbl, FILE *fp_dts, FILE *fp_res, FILE *fp_out,
		int* jndx, float *a, int *na, int *ja, float *b) {
	int l = 0, mm = 0, namm = 0;
// c  read in velocity information, along with the data vector, first
	while (get_namm(&namm, fp_dts, fp_out)) {
		assert(get_size(fp_dts) == TYPE_INT);
		int header = 0, ender = 0;
		fread(&header, sizeof(header), 1, fp_dts);

		if(get_b(b, mm, fp_res) == EOF) {
			printf(" Error: Ran out of data! \n");
			fprintf(fp_out, " Error: Ran out of data! \n");
			assert(0);
		}

// c  test to check for overflow of array
		int l1 = l;
		int l2 = l + namm + 1;
		int l3 = l2;
		if (l3 > SIZEOFA) {
			printf(" FATAL ERROR (makea):  work array a is full\n");
			printf(" Increase size parameter SIZEOFA in source code\n");
			fprintf(fp_out, " FATAL ERROR (makea):  work array a is full\n");
			fprintf(fp_out,
					" Increase size parameter SIZEOFA in source code\n");
			assert(0);
		}

		for(int i = l1; i < l2 - 1; i++) {
			fread(ja + i, sizeof(ja[0]), 1, fp_dts);
			fread(a + i, sizeof(a[0]), 1, fp_dts);
			ja[i]--;
		}
		fread(&ender, sizeof(ender), 1, fp_dts);
		assert(header == ender);
		
		for (int i = l1; i < l2 - 1; i++) {
			if (ja[i] > NMAX) {
				printf(" FATAL ERROR (makea):  grid overflow \n");
				printf(" Datum number %d(row)\n", mm);
				fprintf(fp_out, " FATAL ERROR (makea):  grid overflow \n");
				fprintf(fp_out, " Datum number %d(row)\n", mm);
				assert(0);
			}
		}
		na[mm] = namm;
		l = l2 - 1;
		mm++;

// c  error exit for too many data values (too many rows)
		if (mm > MMAX) {
			printf(" FATAL ERROR (makea):  grid overflow \n");
			printf(" Maximum allowed number of rows = %d(rows)\n", MMAX);
			fprintf(fp_out, " FATAL ERROR (makea):  grid overflow \n");
			fprintf(fp_out, " Maximum allowed number of rows = %d(rows)\n", MMAX);
			assert(0);
		}
	}

// c------now read in the constraint matrix
	printf("%12d  rows of vel, hyp, and data read in \n", mm);
	fprintf(fp_out, "%12d  rows of vel, hyp, and data read in \n", mm);
// c---read in the index array

	fread(nbl, sizeof(*nbl), 1, fp_dts);
	fread(nbl, sizeof(*nbl), 1, fp_dts);
	fread(nbl, sizeof(*nbl), 1, fp_dts);
	get_size(fp_dts);
	int header;
	fread(&header, sizeof(header), 1, fp_dts);
	if (!fread(jndx, sizeof(jndx[0]), *nbl, fp_dts)) {
		printf("error on read jndx: fp_dts\n");
	}
	*m = mm;
	printf("  A total of %12d  rows read in \n", *m);
	printf("  A total of %12d  elements of A read in \n", l);
	printf("  Total number of variables = %14d \n", *nbl);
	fprintf(fp_out, "  A total of %12d  rows read in \n", *m);
	fprintf(fp_out, "  A total of %12d  elements of A read in \n", l);
	fprintf(fp_out, "  Total number of variables = %14d \n", *nbl);
}

int get_size(FILE *fptr) {
	int size;
	if(fread(&size, sizeof(size), 1, fptr) == EOF) {
		assert(0);
		return -1;
	}
	if(size != TYPE_INT) {
		printf("%d\n", size);
		assert(0);
	}
	return TYPE_INT;
}

int get_namm(int *namm, FILE *fptr, FILE *fp_out) {
// c----Abnormal Endings
	*namm = 0;
	assert(get_size(fptr) == TYPE_INT);
	if (!fread(namm, sizeof(*namm), 1, fptr)) {
		printf(" Error: Ran out of velocity info!\n");
		fprintf(fp_out, " Error: Ran out of velocity info!\n");
		assert(0);
	}
	return (*namm != -1);
}

int get_b(float *b, int mm, FILE *fptr) {
	int header;
	if(fread(&header, sizeof(header), 1, fptr) == EOF) {
		return EOF;
	}
	if (mm != 0) {
		if(fread(&header, sizeof(header), 1, fptr) == EOF) {
			return EOF;
		}
	}
	if(fread(&b[mm], sizeof(b[0]), 1, fptr) == EOF) {
		return EOF;
	}
	return 0;
}