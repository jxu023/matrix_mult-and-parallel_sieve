#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

// set these global variables before calling functions
int n,B;
double *a,*b,*c;

// these vars get intialized ..
int i,j,k;
// produces random double between 0 and 1
double rand_double() {
    return ((double)rand())/((double)RAND_MAX);
}

void alg1_ijk() {
	for (i=0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			register double r = c[i*n+j];
			for (k=0; k < n; ++k) {
				r += (a[i*n+k] * b[k*n+j]);
			}
			c[i*n+j] = r;
		}
	}
}
void alg1_ikj() {
	for (i=0; i < n; ++i) {
		for (k = 0; k < n; ++k) {
			register double r = a[i*n+k];
			for (j=0; j < n; ++j) {
				c[i*n+j] += (r * b[k*n+j]);
			}
		}
	}
}
void alg1_jik() {
	for (j=0; j < n; ++j) {
		for (i = 0; i < n; ++i) {
			register double r = c[i*n+j];
			for (k=0; k < n; ++k) {
				r += (a[i*n+k] * b[k*n+j]);
			}
			c[i*n+j] = r;
		}
	}
}
void alg1_jki() {
	for (j=0; j < n; ++j) {
		for (k = 0; k < n; ++k) {
			register double r = b[k*n+j];
			for (i=0; i < n; ++i) {
				c[i*n+j] += (a[i*n+k] * r);
			}
		}
	}
}
void alg1_kij() {
	for (k = 0; k < n; ++k) {
		for (i=0; i < n; ++i) {
			register double r = a[i*n+k];
			for (j=0; j < n; ++j) {
				c[i*n+j] += (r * b[k*n+j]);
			}
		}
	}
}
void alg1_kji() {
	for (k = 0; k < n; ++k) {
		for (j=0; j < n; ++j) {
			register double r = b[k*n+j];
			for (i=0; i < n; ++i) {
				c[i*n+j] += (a[i*n+k] * r);
			}
		}
	}
}
int partial;
void do_remainder() {
	//take care of remainder	
	if (n % B == 0) {
		return;	
	}
	// these are all n^2
	// first do partially completed c
	for (i = 0; i < partial; ++i) {
		for (j = 0; j < partial; ++j) {
			register double r = c[i*n+j];
			for (k = partial; k < n; ++k) {
				r += a[i*n + k] * b[k*n + j];
			}
			c[i*n+j] = r;
		}
	}
	// then do c where still 0
		// corner
	for (i = partial; i < n; ++i) {
		for (j = partial; j < n; ++j) {
			register double r = c[i*n+j];
			for (k = 0; k < n; ++k) {
				r += a[i*n + k] * b[k*n + j];
			}
			c[i*n+j] = r;
		}
	}
		// side
	for (i = partial; i < n; ++i) {
		for (j = 0; j < partial; ++j) {
			register double r = c[i*n+j];
			for (k = 0; k < n; ++k) {
				r += a[i*n + k] * b[k*n + j];
			}
			c[i*n+j] = r;
		}
	}
		// side
	for (i = 0; i < partial; ++i) {
		for (j = partial; j < n; ++j) {
			register double r = c[i*n+j];
			for (k = 0; k < n; ++k) {
				r += a[i*n + k] * b[k*n + j];
			}
			c[i*n+j] = r;
		}
	}
}
void blockedMM_ijk() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (i = 0; i < partial; i+=B)
		for (j = 0; j < partial; j+=B)
			for (k = 0; k < partial; k+=B)
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < i+B; i1++)
					for (j1 = j; j1 < j+B; j1++)
					{
						register double r=c[i1*n+j1];
						for (k1 = k; k1 < k+B; k1++)
							r += a[i1*n + k1]*b[k1*n + j1];
						c[i1*n+j1]=r;
					}
	do_remainder();
}
void blockedMM_ikj() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (i = 0; i < partial; i+=B)
		for (k = 0; k < partial; k+=B)
			for (j = 0; j < partial; j+=B)
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < i+B; i1++)
					for (k1 = k; k1 < k+B; k1++)
					{
						register double r=a[i1*n+k1];
						for (j1 = j; j1 < j+B; j1++)
							c[i1*n+j1] += r*b[k1*n + j1];
					}
	do_remainder();
}
void blockedMM_jik() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (j = 0; j < partial; j+=B)
		for (i = 0; i < partial; i+=B)
			for (k = 0; k < partial; k+=B)
				/* B x B mini matrix multiplications */
				for (j1 = j; j1 < j+B; j1++)
					for (i1 = i; i1 < i+B; i1++)
					{
						register double r=c[i1*n+j1];
						for (k1 = k; k1 < k+B; k1++)
							r += a[i1*n + k1]*b[k1*n + j1];
						c[i1*n+j1]=r;
					}
	do_remainder();
}
void blockedMM_jki() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (j = 0; j < partial; j+=B)
		for (k = 0; k < partial; k+=B)
			for (i = 0; i < partial; i+=B)
				/* B x B mini matrix multiplications */
				for (j1 = j; j1 < j+B; j1++)
					for (k1 = k; k1 < k+B; k1++)
					{
						register double r=b[k1*n+j1];
						for (i1 = i; i1 < i+B; i1++)
							c[i1*n+j1] += a[i1*n + k1]*r;
					}
	do_remainder();
}
void blockedMM_kij() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (k = 0; k < partial; k+=B)
		for (i = 0; i < partial; i+=B)
			for (j = 0; j < partial; j+=B)
				/* B x B mini matrix multiplications */
				for (k1 = k; k1 < k+B; k1++)
					for (i1 = i; i1 < i+B; i1++)
					{
						register double r=a[i1*n+k1];
						for (j1 = j; j1 < j+B; j1++)
							c[i1*n+j1] += r*b[k1*n + j1];
					}
	do_remainder();
}
void blockedMM_kji() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (k = 0; k < partial; k+=B)
		for (j = 0; j < partial; j+=B)
			for (i = 0; i < partial; i+=B)
				/* B x B mini matrix multiplications */
				for (k1 = k; k1 < k+B; k1++)
					for (j1 = j; j1 < j+B; j1++)
					{
						register double r=b[k1*n+j1];
						for (i1 = i; i1 < i+B; i1++)
							c[i1*n+j1] += a[i1*n + k1]*r;
					}
	do_remainder();
}
double * reset_c() {
	for (i = 0; i < n*n; ++i) {
		c[i] = 0;
	}
}
double * init() {
	double * arr = (double *)malloc(n*n*sizeof(double));
	for (i = 0; i < n*n; ++i) {
		arr[i] = rand_double();
	}
	return arr;
}
// code taken from ilearn
void max_diff(double *x, double *y) {
	double diff = fabs(y[0] - x[0]);
	double maxA = fabs(a[0]);
	double maxB = fabs(b[0]);
	for (i = 0; i < n*n;++i) {
		if (fabs(y[i] - x[i]) > diff)
			diff = fabs(y[i]-x[i]);
		if (fabs(a[i] > maxA))
			maxA = fabs(a[i]);
		if (fabs(b[i]) > maxB)
			maxB = fabs(b[i]);
	}
	printf ("maximum difference: %f\n", diff/(maxA*maxB));
}
void badMM() {
	// does not even depend on A or B
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			c[i*n+j] = 13;
		}
	}
}
void timed_run(void (*fn)(),const char * name) {
	reset_c();
	uint64_t t0;
	struct timespec begin, end;
	clock_gettime(CLOCK_MONOTONIC,&begin);
	(*fn)();
	clock_gettime(CLOCK_MONOTONIC,&end);
	t0 = 1000000000L * (end.tv_sec - begin.tv_sec) + end.tv_nsec - begin.tv_nsec;
	printf ("elasped time of %s: %llu nanoseconds\n", name, (long long unsigned int) t0);
	printf("%llu seconds\n",t0/1000000000);
}
// used for testing ...
// i found out if wrong values were more than or less than the correct ones using this
// while testing on low values of n like 4 or 6
void print(double * c) {
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			printf("%f ",c[i*n+j]);
		}
		printf("\n");
	}
	printf("\n");
}
void dgemm3_kij() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (k = 0; k < partial; k+=B)
		for (i = 0; i < partial; i+=B)
			for (j = 0; j < partial; j+=B)
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < i+B; i1+=3)
					for (j1 = j; j1 < j+B; j1+=3)
					{
						register int t = i1*n+j1; register int tt = t+n; register int ttt = t+n+n;
						register double c00 = c[t]; register double c01 = c[t+1]; register double c02 = c[t+2];
						register double c10 = c[tt]; register double c11 = c[tt+1]; register double c12 = c[tt+2];
						register double c20 = c[ttt]; register double c21 = c[ttt+1]; register double c22 = c[ttt+2];
						for (k1 = k; k1 < k+B; k1+=3)
						{
							/* 3 by 3 mini matrix multiplication using registers*/
							register int ta = i1*n+k1; register int tta = ta+n; register int ttta = ta+n+n;
							register int tb = k1*n+j1; register int ttb = tb+n; register int tttb = tb+n+n;
							register double a00 = a[ta]; register double b00 = b[tb];
							register double a10 = a[tta]; register double b01 = b[tb+1];
							register double a20 = a[ttta]; register double b02 = b[tb+2];

							// advantage is given in that it allows the compiler to reorder the instruction
							// somewhere else
							register double extra = a[ta+1];

							c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;

							/*a00 = a[ta+1];*/  b00 = b[ttb];
							a10 = a[tta+1];  b01 = b[ttb+1];
							a20 = a[ttta+1];  b02 = b[ttb+2];

							c00 += extra*b00; c01 += extra*b01; c02 += extra*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;

							a00 = a[ta+2]; b00 = b[tttb];
							a10 = a[tta+2]; b01 = b[tttb+1];
							a20 = a[ttta+2]; b02 = b[tttb+2];

							c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;
						}
						c[t] = c00; c[t+1] = c01; c[t+2] = c02;
						c[tt] = c10; c[tt+1] = c11; c[tt+2] = c12;
						c[ttt] = c20; c[ttt+1] = c21; c[ttt+2] = c22;
					}
	do_remainder();
}
void dgemm3() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (i = 0; i < partial; i+=B)
		for (j = 0; j < partial; j+=B)
			for (k = 0; k < partial; k+=B)
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < i+B; i1+=3)
					for (j1 = j; j1 < j+B; j1+=3)
					{
						register int t = i1*n+j1; register int tt = t+n; register int ttt = t+n+n;
						register double c00 = c[t]; register double c01 = c[t+1]; register double c02 = c[t+2];
						register double c10 = c[tt]; register double c11 = c[tt+1]; register double c12 = c[tt+2];
						register double c20 = c[ttt]; register double c21 = c[ttt+1]; register double c22 = c[ttt+2];
						for (k1 = k; k1 < k+B; k1+=3)
						{
							/* 3 by 3 mini matrix multiplication using registers*/
							register int ta = i1*n+k1; register int tta = ta+n; register int ttta = ta+n+n;
							register int tb = k1*n+j1; register int ttb = tb+n; register int tttb = tb+n+n;
							register double a00 = a[ta]; register double b00 = b[tb];
							register double a10 = a[tta]; register double b01 = b[tb+1];
							register double a20 = a[ttta]; register double b02 = b[tb+2];

							// advantage is given in that it allows the compiler to reorder the instruction
							// somewhere else
							register double extra = a[ta+1];

							c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;

							/*a00 = a[ta+1];*/  b00 = b[ttb];
							a10 = a[tta+1];  b01 = b[ttb+1];
							a20 = a[ttta+1];  b02 = b[ttb+2];

							c00 += extra*b00; c01 += extra*b01; c02 += extra*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;

							a00 = a[ta+2]; b00 = b[tttb];
							a10 = a[tta+2]; b01 = b[tttb+1];
							a20 = a[ttta+2]; b02 = b[tttb+2];

							c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;
						}
						c[t] = c00; c[t+1] = c01; c[t+2] = c02;
						c[tt] = c10; c[tt+1] = c11; c[tt+2] = c12;
						c[ttt] = c20; c[ttt+1] = c21; c[ttt+2] = c22;
					}
	do_remainder();
}


void dgemm2() {
	int i1,j1,k1;
	partial = n - (n % B);
	for (i = 0; i < partial; i+=B)
		for (j = 0; j < partial; j+=B)
			for (k = 0; k < partial; k+=B) {
				int partiali = i + B - (B % 3);
				int partialj = j + B - (B % 3);
				int partialk = k + B - (B % 3);
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < partiali; i1+=3)
					for (j1 = j; j1 < partialj; j1+=3)
					{
						register int t = i1*n+j1; register int tt = t+n; register int ttt = t+n+n;
						register double c00 = c[t]; register double c01 = c[t+1]; register double c02 = c[t+2];
						register double c10 = c[tt]; register double c11 = c[tt+1]; register double c12 = c[tt+2];
						register double c20 = c[ttt]; register double c21 = c[ttt+1]; register double c22 = c[ttt+2];
						for (k1 = k; k1 < partialk; k1+=3)
						{
							/* 3 by 3 mini matrix multiplication using registers*/
							register int ta = i1*n+k1; register int tta = ta+n; register int ttta = ta+n+n;
							register int tb = k1*n+j1; register int ttb = tb+n; register int tttb = tb+n+n;
							register double a00 = a[ta]; register double b00 = b[tb];
							register double a10 = a[tta]; register double b01 = b[tb+1];
							register double a20 = a[ttta]; register double b02 = b[tb+2];

							// advantage is given in that it allows the compiler to reorder the instruction
							// somewhere else
							register double extra = a[ta+1];

							c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;

							/*a00 = a[ta+1];*/  b00 = b[ttb];
							a10 = a[tta+1];  b01 = b[ttb+1];
							a20 = a[ttta+1];  b02 = b[ttb+2];

							c00 += extra*b00; c01 += extra*b01; c02 += extra*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;

							a00 = a[ta+2]; b00 = b[tttb];
							a10 = a[tta+2]; b01 = b[tttb+1];
							a20 = a[ttta+2]; b02 = b[tttb+2];

							c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
							c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
							c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;
						}
						c[t] = c00; c[t+1] = c01; c[t+2] = c02;
						c[tt] = c10; c[tt+1] = c11; c[tt+2] = c12;
						c[ttt] = c20; c[ttt+1] = c21; c[ttt+2] = c22;
					}
				// copy remainder func with fixed indices here and test..
				//take care of remainder	
				if (B % 3 != 0) {
					// these are all n^2
					// first do partially completed c
					for (i1 = i; i1 < partiali; ++i1) {
						for (j1 = j; j1 < partialj; ++j1) {
							register double r = c[i1*n+j1];
							for (k1 = partialk; k1 < k+B; ++k1) {
								r += a[i1*n + k1] * b[k1*n + j1];
							}
							c[i1*n+j1] = r;
						}
					}
					// then do c where sti1ll 0
					// corner
					for (i1 = partiali; i1 < i+B; ++i1) {
						for (j1 = partialj; j1 < j+B; ++j1) {
							register double r = c[i1*n+j1];
							for (k1 = k; k1 < k+B; ++k1) {
								r += a[i1*n + k1] * b[k1*n + j1];
							}
							c[i1*n+j1] = r;
						}
					}
					// si1de
					for (i1 = partiali; i1 < i+B; ++i1) {
						for (j1 = j; j1 < partialj; ++j1) {
							register double r = c[i1*n+j1];
							for (k1 = k; k1 < k+B; ++k1) {
								r += a[i1*n + k1] * b[k1*n + j1];
							}
							c[i1*n+j1] = r;
						}
					}
					// si1de
					for (i1 = i; i1 < partiali; ++i1) {
						for (j1 = partialj; j1 < j+B; ++j1) {
							register double r = c[i1*n+j1];
							for (k1 = k; k1 < k+B; ++k1) {
								r += a[i1*n + k1] * b[k1*n + j1];
							}
							c[i1*n+j1] = r;
						}
					}
				}
			}	
	do_remainder();
}
int main()
{
	srand(13);
	n = 2048;
	a = init();
	b = init();
	double * c1 = init();
	double * c2 = init();

	c = c1;
	timed_run(alg1_kij,"alg1_kij");
	c = c2;
	B = 8;
	timed_run(blockedMM_kij,"blockedMM_kij");
	B = 28;
	timed_run(dgemm2,"dgemm2");
	max_diff(c1,c2);
	B = 33;
	timed_run(dgemm3,"dgemm3");
	max_diff(c1,c2);
	timed_run(dgemm3_kij,"dgemm3_kij");
	max_diff(c1,c2);

	free(a);
	free(b);
	free(c1);
	free(c2);
	return 0;
}
