#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

double *a,*b,*c0,*c1,*c2,*c3;
int n;
// produces random double between 0 and 1
double rand_double() {
    return ((double)rand())/((double)RAND_MAX);
}

//code taken from given 
void dgemm0(double * &c) {
	int i,j,k;
	for (i=0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			for (k=0; k < n; ++k) {
				c[i*n+j] += (a[i*n+k] * b[k*n+j]);
			}
		}
	}
}
// code taken from given
void dgemm1(double * &c) {
	int i,j,k;
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
// code taken from slide
void dgemm2(double * &c) {
	int i,j,k;
	for(i = 0; i < n; i += 2) {
		   for(j = 0; j < n; j += 2)  {
				register int t = i*n+j; register int tt = t+n; 
				register double c00 = c[t]; register double c01 = c[t+1];  register double c10 = c[tt]; register double c11 = c[tt+1];

				for(k = 0; k < n; k += 2) {
					/* 2 by 2 mini matrix multiplication using registers*/
					register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
					register double a00 = a[ta]; register double a01 = a[ta+1]; register double a10 = a[tta]; register double a11 = a[tta+1];
					register double b00 = b[tb]; register double b01 = b[tb+1]; register double b10 = b[ttb]; register double b11 = b[ttb+1];
					c00 += a00*b00 + a01*b10;
					c01 += a00*b01 + a01*b11;
					c10 += a10*b00 + a11*b10;
					c11 += a10*b01 + a11*b11;
				 }

				 c[t] = c00; c[t+1] = c01;
				 c[tt] = c10; c[tt+1] = c11;
			}
	}
}
// maximize use of 16 double registers
// ignore use of int registers
// dgemm2 used 12 double registers
// could reduce registers used by setting a00 = a01 after a00 is of no use as in the next slide
// uses 15 registers
// note : only works for n >= 3
// time to edit to make use 16 registers ..
void dgemm3(double * &c) {
	int i,j,k;
	int r = n % 3;
	int partial = n - r;
	for(i = 0; i < partial; i += 3) {
		   for(j = 0; j < partial; j += 3)  {
				register int t = i*n+j; register int tt = t+n; register int ttt = t+n+n;
				register double c00 = c[t]; register double c01 = c[t+1]; register double c02 = c[t+2];
				register double c10 = c[tt]; register double c11 = c[tt+1]; register double c12 = c[tt+2];
				register double c20 = c[ttt]; register double c21 = c[ttt+1]; register double c22 = c[ttt+2];

				for(k = 0; k < partial; k += 3) {
					/* 3 by 3 mini matrix multiplication using registers*/
					register int ta = i*n+k; register int tta = ta+n; register int ttta = ta+n+n;
					register int tb = k*n+j; register int ttb = tb+n; register int tttb = tb+n+n;
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
	}
	if (n % 3 == 0) {
		return;	
	}
	// these are all n^2
	// first do partially completed c
	for (i = 0; i < partial; ++i) {
		for (j = 0; j < partial; ++j) {
			for (k = partial; k < n; ++k) {
				c[i*n + j] += a[i*n + k] * b[k*n + j];
			}
		}
	}
	// then do c where still 0
		// corner
	for (i = partial; i < n; ++i) {
		for (j = partial; j < n; ++j) {
			for (k = 0; k < n; ++k) {
				c[i*n + j] += a[i*n + k] * b[k*n + j];
			}
		}
	}
		// side
	for (i = partial; i < n; ++i) {
		for (j = 0; j < partial; ++j) {
			for (k = 0; k < n; ++k) {
				c[i*n + j] += a[i*n + k] * b[k*n + j];
			}
		}
	}
		// side
	for (i = 0; i < partial; ++i) {
		for (j = partial; j < n; ++j) {
			for (k = 0; k < n; ++k) {
				c[i*n + j] += a[i*n + k] * b[k*n + j];
			}
		}
	}
}

void init(double *& arr) {
    arr = (double *)malloc(n*n*sizeof(double));
    for (int i = 0; i < n*n; ++i) {
        arr[i] = rand_double();
    }
}
// code taken from ilearn
void max_diff(double *x, double *y) {
	double diff = fabs(y[0] - x[0]);
	double maxA = fabs(a[0]);
	double maxB = fabs(b[0]);
	for (int i = 0; i < n*n;++i) {
		if (fabs(y[i] - x[i]) > diff)
			diff = fabs(y[i]-x[i]);
		if (fabs(a[i] > maxA))
			maxA = fabs(a[i]);
		if (fabs(b[i]) > maxB)
			maxB = fabs(b[i]);
	}
	printf ("maximum difference: %f", diff/(maxA*maxB));
}
void run(int N) {     
	
	printf("\nN is %d\n",N);
	n = N;
	uint64_t t0;
	struct timespec begin, end;

	clock_gettime(CLOCK_MONOTONIC,&begin);
    dgemm0(c0);
	clock_gettime(CLOCK_MONOTONIC,&end);
	t0 = 1000000000L * (end.tv_sec - begin.tv_sec) + end.tv_nsec - begin.tv_nsec;
	printf ("elasped time of dgemm0: %llu nanoseconds\n", (long long unsigned int) t0);
	printf("%llu seconds\n",t0/1000000000);

	clock_gettime(CLOCK_MONOTONIC,&begin);
    dgemm1(c1);
	clock_gettime(CLOCK_MONOTONIC,&end);
	t0 = 1000000000L * (end.tv_sec - begin.tv_sec) + end.tv_nsec - begin.tv_nsec;
	printf ("elasped time of dgemm1: %llu nanoseconds\n", (long long unsigned int) t0);
	printf("%llu seconds\n",t0/1000000000);

	clock_gettime(CLOCK_MONOTONIC,&begin);
    dgemm2(c2);
	clock_gettime(CLOCK_MONOTONIC,&end);
	t0 = 1000000000L * (end.tv_sec - begin.tv_sec) + end.tv_nsec - begin.tv_nsec;
	printf ("elasped time of dgemm2: %llu nanoseconds\n", (long long unsigned int) t0);
	printf("%llu seconds\n",t0/1000000000);

	clock_gettime(CLOCK_MONOTONIC,&begin);
    dgemm3(c3);
	clock_gettime(CLOCK_MONOTONIC,&end);
	t0 = 1000000000L * (end.tv_sec - begin.tv_sec) + end.tv_nsec - begin.tv_nsec;
	printf ("elasped time of dgemm3: %llu nanoseconds\n", (long long unsigned int) t0);
	printf("%llu seconds\n",t0/1000000000);


	max_diff(c3,c2);
	max_diff(c3,c1);
}
// used for testing ...
// i found out if wrong values were more than or less than the correct ones using this
// while testing on low values of n like 4 or 6
void print(double * c) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%f ",c[i*n+j]);
		}
		printf("\n");
	}
	printf("\n");
}
int main()
{
    srand(13);
    n = 2048;
    init(a);
    init(b);
    c3 = (double *) calloc(n*n,sizeof(double));
    c2 = (double *) calloc(n*n,sizeof(double));
	c1 = (double *) calloc(n*n,sizeof(double));
	c0 = (double *) calloc(n*n,sizeof(double));
	int N[6] = {10,10,10,10,10,10};
	for (int i = 0; i < 6; ++i) {
		run(N[i]);
	}
    free(a);
    free(b);
	free(c0);
    free(c1);
    free(c2);
    free(c3);
	return 0;
}
