#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

double *a,*b,*c3;
int n;
// produces random double between 0 and 1
double rand_double() {
    return ((double)rand())/((double)RAND_MAX);
}

// maximize use of 16 double registers
// ignore use of int registers
// dgemm2 used 12 double registers
// could reduce registers used by setting a00 = a01 after a00 is of no use as in the next slide
// uses 15 registers
// note : only works for n >= 3
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

					c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
					c10 += a10*b00; c11 += a10*b01; c12 += a10*b02;
					c20 += a20*b00; c21 += a20*b01; c22 += a20*b02;

					a00 = a[ta+1];  b00 = b[ttb];
					a10 = a[tta+1];  b01 = b[ttb+1];
					a20 = a[ttta+1];  b02 = b[ttb+2];

					c00 += a00*b00; c01 += a00*b01; c02 += a00*b02;
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
void run_part3(int N) {
	
	printf("\nN is %d\n",N);
	n = N;
	uint64_t t0;
	struct timespec begin, end;

	clock_gettime(CLOCK_MONOTONIC,&begin);
    dgemm3(c3);
	clock_gettime(CLOCK_MONOTONIC,&end);
	t0 = 1000000000L * (end.tv_sec - begin.tv_sec) + end.tv_nsec - begin.tv_nsec;
	printf ("elasped time of dgemm3: %llu nanoseconds\n", (long long unsigned int) t0);
	printf("%llu seconds\n",t0/1000000000);
}
int main()
{
    srand(13);
    n = 2048;
    init(a);
    init(b);
    c3 = (double *) calloc(n*n,sizeof(double));
	int N[6] = {64,128,256,512,1024,2048};
	for (int i = 0; i < 6; ++i) {
		run_part3(N[i]);
	}
    free(a);
    free(b);
    free(c3);
	return 0;
}
