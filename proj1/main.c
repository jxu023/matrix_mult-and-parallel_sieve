#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

double *a,*b,*c2,*c3;
int n;
// produces random double between 0 and 1
double rand_double() {
    return ((double)rand())/((double)RAND_MAX);
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
void dgemm3(double * &c) {
	int i,j,k;
	for(i = 0; i < n; i += 3) {
		   for(j = 0; j < n; j += 3)  {
				register int t = i*n+j; register int tt = t+n; register int ttt = t+n+n;
				register double c00 = c[t]; register double c01 = c[t+1]; register double c02 = c[t+2];
				register double c10 = c[tt]; register double c11 = c[tt+1]; register double c12 = c[tt+2];
				register double c20 = c[ttt]; register double c21 = c[ttt+1]; register double c22 = c[ttt+2];

				for(k = 0; k < n; k += 3) {
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
}

void init(double *& arr) {
    arr = (double *)malloc(n*n*sizeof(double));
    for (int i = 0; i < n*n; ++i) {
        arr[i] = rand_double();
    }
}
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
void run_part3(int N) {
	
	printf("\nN is %d\n",N);
	n = N;
	uint64_t t0;
	struct timespec begin, end;

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
}
int main()
{
    srand(13);
    n = 2048;
    init(a);
    init(b);
    c3 = (double *) calloc(n*n,sizeof(double));
    c2 = (double *) calloc(n*n,sizeof(double));
	int N[6] = {64,128,256,512,1024/*,2048*/};
	//int N[6] = {36,72,144,288,576,576};
	for (int i = 0; i < 6; ++i) {
		run_part3(N[i]);
	}
    free(a);
    free(b);
    free(c2);
    free(c3);
	return 0;
}
