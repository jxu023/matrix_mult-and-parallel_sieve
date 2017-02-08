#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

double *a,*b,*c1,*c2;
int n;
// produces random double between 0 and 1
double rand_double() {
    return ((double)rand())/((double)RAND_MAX);
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

				 c[t] = c00;
				 c[t+1] = c01;
				 c[tt] = c10;
				 c[tt+1] = c11;
			}
	}
}

void init(double *& arr) {
    arr = (double *)malloc(n*n*sizeof(double));
    for (int i = 0; i < n*n; ++i) {
        arr[i] = rand_double();
    }
}
void max_diff() {
	double diff = fabs(c2[0] - c1[0]);
	double maxA = fabs(a[0]);
	double maxB = fabs(b[0]);
	for (int i = 0; i < n*n;++i) {
		if (fabs(c2[i] - c1[i]) > diff)
			diff = fabs(c2[i]-c1[i]);
		if (fabs(a[i] > maxA))
			maxA = fabs(a[i]);
		if (fabs(b[i]) > maxB)
			maxB = fabs(b[i]);
	}
	printf ("maximum difference between dgemm1 and dgemm2: %f", diff/(maxA*maxB));
}
void run_part2(int N) {
	
	printf("\nN is %d\n",N);
	n = N;
	uint64_t t0;
	struct timespec begin, end;

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

	max_diff();
}
int main()
{
    srand(13);
    n = 2048;
    init(a);
    init(b);
    c1 = (double *) calloc(n*n,sizeof(double));
    c2 = (double *) calloc(n*n,sizeof(double));
	int N[6] = {64,128,256,512,1024/*,2048*/};
	for (int i = 0; i < 6; ++i) {
		run_part2(N[i]);
	}
    free(a);
    free(b);
    free(c1);
    free(c2);
	return 0;
}
