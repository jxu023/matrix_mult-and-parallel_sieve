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
int main()
{
	srand(13);
	n = 1024;
	a = init();
	b = init();
	double * c1 = init();
	double * c2 = init();
	c = c1;
	timed_run(alg1_ijk,"alg1_ijk");
	c = c2;
	/*timed_run(alg1_ikj,"alg1_ikj");
	max_diff(c1,c2);
	timed_run(alg1_jik,"alg1_jik");
	max_diff(c1,c2);
	timed_run(alg1_jki,"alg1_jki");
	max_diff(c1,c2);
	timed_run(alg1_kij,"alg1_kij");
	max_diff(c1,c2);
	timed_run(alg1_kji,"alg1_kji");
	max_diff(c1,c2);*/
	B = 10;
	for (B = 5; B <= 15; ++B) {
		printf("B is : %d\n",B);
		timed_run(blockedMM_ijk,"blockedMM_ijk");
		max_diff(c1,c2);
		/*timed_run(blockedMM_ikj,"blockedMM_ikj");
		max_diff(c1,c2);
		timed_run(blockedMM_jik,"blockedMM_jik");
		max_diff(c1,c2);
		timed_run(blockedMM_jki,"blockedMM_jki");
		max_diff(c1,c2);
		timed_run(blockedMM_kij,"blockedMM_kij");
		max_diff(c1,c2);
		timed_run(blockedMM_kji,"blockedMM_kji");
		max_diff(c1,c2);*/
	}

	free(a);
	free(b);
	free(c1);
	free(c2);
	return 0;
}
