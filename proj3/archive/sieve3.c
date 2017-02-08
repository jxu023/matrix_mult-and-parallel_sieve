#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(A,B) ((A < B) ? A : B)
#define MAX(A,B) ((A > B) ? A : B)

int main (int argc, char *argv[])
{
	int    count;        /* Local prime count */
	double elapsed_time; /* Parallel execution time */
	int    global_count; /* Global prime count */
	long    high_value;   /* Highest value on this proc */
	int    i;
	int    id;           /* Process ID number */
	int    index;        /* Index of current prime */
	long    low_value;    /* Lowest value on this proc */
	char  *marked;       /* Portion of 2,...,'n' */
	long    n;            /* Sieving from 2, ..., 'n' */
	int    p;            /* Number of processes */
	int    proc0_size;   /* Size of proc 0's subarray */
	int    prime;        /* Current prime */
	int    size;         /* Elements in 'marked' */
	int	   first;

	MPI_Init (&argc, &argv);

	/* Start the timer */

	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();

	if (argc != 3) {
		if (!id) printf ("Command line: %s <m>\n", argv[0]);
		MPI_Finalize();
		exit (1);
	}

	n = atol(argv[1]);
	int block_size = atoi(argv[2]);
	long Llow_value = 2 + (id)*(n-1)/p;
	long Hhigh_value = 1 + (id+1)*(n-1)/p;
	//long Llow_value = 2 + (id)*((n-1)/p);
	//long Hhigh_value = 1 + (id+1)*((n-1)/p);
	int sieve_sz = (int)ceil(sqrt(n));
	proc0_size = (n-1)/p;

	if ((2 + proc0_size) < (int) sqrt((double) n)) {
		if (!id) printf ("Too many processes\n");
		MPI_Finalize();
		exit (1);
	}
	marked = (char *) calloc(MAX(block_size,sieve_sz),sizeof(char));
	if (marked == NULL) {
		printf ("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit (1);
	}

	// use sequential algorithm to get sieving primes
	// could put both sieving primes, and blocked numbers in same array ... improves locality?
	int k = 2;
	while (k*k < sieve_sz) {
		for (i = k*k; i < sieve_sz; i+=k) {
			marked[i] = 1; 
		}
		i = k;
		while(marked[++i]);
		k = i;
	}
	count = 0;
	int * primes = (int *)malloc(sieve_sz);
	int j = 0;
	for (i = 2; i < sieve_sz; ++i) {
		if (!marked[i]) {
			++count;
			primes[j++] = i;
		}
	}
	for (i = 0; i < sieve_sz; ++i) {
		marked[i] = 0;
	}
	sieve_sz = count;
	count = 0;

	int b = 0; // current block
	high_value = 0;
	while (high_value < Hhigh_value) {
		low_value = Llow_value + b*block_size*2;
		high_value = MIN(low_value+block_size*2-1,Hhigh_value);
	//	printf("low_value:%d\n",low_value);
	//	printf("high_value:%d\n",high_value);
		if (high_value % 2 == 1) {
			++high_value;
		}
		if (low_value % 2 == 1) {
			--low_value;
		}
		size = (high_value - low_value)/2;
		int j = 1;
		do {
			prime = primes[j++];
			if (prime * prime > low_value) {
				first = prime * prime - low_value;
			}
			else {
				if (!(low_value % prime)) {
					first = 0;
					if (low_value % 2 == 0) {
						first+=prime;
					}
				}
				else {
					first = prime - (low_value % prime);
					if ((low_value + first) % 2 == 0) {
						first+=prime;
					}
				}
			}
			// this commented piece of code helped me find a bug ...
			// this bug only showed up when running n=10^10 .. tricky tricky..
			/*if (first > high_value) {
				printf("first:%d\n",first);
				printf("size:%d\n",size);
			}*/
			first/=2;

			// now mark the multiples
			for (i = first; i < size; i += prime) {
				marked[i] = 1;
			}
		} while (j < sieve_sz && prime < (int)floor(sqrt(high_value)));
		for (i = 0; i < size; i++)
			if (!marked[i]) {
				++count;
			}
		for (i = 0; i < size; ++i) {
			marked[i] = 0;
		}
		++b;
	}

	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
			0, MPI_COMM_WORLD);
	/*printf("\n%lld\t%lld\t%d\n",Llow_value,Hhigh_value,count);
	fflush(stdout);*/
	/* Stop the timer */

	elapsed_time += MPI_Wtime();

	/* Print the results */

	if (!id) {
		// global_count+1 to account for 2
		if (p == 1) {
			global_count = count;
		}
		printf ("There are %d primes less than %lld\n",
				global_count+1, n);
		printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
	}
	MPI_Finalize ();
	return 0;
}
