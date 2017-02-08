#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
	int    count;        /* Local prime count */
	double elapsed_time; /* Parallel execution time */
	int    first;        /* Index of first multiple */
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

	MPI_Init (&argc, &argv);

	/* Start the timer */

	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();

	if (argc != 2) {
		if (!id) printf ("Command line: %s <m>\n", argv[0]);
		MPI_Finalize();
		exit (1);
	}

	n = atol(argv[1]);

	/* Figure out this process's share of the array, as
	   well as the integers represented by the first and
	   last array elements */

	low_value = 2 + id*(n-1)/p;
	high_value = 1 + (id+1)*(n-1)/p;
	if (high_value % 2 == 1) {
		++high_value;
	}
	if (low_value % 2 == 1) {
		--low_value;
	}
	size = (high_value - low_value)/2;

	/* Bail out if all the primes used for sieving are
	   not all held by process 0 */

	proc0_size = (n-1)/p;

	if ((2 + proc0_size) < (int) sqrt((double) n)) {
		if (!id) printf ("Too many processes\n");
		MPI_Finalize();
		exit (1);
	}

	/* Allocate this process's share of the array. */
	marked = (char *) calloc(size,sizeof(char));
	if (marked == NULL) {
		printf ("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit (1);
	}

	// use sequential algorithm to get sieving primes
	// calculate sieving primes by borrowing the marked array,
	// then transfer primes into int array

	// note i guess process 0 will just reprocess these sieving primes .. o_o
	int sieve_sz = (int)floor(sqrt(n));
	marked[0] = 1;
	marked[1] = 1;
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
	for (i = 0; i < sieve_sz; ++i) {
		if (!marked[i]) {
			++count;
		}
	}
	int * primes = (int *)malloc(sieve_sz);
	int j = 0;
	for (i = 0; i < sieve_sz; ++i) {
		if (!marked[i]) {
			primes[j++] = i;
			//printf("%d ",i);
		}
	}
	for (i = 0; i < sieve_sz; ++i) {
		marked[i] = 0;
	}
	sieve_sz = count;
	if (!id)
		index = 0;
	j = 1;
	do {
		prime = primes[j++];
		// global index(gi) goes from low_value to high value
		// local index(li) goes from 0 to size
		// local index = gi - low_value

		// we set first to be the first thing to be marked
		if (prime * prime > low_value) {
			// everything below k^2 has been marked already
			// global index = k^2
			// local index = k^2 - low value
			first = prime * prime - low_value;
		}
		else {
			if (!(low_value % prime)) {
				// global index = low_value
				// local index = gi - low_value = 0
				first = 0;
				if (low_value % 2 == 0) {
					first+=prime;
				}
			}
			else {
				// add remainder to non-multiple to get multiple
				// there is (prime - low_value % prime) values until the next multiple
				// global index = low_value + prime - (low_value % prime)
				first = prime - (low_value % prime);
				if ((low_value + first) % 2 == 0) {
					first+=prime;
				}
			}
		}
		first/=2;

		// now mark the multiples
		for (i = first; i < size; i += prime) {
			marked[i] = 1;
		}
		//printf("\nlow_value:%d high_value:%d count:%d size:%d first:%d \n",low_value,high_value,count,size,first*2+3+low_value);
	} while (j < sieve_sz);
	count = 0;
	for (i = 0; i < size; i++)
		if (!marked[i]) {
			++count;
			/*if (!id) {
				printf("%d ", i*2+3+low_value-2);
			}
			else {
				if ((i*2+3+low_value-2) % 2 == 0) {
					printf("\nprime:%d low_value:%d high_value:%d count:%d size:%d first:%d \n",i*2+3 + low_value-2,low_value,high_value,count,size,first*2+3+low_value);
				}
			}*/
		}
	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
			0, MPI_COMM_WORLD);

	/* Stop the timer */

	elapsed_time += MPI_Wtime();

	/* Print the results */

	if (!id) {
		// global_count+1 to account for 2
		if (p == 1) {
			global_count = count;
		}
		printf ("There are %d primes less than %d\n",
				global_count+1, n);
		printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
	}
	MPI_Finalize ();
	return 0;
}
