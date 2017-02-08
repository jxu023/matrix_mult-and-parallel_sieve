a: main.c
	mpicc -o sieve3 sieve3.c
run:
	mpiexec -np 128 ./sieve3 10000000000 100000 #10^10 10^5
run124:
	mpiexec -np 32 ./main 10000000000 100000 #10^10 10^5
	mpiexec -np 64 ./main 10000000000 100000 #10^10 10^5
	mpiexec -np 128 ./main 10000000000 100000 #10^10 10^5
comp:
	mpicc -o sieve0 sieve0.c
	mpicc -o sieve1 sieve1.c
	mpicc -o sieve2 sieve2.c
	mpicc -o sieve3 sieve3.c
comp_run:
	mpiexec -np 4 ./sieve0 100000000 #10^8
	mpiexec -np 4 ./sieve1 100000000 #10^8
	mpiexec -np 4 ./sieve2 100000000 #10^8
	mpiexec -np 4 ./sieve3 100000000 100000 #10^8 10^5
