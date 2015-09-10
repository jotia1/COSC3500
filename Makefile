CC=gcc
CFLAGS=-I. -std=c99 

omp: paraOMP.c
	icc -o OMPrun paraOMP.c -pg -O3 $(CFLAGS)

clean: 
	rm OMPrun
