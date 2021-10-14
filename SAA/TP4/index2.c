// vectorisation-loops-check.c

/* 
README : vectorise a non-vectorisable loop via a temporary array with re-indexing

compile/run as gcc -O2 -fopt-info-vec -fopt-info-vec-missed vectorisation-reindexation.c  && time ./a.out  
 */

#include <stdio.h>
#include <assert.h>


int main()
{
	const int dim=100000;
	double A[dim];
	double S=0;
	double X[dim];

	// initialise array
	for (int i=0; i<dim; i++)
		X[i] = (i+1)*1./(i+1) +1./(i+1)*1./(dim-i);

	// calculate required sum
	for (int i=0; i<dim; i++)
		S += X[i];

  	printf("S: %lf",S); // dummy print
}
