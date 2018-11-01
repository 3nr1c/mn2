#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NUM 10
#define MAX_ERROR 1e-14

void fill_matrix(double **A) 
{
	for (int i = 0; i < NUM; i++) {
		A[i] = (double *)malloc(NUM * sizeof(double));

		// fill the matrix however you like
		for (int j = 0; j < NUM; j++) {
			A[i][j] = (double)(i + 2*j + i*i);
		}
	}
}

void print_matrix(double **A)
{
	printf("Matrix A\n");

	for (int i = 0; i < NUM; i++) {
		for (int j = 0; j < NUM; j++) {
			printf("%5.1lf ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

double dot_product(double *a, double *b)
{
	double res = 0;

	for (int i = 0; i < NUM; i++) {
		res += a[i] * b[i];
	}

	return res;
}

void normalize(double *src, double *dest)
{
	double norm = sqrt(dot_product(src, src));

	for (int i = 0; i < NUM; i++) {
		dest[i] = src[i] / norm;
	}
}

void matrix_product(double **A, double *x, double *dest)
{
	for (int i = 0; i < NUM; i++) {
		dest[i] = 0;
		for (int j = 0; j < NUM; j++) {
			dest[i] += A[i][j] * x[j];
		}
	}
}

int main(int argc, char const *argv[])
{
	
	double **A = (double **)malloc(NUM * sizeof(double *));
	double *y = (double *)malloc(NUM * sizeof(double));
	double *z = (double *)malloc(NUM * sizeof(double));
	double eigenvalue = 1e3;
	double prev_eigenvalue = 1e2;
	double error = 1;
	int k = 0;
	char seed[4];
	FILE *fp;

	fp = fopen("/dev/urandom", "r");
	fread(&seed, 1, 4, fp);
	fclose(fp);
	printf("Using seed %d\n", (int)seed);

	fill_matrix(A);
	print_matrix(A);

	srandom((int)seed);

	// fill y and z
	for (int i = 0; i < NUM; i++) {
		z[i] = random();
		y[i] = random();
	}

	// until convergence, actually
	while (error >= MAX_ERROR && k++ < 100) {
		normalize(z, y);
		matrix_product(A, y, z);

		prev_eigenvalue = eigenvalue;
		eigenvalue = dot_product(y, z);

		error = fabs(prev_eigenvalue - eigenvalue);

		printf("Approximation #%3d: %.12lf\n", k, eigenvalue);
	}

	free(y);
	free(z);
	for (int i = 0; i < NUM; i++) {
		free(A[i]);
	}
	free(A);

	return 0;
}