#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM 1000000
#define MAX_ERROR 10e-12

double norma_inf_vect(double *v);
void fillB(double *b);
void swap(double **a, double **b);
int jacobi(double *b, double *x_k, double *x_k_1, double *error_vect);
int gauss_seidel(double *b, double *x_k, double *x_k_1, double *error_vect);
/** BEGIN TEST **/

#include <assert.h>

#undef MAX_ERROR
#define MAX_ERROR 10e-16
#define TEST_ERROR 2


void write_vector(double *v, char *file)
{
	FILE *fp;
	int i;
	
	fp = fopen(file, "w");
	assert(fp != NULL);

	// print start of b
	for (i = 0; i < NUM; i++) {
		fprintf(fp, "%.12lf\n", v[i]);
	}
	fclose(fp);
}

void test_matrix_error(double *b, double *v, double *error_vect)
{
	int i;
	double divider;

	for (i = 0; i < NUM; i++) {
		divider = (double)(!(i % 2) ? 3 : 4);

		error_vect[i] = b[i] - divider * v[i];
		if (i + 2 < NUM)
			error_vect[i] += v[i+2];
		if (i - 2 < NUM)
			error_vect[i] += v[i-2];

		if (i == 0)
			error_vect[i] -= v[NUM - 2];
		if (i == 1)
			error_vect[i] -= v[NUM - 1];

		if (i == NUM - 2)
			error_vect[i] -= v[0];
		if (i == NUM - 1)
			error_vect[i] -= v[1];
	}
	assert(norma_inf_vect(error_vect) < TEST_ERROR);
}
/** END TEST **/

double norma_inf_vect(double *v) 
{
	int i = 0;
	double max = 0;

	for (i = 0; i < NUM; i++) {
		if (fabs(v[i]) > max) max = fabs(v[i]);
	}

	return max;
}

void fillB(double *b) 
{
	int i;
	for (i = 0; i < NUM; i++) {
		if (!(i % 2)) {
			b[i] = (double)(i+2) / (double)NUM;
			b[i+1] = (double)(i+2) / (double)NUM;
		}
	}
}

void swap(double **a, double **b)
{
	double *temp = *a;
	*a = *b;
	*b = temp;
}

int jacobi(double *b, double *x_k, double *x_k_1, double *error_vect)
{
	double divider;
	int i,k;

	// set the norm to 1 to avoid skipping any iterations
	error_vect[0] = 1;

	// iterations
	k = 0;
	while (2 * norma_inf_vect(error_vect) >= MAX_ERROR && ++k < 1000) {
		for (i = 0; i < NUM; i++) {
			divider = (double)(1 + i % 2 ? 3 : 4);

			x_k_1[i] = b[i];
			if (i - 2 >= 0) 
				x_k_1[i] += x_k[i-2];
			if (i + 2 < NUM)
				x_k_1[i] += x_k[i+2];

			if (i == 0)
				x_k_1[i] += -x_k[NUM - 2];
			if (i == 1)
				x_k_1[i] += -x_k[NUM - 1];

			if (i == NUM - 2)
				x_k_1[i] += -x_k[0];
			if (i == NUM - 1)
				x_k_1[i] += -x_k[1];

			x_k_1[i] /= divider;

			error_vect[i] = x_k_1[i] - x_k[i];
		}
		swap(&x_k, &x_k_1);
	}
	/** BEGIN TEST **/
	write_vector(x_k, "jacobi.txt");
	test_matrix_error(b, x_k, error_vect);
	/** END TEST **/

	return k;
}

int gauss_seidel(double *b, double *x_k, double *x_k_1, double *error_vect)
{	
	double divider;
	int i,k;

	// set the norm to 1 to avoid skipping any iterations
	error_vect[0] = 1;

	// actual iteration
	k = 0;
	while (2 * norma_inf_vect(error_vect) >= MAX_ERROR && ++k < 1000) {
		for (i = 0; i < NUM; i++) {
			divider = (double)(1 + i % 2 ? 3 : 4);

			x_k_1[i] = b[i];
			if (i - 2 >= 0) 
				x_k_1[i] += x_k_1[i-2];
			if (i + 2 < NUM)
				x_k_1[i] += x_k[i+2];
			
			if (i == 0)
				x_k_1[i] += -x_k[NUM - 2];
			if (i == 1)
				x_k_1[i] += -x_k[NUM - 1];

			if (i == NUM - 2)
				x_k_1[i] += -x_k_1[0];
			if (i == NUM - 1)
				x_k_1[i] += -x_k_1[1];

			x_k_1[i] /= divider;

			error_vect[i] = x_k_1[i] - x_k[i];
		}
		swap(&x_k, &x_k_1);
	}
	/** BEGIN TEST **/
	write_vector(x_k, "gs.txt");
	test_matrix_error(b, x_k, error_vect);
	/** END TEST **/

	return k;
}

int main()
{
	double *b = (double*)calloc(NUM, sizeof(double));
	double *x_k = (double*)calloc(NUM, sizeof(double));
	double *x_k_1 = (double*)calloc(NUM, sizeof(double));
	double *error_vect = (double*)calloc(NUM, sizeof(double));
	int k_jacobi;
	int k_gs;

	// create b vector
	fillB(b);

	k_jacobi = jacobi(b, x_k, x_k_1, error_vect);
	printf("Jacobi: %2d iterations\n", k_jacobi);

	free(x_k);
	free(x_k_1);
	x_k = (double*)calloc(NUM, sizeof(double));
	x_k_1 = (double*)calloc(NUM, sizeof(double));

	k_gs = gauss_seidel(b, x_k, x_k_1, error_vect);
	printf("Gauss-Seidel: %2d iterations\n", k_gs);

	free(b);
	free(x_k);
	free(x_k_1);
	free(error_vect);

	return 0;
}