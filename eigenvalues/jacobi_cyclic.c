#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define NUM 100
#define TOL 1e-30
#define MAX_ITER 1000
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define sign(x) (((x)) ? (x/fabs(x)) : (1))

void fill_matrix(double **A) 
{
	for (int i = 0; i < NUM; i++) {
		A[i] = (double *)calloc(NUM, sizeof(double));

		// fill the matrix however you like
		for (int j = 0; j < NUM; j++) {
			A[i][j] = (double)((i+1) * (j+1));
			//if (i == j) A[i][j] = 4;
			//if (i == j + 1 || i == j - 1)
			//	A[i][j] = 1;
		}
	}

	/*A[0][0] = 2;
	A[0][1] = 2;
	A[0][2] = 1;
	A[1][0] = 2;
	A[1][1] = 3;
	A[1][2] = 1;
	A[2][0] = 1;
	A[2][1] = 1;
	A[2][2] = 2;*/
}

void print_matrix(double **A)
{
	for (int i = 0; i < NUM; i++) {
		if (i==0) printf("/ ");
		else if (i==NUM - 1) printf("\\ ");
		else printf("| ");

		for (int j = 0; j < NUM; j++) {
			printf("%5.1lf ", A[i][j]);
		}
		if (i==0) printf(" \\");
		else if (i==NUM - 1) printf(" /");
		else printf(" |");
		printf("\n");
	}
	printf("\n");
}

double frobenius_norm(double **A)
{
	double norm = 0;
	for (int i = 0; i < NUM; i++) {
		for (int j = 0; j < NUM; j++) {
			norm += A[i][j] * A[i][j];
		}
	}

	return sqrt(norm);
}

double jacobi_norm(double **A)
{
	double norm = 0;
	for (int i = 0; i < NUM; i++) {
		for (int j = 0; j < NUM; j++) {
			if (i != j)
				norm += A[i][j] * A[i][j];
		}
	}

	return sqrt(norm);
}

void do_rotation(double **A, int p, int q, double **B)
{
	double tau, sign, t, theta, c, s;

	if (A[p][q] != 0) {
		tau = (A[q][q] - A[p][p])/(2 * A[p][q]);
		t = sign(tau) / (fabs(tau) + sqrt(1 + tau * tau));
		theta = atan(t);
		c = cos(theta);
		s = sin(theta);
	} else {
		c = 1;
		s = 0;
	}

	for (int j = 0; j < NUM; j++) {
		if (j == p || j == q) continue;
		
		// printf("(%d,%d) <- %5.3lf * %5.3lf - %5.3lf * %5.3lf\n",
		//	p, j, A[p][j], c, A[q][j], s);
		B[p][j] = A[p][j] * c - A[q][j] * s;
		B[j][p] = B[p][j];

		// printf("(%d,%d) <- %5.3lf * %5.3lf + %5.3lf * %5.3lf\n",
		//	q, j, A[p][j], s, A[q][j], c);
		B[q][j] = A[p][j] * s + A[q][j] * c;
		B[j][q] = B[q][j];
	}

	// printf("(%d,%d) <- %5.3lf - %5.3lf * %5.3lf\n",
	//		p, p, A[p][p], A[p][q], t);
	// printf("(%d,%d) <- %5.3lf + %5.3lf * %5.3lf\n",
	//		q, q, A[q][q], A[p][q], t);
	B[p][p] = A[p][p] - A[p][q] * t;
	B[q][q] = A[q][q] + A[p][q] * t;

	// this should give 0
	B[p][q] = 0;
	B[q][p] = 0;
}

void matrix_copy(double **src, double **dest, int p, int q)
{
	for (int i = 0; i < NUM; i++) {
		dest[i][p] = src[i][p];
		dest[p][i] = src[p][i];
		dest[i][q] = src[i][q];
		dest[q][i] = src[q][i];
	}
}

int main(int argc, char const *argv[])
{
	double **A = (double **)malloc(NUM * sizeof(double *));
	double **B = (double **)malloc(NUM * sizeof(double *));
	double epsilon;
	double error = 100;
	int k = 0;

	fill_matrix(A);
	fill_matrix(B);
	//print_matrix(A);

	epsilon = TOL * jacobi_norm(A);

	while (jacobi_norm(A) >= epsilon && ++k < MAX_ITER) {	
		for (int i = 0; i < NUM; i++) {
			for (int j = i+1; j < NUM; j++) {
				if (fabs(B[i][j]) < TOL) continue;

				do_rotation(A, i, j, B);
				// print_matrix(B);
				matrix_copy(B, A, i, j);
				// printf("(%d, %d)\n", i, j);
				// print_matrix(A);
			}
		}
		printf("iteration %d, N(A) = %10.12lf\n", k, jacobi_norm(A));
	}

	printf("Finished in %d iterations, N(A) = %10.12lf\n", k, jacobi_norm(A));
	//print_matrix(B);
	printf("Eigenvalues:\n");
	for (int i = 0; i < NUM; i++) {
		printf("%+.12lf\n", B[i][i]);
		free(A[i]);
		free(B[i]);
	}
	printf("\n");

	free(A);
	free(B);

	return 0;
}