#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM 1000000
#define MAX_ERROR 10e-12

double norma_inf_vect(double *v) {
	int i = 0;

	double max = 0;

	for (i = 0; i < NUM; i++) {
		if (fabs(v[i]) > max) max = fabs(v[i]);
	}

	return max;
}

void fillB_J(double *b) 
{
	int i;
	for (i = 0; i < NUM; i++) {
		if (!(i % 2)) {
			b[i] = (double)(i+2) / (double)(NUM * (!(i % 2) ? 3 : 4));
			b[i+1] = (double)(i+2) / (double)(NUM * (i % 2 ? 3 : 4));
		}
	}
}

void swap(double **a, double **b)
{
	double *temp = *a;
	*a = *b;
	*b = temp;
}

int main()
{
	// init memory
	double *x_k = calloc(NUM, sizeof(double));
	double *x_k_1 = malloc(NUM * sizeof(double));
	double *b = malloc(NUM * sizeof(double));
	double *error_vect = calloc(NUM, sizeof(double));
	double divider;
	int i,k;

	// printf("Memory allocated\n");

	// create b vector
	fillB_J(b);

	// set the norm to 1 to avoid skipping any iterations
	error_vect[0] = 1;

	// printf("Vector b initialized\n");

	// actual iteration
	k = 0;
	while (2 * norma_inf_vect(error_vect) >= MAX_ERROR && ++k) {
		// printf("Starting iteration %d\n", k+1);

		for (i = 0; i < NUM; i++) {
			divider = (double)(!(i % 2) ? 3 : 4);

			x_k_1[i] = b[i];
			if (i - 2 >= 0) 
				x_k_1[i] += x_k[i-2] / divider;
			if (i + 2 < NUM)
				x_k_1[i] += x_k[i+2] / divider;
			if (i == NUM - 2)
				x_k_1[i] += -x_k[0] / divider;
			if (i == NUM - 1)
				x_k_1[i] += -x_k[1] / divider;
			if (i == 0)
				x_k_1[i] += -x_k[NUM - 2] / divider;
			if (i == 1)
				x_k_1[i] += -x_k[NUM - 1] / divider;

			error_vect[i] = x_k_1[i] - x_k[i];
		}
		//printf("Finished iteration %d with max error %.12lf\n", k, 2 * norma_inf_vect(error_vect));

		swap(&x_k, &x_k_1);
	}

	// print start of b
	for (i = 0; i < NUM; i++) {
		printf("%.12lf\n", x_k[i]);
	}

	// printf("\n");
	
	/*printf("%d:\t\t %.12lf\n", 1, x_k_1[1]);
	printf("%d:\t\t %.12lf\n", 12335, x_k_1[12335]);
	printf("%d:\t\t %.12lf\n", 34987, x_k_1[34987]);
	printf("%d:\t\t %.12lf\n", 98765, x_k_1[98765]);
	printf("%d:\t\t %.12lf\n", 444555, x_k_1[444555]);*/

	free(x_k);
	free(x_k_1);
	free(b);

	return 0;
}