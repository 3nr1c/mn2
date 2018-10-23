#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM 1000000

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
	double divider;
	int i,k;

	// printf("Memory allocated\n");

	// create b vector
	fillB_J(b);

	// printf("Vector b initialized\n");

	// actual iteration
	for (k = 0; k < 200; k++) {
		// printf("Starting iteration %d\n", k+1);

		for (i = 0; i < NUM; i++) {
			divider = (double)(!(i % 2) ? 3 : 4);

			x_k_1[i] = b[i];
			if (i - 2 >= 0) 
				x_k_1[i] += x_k_1[i-2] / divider;
			if (i + 2 < NUM)
				x_k_1[i] += x_k[i+2] / divider;
			if (i == NUM - 2)
				x_k_1[i] += -x_k_1[0] / divider;
			if (i == NUM - 1)
				x_k_1[i] += -x_k_1[1] / divider;
			if (i == 0)
				x_k_1[i] += -x_k[NUM - 2] / divider;
			if (i == 1)
				x_k_1[i] += -x_k[NUM - 1] / divider;

		}
		//printf("Finished iteration %d with increment %.13lf\n", k+1, error);

		swap(&x_k, &x_k_1);
	}

	// print start of b
	for (i = 0; i < NUM; i++) {
		printf("%.12lf\n", x_k[i]);
	}

	free(x_k);
	free(x_k_1);
	free(b);

	return 0;
}