#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define TOL 1e-10

double f(double x, double y)
{
	double x2 = x*x;
	double y2 = y*y;
	return x2*x2 + 2*x*x2*y + 3*x2*y2 + 2*x*y*y2 + 2*y2*y2 - 2.6*x*x2 - 3.2*x2*y
		-2.6*x*y2 - 3.2*y*y2 - 7.22*x2 - 16*x*y - 15.22*y2 + 20.8*x + 25.6*y - 5.94;
}

double df_x(double x, double y)
{
	double x2 = x*x;
	double y2 = y*y;

	return 4*x*x2 + 6*x2*y + 6*x*y2 + 0 + 0 + 7.8*x2 - 6.4*x*y
		-2.6*y2 - 0 - 14.44*x - 16*y - 0 + 20.8 + 0 - 0;
}

double df_y(double x, double y)
{
	double x2 = x*x;
	double y2 = y*y;

	return 0 + 2*x*x2 + 6*x2*y + 6*x*y2 + 8*y*y2 - 0 - 3.2*x2
		-5.2*x*y - 9.6*y2 - 0 - 16*x - 30.44*y + 0 + 25.6 - 0;
}

void grad_f(double x, double y, double *result)
{
	result[0] = df_x(x, y);
	result[1] = df_y(x, y);
}

double first_approx_x(double x0)
{
	int n = 20;
	double x1, error;
	
	do {
		x1 = x0 - f(x0, 0) / df_x(x0, 0);
		error = x1 - x0;
		x0 = x1;
	} while (fabs(error) >= TOL && --n > 0);

	return x1;
}

double first_approx_y(double y0)
{
	int n = 20;
	double y1, error;
	
	do {
		y1 = y0 - f(0, y0) / df_y(0, y0);
		error = y1 - y0;
		y0 = y1;
	} while (fabs(error) >= TOL && --n > 0);

	return y1;
}

void newton_step(double x0, double y0, 
				 double x1, double y1, 
				 double h, double *result)
{
	double H[2];
	double detH;
	double *partial_f = (double*)calloc(2, sizeof(double));
	double x2, y2;
	int n = 1000;
	
	do {
		H[0] = f(x1, y1);
		H[1] = (x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) - h*h;
		grad_f(x1, y1, partial_f);

		detH = 2*(y1 - y0)*partial_f[0] - 2*(x1 - x0)*partial_f[1];

		x2 = x1 - (2*(y1 - y0)*H[0] - partial_f[1]*H[1]) / detH;
		y2 = y1 - (-2*(x1 - x0)*H[0] + partial_f[0]*H[1]) / detH;

		x1 = x2;
		y1 = y2;
	} while (fabs(f(x1,y1)) >= TOL 
		   && fabs((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) - h*h) >= TOL 
		   && --n > 0);

	free(partial_f);

	result[0] = x1;
	result[1] = y1;
}

/** BEGIN TEST **/
void test()
{
	double *temp = (double *)calloc(2, sizeof(double));
	assert(f(0, 0) == -5.94);

	assert(df_x(0, 0) == 20.8);
	assert(df_y(0, 0) == 25.6);

	grad_f(0, 0, temp);
	assert(temp[0] == 20.8);
	assert(temp[1] == 25.6);

	assert(first_approx_x(0) - 0.326344025693 < 1e-10);

	free(temp);
}
/** END TEST **/

int main()
{	
	/** BEGIN TEST **/
	test();
	/** END TEST **/
	double delta = 0.1;
	double firstx = first_approx_x(0);
	double firsty = first_approx_y(-3);
	double *results = (double *)calloc(2, sizeof(double));
	double *yprime = (double *)calloc(2, sizeof(double));
	double *prev = (double *)calloc(2, sizeof(double));
	double norm;
	FILE* output = fopen("results.txt", "w");

	results[0] = firstx;
	results[1] = 0;

	//results[0] = 0;
	//results[1] = firsty;

	//yprime[0] = delta/2;
	//yprime[1] = 0;

	printf("%.12lf %.12lf\n", results[0], results[1]);
	fprintf(stderr, "(%.12lf, %.12lf)\n", results[0], results[1]);
	for (int i = 0; i < 10000; i++) {
		/*yprime[0] = (results[0] - prev[0]);
		yprime[1] = (results[1] - prev[1]);
		norm = sqrt(yprime[0]*yprime[0] + yprime[1]*yprime[1]);
		yprime[0] = delta*yprime[0]/norm;
		yprime[1] = delta*yprime[1]/norm;*/

		yprime[0] = -df_y(results[0], results[1]);
		yprime[1] = df_x(results[0], results[1]);
		norm = sqrt(yprime[0]*yprime[0] + yprime[1]*yprime[1]);
		yprime[0] = delta*yprime[0]/norm;
		yprime[1] = delta*yprime[1]/norm;

		fprintf(stderr, "y' = (%.6lf, %.6lf); ", yprime[0], yprime[1]);

		prev[0] = results[0];
		prev[1] = results[1];

		newton_step(results[0], results[1],
					results[0] + yprime[0], results[1] + yprime[1],
					delta, results);
		fprintf(stderr, "(%.12lf, %.12lf)\n", results[0], results[1]);
		if (results[0] != results[0]
			|| results[1] != results[1]
			|| fabs(results[0])+fabs(results[1]) > 100) {
			fprintf(stderr, "\nHALT AT %d STEPS\n", i);
			break;
		}
		fprintf(output, "%.12lf %.12lf\n", results[0], results[1]);
	}

	fclose(output);
	free(results);

	return 0;
}
