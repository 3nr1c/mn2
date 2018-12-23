#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define TOL 1e-10

double dot_product(double *x, double *y)
{
	return x[0]*y[0] + x[1]*y[1];
}

double norm_inf(double a, double b) 
{
	if (fabs(a) > fabs(b)) return fabs(a);
	else return fabs(b);
}

double f(double x, double y)
{
	return pow(x,4) + 2*pow(x,3)*y + 3*pow(x,2)*pow(y,2) + 2*x*pow(y,3) 
		+ 2*pow(y,4) - 2.6*pow(x,3) - 3.2*pow(x,2)*y -2.6*x*pow(y,2)
		- 3.2*pow(y,3) - 7.22*pow(x,2) - 16*x*y - 15.22*pow(y,2) + 20.8*x + 25.6*y - 5.94;
}

double df_x(double x, double y)
{
	return 4*pow(x,3) + 6*pow(x,2)*y + 6*x*pow(y,2) + 2*pow(y,3) - 7.8*pow(x,2) - 6.4*x*y
		-2.6*pow(y,2) - 14.44*x - 16*y + 20.8;
}

double df_y(double x, double y)
{
	return 2*pow(x,3) + 6*pow(x,2)*y + 6*x*pow(y,2) + 8*pow(y,3) - 3.2*pow(x,2)
		-5.2*x*y - 9.6*pow(y,2) - 16*x - 30.44*y + 25.6;
}

void grad_f(double x, double y, double *result)
{
	result[0] = df_x(x, y);
	result[1] = df_y(x, y);
}

double first_approx_x(double x0)
{
	int n = 20;
	double x1;
	
	do {
		x1 = x0 - f(x0, 0) / df_x(x0, 0);
		x0 = x1;
	} while (fabs(f(x0, 0)) >= TOL && --n > 0);

	return x1;
}

double first_approx_y(double y0)
{
	int n = 20;
	double y1;
	
	do {
		y1 = y0 - f(0, y0) / df_y(0, y0);
		y0 = y1;
	} while (fabs(f(0,y0)) >= TOL && --n > 0);

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
	int n = 100;
	double incr_x, incr_y;
	
	while (fabs(f(x1,y1)) >= TOL 
		   //&& fabs((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) - h*h) > TOL 
		   && --n > 0) {
		H[0] = f(x1, y1);
		H[1] = pow(x1 - x0, 2) + pow(y1 - y0, 2) - pow(h,2);
	
		grad_f(x1, y1, partial_f);
		detH = 2*(y1 - y0)*partial_f[0] - 2*(x1 - x0)*partial_f[1];

		incr_x = (2*(y1 - y0)*H[0] - partial_f[1]*H[1]) / detH;
		incr_y = (-2*(x1 - x0)*H[0] + partial_f[0]*H[1]) / detH;

		x2 = x1 - incr_x;
		y2 = y1 - incr_y;

		x1 = x2;
		y1 = y2;

		//fprintf(stderr, "it%d (%.6lf, %.6lf)\n", n, x1, y1);
	} 

	free(partial_f);

	result[0] = x1;
	result[1] = y1;
}

void tangent(double *point, double *vect, double sign)
{
	double norm;

	vect[0] = -sign*df_y(point[0], point[1]);
	vect[1] = sign*df_x(point[0], point[1]);
	norm = sqrt(pow(vect[0],2) + pow(vect[1],2));
	vect[0] = vect[0]/norm;
	vect[1] = vect[1]/norm;
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
	double delta = 0.01;

	// first points
	double firstx = first_approx_x(0);
	double firsty = first_approx_y(0);

	// temp storage for the current point
	double *results = (double *)calloc(2, sizeof(double));
	double *prev = (double *)calloc(2, sizeof(double));
	double *yprime = (double *)calloc(2, sizeof(double));

	// store last iteration derivative to check we're going
	// in the right direction
	double *yprime0 = (double *)calloc(2, sizeof(double));
	// this is to normalize the derivative vectors
	double norm;


	FILE* output = fopen("results.txt", "w");

	results[0] = firstx;
	results[1] = 0;

	//results[0] = 0;
	//results[1] = firsty;

	fprintf(output, "%.12lf %.12lf\n", results[0], results[1]);
	fprintf(stderr, "(%.12lf, %.12lf); ", results[0], results[1]);


	
	for (int i = 0; i < 10000; i++) {

		/*if (i > 1) {
			yprime[0] = results[0] - prev[0];
			yprime[1] = results[1] - prev[1];
			norm = sqrt(dot_product(yprime, yprime));
			yprime[0] /= norm;
			yprime[1] /= norm;
		}*/

		tangent(results, yprime, 1);
		if (i > 0 && dot_product(yprime, yprime0) < 0) {
			yprime[0] *= -1;
			yprime[1] *= -1;
		}
		yprime0[0] = yprime[0];
		yprime0[1] = yprime[1];

		prev[0] = results[0];
		prev[1] = results[1];

		fprintf(stderr, "gradf(x,y) = (%.6lf, %.6lf); ", 
				df_x(results[0], results[1]),
				df_y(results[0], results[1]));
		fprintf(stderr, "y' = (%.6lf, %.6lf)\n", yprime[0], yprime[1]);

		
		newton_step(results[0], results[1],
					results[0] + delta * yprime[0], results[1] + delta * yprime[1],
					delta, results);
		
		fprintf(stderr, "(%.12lf, %.12lf); ", results[0], results[1]);

		if (results[0] != results[0]
			|| results[1] != results[1]
			|| fabs(results[0]) + fabs(results[1]) > 10
			|| fabs(f(results[0], results[1])) >= TOL) {
			fprintf(stderr, "\nHALT AT %d STEPS\n", i);
			break;
		}
		fprintf(output, "%.12lf %.12lf\n", results[0], results[1]);
	}

	fclose(output);
	free(results);
	free(prev);
	free(yprime);
	free(yprime0);

	return 0;
}
