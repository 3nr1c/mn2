#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/** BEGIN TEST **/
#include <assert.h>
/** BEGIN TEST **/

#define TOL 1e-10
#define MAX_ITER 1000
#define ERROR 1e6

/*
 * Computes the dot product (scalar product)
 * of two vectors of $\mathbb{R}^2$
 */
double dot_product(double *x, double *y)
{
	return x[0]*y[0] + x[1]*y[1];
}

/*
 * Compute the value of f(x,y)
 */
double f(double x, double y)
{
	return pow(x,4) + 2*pow(x,3)*y + 3*pow(x,2)*pow(y,2) + 2*x*pow(y,3) 
		+ 2*pow(y,4) - 2.6*pow(x,3) - 3.2*pow(x,2)*y -2.6*x*pow(y,2)
		- 3.2*pow(y,3) - 7.22*pow(x,2) - 16*x*y - 15.22*pow(y,2) + 20.8*x + 25.6*y - 5.94;
}

/*
 * Compute the partial derivative
 * $\frac{\partial f}{\partial x}$ at (x,y)
 */
double df_x(double x, double y)
{
	return 4*pow(x,3) + 6*pow(x,2)*y + 6*x*pow(y,2) + 2*pow(y,3) - 7.8*pow(x,2) - 6.4*x*y
		-2.6*pow(y,2) - 14.44*x - 16*y + 20.8;
}

/*
 * Compute the partial derivative
 * $\frac{\partial f}{\partial y}$ at (x,y)
 */
double df_y(double x, double y)
{
	return 2*pow(x,3) + 6*pow(x,2)*y + 6*x*pow(y,2) + 8*pow(y,3) - 3.2*pow(x,2)
		-5.2*x*y - 9.6*pow(y,2) - 16*x - 30.44*y + 25.6;
}

/*
 * Computes the gradient of f at (x,y)
 */
void grad_f(double x, double y, double *result)
{
	result[0] = df_x(x, y);
	result[1] = df_y(x, y);
}

/*
 * Uses Newton's method to find a first
 * approximation (x,0) near (x0,0)
 */
double first_approx_x(double x0)
{
	int n = MAX_ITER;
	double x1;
	
	do {
		if (df_x(x0, 0) == 0) return ERROR;
		x1 = x0 - f(x0, 0) / df_x(x0, 0);
		x0 = x1;
	} while (fabs(f(x0, 0)) >= TOL && --n > 0);

	return x1;
}

/*
 * Uses Newton's method to find a first
 * approximation (0,y) near (0,y0)
 */
double first_approx_y(double y0)
{
	int n = MAX_ITER;
	double y1;
	
	do {
		y1 = y0 - f(0, y0) / df_y(0, y0);
		y0 = y1;
	} while (fabs(f(0,y0)) >= TOL && --n > 0);

	return y1;
}

/*
 * Given a point (x0,y0) on the curve f(x,y)=0
 * and a first approximation (x1,y1) at a distance
 * about h from (x0,y0), uses Newton's method
 * to find another point on the curve
 */
void newton_next_point(double x0, double y0, 
				 double x1, double y1, 
				 double h, double *result)
{
	double *H = (double*) calloc(2, sizeof(double));
	double *partial_f = (double*)calloc(2, sizeof(double));
	double detH, x2, y2;
	double incr_x, incr_y;
	int n = MAX_ITER;
	
	while ((fabs(f(x1,y1)) >= TOL 
		   || fabs(pow(x1 - x0,2) + pow(y1 - y0,2) - pow(h,2)) >= TOL)
		   && --n > 0) {
		H[0] = f(x1, y1);
		H[1] = pow(x1 - x0, 2) + pow(y1 - y0, 2) - pow(h,2);
	
		grad_f(x1, y1, partial_f);
		detH = 2*(y1 - y0)*partial_f[0] - 2*(x1 - x0)*partial_f[1];

		if (detH == 0) break; // ERROR

		incr_x = (2*(y1 - y0)*H[0] - partial_f[1]*H[1]) / detH;
		incr_y = (-2*(x1 - x0)*H[0] + partial_f[0]*H[1]) / detH;

		x2 = x1 - incr_x;
		y2 = y1 - incr_y;

		x1 = x2;
		y1 = y2;
	} 

	free(H);
	free(partial_f);

	result[0] = x1;
	result[1] = y1;
}

/*
 * Given a point on the curve f(x,y)=0
 * finds a vector tangent to the curve
 * with norm 1
 */
void tangent(double *point, double *vect)
{
	double norm;

	vect[0] = -df_y(point[0], point[1]);
	vect[1] = df_x(point[0], point[1]);
	norm = sqrt(pow(vect[0],2) + pow(vect[1],2));
	if (norm == 0) return; // ERRROR
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
	int i;

	// temp storage for the current point
	double *point = (double *)calloc(2, sizeof(double));
	double *tgt = (double *)calloc(2, sizeof(double));

	// store last iteration derivative to check we're going
	// in the right direction
	double *tgt0 = (double *)calloc(2, sizeof(double));

	FILE* output = fopen("results.txt", "w");

	// first points
	double firstx = first_approx_x(0);

	point[0] = firstx;
	point[1] = 0;

	// Save the first point
	fprintf(output, "%.12lf %.12lf\n", point[0], point[1]);

	// Print the first point (debugging)
	fprintf(stderr, "(%.12lf, %.12lf); ", point[0], point[1]);
	
	for (i = 0; i < 10000; i++) {
		// Compute tangent vector and change
		// sign if necessary
		tangent(point, tgt);
		if (i > 0 && dot_product(tgt, tgt0) < 0) {
			tgt[0] *= -1;
			tgt[1] *= -1;
		}
		tgt0[0] = tgt[0];
		tgt0[1] = tgt[1];

		// Print the gradient and tangent vector
		// for debugging purposes
		fprintf(stderr, "gradf(x,y) = (%.6lf, %.6lf); ", 
				df_x(point[0], point[1]),
				df_y(point[0], point[1]));
		fprintf(stderr, "y' = (%.6lf, %.6lf)\n", tgt[0], tgt[1]);

		/**
		 * Do Newton's method to find the next point
		 */
		newton_next_point(point[0], point[1],
					point[0] + delta * tgt[0], point[1] + delta * tgt[1],
					delta, point);
		
		// Print the next point (debugging)
		fprintf(stderr, "(%.12lf, %.12lf); ", point[0], point[1]);

		// Check if the point is too large, Not a Number,
		// or the value of f(x,y) is not zero
		if (isnan(point[0])
			|| isnan(point[1])
			|| fabs(point[0]) + fabs(point[1]) > 10
			|| fabs(f(point[0], point[1])) >= TOL) {
			fprintf(stderr, "\nHALT AT %d STEPS\n", i);
			break;
		}

		// Save the point to a file
		fprintf(output, "%.12lf %.12lf\n", point[0], point[1]);
	}
	fprintf(stderr, "\n");

	// Free memory
	fclose(output);
	free(point);
	free(tgt);
	free(tgt0);

	return 0;
}
