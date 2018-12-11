#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double f(double x, double y)
{
	double x2 = x*x;
	double y2 = y*y;
	return x2*x2 + 2*x*x2*y + 3*x2*y2 + 2*x*y*y2 + 2*y2*y2 - 2.6*x*x2 - 3.2*x2*y
		-2.6*x*y2 - 3.2*y*y2 - 7.22*x2 - 16*x*y - 15.22*y2 + 20.8*x + 25.6*y - 5.94;
}

void grad_f(double x, double y, double *result)
{
	double x2 = x*x;
	double y2 = y*y;

	double partial_x = 4*x*x2 + 6*x2*y + 6*x*y2 + 0 + 0 + 7.8*x2 - 6.4*x*y
		-2.6*y2 - 0 - 14.44*x - 16*y - 0 + 20.8 + 0 - 0;
	double partial_y = 0 + 2*x*x2 + 6*x2*y + 6*x*y2 + 8*y*y2 - 0 - 3.2*x2
		-5.2*x*y - 9.6*y2 - 0 - 16*x - 30.44*y + 0 + 25.6 - 0;

	result[0] = partial_x;
	result[1] = partial_y;
}

/** BEGIN TEST **/
void test()
{
	double *temp = (double *)calloc(2, sizeof(double));
	assert(f(0, 0) == -5.94);
	grad_f(0, 0, temp);
	assert(temp[0] == 20.8);
	assert(temp[1] == 25.6);
}
/** END TEST **/

int main()
{	
	/** BEGIN TEST **/
	test();
	/** END TEST **/

	return 0;
}
