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

int main()
{
	assert(f(0,0) == -5.94);

	return 0;
}
