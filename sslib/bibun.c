/*	bibun.c	*/
#include <stdio.h>
#include "sslib.h"

double diff2(double x, double h)
{
	double f1, f2;

	f1 = _f(x + h) - _f(x);
	f2 = _f(x + 2. * h) - _f(x);
	return (2. * f1 - 0.5 * f2) / h;
}

double diff3(double x, double h)
{
	double f1, f2, f3;

	f1 = _f(x + h) - _f(x);
	f2 = _f(x + 2. * h) - _f(x);
	f3 = _f(x + 3. * h) - _f(x);
	return (3. * f1 - 1.5 * f2 + f3 / 3.) / h;
}

double diff4(double x, double h)
{
	double f1, f2, f3, f4;

	f1 = _f(x + h) - _f(x);
	f2 = _f(x + 2. * h) - _f(x);
	f3 = _f(x + 3. * h) - _f(x);
	f4 = _f(x + 4. * h) - _f(x);
	return (4. * f1 - 3. * f2 + 4. * f3 / 3. - 0.25 * f4) / h;
}

double diff5(double x, double h)
{
	double f1, f2, f3, f4, f5;

	f1 = _f(x + h) - _f(x);
	f2 = _f(x + 2. * h) - _f(x);
	f3 = _f(x + 3. * h) - _f(x);
	f4 = _f(x + 4. * h) - _f(x);
	f5 = _f(x + 5. * h) - _f(x);
	return (5. * (f1 - f2) + 10. * f3 / 3. - 1.25 * f4 + 0.2 * f5) / h;
}

double difb2(double x, double h)
{
	double f1, f2;

	f1 = _f(x - h) - _f(x);
	f2 = _f(x - 2. * h) - _f(x);
	return - (2. * f1 - 0.5 * f2) / h;
}

double difb3(double x, double h)
{
	double f1, f2, f3;

	f1 = _f(x - h) - _f(x);
	f2 = _f(x - 2. * h) - _f(x);
	f3 = _f(x - 3. * h) - _f(x);
	return - (3. * f1 - 1.5 * f2 + f3 / 3.) / h;
}

double difb4(double x, double h)
{
	double f1, f2, f3, f4;

	f1 = _f(x - h) - _f(x);
	f2 = _f(x - 2. * h) - _f(x);
	f3 = _f(x - 3. * h) - _f(x);
	f4 = _f(x - 4. * h) - _f(x);
	return - (4. * f1 - 3. * f2 + 4. * f3 / 3. - 0.25 * f4) / h;
}

double difb5(double x, double h)
{
	double f1, f2, f3, f4, f5;

	f1 = _f(x - h) - _f(x);
	f2 = _f(x - 2. * h) - _f(x);
	f3 = _f(x - 3. * h) - _f(x);
	f4 = _f(x - 4. * h) - _f(x);
	f5 = _f(x - 5. * h) - _f(x);
	return - (5. * (f1 - f2) + 10. * f3 / 3. - 1.25 * f4 + 0.2 * f5) / h;
}
