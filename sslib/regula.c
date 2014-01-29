/*		regula.c		*/
#include <stdio.h>
#include "sslib.h"

double regula(double xs, double xe, double h, double eps, int iter)
{
	int ic;
	double x0, x1, x2, y0, y1, y2;

	if(xs == xe || h <= 0. || eps <= 0. || iter < 1)
	{
		fprintf(stderr, "Error : Illegal parameter in regula()\n");
		return 0.;
	}
	if(xs > xe)	swapd(&xs, &xe);
	x0 = xs;
	x1 = xe;
	y1 = _f(x1);
	do
	{
		y0 = _f(x0);
		if(y0 * y1 > 0.)	x0 += h;
		else
		{
			for(ic = 0; ic < iter; ic++)
			{
				x2 = x0 - y0 * (x1 - x0) / (y1 - y0);
				y2 = _f(x2);
				if(y0 * y2 >= 0.)
				{
					x0 = x2;
					y0 = y2;
				}
				else
				{
					x1 = x2;
					y1 = y2;
				}
				if(fabs(y2) <= eps)	return x2;
			}
			fprintf(stderr, "Error : No convergence in regula()\n");
			return x2;
		}
	} while(x0 <= xe);
	fprintf(stderr, "Error : No root in range(%e - %e) on regula()\n", xs, xe);
	return 0.;
}
