/*		rkg.c		*/
#include <stdio.h>
#include "sslib.h"

extern double _fxy(double x, double y);
extern void _fmxy(double x, double y[], double dif[], int multi);


int rngkg(double x0, double y0, int n, double h, double y[])
{
	static double cq[4] = { 2.0, 1.0, 1.0, 2.0 };
	static double ckq[4] = { 0.5, 0.29289321881345248, 1.70710678118654752, 0.16666666666666667 };
	static double ck[4] = { 0.5, 0.29289321881345248, 1.70710678118654752, 0.5 };
	double k, q, r, xx, yy;
	int i, j;

	if(h <= 0.0)	return 999;
	q = 0.0;
	xx = x0;
	yy = y0;
	for(i = 1; i <= n; i++)
	{
		for(j = 0; j <= 3; j++)
		{
			k = h * _fxy(xx, yy);
			r = (k - cq[j] * q) * ckq[j];
			yy += r;
			q += 3.0 * r - ck[j] * k;
			xx += h * ck[j];
			if(j == 3)	y[i] = yy;
		}
	}
	return 0;
}

int hamng(double x0, double y0, int n, double h, double y[])
{
	double f[100];
	double x, xp, yi, yc, ym, yp;
	static double yim3 = 0.0, ycm = 0.0, ypm = 0.0;
	int i;

	if(h <= 0.0 || n > 99)	return 999;
	if(rngkg(x0, y0, n, h, y) != 0)	return 999;
	f[1] = _fxy(x0 + h, y[1]);
	f[2] = _fxy(x0 + 2.0 * h, y[2]);
	for(i = 3; i < n - 1; i++)
	{
		x = x0 + i * h;
		xp = x + h;
		f[i] = _fxy(x, y[i]);
		if(i > 3)	yim3 = y[i - 3];
		yp = yim3 + h * (2.0 * (f[i] + f[i - 2]) - f[i - 1]) * 1.333333333333333333;
		ym = yp - (ypm - ycm) * 112.0 / 121.0;
		yc = (9.0 * y[i] - y[i - 2] + 3.0 * h * (_fxy(xp, ym) + 2.0 * f[i] - f[i - 1])) * 0.125;
		y[i + 1] = yc + (yp - yc) * 9.0 / 121.0;
		ypm = yp;
		ycm = yc;
	}
	x += h;
	return 0;
}

int rngkgm(double x, double y[], double h, int multi, int n)
{
	static double cq[4] = { 2.0, 1.0, 1.0, 2.0 };
	static double ch[4] = { 0.5, 0.0, 0.5, 0.0 };
	static double ckq[4] = { 0.5, 0.29289321881345248, 1.70710678118654752, 0.16666666666666667 };
	static double ck[4] = { 0.5, 0.29289321881345248, 1.70710678118654752, 0.5 };
	double dif[10], q[10], ak, r, x0;
	int i, j, k;

	if(h <= 0.0 || multi > 9)	return 999;
	for(i = 0; i < multi; i++)	q[i] = 0.0;
	x0 = x;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < 4; j++)
		{
			_fmxy(x, y, dif, multi);
			for(k = 0; k < multi; k++)
			{
				ak = h * dif[k];
				r = (ak - cq[j] * q[k]) * ckq[j];
				y[k] += r;
				q[k] += 3.0 * r - ck[j] * ak;
			}
			x = x0 + h * ch[j];
		}
	}
	return 0;
}
