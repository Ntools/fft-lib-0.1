/*		corr.c				*/
#include <stdio.h>
#include "sslib.h"

/* correlation coefficient (‘ŠŠÖŒW”) */

int corr(double x[], double y[], int n, double xbar, double ybar,
	double *uxy, double *r0)
{
	int i;
	double w1, w2, w3, wx, wy;

	if(n < 2)	return 999;
	w1 = w2 = w3 = 0.0;
	n -= 1;
	for(i = 0; i <= n; i++)
	{
		wx = x[i] - xbar;
		wy = y[i] - ybar;
		w1 += (wx * wy);
		w2 += (wx * wx);
		w3 += (wy * wy);
	}
	*uxy = w1 / n;
	*r0 = w1 / sqrt(w2 * w3);
	return 0;
}

/* simple regression (’P‰ñ‹A•ªÍ) */

int sreg(double x[], double y[], int n, double sx, double g,
	double *a, double *b, double *dyx, double *va, double *vb, double *al,
	double *au, double *bl, double *bu)
{
	int nu, i;
	double c[2];
	double fn, fnu, qt, syx, xx, wx, yc, dy, t;
	double wa, wb, wdyx, wva, wvb, wat, wbt;

	if(n < 3 || g <= 0.0 || g >= 1.0)	return 999;
	nu = n - 2;
	fn = n;
	fnu = nu;
	qt = 1.0 - g;
	syx = xx = 0.0;

	lstsq(x, y, n, 1, c);

	wa = c[0];
	wb = c[1];
	for(i = 0; i < n; i++)
	{
		wx = x[i];
		yc = wa + wb * wx;
		dy = y[i] - yc;
		syx += (dy * dy);
		xx += (wx * wx);
	}
	wdyx = syx / fnu;
	wva = xx * wdyx / (fn * sx);
	wvb = wdyx / sx;

	t = pt(qt, nu);

	wat = t * sqrt(wva);
	wbt = t * sqrt(wvb);
	*al = wa - wat;
	*au = wa + wat;
	*bl = wb - wbt;
	*bu = wb + wbt;
	*a = wa;
	*b = wb;
	*dyx = wdyx;
	*va = wva;
	*vb = wvb;
	return 0;
}
