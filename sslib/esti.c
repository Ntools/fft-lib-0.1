/*	esti.c	*/
#include <stdio.h>
#include "sslib.h"

int mest1(int n, double xbar, double pv, double g, double *xl, double *xu)
{
	double w;

	if(n < 1 || g <= 0.0 || g >= 1.0)	return 999;
	w = pnorm((1.0 + g) * 0.5) * sqrt(pv / (double)n);
	*xl = xbar - w;
	*xu = xbar + w;
	return 0;
}

int mest2(int n, double xbar, double squv, double g, double *xl, double *xu)
{
	double w;

	if(n < 1 || g <= 0.0 || g >= 1.0)	return 999;
	w = pt(1.0 - g, n - 1) * squv / sqrt((double)n);
	*xl = xbar - w;
	*xu = xbar + w;
	return 0;
}

int mdest1(int n1, int n2, double xbar1, double xbar2, double pv1, double pv2,
	double g, double *xl, double *xu)
{
	double xb, w;

	if(n1 < 1 || n2 < 1 || g <= 0.0 || g >= 1.0)	return 999;
	xb = fabs(xbar1 - xbar2);
	w = pnorm((1.0 + g) * 0.5) * sqrt(pv1 / (double)n1 + pv2 / (double)n2);
	*xl = xb - w;
	*xu = xb + w;
	return 0;
}

int mdest2(int n1, int n2, double xbar1, double xbar2, double v1, double v2,
	double g, double *xl, double *xu)
{
	double xb, w;

	if(n1 < 1 || n2 < 1 || g <= 0.0 || g >= 1.0)	return 999;
	xb = fabs(xbar1 - xbar2);
	w = pt(1.0 - g, n1 + n2 - 2)
	  * sqrt(1.0 / (double)n1 + 1.0 / (double)n2)
	  * sqrt(((double)n1 * v1 + (double)n2 * v2) / (double)(n1 + n2 - 2));
	*xl = xb - w;
	*xu = xb + w;
	return 0;
}

int vest(int sw, int n, double s, double g, double *xl, double *xu)
{
	double q1;

	if(n < 2 || g <= 0.0 || g >= 1.0)	return 999;
	if(sw < 0)	n--;
	q1 = (1.0 + g) * 0.5;
	*xl = s / pchi(1.0 - q1, n);
	*xu = s / pchi(q1, n);
	return 0;
}

int vpest(int n1, int n2, double uv1, double uv2, double g, double *xl, double *xu)
{
	double f, w;

	if(n1 < 2 || n2 < 2 || g <= 0.0 || g >= 1.0)	return 999;
	w = uv1 / uv2;
	f = (1.0 - g) * 0.5;
	*xl = w / pf(f, n1 - 1, n2 - 1);
	*xu = w * pf(f, n2 - 1, n1 - 1);
	return 0;
}

int bpest(int n, double ps, double g, double *xl, double *xu)
{
	double w;

	if(n < 1 || ps <= 0.0 || ps >= 1.0 || g <= 0.0 || g >= 1.0)	return 999;
	w = pnorm((1.0 + g) * 0.5) * sqrt(ps * (1.0 - ps) / (double)n);
	*xl = ps - w;
	*xu = ps + w;
	return 0;
}

int bpdest(int n1, int n2, double ps1, double ps2, double g, double *xl, double *xu)
{
	double p, w;

	if(n1 < 1 || n2 < 1 || g <= 0.0 || g >= 1.0 || ps1 <= 0.0 || ps1 >= 1.0 ||
		ps2 <= 0.0 || ps2 >= 1.0)	return 999;
	p = fabs(ps1 - ps2);
	w = pnorm((1.0 + g) * 0.5)
	  * sqrt(ps1 * (1.0 - ps1) / (double)n1 + ps2 * (1.0 - ps2) / (double)n2);
	*xl = p - w;
	*xu = p + w;
	return 0;
}

int rest(int n, double r0, double g, double *xl, double *xu)
{
	double qn, w, w1, w2, w3;

	if(n < 4 || fabs(r0) >= 1.0 || g <= 0.0 || g >= 1.0)	return 999;
	w = log((1.0 + r0) / (1.0 - r0)) * 0.5;
	w1 = pnorm((1.0 + g) * 0.5) / sqrt((double)(n - 3));
	w2 = exp(2.0 * (w - w1));
	w3 = exp(2.0 * (w + w1));
	*xl = (w2 - 1.0) / (w2 + 1.0);
	*xu = (w3 - 1.0) / (w3 + 1.0);
	return 0;
}
