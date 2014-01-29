/*	test.c	*/
#include <stdio.h>
#include "sslib.h"

int mtst1(int n, double xbar, double pm, double pv, double al,
	int *sw, double *u0, double *u)
{
	double qn;

	if(n < 1 || al <= 0.0 || al >= 1.0)	return 999;
	*u0 = (xbar - pm) * sqrt(n / pv);
	if(*sw < 0)	qn = 1.0 - al;
	else		qn = 1.0 - 0.5 * al;
	*u = pnorm(qn);
	if(fabs(*u0) >= *u)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int mtst2(int n, double xbar, double pm, double uv, double al,
	int *sw, double *t0, double *t)
{
	double qt;

	if(n < 1 || al <= 0.0 || al >= 0.5)	return 999;
	*t0 = (xbar - pm) * sqrt(n / uv);
	if(*sw < 0)	qt = 2.0 * al;
	else		qt = al;
	n--;
	*t = pt(qt, n);
	if(fabs(*t0) >= *t)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int mdtst1(int n1, int n2, double xbar1, double xbar2, double pv1, double pv2,
	double al, int *sw, double *u0, double *u)
{
	double w1, w2, qn;

	if(n1 < 1 || n2 < 1 || al <= 0.0 || al >= 1.0)	return 999;
	w1 = xbar1 - xbar2;
	w2 = sqrt(pv1 / n1 + pv2 / n2);
	*u0 = w1 / w2;
	if(*sw < 0)	qn = 1.0 - al;
	else		qn = 1.0 - 0.5 * al;
	*u = pnorm(qn);
	if(fabs(*u0) >= *u)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int mdtst2(int n1, int n2, double xbar1, double xbar2, double s1, double s2,
	double al, int *sw, double *t0, double *t)
{
	int nu;
	double w1, w2, w3, qt;

	if(n1 < 2 || n2 < 2 || al <= 0.0 || al >= 0.5)	return 999;
	nu = n1 + n2 - 2;
	w1 = xbar1 - xbar2;
	w2 = sqrt((s1 + s2) / nu);
	w3 = sqrt(1.0 / n1 + 1.0 / n2);
	*t0 = w1 / w2 / w3;
	if(*sw < 0)	qt = 2.0 * al;
	else		qt = al;
	*t = pt(qt, nu);
	if(fabs(*t0) >= *t)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int mdtst3(int n1, int n2, double xbar1, double xbar2, double uv1, double uv2,
	double al, int *sw, double *t0, double *t)
{
	int nu;
	double w, w1, w2, w3, qt;

	if(n1 < 2 || n2 < 2 || al <= 0.0 || al >= 0.5)	return 999;
	w1 = xbar1 - xbar2;
	w2 = uv1 / n1;
	w3 = uv2 / n2;
	w = w2 + w3;
	*t0 = w1 / sqrt(w);
	w1 = w2 / w;
	w2 = w1 * w2 / (double)(n1 - 1);
	w = 1.0 - w1;
	w3 = w * w /(double)(n2 - 1);
	nu = 1.0 / (w2 + w3);
	if(*sw < 0)	qt = 2.0 * al;
	else		qt = al;
	*t = pt(qt, nu);
	if(fabs(*t0) >= *t)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int nbptst(int n, double ps, double pp, double al, int *sw, double *u0, double *u)
{
	double qn;

	if(n < 1 || ps <= 0.0, ps >= 1.0 || pp <= 0. || pp >= 1.0 || al <= 0.0 ||
	   al >= 1.0)	return 999;
	*u0 = (ps - pp) / sqrt(pp * (1.0 - pp) / n);
	if(*sw < 0)	qn = 1.0 - al;
	else		qn = 1.0 - 0.5 * al;
	*u = pnorm(qn);
	if(abs(*u0) >= *u)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int fbptst(int n, double ps, double pp, double al, int *sw, double *f0, double *f)
{
	int nu1, nu2;
	double w;

	if(n < 1 || ps <= 0.0 || ps >= 1.0 || pp <= 0.0 || pp >= 1.0 ||
		al <= 0.0 || al >= 1.0 || ps == pp)	return 999;
	w = ps * n;
	if(ps < pp)
	{
		nu1 = 2.0 * (w + 1.0);
		nu2 = 2.0 * ((double)n - w);
		*f0 = nu2 * pp / nu1 / (1.0 - pp);
	}
	else
	{
		nu1 = 2.0 * ((double)n - w + 1.0);
		nu2 = 2.0 * w;
		*f0 = nu2 * (1.0 - pp) / nu1 / pp;
	}
	*f = pf(al, nu1, nu2);
	if(fabs(*f0) >= *f)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int bpdtst(int n1, int n2, double ps1, double ps2, double al,
	int *sw, double *u0, double *u)
{
	double w, w1, w2, qn;

	if(n1 < 1 || n2 < 1 || ps1 <= 0.0 || ps1 >= 1.0 || ps2 <= 0.0 ||
		ps2 >= 1.0 || al <= 0.0 || al >= 1.0)	return 999;
	w = (ps1 * (double)n1 + ps2 * (double)n2) / (double)(n1 + n2);
	w1 = ps1 - ps2;
	w2 = w * (1.0 - w) * (1.0 / (double)n1 + 1.0 / (double)n2);
	*u0 = w1 / sqrt(w2);
	if(*sw < 0)	qn = 1.0 - al;
	else		qn = 1.0 - 0.5 * al;
	*u = pnorm(qn);
	if(fabs(*u0) >= *u)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int vtst(int sw1, int n, double s, double v, double al,
	int *sw2, double *x20, double *x21, double *x22)
{
	double qc1;

	if(n < 2 || al <= 0.0 || al >= 1.0)	return 999;
	if(sw1 < 0)	n--;
	*x20 = s / v;
	if(*sw2 < 0)	qc1 = al;
	else			qc1 = 0.5 * al;
	*x21 = pchi(qc1, n);
	*x22 = pchi(1.0 - qc1, n);
	if(*x20 < *x21 && *x20 > *x22)	*sw2 = 0;
	else							*sw2 = 1;
	return 0;
}

int vdtst(int n1, int n2, double uv1, double uv2, double al,
	int *sw, double *f0, double *f)
{
	if(n1 < 2 || n2 < 2 || al <= 0.0 || al >= 1.0)	return 999;
	n1--;
	n2--;
	al *= 0.5;
	if(uv1 >= uv2)
	{
		*f0 = uv1 / uv2;
		*f = pf(al, n1, n2);
	}
	else
	{
		*f0 = uv2 / uv1;
		*f = pf(al, n2, n1);
	}
	if(*f0 >= *f)	*sw = 1;
	else			*sw = 0;
	return 0;
}

int cont22(int *sw, double a, double b, double c, double d, double al,
	double *x20, double *x2)
{
	double w1, w2, w3, w21, w22, w23, w24;

	if(al <= 0.0 || al >= 1.0)	return 999;
	w1 = a + b + c + d;
	w3 = (a + b) * (c + d) * (a + c) * (b + d);
	w21 = a * d - b * c;
	if(*sw >= 0)	w2 = w21 * w21;
	else
	{
		w22 = w1 * 0.5;
		w23 = w21 + w22;
		w24 = w21 - w22;
		if(fabs(w23) > fabs(w24))	w2 = w24 * w24;
		else						w2 = w23 * w23;
	}
	*x20 = w1 * w2 / w3;
	*x2 = pchi(al, 1);
	if(*x20 >= *x2)	*sw = 1;
	else			*sw = 0;
	return 0;
}

int contlm(int l, int m, double a[], double al, double ac[],
	double ar[], double *at, double ef[], double *x20, double *x2, int *sw)
{
	int i, j, k, n, nu;
	double w;

	if(l < 2 || m < 2 || al <= 0.0 || al >= 1.0)	return 999;
	for(i = 0; i < m; i++)
	{
		w = 0.0;
		for(j = 0; j < l; j++)
		{
			k = j * m + i;
			w += a[k];
		}
		ac[i] = w;
	}
	for(i = 0; i < l; i++)
	{
		n = i * m;
		w = 0.0;
		for(j = 0; j < m; j++)
		{
			k = n + j;
			w += a[k];
		}
		ar[i] = w;
	}
	*at = 0.0;
	for(i = 0; i < l; i++)	*at += ar[i];
	for(i = 0; i < l; i++)
	{
		n = i * m;
		for(j = 0; j < m; j++)
		{
			k = n + j;
			ef[k] = ar[i] * ac[j] / *at;
		}
	}
	*x20 = 0.0;
	for(i = 0; i < l; i++)
	{
		n = i * m;
		for(j = 0; j < m; j++)
		{
			k = n + j;
			w = a[k] - ef[k];
			*x20 += (w * w / ef[k]);
		}
	}
	nu = (l - 1) * (m - 1);
	*x2 = pchi(al, nu);
	if(*x20 >= *x2)	*sw = 1;
	else			*sw = 0;
	return 0;
}

int thomp(int n, double xk, double xbar, double v, double al,
	int *sw, double *t0, double *t)
{
	double w;

	if(n < 3 || al <= 0.0 || al >= 1.0)	return 999;
	w = (xk - xbar) / sqrt(v);
	n -= 2;
	*t0 = sqrt((double)n) * w / sqrt((double)n + 1.0 - w * w);
	*t = pt(al, n);
	if(fabs(*t0) >= *t)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int rttst(int n, double r0, double al, int *sw, double *t0, double *t)
{
	if(n < 2 || al <= 0.0 || al >= 1.0 || fabs(r0) >= 1.0)	return 999;
	n -= 2;
	*t0 = r0 * sqrt((double)n) / sqrt(1.0 - r0 * r0);
	*t = pt(al, n);
	if(fabs(*t0) >= *t)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int rptst(int n, double r0, double rp, double al, int *sw, double *u0, double *u)
{
	double w1, w2;

	if(n < 3 || al <= 0.0 || al >= 1.0 || fabs(r0) >= 1.0 || fabs(rp) >= 1.0)
		return 999;
	w1 = 0.5 * log((1.0 + r0) / (1.0 - r0));
	w2 = 0.5 * log((1.0 + rp) / (1.0 - rp));
	*u0 = (w1 - w2) * sqrt((double)(n - 3));
	*u = pnorm(1.0 - 0.5 * al);
	if(fabs(*u0) >= *u)	*sw = 1;
	else				*sw = 0;
	return 0;
}

int rptst2(int n1, int n2, double r01, double r02, double al,
	int *sw, double *u0, double *u)
{
	double w1, w2;

	if(n1 < 4 || n2 < 4 || al <= 0.0 || al >= 1.0 ||
		fabs(r01) >= 1.0 || fabs(r02) >= 1.0)	return 999;
	w1 = 0.5 * log((1.0 + r01) / (1.0 - r01));
	w2 = 0.5 * log((1.0 + r02) / (1.0 - r02));
	*u0 = (w1 - w2) / sqrt(1.0 / (double)(n1 - 3) + 1.0 / (double)(n2 - 3));
	*u = pnorm(1.0 - 0.5 * al);
	if(fabs(*u0) >= *u)	*sw = 1;
	else				*sw = 0;
	return 0;
}
