/*		vari.c				*/
#include <stdio.h>
#include <stdlib.h>
#include "sslib.h"

/* analysis of variance (one-way) */

int aov1(double x[], int n[], int l, int m, double al, double *sa, double *se,
	double *st, int *nua, int *nue, int *nut, double *va, double *ve,
	double *f0, double *f)
{
	int i, j, k, t, it, in;
	double xt, xi, xtb, w;
	double *row;

	if(l < 1 || m < 1 || al <= 0.0 || al >= 1.0)	return 999;
	row = (double *)malloc(m * sizeof(double));
	if(row == NULL)
	{
		fprintf(stderr, "Error : Out of memory in aov1()\n");
		return -1;
	}
	it = 0;
	xt = 0.0;
	for(i = 0; i < m; i++)
	{
		k = i * l;
		xi = 0.0;
		in = n[i];
		if(in < 1 || in > l)	return 999;
		it += in;
		for(j = 0; j < in; j++)	xi += x[k + j];
		row[i] = xi / in;
		xt += row[i];
	}
	xtb = xt / m;
	*se = *st = 0.0;
	for(i = 0; i < m; i++)
	{
		k = i * l;
		in = n[i];
		for(j = 0; j < in; j++)
		{
			t = k + j;
			w = x[t] - row[i];
			*se += (w * w);
			w = x[t] - xtb;
			*st += (w * w);
		}
	}
	*sa = *st - *se;
	*nua = m - 1;
	*nut = it - 1;
	*nue = *nut - *nua;
	*va = *sa / *nua;
	*ve = *se / *nue;
	*f0 = *va / *ve;

	*f = pf(al, *nua, *nue);
	free(row);
	return 0;
}

/* analysis of variance (two-way) */

int aov2(double x[], int l, int m, int n, double al, double *sa, double *sb,
	double *se, double *st, int *nua, int *nub, int *nue, int *nut, double *va,
	double *vb, double *ve, double *f0a, double *f0b, double *fa, double *fb)
{
	int i, j, k;
	double xt, xa, xb, xtb, w;
	double *xat, *xbt, *xab, *xbb;

	if(l < 1 || m < 1 || n < 1 || al <= 0.0 || al >= 1.0)	return 999;
	xat = (double *)malloc(m * sizeof(double));
	if(xat == NULL)
	{
		fprintf(stderr, "Error : Out of memory in aov2()\n");
		return -1;
	}
	xbt = (double *)malloc(n * sizeof(double));
	if(xbt == NULL)
	{
		fprintf(stderr, "Error : Out of memory in aov2()\n");
		free(xat);
		return -1;
	}
	xab = (double *)malloc(m * sizeof(double));
	if(xab == NULL)
	{
		fprintf(stderr, "Error : Out of memory in aov2()\n");
		free(xbt);
		free(xat);
		return -1;
	}
	xbb = (double *)malloc(n * sizeof(double));
	if(xbb == NULL)
	{
		fprintf(stderr, "Error : Out of memory in aov2()\n");
		free(xab);
		free(xbt);
		free(xat);
		return -1;
	}
	xt = 0.0;
	for(i = 0; i < m; i++)
	{
		k = i * l;
		xa = 0.0;
		for(j = 0; j < n; j++)	xa += x[k + j];
		xt += xa;
		xat[i] = xa;
	}
	for(i = 0; i < n; i++)
	{
		xb = 0.0;
		for(j = 0; j < m; j++)	xb += x[i + j * l];
		xbt[i] = xb;
	}
	for(i = 0; i < m; i++)	xab[i] = xat[i] / (double)n;
	for(i = 0; i < n; i++)	xbb[i] = xbt[i] / (double)m;
	xtb = xt / (double)m / (double)n;
	*sa = 0.0;
	for(i = 0; i < m; i++)
	{
		w = xab[i] - xtb;
		*sa += (w * w);
	}
	*sa *= (double)n;
	*sb = 0.0;
	for(i = 0; i < n; i++)
	{
		w = xbb[i] - xtb;
		*sb += (w * w);
	}
	*sb *= (double)m;
	*st = 0.0;
	for(i = 0; i < m; i++)
	{
		k = i * l;
		for(j = 0; j < n; j++)
		{
			w = x[k + j] - xtb;
			*st += (w * w);
		}
	}
	*se = *st - *sa - *sb;
	*nua = m - 1;
	*nub = n - 1;
	*nue = *nua * *nub;
	*nut = *nua + *nub + *nue;
	*va = *sa / *nua;
	*vb = *sb / *nub;
	*ve = *se / *nue;
	*f0a = *va / *ve;
	*f0b = *vb / *ve;

	*fa = pf(al, *nua, *nue);
	*fb = pf(al, *nub, *nue);

	free(xbb);
	free(xab);
	free(xbt);
	free(xat);
	return 0;
}
