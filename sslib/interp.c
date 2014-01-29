/*		interp.c		*/
#include <stdio.h>
#include <stdlib.h>
#include "sslib.h"

void lstsq(double x[], double y[], int n, int m, double c[])
{
	int i, j, k, m2, mp1, mp2;
	double *a, aik, pivot, *w, w1, w2, w3;

	if(m >= n || m < 1)
	{
		fprintf(stderr, "Error : Illegal parameter in lstsq()\n");
		return;
	}
	mp1 = m + 1;
	mp2 = m + 2;
	m2 = 2 * m;
	a = (double *)malloc(mp1 * mp2 * sizeof(double));
	if(a == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in lstsq()\n");
		return;
	}
	w = (double *)malloc(mp1 * 2 * sizeof(double));
	if(w == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in lstsq()\n");
		free(a);
		return;
	}
	for(i = 0; i < m2; i++)
	{
		w1 = 0.;
		for(j = 0; j < n; j++)
		{
			w2 = w3 = x[j];
			for(k = 0; k < i; k++)	w2 *= w3;
			w1 += w2;
		}
		w[i] = w1;
	}
	a[0] = n;
	for(i = 0; i < mp1; i++)
		for(j = 0; j < mp1; j++)	if(i || j)	a[i * mp2 + j] = w[i + j - 1];

	w1 = 0.;
	for(i = 0; i < n; i++)	w1 += y[i];
	a[mp1] = w1;
	for(i = 0; i < m; i++)
	{
		w1 = 0.;
		for(j = 0; j < n; j++)
		{
			w2 = w3 = x[j];
			for(k = 0; k < i; k++)	w2 *= w3;
			w1 += y[j] * w2;
		}
		a[mp2 * (i + 1) + mp1] = w1;
	}

	for(k = 0; k < mp1; k++)
	{
		pivot = a[mp2 * k + k];
		a[mp2 * k + k] = 1.0;
		for(j = k + 1; j < mp2; j++)	a[mp2 * k + j] /= pivot;
		for(i = 0; i < mp1; i++)
		{
			if(i != k)
			{
				aik = a[mp2 * i + k];
				for(j = k; j < mp2; j++)
					a[mp2 * i + j] -= aik * a[mp2 * k + j];
			}
		}
	}
	for(i = 0; i < mp1; i++)	c[i] = a[mp2 * i + mp1];
	free(w);
	free(a);
	return;
}

double lagra(double xx[], double yy[], int n, double xi)
{
	int flag, i, j, *jun, *k;
	double *p, *q, *r, *s, w, w1, w2, w3, *x, *y;

	x = (double *)malloc(n * sizeof(double));
	if(x == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in lagra()\n");
		return 0.;
	}
	y = (double *)malloc(n * sizeof(double));
	if(y == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in lagra()\n");
		free((char *)x);
		return 0.;
	}
	for(i = 0, p = xx, flag = 0; i < n - 1; i++, p++)
	{
		if(*p >= *(p + 1))
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		jun = (int *)malloc(n * sizeof(int));
		if(jun == NULL)
		{
			fprintf(stderr, "Error : Out of memory  in lagra()\n");
			free((char *)x);
			free((char *)y);
			return 0.;
		}
		sortdi1(xx, n, jun);
		for(i = 0, p = x, q = y, k = jun; i < n; i++)
		{
			*p++ = *(xx + *k);
			*q++ = *(yy + *k++);
		}
		free((char *)jun);
	}
	else
	{
		for(i = 0, p = x, q = y, r = xx, s = yy; i < n; i++)
		{
			*p++ = *r++;
			*q++ = *s++;
		}
	}

	if(n < 2 || xi < *x || xi > *(x + n - 1))
	{
		fprintf(stderr, "Error : Illegal parameter  in lagra()\n");
		free((char *)x);
		free((char *)y);
		return 0.;
	}
	for(i = 0, p = x; i < n; i++)
	{
		if(*p++ == xi)
		{
			free((char *)x);
			free((char *)y);
			return *(y + i);
		}
	}
	w = 0.;
	for(i = 0, p = x, q = y; i < n; i++)
	{
		w1 = 1.;
		w2 = *p++;
		for(j = 0, r = x; j < n; j++)
		{
			w3 = *r++;
			if(i != j)	w1 *= (xi - w3) / (w2 - w3);
		}
		w += w1 * *q++;
	}
	free((char *)x);
	free((char *)y);
	return w;
}

double splint(double xx[], double yy[], int n, double xi)
{
	int flag, i, *jun, *k;
	double dxp, dxm, *h, hi, hi2, *p, *q, *r, *s, sm, si, *sp, *x, *y, yi;

	x = (double *)malloc(n * sizeof(double));
	if(x == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		return 0.;
	}
	y = (double *)malloc(n * sizeof(double));
	if(y == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		free((char *)x);
		return 0.;
	}
	for(i = 0, p = xx, flag = 0; i < n - 1; i++, p++)
	{
		if(*p >= *(p + 1))
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		jun = (int *)malloc(n * sizeof(int));
		if(jun == NULL)
		{
			fprintf(stderr, "Error : Out of memory  in splint()\n");
			free((char *)x);
			free((char *)y);
			return 0.;
		}
		sortdi1(xx, n, jun);
		for(i = 0, p = x, q = y, k = jun; i < n; i++)
		{
			*p++ = *(xx + *k);
			*q++ = *(yy + *k++);
		}
		free((char *)jun);
	}
	else
	{
		for(i = 0, p = x, q = y, r = xx, s = yy; i < n; i++)
		{
			*p++ = *r++;
			*q++ = *s++;
		}
	}

	if(n < 2 || xi < *x || xi > *(x + n - 1))
	{
		fprintf(stderr, "Error : Illegal parameter  in splint()\n");
		free((char *)x);
		free((char *)y);
		return 0.;
	}
	for(i = 0, p = x; i < n; i++)
	{
		if(*p++ == xi)
		{
			free((char *)x);
			free((char *)y);
			return *(y + i);
		}
	}
	h = (double *)malloc(n * sizeof(double));
	if(h == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		free((char *)x);
		free((char *)y);
		return 0.;
	}
	sp = (double *)malloc(n * sizeof(double));
	if(sp == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		free((char *)x);
		free((char *)y);
		free((char *)h);
		return 0.;
	}

	subspl(x, y, n - 1, h, sp);
	for(i = 1, p = x + 1; i <= n; i++, p++)
	{
		if(*(p - 1) <= xi && xi < *p)
		{
			sm = *(sp + i - 1);
			si = *(sp + i);
			hi = *(h + i);
			hi2 = hi * hi;
			dxp = *p - xi;
			dxm = xi - *(p - 1);
			yi = (  sm * dxp * dxp * dxp + si * dxm * dxm * dxm
				  + (6. * *(y + i - 1) - hi2 * sm) * dxp
				  + (6. * *(y + i)     - hi2 * si) * dxm) / hi / 6.;
			free((char *)x);
			free((char *)y);
			free((char *)h);
			free((char *)sp);
			return yi;
		}
	}
	return 0;
}

double polynomial(double a[], double x, int n)
{
	int i;
	double w;

	w = a[n];
	for(i = n - 1; i >= 0; i--)	w = w * x + a[i];
	return w;
}

void chebyshev(double min, double max, int n, double a[])
{						/* min, max : range of approximation */
	double **T;
	double *c, *func, *x;
	double k1, k2, n1, w;
	int i, j;

	if(min == max || n <= 0)
	{
		fprintf(stderr, "Error : illegal parameter input in chebyshev()\n");
		return;
	}

	T = (double **)malloc((n + 1) * sizeof(double *));
	if(T == NULL)
	{
		fprintf(stderr, "Error : out of memory in chebyshev().\n");
		return;
	}
	for(i = 0; i <= n; i++)
	{
		T[i] = (double *)calloc(i + 1, sizeof(double));
		if(T[i] == NULL)
		{
			fprintf(stderr, "Error : out of memory in chebyshev().\n");
			return;
		}
	}

	c = (double *)calloc(n + 1, sizeof(double));
	func = (double *)calloc(n + 1, sizeof(double));
	x = (double *)calloc(n + 1, sizeof(double));
	if(c == NULL || func == NULL || x == NULL)
	{
		fprintf(stderr, "Error : out of memory in chebyshev().\n");
		return;
	}

	/*	Chebyshev Polynomial	*/
	T[0][0] = 1;
	T[1][0] = 0;
	T[1][1] = 1;
	for(i = 2; i <= n; i++)
	{
		for(j = 1; j <= i; j++)	T[i][j] = 2. * T[i - 1][j - 1];
		for(j = 0; j < i - 1; j++)	T[i][j] -= T[i - 2][j];
	}

	k1 = (max - min) / 2.;
	k2 = (min + max) / 2.;
	n1 = n + 1;

	w = 0.5 * M_PI / n1;
	for(i = 0; i <= n; i++)
	{
		x[i] = cos((2 * i + 1) * w);
		func[i] = _f(k1 * x[i] + k2);
	}

	w = 0.;
	for(j = 0; j <= n; j++)	w += func[j];
	c[0] = w / n1;

	for(i = 1; i <= n; i++)
	{
		w = 0.;
		for(j = 0; j <= n; j++)	w += func[j] * polynomial(T[i], x[j], i);
		c[i] = 2. * w / n1;
	}

	for(i = 0; i <= n; i++)	a[i] = 0;
	for(i = 0; i <= n; i++)
	{
		w = c[i];
		for(j = 0; j <= i; j++)	a[j] += w * T[i][j];
	}

	if(min != -1. || max != 1.)
	{
		for(i = 1; i <= n; i++)	for(j = 1; j <= i; j++)	a[i] /= k1;

		if(k2 != 0.)

		T[0][0] = 1;				/* Pascal's Triangle */
		for(i = 1; i <= n; i++)
		{
			T[i][0] = 1;
			T[i][i] = 1;
			for(j = 1; j < i; j++)	T[i][j] = T[i - 1][j] + T[i - 1][j - 1];
		}

		for(i = 0; i <= n; i++)
		{
			c[i] = a[i];
			w = -k2 * a[i];
			for(j = i - 1; j >= 0; j--)
			{
				c[j] += T[i][j] * w;
				w *= -k2;
			}
		}
		for(i = 0; i <= n; i++)	a[i] = c[i];
	}
	free((char *)c);
	free((char *)func);
	free((char *)x);
	for(i = 0; i <= n; i++)	free((char *)T[i]);
	free((char *)T);
	return;
}
