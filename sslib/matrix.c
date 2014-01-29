/*		matrix.c		*/
#include <stdio.h>
#include <stdlib.h>
#include "sslib.h"

void madd(double a[], double b[], double c[], int la, int lb, int lc, int m, int n)
{
	int i, j, ka, kb, kc;
	double *p, *q, *r;

	if(m <= 0 || n <= 0 || la < n || lb < n || lc < n )
	{
		fprintf(stderr, "Error : Illegal parameter  in madd()\n");
		return;
	}
	for(i = 0, ka = kb = kc = 0; i < m; i++, ka += la, kb += lb, kc += lc)
		for(j = 0, p = a + ka, q = b + kb, r = c + kc; j < n; j++)
			*r++ = *p++ + *q++;
	return;
}

void msub(double a[], double b[], double c[], int la, int lb, int lc, int m, int n)
{
	int i, j, ka, kb, kc;
	double *p, *q, *r;

	if(m <= 0 || n <= 0 || la < n || lb < n || lc < n )
	{
		fprintf(stderr, "Error : Illegal parameter  in msub()\n");
		return;
	}
	for(i = 0, ka = kb = kc = 0; i < m; i++, ka += la, kb += lb, kc += lc)
		for (j = 0, p = a + ka, q = b + kb, r = c + kc; j < n; j++)
			*r++ = *p++ - *q++;
	return;
}

void mmul1(double a[], double b[], int la, int lb, int m)
{
	int i, j, ka, kc, l;
	double *p, *q, *r, w, *work;

	if(m <= 0 || la < m || lb < m)
	{
		fprintf(stderr, "Error : Illegal parameter  in mmul1()\n");
		return;
	}
	work = (double *)malloc(m * m * sizeof(double));
	if(work == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in mmul1()\n");
		return;
	}
	for(i = 0, ka = 0, r = work; i < m; i++, ka += la)
	{
		for(l = 0; l < m; l++)
		{
			w = 0.;
			for(j = 0, p = a + ka, q = b + l; j < m; j++)
			{
				w += *p++ * *q;
				q += lb;
			}
			*r++ = w;
		}
	}
	mmove(work, a, m, la, m, m);
	free((char *)work);
	return;
}

void mmul2(double a[], double b[], double c[], int la, int lb, int lc, int m, int n, int k)
{
	int i, j, ka, kc, l;
	double *p, *q, *r, w;

	if(m <= 0 || n <= 0 || k <= 0 || la < n || lb < k || lc < k )
	{
		fprintf(stderr, "Error : Illegal parameter  in mmul2()\n");
		return;
	}
	for(i = 0, ka = kc = 0; i < m; i++, ka += la, kc += lc)
	{
		for(l = 0, r = c + kc; l < k; l++)
		{
			w = 0.;
			for(j = 0, p = a + ka, q = b + l; j < n; j++)
			{
				w += *p++ * *q;
				q += lb;
			}
			*r++ = w;
		}
	}
}

void mtra1(double a[], int l, int m, int n)
{
	int i, j, ka, mn;
	double *p, *q;

	if(m <= 0 || n <= 0 || l < m || l < n || m < n)
	{
		fprintf(stderr, "Error : Illegal parameter  in mtra1()\n");
		return;
	}
	mn = m;
	if(mn < n)	mn = n;
	for(i = 0, ka = 0; i < mn - 1; i++, ka += l)
		for(j = i + 1, p = a + ka + j, q = p + l - 1; j < mn; j++, q += l)
			swapd(p++, q);
}

void mtra2(double a[], double b[], int la, int lb, int m, int n)
{
	int i, j, ka;
	double *p, *q;

	if(m <= 0 || n <= 0 || la < n || lb < m)
	{
		fprintf(stderr, "Error : Illegal parameter  in mtra2()\n");
		return;
	}
	for(i = 0, ka = 0; i < m; i++, ka += la)
		for(j = 0, p = a + ka, q = b + i; j < n; j++, q += lb)	*q = *p++;
}

double minver(double a[], int l, int m, double eps)
{
	int i, iw, j, k, *p, r, s, t, u, v, *work;
	double api, pivot, *q, *q1, w, w1, wmax;

	if(m < 2 || l < m || eps <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter  in minver()\n");
		return 0.;
	}
	work = (int *)malloc(m * sizeof(int));
	if(work == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in minver()\n");
		return 0.;
	}
	w1 = 1.;
	for(i = 0, p = work; i < m; i++)	*p++ = i;
	for(k = 0, u = 0; k < m; k++, u += l)
	{
		wmax = 0.;
		for(i = k; i < m; i++)
		{
			w = fabs(*(a + i * l + k));
			if(w > wmax)
			{
				wmax = w;
				r = i;
			}
		}
		api = fabs(pivot = *(a + r * l + k));
		if(api < eps)
		{
			fprintf(stderr, "Error : api < eps  in minver()\n");
			free((char *)work);
			return w1;
		}
		w1 *= pivot;
		v = r * l;
		if(r != k)
		{
			w1 = -w1;
			swapi(work + k, work + r);
			for(j = 0, q = a + u, q1 = a + v; j < m; j++)	swapd(q++, q1++);
		}
		for(i = 0, q = a + u; i < m; i++)	*q++ /= pivot;
		for(i = 0, v = 0; i < m; i++, v += l)
		{
			if(i != k)
			{
				s = v + k;
				w = *(a + s);
				if(w != 0.)
				{
					for(j = 0, q = a + u, q1 = a + v; j < m; j++, q++, q1++)
						if (j != k)	*q1 -= w * *q;
					*(a + s) = - w / pivot;
				}
			}
		}
		*(a + u + k) = 1. / pivot;
	}
	for(i = 0; i < m; i++)
	{
		for(;;)
		{
			k = *(work + i);
			if(k == i)	break;
			swapi(work + k, work + i);
			for(j = 0, u = 0; j < m; j++, u += l)	swapd(a + u + i, a + u + k);
		}
	}
	free((char *)work);
	return w1;
}

void mmove(double a[], double b[], int la, int lb, int m, int n)
{
	int i, j, ka, kb;
	double *p, *q;

	if(m <= 0 || n <= 0 || la < n || lb < n)
	{
		fprintf(stderr, "Error : Illegal parameter  in mmove()\n");
		return;
	}
	for(i = ka = kb = 0; i < m; i++, ka += la, kb += lb)
		for(j = 0, p = a + ka, q = b + kb; j < n; j++)	*q++ = *p++;
	return;
}

void mswap(double a[], double b[], int la, int lb, int m, int n)
{
	int i, j, ka, kb;
	double *p, *q;

	if(m <= 0 || n <= 0 || la < n || lb < n)
	{
		fprintf(stderr, "Error : Illegal parameter  in mswap()\n");
		return;
	}
	for(i = ka = kb = 0; i < m; i++, ka += la, kb += lb)
		for(j = 0, p = a + ka, q = b + kb; j < n; j++)	swapd(p++, q++);
	return;
}

int jacobi(double a[], double v[], int l, int m, int *nr, double eps)
{
	int n, i, j, r, c, k1, k2, s1, s2, s3, s5, s6;
	double wmax, w, a1, a2, a3, t1, ta, si, co;

	if(m < 2 || *nr < 1 || eps <= 0.0)	return 999;
	n = 0;
	for(i = 0; i < m - 1; i++)
	{
		k1 = i * l;
		for(j = i + 1; j < m; j++)	if(a[k1 + j] != a[j * l + i])	return 999;
	}
	for(i = 0; i < m; i++)
	{
		k1 = i * l;
		for(j = 0; j < m; j++)	if(i != j)	v[k1 + j] = 0.0;
		v[k1 + i] = 1.0;
	}
	while(1)
	{
		wmax = 0.0;
		for(i = 0; i < m; i++)
		{
			k1 = i * l;
			for(j = i + 1; j < m; j++)
			{
				w = fabs(a[k1 + j]);
				if(w > wmax)
				{
					wmax = w;
					r = i;
					c = j;
				}
			}
		}
		if(wmax <= eps)
		{
			*nr = n;
			return 0;
		}
		if(n >= *nr)
		{
			*nr = n;
			return 1;
		}
		k1 = r * l;
		k2 = c * l;
		s1 = k1 + r;
		s2 = k2 + c;
		s3 = k1 + c;
		n++;
		a1 = a[s1];
		a2 = a[s2];
		a3 = a[s3];
		t1 = fabs(a1 - a2);
		ta = 2.0 * a3 / (t1 + sqrt(t1 * t1 + 4.0 * a3 * a3));
		if(a1 < a2)	ta = - ta;
		co = sqrt(1.0 / (ta * ta + 1.0));
		si = ta * co;
		for(i = 0; i < m; i++)
		{
			s5 = i * l + r;
			s6 = i * l + c;
			w = v[s5];
			v[s5] = w * co + v[s6] * si;
			v[s6] = - w * si + v[s6] * co;
			if(i != r && i != c)
			{
				w = a[s5];
				a[s5] = w * co + a[s6] *si;
				a[s6] = - w * si + a[s6] * co;
				a[k1 + i] = a[s5];
				a[k2 + i] = a[s6];
			}
		}
		a[s1] = a1 * co * co + a2 * si * si + 2.0 * a3 * co * si;
		a[s2] = a1 + a2 - a[s1];
		a[s3] = 0.0;
		a[k2 + r] = 0.0;
	}
}
