/*		linear.c		*/
#include <stdio.h>
#include <stdlib.h>
#include "sslib.h"

void gausei(double a[], int l, int m, int iter, double eps, double x[])
{
	int i, j, k;
	double ay, def, *p, *q, *r, sum, w, y;

	if(l < m || m < 2 || iter < 1 || eps <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter in gausei()\n");
		return;
	}
	for(i = 0, p = x; i < m; i++)	*p++ = 0.;
	k = 1;
	do
	{
		def = 0.;
		for(i = 0, p = x; i < m; i++)
		{
			w = *(a + i * l + i);
			sum = 0.;
			for(j = 0, q = x, r = a + i * l; j < m; j++, r++, q++)
				if (i != j)	sum += *r * *q;
			y = (*(a + i * l + m) - sum) / w;
			ay = fabs(y - *p);
			if(ay > def)	def = ay;
			*p++ = y;
		}
		if(def <= eps)	return;
	} while(k++ < iter);
	fprintf(stderr, "Error : No convergence in gausei()\n");
	return;
}

void gauss(double a[], int l, int m, int n, double eps)
{
	int i, iw, j, k, kp1, mm1, *p, r, s, t, u, v, *work;
	double api, pivot, *q, *q1, w, wmax;

	if(l < n || m < 2 || m >= n || eps <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter in gauss()\n");
		return;
	}
	work = (int *)malloc(m * sizeof(int));
	if(work == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in gauss()\n");
		return;
	}

	for(i = 0, p = work; i < m; i++)	*p++ = i;
	for(k = 0; k < m; k++)
	{
		wmax = 0.;
		for(i = k, q = a + i * l + k; i < m; i++)
		{
			w = fabs(*q++);
			if(w > wmax)
			{
				wmax = w;
				r = i;
			}
		}
		u = r * l;
		v = k * l;
		s = u + k;
		pivot = *(a + r * l + k);
		api = fabs(pivot);
		if(api < eps)
		{
			fprintf(stderr, "Error : pivot < eps in gauss()\n");
			free((char *)work);
			return;
		}
		if(r != k)
		{
			iw = *(work + k);
			*(work + k) = *(work + r);
			*(work + r) = iw;
			for(j = 0; j < n; j++)
			{
				s = k * l + j;
				t = r * l + j;
				w = *(a + s);
				*(a + s) = *(a + t);
				*(a + t) = w;
			}
		}
		for(i = k, q = a + v + k; i < n; i++)	*q++ /= pivot;
		kp1 = k + 1;
		if(kp1 > m)	break;
		for(i = kp1; i < m; i++)
		{
			w = *(a + i * l + k);
			if(w != 0.)
			{
				q = a + i * l + kp1;
				q1 = a + k * l + kp1;
				for(j = kp1; j < n; j++)	*q++ -= w * *q1++;
			}
		}
	}
	mm1 = m - 1;
	for(k = m; k < n; k++)
	{
		for(i = 0, iw = m - 2; i < mm1; i++, iw--)
		{
			w = 0.;
			for(j = m - 1, q = a + iw * l + j; j > iw; j--)
				w += *q-- * *(a + j * l + k);
			*(a + iw * l + k) -= w;
		}
	}
	free((char *)work);
	return;
}

void gaujor(double a[], int l, int m, int n, double eps)
{
	int c, i, iw, j, k, kp1, mm1, r, s, t, u, v, *work, x;
	double api, pivot, w, wmax;

	if(l < n || m < 2 || m >= n || eps <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter in gaujor()\n");
		return;
	}
	work = (int *)malloc(m * sizeof(int));
	if(work == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in gaujor()\n");
		return;
	}

	for(i = 0; i < m; i++)	work[i] = i;
	for(k = 0; k < m; k++)
	{
		x = k * l;
		wmax = 0.;
		for(i = 0; i < m; i++)
		{
			s = i * l;
			for(j = k; j < m; j++)
			{
				w = fabs(a[s + j]);
				if(w > wmax)
				{
					wmax = w;
					r = i;
					c = j;
				}
			}
		}
		s = r * l;
		pivot = a[s + c];
		api = fabs(pivot);
		if(api < eps)
		{
			fprintf(stderr, "Error : api < eps  in gaujor()\n");
			free(work);
			return;
		}
		if(c != k)
		{
			iw = work[k];
			work[k] = work[c];
			work[c] = iw;
			for(i = 0; i < m; i++)
			{
				t = i * l;
				u = t + k;
				v = t + c;
				w = a[u];
				a[u] = a[v];
				a[v] = w;
			}
		}
		if(r != k)
		{
			for(j = k; j < n; j++)
			{
				t = s + j;
				u = x + j;
				w = a[t];
				a[t] = a[u];
				a[u] = w;
			}
		}
		kp1 = k + 1;
		for(i = kp1; i < n; i++)	a[x + i] /= pivot;
		for(i = 0; i < m; i++)
		{
			t = i * l;
			if(i != k)
			{
				w = a[t + k];
				for(j = kp1; j < n; j++)	a[t + j] -= w * a[x + j];
			}
		}
	}
	mm1 = m - 1;
	for(i = 0; i < mm1; i++)
	{
		for(;;)
		{
			k = work[i];
			if(k == i)	break;
			iw = work[k];
			work[k] = work[i];
			work[i] = iw;
			s = k * l;
			t = i * l;
			for(j = m; j < n; j++)
			{
				u = s + j;
				v = t + j;
				w = a[u];
				a[u] = a[v];
				a[v] = w;
			}
		}
	}
	free(work);
	return;
}

int ludcmp(double aa[], double b[], double x[], double eps, int n)
{
	int i, j, k, n1, n2;
	double w, *y, *a;

	if(n < 2 || eps <= 0.)
	{
		fprintf(stderr, "Error : eps <= 0  in ludcmp()\n");
		return -1;
	}
	n1 = n + 1;
	n2 = n1 * n1;
	a = (double *)malloc(n2 * sizeof(double));
	if(a == NULL)
	{
		fprintf(stderr, "Error : Out of Memory in ludcmp()\n");
		return 999;
	}
	y = (double *)malloc(n1 * sizeof(double));
	if(y == NULL)
	{
		fprintf(stderr, "Error : Out of Memory in ludcmp()\n");
		free((char *)a);
		return 999;
	}
	for(i = 0; i < n2; i++)	*(a + i) = *(aa + i);
	for(i = 0; i < n; i++)
	{
		if(fabs(*(a + i + i * n1)) < eps)
		{
			fprintf(stderr, "Error : a[%d][%d] < eps in ludcmp()\n", i, i);
			free((char *)a);
			free((char *)y);
			return 1;
		}
		for(j = i + 1; j <= n; j++)
		{
			w = *(a + i + j * n1);
			if(i != 0)
			{
				for(k = 0; k < i; k++)
					w -= (*(a + i + k * n1) * *(a + k + j * n1));
			}
			*(a + i + j * n1) = w / *(a + i + i * n1);
		}
		for(j = i + 1; j <= n; j++)
		{
			w = *(a + j + (i + 1) * n1);
			for(k = 0; k <= i; k++)
				w -= (*(a + k + (i + 1) * n1) * *(a + j + k * n1));
			*(a + j + (i + 1) * n1) = w;
		}
	}

	*y = *b;
	for(i = 1; i <= n; i++)
	{
		w = *(b + i);
		for(j = 0; j < i; j++)	w -= (*(a + j + i * n1) * *(y + j));
		*(y + i) = w;
	}
	*(x + n) = *(y + n) / *(a + n + n * n1);
	for(i = n - 1; i >= 0; i--)
	{
		w = *(y + i);
		for(j = i + 1; j <= n; j++)
			w -= (*(a + j + i * n1) * *(x + j));
		*(x + i) = w / *(a + i + i * n1);
	}
	free((char *)a);
	free((char *)y);
	return(0);
}

/*		conjugate gradient method		*/
int cgm(double a[], double c[], double x[], int n)
{
	double *r1, *r2, *p1, *p2, *w2;
	double ai, atr1, bi, ww;
	int i, j, k = 1;

	if(n < 2)
	{
		fprintf(stderr, "Error : illegal parameter input in cgm()\n");
		return -1;
	}
	r1 = (double *)calloc(n, sizeof(double));
	r2 = (double *)calloc(n, sizeof(double));
	p1 = (double *)calloc(n, sizeof(double));
	p2 = (double *)calloc(n, sizeof(double));
	w2 = (double *)calloc(n, sizeof(double));
	if(r1 == NULL || r2 == NULL || p1 == NULL || p2 == NULL || w2 == NULL)
	{
		fprintf(stderr, "Error : Out of Memory in cgm()\n");
		return 999;
	}

	for(i = 0; i < n; i++)
	{
		x[i] = 0;
		r1[i] = c[i];
		ww = 0;
		for(j = 0; j < n; j++)	ww += a[j * n + i] * c[j];
		p1[i] = ww;
	}

	for(;;)
	{
		atr1 = 0;
		for(i = 0; i < n; i++)
		{
			ww = 0;
			for(j = 0; j < n; j++)	ww += a[j * n + i] * r1[j];
			atr1 += ww * ww;
		}
		bi = 0;
		for(i = 0; i < n; i++)
		{
			ww = 0;
			for(j = 0; j < n; j++)	ww += a[i * n + j] * p1[j];
			w2[i] = ww;
			bi += ww * ww;
		}
		ai = atr1 / bi;
		for(i = 0; i < n; i++)	x[i] += ai * p1[i];
		if(k == n)	break;

		for(i = 0; i < n; i++)	r2[i] = r1[i] - ai * w2[i];
		bi = 0;
		for(i = 0; i < n; i++)
		{
			ww = 0;
			for(j = 0; j < n; j++)	ww += a[j * n + i] * r2[j];
			p2[i] = ww;
			bi += ww * ww;
		}
		bi /= atr1;
		for(i = 0; i < n; i++)
		{
			p1[i] = p2[i] + bi * p1[i];
			r1[i] = r2[i];
		}
		k++;
	}

	free((char *)r1);
	free((char *)r2);
	free((char *)p1);
	free((char *)p2);
	free((char *)w2);
	return 0;
}

int trdiam(double al[], double am[], double au[], double b[], int n, double x[])
{
	double w;
	int i;

	if(n < 1 || am[0] == 0.0)	return 999;
	au[0] /= am[0];
	x[0] = b[0] / am[0];
	for(i = 1; i <= n; i++)
	{
		w = 1.0 / (am[i] - au[i - 1] * al[i]);
		au[i] *= w;
		x[i] = (b[i] - x[i - 1] * al[i]) * w;
	}
	for(i = n - 1; i >= 0; i--)	x[i] -= x[i + 1] * au[i];
	return 0;
}

void cgauj(Complex a[], int l, int m, int n, double eps)
{
	int c, i, iw, j, k, mm1, *p, r, *work;
	double api, w, wmax;
	Complex pivot, *q, *q1, wk, wk1;

	if(l < n || m < 2 || m >= n || eps <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter in cgauj()\n");
		return;
	}
	work = (int *)malloc(m * sizeof(int));
	if(work == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in cgauj()\n");
		return;
	}

	for(i = 0, p = work; i < m; i++)	*p++ = i;
	for(k = 0; k < m; k++)
	{
		wmax = 0.;
		for(i = k; i < m; i++)
		{
			for(j = k, q = a + i * l + k; j < m; j++)
			{
				w = cabslt(*q++);
				if(w > wmax)
				{
					wmax = w;
					r = i;
					c = j;
				}
			}
		}
		pivot = *(a + r * l + c);
		api = cabslt(pivot);
		if(api < eps)
		{
			fprintf(stderr, "Error : api < eps  in cgauj()\n");
			free((char *)work);
			return;
		}
		if(c != k)
		{
			iw = *(work + k);
			*(work + k) = *(work + c);
			*(work + c) = iw;
			for(i = 0, q = a + i * l; i < m; i++, q += l)
			{
				wk = *(q + k);
				*(q + k) = *(q + c);
				*(q + c) = wk;
			}
		}
		if(r != k)
		{
			for(j = k; j < n; j++)
			{
				wk = *(a + r * l + j);
				*(a + r * l + j) = *(a + k * l + j);
				*(a + k * l + j) = wk;
			}
		}

		*(a + k * l + k) = tocomplex(1., 0.);
		for(i = k + 1, q = a + k * l + k + 1; i < n; i++, q++)
			*q = cdiv(*q, pivot);
		for(i = 0; i < m; i++)
		{
			if(i != k)
			{
				wk = *(a + i * l + k);
				for(j = k; j < n; j++)
				{
					wk1 = cmul(wk, *(a + k * l + j));
					*(a + i * l + j) = csub(*(a + i * l + j), wk1);
				}
			}
		}
	}
	mm1 = m - 1;
	for(i = 0; i < mm1; i++)
	{
		for(;;)
		{
			k = *(work + i);
			if(k == i)	break;
			iw = *(work + k);
			*(work + k) = *(work + i);
			*(work + i) = iw;
			for(j = m, q = a + k * l + m, q1 = a + i * l + m; j < n; j++)
			{
				wk = *q;
				*q++ = *q1;
				*q1++ = wk;
			}
		}
	}
	free((char *)work);
	return;
}
