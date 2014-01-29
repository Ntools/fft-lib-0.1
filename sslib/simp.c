/*	simp.c	*/
#include <stdio.h>
#include <stdlib.h>
#include "sslib.h"

double *h, *sp, *x, *y;
int *jun;

int h_malloc(int n, char *proname);
int sp_malloc(int n, char *proname);
int x_malloc(int n, char *proname);
int y_malloc(int n, char *proname);
int jun_malloc(int n, char *proname);

int h_malloc(int n, char *proname)
{
	h = (double *)calloc(n, sizeof(double));
	if(h == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in %s()\n", proname);
		return -1;
	}
	return 0;
}

int sp_malloc(int n, char *proname)
{
	sp = (double *)calloc(n, sizeof(double));
	if(sp == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in %s()\n", proname);
		return -1;
	}
	return 0;
}

int x_malloc(int n, char *proname)
{
	x = (double *)calloc(n, sizeof(double));
	if(x == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in %s()\n", proname);
		return -1;
	}
	return 0;
}

int y_malloc(int n, char *proname)
{
	y = (double *)calloc(n, sizeof(double));
	if(y == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in %s()\n", proname);
		return -1;
	}
	return 0;
}

int jun_malloc(int n, char *proname)
{
	jun = (int *)calloc(n, sizeof(int));
	if(jun == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in %s()\n", proname);
		return -1;
	}
	return 0;
}

double lagdif(double xd[], double yd[], int n, double xx)
{
	int flag, i, j, k;
	double sm0, sm1, sm2, w, xi;

	if(x_malloc(n, "lagdif"))	return 0.;
	if(y_malloc(n, "lagdif"))
	{
		free(x);
		return 0.;
	}
	for(i = 0, flag = 0; i < n - 1; i++)
	{
		if(xd[i] >= xd[i + 1])
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		if(jun_malloc(n, "lagdif"))
		{
			free(y);
			free(x);
			return 0.;
		}
		sortdi1(xd, n, jun);
		for(i = 0; i < n; i++)
		{
			x[i] = xd[jun[i]];
			y[i] = yd[jun[i]];
		}
		free(jun);
	}
	else
	{
		for(i = 0; i < n; i++)
		{
			x[i] = xd[i];
			y[i] = yd[i];
		}
	}

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in lagdif()\n");
		free(y);
		free(x);
		return 0.;
	}
	if(xx < x[0] || xx > x[n - 1])
	{
		fprintf(stderr, "Error : xx is out of range in lagdif()\n");
		free(y);
		free(x);
		return 0.;
	}

	w = 0.;
	for(i = 0; i < n; i++)
	{
		sm0 = 0.;
		sm2 = 1.;
		xi = x[i];
		for(j = 0; j < n; j++)
		{
			if(i != j)
			{
				sm2 *= (xi - x[j]);
				for(k = 0, sm1 = 1.; k < n; k++)
					if((k != i) && (k != j))	sm1 *= (xx - x[k]);
				sm0 += sm1;
			}
		}
		w += (y[i] * sm0 / sm2);
	}
	free(y);
	free(x);
	return w;
}

double spldif(double xd[], double yd[], int n, double xx)
{
	int flag, i;
	double diff, dxp2, dxm2, hi, hi2;

	if(h_malloc(n, "spldif"))	return 0;
	if(sp_malloc(n, "spldif"))
	{
		free(h);
		return 0.;
	}
	if(x_malloc(n, "spldif"))
	{
		free(sp);
		free(h);
		return 0.;
	}
	if(y_malloc(n, "spldif"))
	{
		free(x);
		free(sp);
		free(h);
		return 0.;
	}
	for(i = 0, flag = 0; i < n - 1; i++)
	{
		if(xd[i] >= xd[i + 1])
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		if(jun_malloc(n, "spldif"))
		{
			free(y);
			free(x);
			free(sp);
			free(h);
			return 0.;
		}
		sortdi1(xd, n, jun);
		for(i = 0; i < n; i++)
		{
			x[i] = xd[jun[i]];
			y[i] = yd[jun[i]];
		}
		free(jun);
	}
	else
	{
		for(i = 0; i < n; i++)
		{
			x[i] = xd[i];
			y[i] = yd[i];
		}
	}

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in spldif()\n");
		free(y);
		free(x);
		return 0.;
	}
	if(xx < x[0] || xx > x[n - 1])
	{
		fprintf(stderr, "Error : xx is out of range in spldif()\n");
		free(y);
		free(x);
		return 0.;
	}

	subspl(x, y, n - 1, h, sp);
	for(i = 1; i < n; i++)
	{
		if(x[i - 1] <= xx && xx <= x[i])
		{
			hi = h[i];
			hi2 = hi * hi;
			dxp2 = (x[i] - xx);
			dxp2 *= dxp2;
			dxm2 = (xx - x[i - 1]);
			dxm2 *= dxm2;
			diff = (sp[i - 1] * (hi2 - 3. * dxp2)
				 +  sp[i]     * (3. * dxm2 - hi2)
				 +  6. * (y[i] - y[i - 1])) / 6. / hi;
			free(y);
			free(x);
			free(sp);
			free(h);
			return diff;
		}
	}
	return 0;
}

void subspl(double x[], double y[], int n, double h[], double sp[])
{
	int i;
	double *dc, *du, g, hip1;

	dc = (double *)calloc(n + 1, sizeof(double));
	du = (double *)calloc(n + 1, sizeof(double));
	if(dc == NULL || du == NULL)
	{
		fprintf(stderr, "Error : out of memory in subspl().\n");
		return;
	}
	for(i = 1; i <= n; i++)	h[i] = x[i] - x[i - 1];
	for(i = 1; i < n; i++)
	{
		hip1 = x[i + 1] - x[i];
		sp[i] = 6.* ((y[i + 1] - y[i]) / h[i] / hip1 - (y[i] - y[i - 1]) / h[i] / h[i]);
		du[i] = h[i + 1] / h[i];
		dc[i] = 2. * (1. + du[i]);
	}
	dc[n] += 1.;
	dc[n - 1] += (h[n] / h[n - 1]);
	du[1] /= dc[1];
	sp[1] /= dc[1];
	for(i = 2; i <= n; i++)
	{
		g = 1. / (dc[i] - du[i - 1]);
		du[i] *= g;
		sp[i] = (sp[i] - sp[i - 1]) * g;
	}
	for(i = n - 1; i >= 1; i--)	sp[i] -= sp[i + 1] * du[i];
	sp[0] = sp[1];
	sp[n] = sp[n - 1];
	free(du);
	free(dc);
}

double trap(double xx[], double yy[], int n)
{
	int flag, i;
	double w;

	if(x_malloc(n, "trap"))	return 0.;
	if(y_malloc(n, "trap"))
	{
		free(x);
		return 0.;
	}
	for(i = 0, flag = 0; i < n - 1; i++)
	{
		if(xx[i] >= xx[i + 1])
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		if(jun_malloc(n, "trap"))
		{
			free(y);
			free(x);
			return 0.;
		}
		sortdi1(xx, n, jun);
		for(i = 0; i < n; i++)
		{
			x[i] = xx[jun[i]];
			y[i] = yy[jun[i]];
		}
		free(jun);
	}
	else
	{
		for(i = 0; i < n; i++)
		{
			x[i] = xx[i];
			y[i] = yy[i];
		}
	}

	if(n <= 1)
	{
		fprintf(stderr, "Error : n <= 1  in trap()\n");
		free(y);
		free(x);
		return 0.;
	}
	for(i = 1, w = 0.; i < n; i++)
		w += (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
	free(y);
	free(x);
	return w / 2.;
}

double simpei(double y[], int n, double h)
{
	int i;
	double w;

	if(n <= 2 || n % 2 == 0 || h <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter  in simpei()\n");
		return 0.;
	}
	w = - y[0] + y[n - 1];
	for(i = 0; i < n - 1; i += 2)	w += (2. * y[i] + 4. * y[i + 1]);
	return w * h / 3.;
}

double simpui(double xx[], double yy[], int n)
{
	int flag, i;
	double d, h, hpd, hmd, w, ww;

	if(x_malloc(n, "simpui"))	return 0.;
	if(y_malloc(n, "simpui"))
	{
		free(x);
		return 0.;
	}
	for(i = 0, flag = 0; i < n - 1; i++)
	{
		if(xx[i] >= xx[i + 1])
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		if(jun_malloc(n, "simpui"))
		{
			free(y);
			free(x);
			return 0.;
		}
		sortdi1(xx, n, jun);
		for(i = 0; i < n; i++)
		{
			x[i] = xx[jun[i]];
			y[i] = yy[jun[i]];
		}
		free(jun);
	}
	else
	{
		for(i = 0; i < n; i++)
		{
			x[i] = xx[i];
			y[i] = yy[i];
		}
	}

	if(n <= 2 || n % 2 == 0)
	{
		fprintf(stderr, "Error : Illegal parameter  in simpui()\n");
		free(y);
		free(x);
		return 0.;
	}
	for(i = 2, w = 0.; i < n; i += 2)
	{
		h = (x[i] - x[i - 2]) / 2.;
		hpd = 1. / (x[i - 1] - x[i - 2]);
		d = x[i - 1] - x[i - 2] - h;
		hmd = 1. / (x[i] - x[i - 1]);
		ww = (1. + 2. * d * hpd) * y[i - 2];
		ww += 2. * h * (hpd + hmd) * y[i - 1];
		ww += (1. - 2. * d * hmd) * y[i];
		w += ww * h;
	}
	free(y);
	free(x);
	return w / 3.;
}

double splitg(double xx[], double yy[], int n)
{
	int flag, i;
	double hi, hi2, w;

	if(h_malloc(n, "splitg"))	return 0.;
	if(sp_malloc(n, "splitg"))
	{
		free(h);
		return 0.;
	}
	if(x_malloc(n, "splitg"))
	{
		free(sp);
		free(h);
		return 0.;
	}
	if(y_malloc(n, "splitg"))
	{
		free(x);
		free(sp);
		free(h);
		return 0.;
	}
	for(i = 0, flag = 0; i < n - 1; i++)
	{
		if(xx[i] >= xx[i + 1])
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		if(jun_malloc(n, "splitg"))
		{
			free(y);
			free(x);
			free(sp);
			free(h);
			return 0.;
		}
		sortdi1(xx, n, jun);
		for(i = 0; i < n; i++)
		{
			x[i] = xx[jun[i]];
			y[i] = yy[jun[i]];
		}
		free(jun);
	}
	else
	{
		for(i = 0; i < n; i++)
		{
			x[i] = xx[i];
			y[i] = yy[i];
		}
	}

	if(n <= 2 || n % 2 == 0)
	{
		fprintf(stderr, "Error : Illegal parameter  in simpui()\n");
		free(y);
		free(x);
		free(sp);
		free(h);
		return 0.;
	}
	subspl(x, y, n - 1, h, sp);

	for(i = 1, w = 0.; i < n; i++)
	{
		hi = h[i];
		hi2 = hi * hi * 0.0833333333333333333;
		w += hi * (y[i] + y[i - 1] - hi2 * (sp[i - 1] + sp[i]));
	}
	free(y);
	free(x);
	free(sp);
	free(h);
	return w / 2.;
}

double simpe2(double y[], int m, int n, double h1, double h2)
{
	int i, ij, j, m1, n1;
	double v;

	if(x_malloc(m, "simpe2"))	return 0.;
	if(m <= 2 || n <= 2 || m % 2 == 0 || n % 2 == 0)
	{
		fprintf(stderr, "Error : Illegal parameter in simpe2()\n");
		free(x);
		return 0.;
	}

	m1 = m - 1;
	n1 = n - 1;
	for(i = 0; i < m; i++)
	{
		v = - y[i * n] + y[i * n + n1];
		for(j = 0, ij = i * n; j < n1 - 1; j += 2, ij += 2)
			v += (2 * y[ij] + 4 * y[ij + 1]);
		x[i] = v;
	}
	v = - x[0] + x[m1];
	for(i = 0; i < m1 - 1; i += 2)	v += (2 * x[i] + 4 * x[i + 1]);
	free(x);
	return v * h1 * h2 / 9.;
}
