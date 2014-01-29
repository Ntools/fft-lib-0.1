/*		etc.c		*/
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include "sslib.h"

#define		uchar		unsigned char

int strcomp(uchar *a, uchar *b);

double normal(double av, double st)
/* av:平均値 , st:標準偏差 */
{
	float r1, r2, z1, z2;

	r1 = rnd();
	r2 = rnd();
	z1 = sqrt(-2. * log(r1));
	z2 = sin(6.283185307179586 * r2);
	return st * z1 * z2 + av;
}

double rnd(void)	/*  0  <=  rnd()  <  1.0  */
{
	return (double)rand() / RAND_MAX;
}

int poison(double av)
/* av:平均値	*/
{
	double c, x, w;
	int rn;

	if(av <= 0.)
	{
		fprintf(stderr, "Error : av <= 0.  in poison()\n");
		return 0.;
	}
	if(av >= 10.)
	{
		fprintf(stderr, "Warning : av >= 10.  in poison(). Use normal Random number!\n");
		return 0.;
	}
	c = exp(-av);
	x = 1.;
	rn = -1;
	while(x > c)
	{
		w = rnd();
		x *= w;
		rn++;
	}
	return rn;
}

double bino(int m, int n)
{
	double a[3], w1, w2;
	int k[3], i, j, l;

	if(m < n || n < 0)
	{
		fprintf(stderr, "Error : Illegal parameter  in bino()\n");
		return 0.;
	}
	if(m == n || n == 0)	return 1.;
	k[0] = m;
	k[1] = n;
	k[2] = m - n;
	for(i = 0; i < 3; i++)
	{
		w1 = 0.;
		l = k[i];
		for(j = 2; j <= l; j++)
		{
			w2 = l + 2 - j;
			w1 += log(w2);
		}
		a[i] = w1;
	}
	return exp(a[0] - a[1] - a[2]);
}

/*					逆ラプラス変換					*/

void filt(double ts, double te, double *f, int nt, int n, int m, double *a1, double *b1, double a, double k1, double k2, int p)
{
	double *am, *ap, expa, *ff, pit, *q, *r, sg, *t, td, w, x, y;
	double sr, si, z1r, z1i, z2r, z2i;
	int i, i1, j, k, l;

	if(ts >= te || nt < 1 || n < 0 || m < 0)
	{
		fprintf(stderr, "Error : Illegal parameter  in filt()\n");
		return;
	}
	t = (double *)malloc(nt * sizeof(double));
	if(t == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in filt()\n");
		return;
	}
	ap = (double *)malloc((p + 1) * sizeof(double));
	if(ap == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in filt()\n");
		free((char *)t);
		return;
	}

/* calc Euler coeffs */
	x = 1.;
	r = ap + p;
	*r-- = x;
	for(y = (double)p, w = 1.; y > 0; y--, w++, r--)
		*r = *(r + 1) + (x *= (y / w));

	expa = exp(a) / *ap;

	td = (te - ts) / (double)nt;
	for(i = 0, q = t; i < nt; i++)	*q++ = ts + td * (double)i;
	q = t;
	ff = f;
	if(*q == 0.)
	{
		*ff++ = 0.;
		i1 = 1;
		q++;
	}
	else	i1 = 0;
	for(i = i1; i < nt; i++)
	{
		k = k1 + *q * k2;
		pit = M_PI / *q;
		sr = a / *q;
		si = pit / 2.;
		sg = -1.;
		w = 0.;

		for(j = 1; j <= k; j++)
		{
			am = a1 + n - 1;
			z1r = *am--;
			z1i = 0.;
			for(l = n; l > 1; l--)
			{
				x   = sr * z1r - si * z1i;
				z1i = sr * z1i + si * z1r;
				z1r = x + *am--;
			}
			am = b1 + m - 1;
			z2r = *am--;
			z2i = 0.;
			for(l = m; l > 1; l--)
			{
				x   = sr * z2r - si * z2i;
				z2i = sr * z2i + si * z2r;
				z2r = x + *am--;
			}
			w += sg * (z1i * z2r - z1r * z2i) / (z2r * z2r + z2i * z2i);
			sg = - sg;
			si += pit;
		}

/* Euler translation */
		r = ap;
		w *= *r++;
		for(j = 1; j <= p; j++)
		{
			am = a1 + n - 1;
			z1r = *am--;
			z1i = 0.;
			for(l = n; l > 1; l--)
			{
				x   = sr * z1r - si * z1i;
				z1i = sr * z1i + si * z1r;
				z1r = x + *am--;
			}
			am = b1 + m - 1;
			z2r = *am--;
			z2i = 0.;
			for(l = m; l > 1; l--)
			{
				x   = sr * z2r - si * z2i;
				z2i = sr * z2i + si * z2r;
				z2r = x + *am--;
			}
			w += sg * (z1i * z2r - z1r * z2i) / (z2r * z2r + z2i * z2i) * *r++;
			sg = - sg;
			si += pit;
		}
		*ff++ = w * expa / *q++;
	}
	free((char *)t);
	free((char *)ap);
	return;
}

double kaijo1(int n)
{
	int i;
	double k;

	if(n < 0)
	{
		fprintf(stderr, "Error : n < 0  in kaijo1()\n");
		return 0.;
	}
	if(n > 170)
	{
		fprintf(stderr, "Error : Too large n (n > 170)  in kaijo1()\n");
		return MAXDOUBLE;
	}
	k = 1.;
	if(n <= 1)	return k;
	for(i = 2; i <= n; i++)	k *= (double)i;
	return k;
}

double kaijo2(int n)
{
	if(n < 0)
	{
		fprintf(stderr, "Error : n < 0  in kaijo2()\n");
		return 0.;
	}
	if(n > 170)
	{
		fprintf(stderr, "Error : Too large n (n > 170)  in kaijo2()\n");
		return MAXDOUBLE;
	}
	return pow(10., log_kai2(n));
}

double log_kai1(int n)
{
	int i;
	double k;

	if(n < 0)
	{
		fprintf(stderr, "Error : n < 0  in log_kai1()\n");
		return 0.;
	}
	k = 0.;
	if(n <= 1)	return k;
	for(i = 2; i <= n; i++)	k += log10((double)i);
	return k;
}

double log_kai2(int n)
{
	static double lp = 0.39908993417905751;
	double le, lf, lnn, nn;

	if(n < 0)
	{
		fprintf(stderr, "Error : n < 0  in log_kai2()\n");
		return 0.;
	}
	nn = (double)n;
	le = nn * M_LOG10E;
	lf = log10((((-571. / 2488320. / nn - 139. / 51840.) / nn + 1. / 288.) / nn
				 + 1. / 12.) / nn + 1.);
	lnn = (nn + 0.5) * log10(nn);
	return lnn - le + lp + lf;
}

void minmax_i(int a[], int n, int *min, int *max)
{
	int i;

	if(n == 1)
	{
		*min = *max = a[0];
		return;
	}
	if(a[0] == a[1])	*min = *max = a[0];
	else if(a[0] < a[1])
	{
		*min = a[0];
		*max = a[1];
	}
	else
	{
		*min = a[1];
		*max = a[0];
	}
	for(i = n - 1; i > 2; i -= 2)
	{
		if(a[i] >= a[i - 1])
		{
			if(*max < a[i])		*max = a[i];
			if(*min > a[i - 1])	*min = a[i - 1];
		}
		else
		{
			if(*max < a[i - 1])	*max = a[i - 1];
			if(*min > a[i])		*min = a[i];
		}
	}
	if(n % 2)
		if(*max < a[2])	*max = a[2];
		if(*min > a[2])	*min = a[2];
}

void minmax_d(double a[], int n, double *min, double *max)
{
	int i;

	if(n == 1)
	{
		*min = *max = a[0];
		return;
	}
	if(a[0] == a[1])	*min = *max = a[0];
	else if(a[0] < a[1])
	{
		*min = a[0];
		*max = a[1];
	}
	else
	{
		*min = a[1];
		*max = a[0];
	}
	for(i = n - 1; i > 2; i -= 2)
	{
		if(a[i] >= a[i - 1])
		{
			if(*max < a[i])		*max = a[i];
			if(*min > a[i - 1])	*min = a[i - 1];
		}
		else
		{
			if(*max < a[i - 1])	*max = a[i - 1];
			if(*min > a[i])		*min = a[i];
		}
	}
	if(n % 2)
		if(*max < a[2])	*max = a[2];
		if(*min > a[2])	*min = a[2];
}

void minmax_c(uchar *a[], int n, uchar **min, uchar **max)
{
	int i, c;

	if(n == 1)
	{
		*min = *max = a[0];
		return;
	}
	c = strcomp(a[0], a[1]);
	if(c == 0)	*min = *max = a[0];
	else if(c < 0)
	{
		*min = a[0];
		*max = a[1];
	}
	else
	{
		*min = a[1];
		*max = a[0];
	}
	for(i = n - 1; i > 2; i -= 2)
	{
		if(strcomp(a[i], a[i - 1]) >= 0)
		{
			if(strcomp(*max, a[i]))		*max = a[i];
			if(strcomp(*min, a[i - 1]))	*min = a[i - 1];
		}
		else
		{
			if(strcomp(*max, a[i - 1]))	*max = a[i - 1];
			if(strcomp(*min, a[i]))		*min = a[i];
		}
	}
	if(n % 2)
		if(strcomp(*max, a[2]))	*max = a[2];
		if(strcomp(*min, a[2]))	*min = a[2];
}

int strcomp(uchar *a, uchar *b)
{
	uchar *p, *q;

	p = a;
	q = b;
	while(*p)
	{
		if(*p > *q)	return 1;
		if(*p++ < *q++)	return -1;
	}
	if(*q)	return -1;
	return 0;
}
