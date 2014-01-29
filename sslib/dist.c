/*		dist.c		*/
#include <stdio.h>
#include "sslib.h"

double qnorm(double u)
{
	static double a[9] = {	 1.24818987e-4, -1.075204047e-3, 5.198775019e-3,
							-0.019198292004, 0.059054035642,-0.151968751364,
							 0.319152932694,-0.5319230073,   0.797884560593};
	static double b[15] = {	-4.5255659e-5,   1.5252929e-4,  -1.9538132e-5,
							-6.76904986e-4,  1.390604284e-3,-7.9462082e-4,
							-2.034254874e-3, 6.549791214e-3,-0.010557625006,
							 0.011630447319,-9.279453341e-3, 5.353579108e-3,
							-2.141268741e-3, 5.35310549e-4,  0.999936657524};
	double w, y, z;
	int i;

	if(u == 0.)	return 0.5;
	y = u / 2.;
	if(y < -6.)	return 0.;
	if(y > 6.)		return 1.;
	if(y < 0.)		y = - y;
	if(y < 1.)
	{
		w = y * y;
		z = a[0];
		for(i = 1; i < 9; i++)		z = z * w + a[i];
		z *= (y * 2.);
	}
	else
	{
		y -= 2.;
		z = b[0];
		for(i = 1; i < 15; i++)	z = z * y + b[i];
	}

	if(u < 0.)	return (1. - z) / 2.;
	return (1. + z) / 2.;
}

double pnorm(double qn)
{
	static double b[11] = {	 1.570796288,     0.03706987906,  -0.8364353589e-3,
							-0.2250947176e-3, 0.6841218299e-5, 0.5824238515e-5,
							-0.104527497e-5,  0.8360937017e-7,-0.3231081277e-8,
							 0.3657763036e-10,0.6936233982e-12};
	double w1, w3;
	int i;

	if(qn < 0. || 1. < qn)
	{
		fprintf(stderr, "Error : qn <= 0 or qn >= 1  in pnorm()!\n");
		return 0.;
	}
	if(qn == 0.5)	return 0.;

	w1 = qn;
	if(qn > 0.5)	w1 = 1. - w1;
	w3 = -log(4. * w1 * (1. - w1));
	w1 = b[0];
	for(i = 1; i < 11; i++)	w1 += (b[i] * pow(w3, (double)i));
	if(qn > 0.5)	return sqrt(w1 * w3);
	return -sqrt(w1 * w3);
}
 
double qchi(double x2, int n)
{
	double w, pw, x, qc;
	int i, i1;

	if(n < 1)
	{
		fprintf(stderr,"Error : Ž©—R“x < 1 in qchi()!\n");
		return 0.;
	}
	if(x2 <= 0.)	return 1.;
	if(x2 > 400.)	return 0.;
	if(n > 10)
	{
		w = 2. / (9. * (double)n);
		pw = pow(x2 / (double)n, 1. / 3.);
		return 1. - qnorm((pw - 1. + w) / sqrt(w));
	}
	w = exp(-x2 / 2.);
	if(n == 2)	return w;
	x = sqrt(x2);
	if(n == 1)	return 2. * (1. - qnorm(x));
	if((n % 2) == 0)
	{
		i1 = 2;
		qc = w;
	}
	else
	{
 		i1 = 1;
		qc = 2. * (1. - qnorm(x));
		w *= (0.797884560750774 / x);
	}
	for(i = i1; i <= n - 2; i += 2)
	{
		w *= (x2 / (double)i);
		qc += w;
	}
	return qc;
}

double pchi(double qc, int n)
{
	double c1, c2, gam, x, w, wx;
	int i, j;

	if(qc <= 0. || qc >= 1. || n < 1)
	{
		fprintf(stderr,"Error : Illigal parameter in pchi()!\n");
		return 0.;
	}
	if(n == 1)
	{
		w = pnorm(qc / 2.);
		return w * w;
	}
	if(n == 2)	return -2. * log(qc);

	x = -pnorm(qc);
	if(n > 10)
	{
		w = x * x;
		wx = sqrt(2. * (double)n);
		c1 = (w - 7.) * x / 9. / wx;
		c2 = ((3. * w + 7.) * w - 16.) * 2. / 405. / (double)n;
		wx = (double)n + wx * x + 2. * (w - 1.) / 3. + c1 - c2;
		if(wx < 0.)	return 0.;
		return wx;
	}

	w = 2. / 9. / (double)n;
	w = 1. - w + x * sqrt(w);
	wx = (double)n * w * w * w;
	if((n % 2) == 0)	gam = 1.;
	else				gam = 1.772453850905516;
	j = (n + 1) / 2 - 1;
	w = (double)n / 2.;
	for(i = 1; i <= j; i++)	gam *= (w - (double)i);
	x = wx / 2.;
	c1 = pow(x, w - 1.);
	c2 = exp(-x) / 2.;
	return wx + (qchi(wx, n) - qc) * gam / c1 / c2;
}

double qt(double t, int n)
{
	double t1, t2, w1, w2, w3, wq;

	if(n < 1)
	{
		fprintf(stderr,"Error : n < 1  in qt()!\n");
		return 0.;
	}

	w1 = 0.636619772284456;
	if(t < 0.)	t = - t;
	t1 = t / sqrt((double)n);
	t2 = 1. / (1. + t1 * t1);
	if((n % 2) != 0)
	{
		wq = 1. - w1 * atan(t1);
		if(n != 1)
		{
			wq -= (w2 = w1 * t1 * t2);
			if (n != 3)	qtsub(&wq, n, w2, 0., t2);
		}
		if(wq > 0.)	return wq;
		return 0.;
	}
	wq = 1. - (w2 = t1 * sqrt(t2));
	if(n != 2)		qtsub(&wq, n, w2, 1., t2);
	if(wq > 0.)	return wq;
	return 0.;
}

void qtsub(double *q, int n, double w2, double w3, double t2)
{
	int i, j;
	double w;

	j = (n - 2) / 2;
	for(i = 1; i <= j; i++)
	{
		w = 2. * (double)i - w3;
		*q -= (w2 *= (t2 * w / (1. + w)));
	}
}
   
double pt(double q, int n)
{
	double f, f1, f2, f3, f4, f5, u, u2, w, w0, w1, w2, w3, w4, x;
	int i, j, k;

	if(q < 1.e-5 || q > 1. || n < 1)
	{
		fprintf(stderr,"Error : Illigal parameter  in pt()!\n");
		return 0.;
	}

	if(n <= 5)	return ptsub(q, n);

	if(q <= 5.e-3 && n <= 13)	return ptsub(q, n);

	f1 = 4. * (f = (double)n);
	f5 = (f4 = (f3 = (f2 = f * f) * f) * f) * f;
	f2 *= 96.;
	f3 *= 384.;
	f4 *= 92160.;
	f5 *= 368640.;
	u = pnorm(1. - q / 2.);

	w0 = (u2 = u * u) * u;
	w1 = w0 * u2;
	w2 = w1 * u2;
	w3 = w2 * u2;
	w4 = w3 * u2;
	w = ((w0 + u) / f1);
	w += ((5. * w1 + 16. * w0 + 3. * u) / f2);
	w += ((3. * w2 + 19. * w1 + 17. * w0 - 15. * u) / f3);
	w += ((79. * w3 + 776. * w2 + 1482. * w1 - 1920. * w0 - 945. * u) / f4);
	w += ((27. * w4 + 339. * w3 + 930. * w2 - 1782. * w1 - 765. * w0
		+ 17955. * u) / f5);
	return u + w;
}

double ptsub(double q, int n)
{
	double eps, qe, s, w;

	if(n == 1 && 0.001 <= q && q < 0.01)	eps = 1.e-4;
	else if (n == 2 && q < 0.0001)			eps = 1.e-4;
	else if (n == 1 && q < 0.001)			eps = 1.e-2;
	else									eps = 1.e-5;
	s = 10000.;
	w = 0.;
	for(;;)
	{
		w += s;
		if(s <= eps)	return w;
		if((qe = qt(w,n) - q) == 0.)	return w;
		if(qe < 0.)
		{
			w -= s;
			s /= 10.;
		}
	}
} 

double qf(double f, int n1, int n2)
{
	double d, p, w, y, z;
	int i, j, k, m, mn, n;

	if(f < 0. || n1 < 1 || n2 < 1)
	{
		fprintf(stderr,"Error : Illegal parameter  in qf()!\n");
		return 0.;
	}

	w = f * (double)n1 / (double)n2;
	z = 1. / (1. + w);
	m = 2 * (n1 / 2) - n1 + 2;
	n = 2 * (n2 / 2) - n2 + 2;
	if((mn = 4 - 2 * (n1 % 2) - (n2 % 2)) == 1)
	{
		y = 0.318309886142228;
		p = sqrt(w);
		d = y * z / p;
		p = 2. * y * atan(p);
	}
	else if(mn == 2)
	{
		p = sqrt(z * w);
		d = p * z / w / 2.;
	}
	else if(mn == 3)
	{
		p = sqrt(z);
		d = z * p / 2.;
		p = 1. - p;
	}
	else
	{
		d = z * z;
		p = z * w;
	}
	if(n1 <= 2 && n2 <= 2)	return 1. - p;

	y = 2. * w / z;
	if((k = n + 2) <= n2)
	{
		for( ; k <= n2; k += 2)
		{
			d *= ((1. + (double)m / (double)(k - 2)) * z);
			if(m != 1)	p = (p + w) * z;
			else		p += (d * y / (double)(k - 1));
		}
	}
	k = m + 2;
	if(k > n1)	return 1. - p;

	y = z * w;
	z = 2. / z;
	j = n2 - 2;
	for( ; k <= n1; k += 2)
	{
		d *= (y * (double)(k + j) / (double)(k - 2));
		p -= (z * d / (double)(k + j));
	}
	return 1. - p;
}

double pf(double q, int n1, int n2)
{
	double a, b, c, d, eps, fw, f1, f2, qe, qn, s, u, u2, w1, w2, w3, w4;

	if(q < 0. || q > 1. || n1 < 1 || n2 < 1)
	{
		fprintf(stderr,"Error : Illegal parameter  in pf()!\n");
		return 0.;
	}

	if(n1 <= 240 || n2 <= 240)
	{
		eps = 1.e-5;
		if(n2 == 1)	eps = 1.e-4;
		s = 1000.;
		fw = 0.;
		for(;;)
		{
			fw += s;
			if(s <= eps)	return fw;
			if((qe = qf(fw,n1,n2) - q) == 0.)	return fw;
			if(qe < 0.)
			{
				fw -= s;
				s /= 10.;
			}
		}
	}

	eps = 1.e-6;
	qn = q;
	if(q < 0.5)	qn = 1. - q;
	u = pnorm(qn);
	w1 = 2. / (double)n1 / 9.;
	w2 = 2. / (double)n2 / 9.;
	w3 = 1. - w1;
	w4 = 1. - w2;
	u2 = u * u;
	a = w4 * w4 - u2 * w2;
	b = -2. * w3 * w4;
	c = w3 * w3 - u2 * w1;
	d = b * b - 4 * a * c;
	if(d < 0.)	fw = pfsub(a, b, 0.);
	else
	{
		if(fabs(a) > eps)	fw = pfsub(a, b, d);
		else
		{
			if(fabs(b) > eps)	return -c / b;
			fw = pfsub(a, b, 0.);
		}
	}
	return fw * fw * fw;
}

double pfsub(double x, double y, double z)
{
	return (sqrt(z) - y) / x / 2.;
}
