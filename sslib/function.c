/*		function.c		*/
#include <stdio.h>
#include "sslib.h"

double cbrt(double a)
{
	double	w;

	if(a == 0.)	return 0.;
	w = pow(fabs(a), 1./ 3.);
	if(a < 0.)	w = -w;
	return (w + 3. * a / (2. * w * w + a / w)) / 2.;
}

double besj0(double x)
{
	int i;
	double w, wp, wq, wx4;
	static double a[7] = {	 1.0,          -3.9999998721, 3.9999973021,
							-1.7777560599, 0.4443584263, -0.0709253492,
							 0.0076771853};
	static double c[5] = {	 0.3989422793, -0.00175302,   0.00017343,
							-0.0000487613,  0.0000173565};
	static double k[5] = {	-0.0124669441,  0.0004564324, -0.0000869791,
							 0.0000342468, -0.0000142078};

	if(x < 0.)
	{
		fprintf(stderr, "Error : x < 0 in besj0()\n");
		return 0.;
	}
	wx4 = x * x / 16.;
	if(x <= 4.)
	{
		w = -0.0005014415;
		for(i = 6; i >= 0; i--)	w = w * wx4 + a[i];
		return w;
	}
	wx4 = 1. / wx4;
	wp = -0.0000037043;
	wq =  0.0000032312;
	for(i = 4; i >= 0; i--)
	{
		wp = wp * wx4 + c[i];
		wq = wq * wx4 + k[i];
	}
	w = x - 0.78539816339744831;
	wp *= cos(w);
	wq *= sin(w);
	return 2. / sqrt(x) * (wp - 4. * wq / x);
}

double besj1(double x)
{
	int i;
	double w, wp, wq, wx4;
	static double a[7] = {	 1.9999999998, -3.999999971,   2.6666660544,
							-0.8888839649,  0.1777582922, -0.0236616773,
							 0.0022069155};
	static double c[5] = {	 0.3989422819,  0.0029218256, -0.000223203,
							 0.0000580759, -0.000020092};
	static double k[5] = {	 0.0374008364, -0.00063904,    0.0001064741,
							-0.0000398708,  0.00001622};

	if(x < 0.)
	{
		fprintf(stderr, "Error : x < 0 in besj1()\n");
		return 0.;
	}
	wx4 = x * x / 16.;
	if(x <= 4.)
	{
		w = -0.0001289769;
		for(i = 6; i >= 0; i--)	w = w * wx4 + a[i];
		return w * x / 4.;
	}
	wx4 = 1. / wx4;
	wp =  0.0000042414;
	wq = -0.0000036594;
	for(i = 4; i >= 0; i--)
	{
		wp = wp * wx4 + c[i];
		wq = wq * wx4 + k[i];
	}
	w = x - 2.356194490192345;
	wp *= cos(w);
	wq *= sin(w);
	return 2. / sqrt(x) * (wp - 4. * wq / x);
}

double besy0(double x)
{
	int i;
	double bj0, w, wf, wp, wx3;
	static double a[6] = {	 0.36746691,  0.60559366, -0.74350384,
							 0.25300117, -0.04261214,  0.00427916};
	static double b[6] = {	 1.0,        -2.2499997,   1.2656208,
							-0.3163866,   0.0444479,  -0.0039444};
	static double c[6] = {	 0.79788456, -0.00000077, -0.0055274,
							-0.00009512,  0.00137237, -0.00072805};
	static double k[6] = {	 0.78539816,  0.04166397,  0.00003954,
							-0.00262573,  0.00054125,  0.00029333};

	if(x <= 0.)
	{
		fprintf(stderr, "Error : x <= 0  in besy0()\n");
		return 0.;
	}
	if(x <= 3.)
	{
		w = -0.00024846;
		bj0 = 0.00021;
		wx3 = x * x / 9.;
		for(i = 5; i >= 0; i--)
		{
			w = w * wx3 + a[i];
			bj0 = bj0 * wx3 + b[i];
		}
		return w + 0.636619772367581 * log(x / 2.) * bj0;
	}
	wx3 = 3. / x;
	wf = 0.00014476;
	wp = -0.00013558;
	for(i = 5; i >= 0; i--)
	{
		wf = wf * wx3 + c[i];
		wp = wp * wx3 + k[i];
	}
	return wf * sin(x - wp) / sqrt(x);
}

double besy1(double x)
{
	int i;
	double bj1, w, wf, wp, wx3;
	static double a[6] = {	-0.6366198,  0.2212091,  2.1682709,
							-1.3164827,  0.3123951, -0.0400976};
	static double b[6] = {	 0.5,        -0.56249985,  0.21093573,
							-0.03954289,  0.00443319, -0.00031761};
	static double c[6] = {	 0.79788456,  0.00000156,  0.01659667,
							 0.00017105, -0.00249511,  0.00113653};
	static double k[6] = {	 0.78539816, -0.12499612, -0.0000565,
							 0.00637879, -0.00074348, -0.00079824};

	if(x <= 0.)
	{
		fprintf(stderr, "Error : x <= 0  in besy1()\n");
		return 0.;
	}
	if(x <= 3.)
	{
		w = 0.0027873;
		bj1 = 0.00001109;
		wx3 = x * x / 9.;
		for(i = 5; i >= 0; i--)
		{
			w = w * wx3 + a[i];
			bj1 = bj1 * wx3 + b[i];
		}
		return w / x + M_2_PI * log(x / 2.) * bj1 * x;
	}
	wx3 = 3. / x;
	wf = -0.00020033;
	wp = 0.00029166;
	for(i = 5; i >= 0; i--)
	{
		wf = wf * wx3 + c[i];
		wp = wp * wx3 + k[i];
	}
	return - wf * cos(x - wp) / sqrt(x);
}

double besi0(double x)
{
	int i;
	double w, wx375;
	static double a[6] = {	 1.0,        3.5156229,  3.0899424,
							 1.2067492,  0.2659732,  0.0360768};
	static double b[8] = {	 0.39894228,   0.013285917,  0.002253187,
							-0.001575649,  0.009162808, -0.020577063,
							 0.026355372, -0.016476329};

	if(x < 0.)
	{
		fprintf(stderr, "Error : x < 0  in besi0()\n");
		return 0.;
	}
	if(x <= 3.75)
	{
		wx375 = x * x / 14.0625;
		w = 0.0045813;
		for(i = 5; i >= 0; i--)	w = w * wx375 + a[i];
		return w;
	}
	wx375 = 3.75 / x;
	w = 0.003923767;
	for(i = 7; i >= 0; i--)	w = w * wx375 + b[i];
	return w / sqrt(x) * exp(x);
}

double besi1(double x)
{
	int i;
	double w, wx375;
	static double a[6] = {	 0.5,        0.87890594,  0.51498869,
							 0.15084934, 0.02658733,  0.00301532};
	static double b[8] = {	 0.39894228,  -0.039880242, -0.003620183,
							 0.001638014, -0.01031555,   0.022829673,
							-0.028953121,  0.017876535};

	if(x < 0.)
	{
		fprintf(stderr, "Error : x < 0  in besi1()\n");
		return 0.;
	}
	if(x <= 3.75)
	{
		wx375 = x * x / 14.0625;
		w = 0.00032411;
		for(i = 5; i >= 0; i--)	w = w * wx375 + a[i];
		return w * x;
	}
	wx375 = 3.75 / x;
	w = -0.004200587;
	for(i = 7; i >= 0; i--)	w = w * wx375 + b[i];
	return w / sqrt(x) * exp(x);
}

double besk0(double x)
{
	int i;
	double bi0, w, wx2, wx375;
	static double a[6] = {	-0.57721566,  0.42278442,  0.23069756,
							 0.0348859,   0.00262698,  0.0001075};
	static double b[6] = {	 1.0,        3.5156229,  3.0899424,
							 1.2067492,  0.2659732,  0.0360768};
	static double c[6] = {	 1.25331414, -0.07832358,  0.02189568,
							-0.01062446,  0.00587872, -0.0025154};

	if(x <= 0.)
	{
		fprintf(stderr, "Error : x <= 0  in besk0()\n");
		return 0.;
	}
	if(x <= 2.)
	{
		wx2 = x * x / 4.;
		wx375 = x * x / 14.0625;
		w = 0.0000074;
		bi0 = 0.0045813;
		for(i = 5; i >= 0; i--)
		{
			w = w * wx2 + a[i];
			bi0 = bi0 * wx375 + b[i];
		}
		return w - log(x / 2.) * bi0;
	}
	wx2 = 2. / x;
	w = 0.00053208;
	for(i = 5; i >= 0; i--)	w = w * wx2 + c[i];
	return w / sqrt(x) / exp(x);
}

double besk1(double x)
{
	int i;
	double bi1, w, wx2, wx375;
	static double a[6] = {	 1.0,         0.15443144, -0.67278579,
							-0.18156897, -0.01919402, -0.00110404};
	static double b[6] = {	 0.5,         0.87890594,  0.51498869,
							 0.15084934,  0.02658733,  0.003011532};
	static double c[6] = {	 1.25331414,  0.23498619, -0.0365562,
							 0.01504268, -0.00780353,  0.00325614};

	if(x <= 0.)
	{
		fprintf(stderr, "Error : x <= 0  in besk1()\n");
		return 0.;
	}
	if(x <= 2.)
	{
		wx2 = x * x / 4.;
		wx375 = x * x / 14.0625;
		w = -0.00004686;
		bi1 = 0.00032411;
		for(i = 5; i >= 0; i--)
		{
			w = w * wx2 + a[i];
			bi1 = bi1 * wx375 + b[i];
		}
		return w / x + log(x / 2.) * bi1 * x;
	}
	wx2 = 2. / x;
	w = -0.00068245;
	for(i = 5; i >= 0; i--)	w = w * wx2 + c[i];
	return w / sqrt(x) / exp(x);
}

double erfnc(double x)
{
	int i;
	double w;
	static double a[6] = {	 1.0,           0.0705230784,  0.0422820123,
							 0.0092705272,	0.0001520143,  0.0002765672};

	if(x < 0.)
	{
		fprintf(stderr, "Error : x < 0  in erfnc()\n");
		return 0.;
	}
	if(x == 0.)	return 0.;
	if(x >= 10.)	return 1.;
	w = 0.0000430638;
	for(i = 5; i >= 0; i--)	w = w * x + a[i];
	for(i = 0; i < 4; i++)	w *= w;
	return 1. - 1. / w;
}

double gammaf(double x)
{
	int i;
	double w, w1;
	static double a[8] = {	 1.0,         -0.577191652,  0.988205891,
							-0.897056937,  0.918206857, -0.756704078,
							 0.482199394, -0.193527818};

	if(x <= 0.)
	{
		fprintf(stderr, "Error : x <= 0  in gammaf()\n");
		return 0.;
	}
	if(x >= 34.)
	{
		fprintf(stderr,	"Error : x >= 34  (Overflow) in gammaf()\n");
		return 0.;
	}
	w1 = 1.;
	if(x > 1.)
	{
		do
		{
			x -= 1.;
			w1 *= x;
		} while(x > 1.);
	}
	w = 0.035868343;
	for(i = 7; i >= 0; i--)	w = w * x + a[i];
	return w1 * w / x;
}

double legend(double x, int n)
{
	double p0, p1, p2, w, wn;

	if(n < 0)
	{
		fprintf(stderr, "Error : n < 0  in legend()\n");
		return 0.;
	}
	if(n == 0)	return 1.;
	if(n == 1)	return x;
	wn = (double)n;
	p0 = 1.;
	p1 = x;
	do
	{
		wn -= 1.;
		w = wn / (wn + 1.);
		p2 = x * p1 + w * (x * p1 - p0);
		p0 = p1;
		p1 = p2;
	} while(wn >= 2.);
	return p2;
}

double celi1(double k, double eps)
{
	double a0, a1, b0, b1;

	if(eps <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter  in celi1()\n");
		return 0.;
	}
	if(fabs(k) > 1.)
	{
		fprintf(stderr, "Error : |k| > 1  in celi1()\n");
		return 0.;
	}
	a0 = 1.;
	b0 = sqrt(1. - k * k);
	for(;;)
	{
		a1 = (a0 + b0) / 2.;
		b1 = sqrt(a0 * b0);
		if(fabs(a1 - b1) <= eps)	return M_PI_2 / a1;
		a0 = a1;
		b0 = b1;
	}
}

double celi2(double k, double eps)
{
	int i, j;
	double a0, a1, b0, b1, c0, c1, cn, w;

	if(eps <= 0.)
	{
		fprintf(stderr, "Error : Illegal parameter  in celi2()\n");
		return 0.;
	}
	if(fabs(k) > 1.)
	{
		fprintf(stderr, "Error : |k| > 1  in celi2()\n");
		return 0.;
	}
	a0 = 1.;
	c0 = k * k;
	b0 = sqrt(1. - c0);
	cn = 0.;
	i = 0;
	for(;;)
	{
		a1 = (a0 + b0) / 2.;
		b1 = sqrt(a0 * b0);
		c1 = (a0 - b0) / 2.;
		w = 1.;
		for(j = 0; j < i; j++)	w *= 2.;
		cn += w * c1 * c1;
		i++;
		if(fabs(a1 - b1) < eps)	return M_PI_2 / a1 * (1. - (c0 / 2. + cn));
		a0 = a1;
		b0 = b1;
	}
}
