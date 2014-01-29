/*		difint.c		*/
#include <stdio.h>
#include <stdlib.h>
#include "sslib.h"

double cheb3(double a, double b)
{
	static double g[2] = { 0.0, M_SQRT_2};
	double w;
	int mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		w = a;
		a = b;
		b = w;
	}
	w =  _f(_Normal(a, b,  g[0]));
	w += _f(_Normal(a, b,  g[1]));
	w += _f(_Normal(a, b, -g[1]));
	w *= ((b - a) / 3.);
	if(mflag)	return -w;
	return w;
}

double cheb4(double a, double b)
{
	static double g[2] = { 0.7946544722917661, 0.1875924740850799};
	double w;
	int mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		w = a;
		a = b;
		b = w;
	}
	w =  _f(_Normal(a, b,  g[0]));
	w += _f(_Normal(a, b, -g[0]));
	w += _f(_Normal(a, b,  g[1]));
	w += _f(_Normal(a, b, -g[1]));
	w *= ((b - a) / 4.);
	if(mflag)	return -w;
	return w;
}

double cheb6(double a, double b)
{
	static double g[3] = 
		{ 0.2666354015167047, 0.4225186537611115, 0.8662468181078206};
	double w;
	int mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		w = a;
		a = b;
		b = w;
	}
	w =  _f(_Normal(a, b,  g[0]));
	w += _f(_Normal(a, b, -g[0]));
	w += _f(_Normal(a, b,  g[1]));
	w += _f(_Normal(a, b, -g[1]));
	w += _f(_Normal(a, b,  g[2]));
	w += _f(_Normal(a, b, -g[2]));
	w *= ((b - a) / 6.);
	if(mflag)	return -w;
	return w;
}

double dgl3(double a, double b)
{
	static double g = 0.774596669241483;
	static double w1 = 0.555555555555556;
	static double w0 = 0.88888888888889;
	double w;
	int mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		w = a;
		a = b;
		b = w;
	}
	w =  _f(_Normal(a, b,  0.)) * w0;
	w += (_f(_Normal(a, b, g)) + _f(_Normal(a, b, -g))) * w1;
	w *= ((b - a) / 2.);
	if(mflag)	return -w;
	return w;
}

double dgl10(double a, double b)
{
	static double g[5] = {	0.148874338981631, 0.433395394129247,
							0.679409568299024, 0.865063366688985,
							0.973906528517172};
	static double w[5] = {	0.295524224714753, 0.269266719309996,
							0.219086362515982, 0.149451349150581,
							0.066671344308688};
	double s;
	int i, mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		s = a;
		a = b;
		b = s;
	}
	for(i = 0, s = 0.; i < 5; i++)	s += w[i] *(_f(_Normal(a, b, g[i])) + _f(_Normal(a, b, -g[i])));
	s *= ((b - a) / 2.);
	if(mflag)	return -s;
	return s;
}

double dgl20(double a, double b)
{
	static double g[10] = {	0.076526521133497333, 0.227785851141645078,
							0.373706088715419561, 0.510867001950827098,
							0.636053680726515025, 0.746331906460150793,
							0.839116971822218823, 0.912234428251325906,
							0.963971927277913791, 0.993128599185094925};
	static double w[10] = {	0.152753387130725851, 0.149172986472603747,
							0.142096109318382051, 0.131688638449176627,
							0.118194531961518417, 0.101930119817240435,
							0.083276741576704749, 0.062672048334109064,
							0.040601429800386941, 0.017614007139152118};
	double s;
	int i, mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		s = a;
		a = b;
		b = s;
	}
	for(i = 0, s = 0.; i < 10; i++)	s += w[i] *(_f(_Normal(a, b, g[i])) + _f(_Normal(a, b, -g[i])));
	s *= ((b - a) / 2.);
	if(mflag)	return -s;
	return s;
}

double dgl32(double a, double b)
{
	static double g[16] = {	0.048307665687738316, 0.144471961582796494,
							0.239287362252137075, 0.331868602282127650,
							0.421351276130635345, 0.506899908932229390,
							0.587715757240762329, 0.663044266930215201,
							0.732182118740289680, 0.794483795967942407,
							0.849367613732569970, 0.896321155766052124,
							0.934906075937739690, 0.964762255587506431,
							0.985611511545268335, 0.997263861849481564};
	static double w[16] = {	0.096540088514727801, 0.095638720079274859,
							0.093844399080804566, 0.091173878695763885,
							0.087652093004403811, 0.083311924226946755,
							0.078193895787070306, 0.072345794108848506,
							0.065822222776361847, 0.058684093478535547,
							0.050998059262376176, 0.042835898022226681,
							0.034273862913021433, 0.025392065309262059,
							0.016274394730905671, 0.007018610009470097};
	double s;
	int i, mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		s = a;
		a = b;
		b = s;
	}
	for(i = 0, s = 0.; i < 16; i++)	s += w[i] *(_f(_Normal(a, b, g[i])) + _f(_Normal(a, b, -g[i])));
	s *= ((b - a) / 2.);
	if(mflag)	return -s;
	return s;
}

double dgl48(double a, double b)
{
	static double g[24] = {	0.032380170962869362, 0.097004699209462699,
							0.161222356068891718, 0.224763790394689061,
							0.287362487355455577, 0.348755886292160738,
							0.408686481990716730, 0.466902904750958405,
							0.523160974722233034, 0.577224726083972704,
							0.628867396776513624, 0.677872379632663905,
							0.724034130923814655, 0.767159032515740339,
							0.807066204029442627, 0.843588261624393531,
							0.876572020274247886, 0.905879136715569673,
							0.931386690706554333, 0.952987703160430861,
							0.970591592546247250, 0.984124583722826858,
							0.993530172266350758, 0.998771007252426119};
	static double w[24] = {	0.064737696812683923, 0.064466164435950082,
							0.063924238584648187, 0.063114192286254026,
							0.062039423159892664, 0.060704439165893880,
							0.059114839698395636, 0.057277292100403216,
							0.055199503699984163, 0.052890189485193667,
							0.050359035553854475, 0.047616658492490475,
							0.044674560856694280, 0.041545082943464749,
							0.038241351065830706, 0.034777222564770439,
							0.031167227832798089, 0.027426509708356948,
							0.023570760839324379, 0.019616160457355528,
							0.015579315722943849, 0.011477234579234539,
							0.007327553901276262, 0.003153346052305839};
	double s;
	int i, mflag;

	if(a == b)	return 0.;
	mflag = 0;
	if(a > b)
	{
		mflag = 1;
		s = a;
		a = b;
		b = s;
	}
	for(i = 0, s = 0.; i < 24; i++)	s += w[i] *(_f(_Normal(a, b, g[i])) + _f(_Normal(a, b, -g[i])));
	s *= ((b - a) / 2.);
	if(mflag)	return -s;
	return s;
}

double dglg3(void)
{
	int i;
	static double g[3] =
		{ 0.41577455678348, 2.294280360279, 6.2899450829375};
	static double w[3] =
		{ 0.711093009929, 0.27851773356924, 0.010389256501586};
	double s;

	for(i = 0, s = 0.; i < 3; i++)	s += _f(g[i]) * w[i];
	return s;
}

double dglg5(void)
{
	int i;
	static double g[5] =
		{ 12.64080084427578,    7.085810005858838, 3.596425771040722,
		   1.41340305910651679, 0.2635603197181409};
	static double w[5] =
		{ 0.0000233699723857762, 0.00361175867992205, 0.0759424496817060,
		  0.398666811083176,     0.521755610582809};
	double s;

	for(i = 0, s = 0.; i < 5; i++)	s += _f(g[i]) * w[i];
	return s;
}

double dglg10(void)
{
	int i;
	static double g[10] =
		{ 29.92069701227389, 21.99658581198076,  16.27925783137810,
		  11.84378583790007,  8.330152746764497,  5.552496140063804,
		   3.401433697854900, 1.808342901740316,  0.7294545495031705,
		   0.1377934705404924};
	static double w[10] =
		{ 9.9118272196090e-13,   1.839564823979631e-9,  4.249313984962686e-7,
		  2.8259233495995656e-5, 7.530083885875388e-4,  9.5015169751811006e-3,
		  6.208745609867775e-2,  2.1806828761180942e-1, 4.011199291552736e-1,
		  3.0844111576502014e-1};
	double s;

	for(i = 0, s = 0.; i < 10; i++)	s += _f(g[i]) * w[i];
	return s;
}

double dgh10(void)
{
	int i;
	static double g[5] =
		{ 0.3429013272237046, 1.036610829789514, 1.756683649299882,
		  2.532731674232790,  3.436159118837738};
	static double w[5] =
		{ 0.6108626337353258, 0.2401386110823147, 0.03387439445548106,
		  0.001343645746781233, 7.640432855232621e-6};
	double s;

	for(i = 0, s = 0.; i < 5; i++)	s += (_f(g[i]) + _f(-g[i])) * w[i];
	return s;
}

double dgh15(void)
{
	int i;
	static double g[8] =
		{ 0.0,               0.56506958325557575, 1.136115585210921,
		  1.719992575186489, 2.325732486173858,   2.967166927905603,
		  3.669950373404453, 4.499990707309392};
	static double w[8] =
		{ 0.5641003087264175,   0.4120286874988986,   0.1584889157959357,
		  0.03078003387254608,  0.002778068842912276, 0.0001000044412324999,
		  1.059115547711067e-6, 1.522475804253517e-9};
	double s;

	s = _f(g[0]) * w[0];
	for(i = 1; i < 8; i++)	s += (_f(g[i]) + _f(-g[i])) * w[i];
	return s;
}

double hardy(double xmin, double xmax, int n)
{
	double h, hh, h3, h5, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in hardy()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - 28. * _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 6.;
	h3 = h / 2.;
	h5 = h / 1.2;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (56. * _f(x) + 162. * (_f(x - h5) + _f(x - hh))
						  + 220. *  _f(x - h3));
	}
	return s * h / 600.;
}

double lomberg(double xmin, double xmax, double eps)
{
	int idx, ierr, j, k;
	double t[434];
	double fxmax, fxmin, fw, hi, hk, xw, w, w1;

	if(eps <= 0.)
	{
		fprintf(stderr, "Error : eps <= 0  in lomberg()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	w = t[0] = ((fxmin = _f(xmin)) + (fxmax = _f(xmax)))
			 * (xw = xmax - xmin) / 2.;
	fw = (fxmin - fxmax) / 2.;
	idx = 0;
	for(k = 2, hk = 2.; k < 30; k++, hk *= 2.)
	{
		w = xw / hk;
		for(hi = 1., w1 = fw; hi <= hk; hi += 1.)	w1 += _f(xmin + w * hi);
		w = (t[++idx] = w1 * w);
		for(j = 1, w1 = 4.; j < k; j++, w1 *= 4.)
		{
			idx++;
			w = (t[idx] = (w1 * w - t[idx - k]) / (w1 - 1.));
		}
		if(fabs(w - t[idx - k]) < fabs(eps * t[idx - k])) return w;
	}
	fprintf(stderr, "Error : No convergence in lomberg()\n");
	return w;
}

double nc1(double xmin, double xmax, int n)
{
	double h, hi, s;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc1()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = (_f(xmin) + _f(xmax)) / 2.;
	h = (xmax - xmin) / (double)n;
	for(i = 1, hi = 1.; i < n; i++, hi += 1.)	s += _f(xmin + h * hi);
	return s * h;
}

double nc2(double xmin, double xmax, int n)
{
	double h, hh, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc2()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 2.;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (4. * _f(x - hh) + 2. * _f(x));
	}
	return s * h / 6.;
}

double nc3(double xmin, double xmax, int n)
{
	double h, hh, h2, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc3()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 3.;
	h2 = h / 1.5;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (3. * (_f(x - h2) + _f(x - hh)) + 2. * _f(x));
	}
	return s * h / 8.;
}

double nc4(double xmin, double xmax, int n)
{
	double h, hh, h2, h3, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc4()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - 7. * _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 4.;
	h2 = h / 2.;
	h3 = h * .75;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (14. * _f(x) + 32. * (_f(x - h3) + _f(x - hh)) + 12. * _f(x - h2));
	}
	return s * h / 90.;
}

double nc5(double xmin, double xmax, int n)
{
	double h, hh, h2, h3, h4, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc5()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - 19. * _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 5.;
	h2 = h * .4;
	h3 = h * .6;
	h4 = h * .8;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (38. * _f(x) + 75. * (_f(x - h4) + _f(x - hh))
						  + 50. * (_f(x - h3) + _f(x - h2)));
	}
	return s * h / 288.;
}

double nc6(double xmin, double xmax, int n)
{
	double h, hh, h2, h3, h4, h5, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc6()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - 41. * _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 6.;
	h2 = h / 3.;
	h3 = h / 2.;
	h4 = h / 1.5;
	h5 = h / 1.2;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (82. * _f(x) + 216. * (_f(x - h5) + _f(x - hh))
						  +  27. * (_f(x - h4) + _f(x - h2))
						  + 272. *  _f(x - h3));
	}
	return s * h / 840.;
}

double nc7(double xmin, double xmax, int n)
{
	double h, hh, h2, h3, h4, h5, h6, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc7()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - 751. * _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 7.;
	h2 = h / 3.5;
	h3 = hh * 3.;
	h4 = h / 1.75;
	h5 = h / 1.4;
	h6 = hh * 6.;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (1502. * _f(x) + 3577. * (_f(x - h6) + _f(x - hh))
							+ 1323. * (_f(x - h5) + _f(x - h2))
							+ 2989. * (_f(x - h4) + _f(x - h3)));
	}
	return s * h / 17280.;
}

double nc8(double xmin, double xmax, int n)
{
	double h, hh, h2, h3, h4, h5, h6, h7, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in nc8()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - 989. * _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 8.;
	h2 = h / 4.;
	h3 = h * .375;
	h4 = h / 2.;
	h5 = h / 1.6;
	h6 = h * .75;
	h7 = h * .875;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (1978. * _f(x) +  5888. * (_f(x - h7) + _f(x - hh))
							-   928. * (_f(x - h6) + _f(x - h2))
							+ 10496. * (_f(x - h5) + _f(x - h3))
							-  4540. *  _f(x - h4));
	}
	return s * h / 28350.;
}

double weddle(double xmin, double xmax, int n)
{
	double h, hh, h2, h3, h4, h5, hi, s, x;
	int i;

	if(n < 2)
	{
		fprintf(stderr, "Error : n < 2  in weddle()\n");
		return 0.;
	}
	if(xmin == xmax)	return 0.;
	s = _f(xmin) - _f(xmax);
	h = (xmax - xmin) / (double)n;
	hh = h / 6.;
	h2 = h / 3.;
	h3 = h / 2.;
	h4 = h / 1.5;
	h5 = h / 1.2;
	for(i = 0, hi = 1.; i < n; i++, hi += 1.)
	{
		x = xmin + h * hi;
		s += (2. * _f(x) + 5. * (_f(x - h5) + _f(x - hh))
						 +       _f(x - h4) + _f(x - h2)
						 + 6. *  _f(x - h3));
	}
	return s * h / 20.;
}

double _Normal(double a, double b, double x)
{
	return ((b - a) * x + a + b) / 2.;
}

double difm1(double x, double h)
{
	double f1;

	f1 = _f(x + h) - _f(x - h);
	return 0.5 * f1 / h;
}

double difm2(double x, double h)
{
	double f1, f2;

	f1 = _f(x + h) - _f(x - h);
	f2 = _f(x + 2. * h) - _f(x - 2. * h);
	return (2. * f1 / 3. - f2 / 12.) / h;
}

double difm3(double x, double h)
{
	double f1, f2, f3;

	f1 = _f(x + h) - _f(x - h);
	f2 = _f(x + 2. * h) - _f(x - 2. * h);
	f3 = _f(x + 3. * h) - _f(x - 3. * h);
	return (0.75 * f1 - 0.15 * f2 + f3 / 60.) / h;
}

double difm4(double x, double h)
{
	double f1, f2, f3, f4;

	f1 = _f(x + h) - _f(x - h);
	f2 = _f(x + 2. * h) - _f(x - 2. * h);
	f3 = _f(x + 3. * h) - _f(x - 3. * h);
	f4 = _f(x + 4. * h) - _f(x - 4. * h);
	return (0.8 * f1 - 0.2 * f2 + 4. * f3 / 105. - f4 / 280.) / h;
}

double difm5(double x, double h)
{
	double f1, f2, f3, f4, f5;

	f1 = _f(x + h) - _f(x - h);
	f2 = _f(x + 2. * h) - _f(x - 2. * h);
	f3 = _f(x + 3. * h) - _f(x - 3. * h);
	f4 = _f(x + 4. * h) - _f(x - 4. * h);
	f5 = _f(x + 5. * h) - _f(x - 5. * h);
	return (2100. * f1 - 600. * f2 + 150. * f3 - 25. * f4 + 2. * f5) / 2520./ h;
}

double difm6(double x, double h)
{
	double f1, f2, f3, f4, f5, f6;

	f1 = _f(x + h) - _f(x - h);
	f2 = _f(x + 2. * h) - _f(x - 2. * h);
	f3 = _f(x + 3. * h) - _f(x - 3. * h);
	f4 = _f(x + 4. * h) - _f(x - 4. * h);
	f5 = _f(x + 5. * h) - _f(x - 5. * h);
	f6 = _f(x + 6. * h) - _f(x - 6. * h);
	return (23760. * f1 - 7425. * f2 + 2200. * f3 - 495. * f4 + 72. * f5 - 5. * f6) / 27720. / h;
}

double difm7(double x, double h)
{
	double f1, f2, f3, f4, f5, f6, f7;

	f1 = _f(x + h) - _f(x - h);
	f2 = _f(x + 2. * h) - _f(x - 2. * h);
	f3 = _f(x + 3. * h) - _f(x - 3. * h);
	f4 = _f(x + 4. * h) - _f(x - 4. * h);
	f5 = _f(x + 5. * h) - _f(x - 5. * h);
	f6 = _f(x + 6. * h) - _f(x - 6. * h);
	f7 = _f(x + 7. * h) - _f(x - 7. * h);
	return (315315. * f1 - 105105. * f2 + 35035. * f3 - 9555. * f4 + 1911. * f5 - 245. * f6 + 15. * f7) / 360360. / h;
}

double difm8(double x, double h)
{
	double f1, f2, f3, f4, f5, f6, f7, f8;

	f1 = _f(x + h) - _f(x - h);
	f2 = _f(x + 2. * h) - _f(x - 2. * h);
	f3 = _f(x + 3. * h) - _f(x - 3. * h);
	f4 = _f(x + 4. * h) - _f(x - 4. * h);
	f5 = _f(x + 5. * h) - _f(x - 5. * h);
	f6 = _f(x + 6. * h) - _f(x - 6. * h);
	f7 = _f(x + 7. * h) - _f(x - 7. * h);
	f8 = _f(x + 8. * h) - _f(x - 8. * h);
	return (1921920. * f1 - 672672. * f2 + 244608. * f3 - 76440. * f4 + 18816. * f5 - 3360. * f6 + 384. * f7 - 21. * f8) / 2162160. / h;
}
