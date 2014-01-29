/*		axis.c			*/
#include "sslib.h"

/* Degree to Radian Transformation */

int dtor(double deg, double *rad)
{
	static double dtrsr = 1.7463292519943296e-2;

	*rad = deg * dtrsr;
	return 0;
}

/* Radian to Degree Transformation */

int rtod(double rad, double *deg)
{
	static double rtrsd = 57.29577951308232;

	*deg = rad * rtrsd;
	return 0;
}

/* Cartesian to Polar Transformation */

int dtop2(double xo, double yo, double *r, double *t)
{
	static double rtod = 57.29577951308232;

	*r = sqrt(xo * xo + yo * yo);
	if(xo != 0.0)			*t = atan2(yo, xo) * rtod;
	else
	{
		if(yo < 0.0)		*t = - 90.0;
		else if(yo > 0.0)	*t = 90.0;
		else				*t = 0.0;
	}
	return 0;
}

/* Polar to Cartesian Transformation */

int ptod2(double r, double t, double *x, double *y)
{
	double tr;
	static double dtor = 1.7453292519943296e-2;

	tr = t * dtor;
	*x = r * cos(tr);
	*y = r * sin(tr);
	return 0;
}

/* Parallel move (Cartesian coordinate) */

int dmov2(double xo, double yo, double xm, double ym, double *x, double *y)
{
	*x = xo + xm;
	*y = yo + ym;
	return 0;
}

/* Parallel move (Polar coordinate) */

int pmov2(double ro, double to, double xm, double ym, double *r, double *t)
{
	double dx, dy, ddx, ddy;

	ptod2(ro, to, &dx, &dy);
	dmov2(dx, dy, xm, ym, &ddx, &ddy);
	dtop2(ddx, ddy, r, t);
	return 0;
}

/* Rotation (Cartesian coordinate) */

int drot2(double xo, double yo, double t, double *x, double *y)
{
	double r, thet, to;
	static double dtor = 1.7453292519943296e-2;

	dtop2(xo, yo, &r, &thet);
	to = (t + thet) * dtor;
	*x = r * cos(to);
	*y = r * sin(to);
	return 0;
}

/* Rotation (Polar coordinate) */

int prot2(double ro, double to, double trot, double *r, double *t)
{
	*t = to + trot;
	*r = ro;
	return 0;
}

/* Cartesian to Cylindorical ( 3-dimansion ) */

int dtoc3(double dx, double dy, double dz, double *cr, double *ct, double *cz)
{
	dtop2(dx, dy, cr, ct);
	*cz = dz;
	return 0;
}

/* Cylindorical to Cartesian ( 3-dimansion ) */

int ctod3(double cr, double ct, double cz, double *dx, double *dy, double *dz)
{
	ptod2(cr, ct, dx, dy);
	*dz = cz;
	return 0;
}

/* Cartesian to Polar ( 3-dimension ) */

int dtop3(double dx, double dy, double dz, double *pr, double *pe, double *pp)
{
	double w;
	static double rtod = 57.295779513082321;

	*pr = sqrt(dx * dx + dy * dy + dz * dz);
	if(dx != 0.0)			*pe = atan2(dy, dx) * rtod;
	else
	{
		if(dy < 0.0)		*pe = -90.0;
		else if(dy > 0.0)	*pe = 90.0;
		else				*pe = 0.0;
	}
	w = dx * dx + dy * dy;
	if(dz != 0.0)			*pp = atan2(sqrt(w), dz) * rtod;
	else
	{
		if(w < 0.0)			*pp = -90.0;
		else if(w > 0.0)	*pp = 90.0;
		else				*pp = 0.0;
	}
	return 0;
}

/* Polar to Cartesian ( 3-dimension ) */

int ptod3(double pr, double pe, double pp, double *dx, double *dy, double *dz)
{
	double rpe, rpp;
	static double dtor = 1.74532925199432958e-2;

	rpe = dtor * pe;
	rpp = dtor * pp;
	*dx = pr * sin(rpp) * cos(rpe);
	*dy = pr * sin(rpp) * sin(rpe);
	*dz = pr * cos(rpp);
	return 0;
}

/* Exchange Coordinate ( 3-dimension ) */

int ptoc3(double pr, double pe, double pp, double *cr, double *ct, double *cz)
{
	double rcp;
	static double dtor = 1.74532925199432958e-2;

	rcp = dtor * pp;
	*cr = pr * sin(rcp);
	*ct = pe;
	*cz = pr * cos(rcp);
	return 0;
}

/* Exchange Coordinate ( 3-dimension ) */

int ctop3(double cr, double ct, double cz, double *pr, double *pe, double *pp)
{
	static double rtod = 57.29577951308232088;

	*pr = sqrt(cr * cr + cz * cz);
	*pe = ct;
	if(cz != 0)				*pp = atan2(cr, cz) * rtod;
	else
	{
		if(cr < 0.0)		*pp = -90.0;
		else if(cr > 0.0)	*pp = 90.0;
		else				*pp = 0.0;
	}
	return 0;
}

/* Cartesian coordinate */

int dpmov3(double xo, double yo, double zo, double dx, double dy, double dz,
	double *x, double *y, double *z)
{
	*x = xo + dx;
	*y = yo + dy;
	*z = zo + dz;
	return 0;
}

/* Cylindorical coordinate */

int cpmov3(double cro, double cto, double czo, double dx, double dy, double dz,
	double *cr, double *ct, double *cz)
{
	double dxo, dyo, dzo, x, y, z;

	ctod3(cro, cto, czo, &dxo, &dyo, &dzo);
	dpmov3(dxo, dyo, dzo, dx, dy, dz, &x, &y, &z);
	dtoc3(x, y, z, cr, ct, cz);
	return 0;
}

/*	Polar coordinate */

int ppmov3(double pro, double peo, double ppo, double dx, double dy, double dz,
	double *pr, double *pe, double *pp)
{
	double dxo, dyo, dzo, x, y, z;

	ptod3(pro, peo, ppo, &dxo, &dyo, &dzo);
	dpmov3(dxo, dyo, dzo, dx, dy, dz, &x, &y, &z);
	dtop3(x, y, z, pr, pe, pp);
	return 0;
}

/* Rotate Move in 3-Dimensional Coordinate ( Cartesian Coordinate ) */

int drot3(double dxo, double dyo, double dzo, double phai, double thet,
	double psai, double *dx, double *dy, double *dz)
{
	double rph, rth, rps, cph, sph, cth, sth, cps, sps;
	static double dtor = 1.74532925199432958e-2;

	rph = dtor * phai;
	rth = dtor * thet;
	rps = dtor * psai;
	cph = cos(rph);
	sph = sin(rph);
	cth = cos(rth);
	sth = sin(rth);
	cps = cos(rps);
	sps = sin(rps);
	*dx = (cph * cps - sph * cth * sps) * dxo
		+ (sph * cps + cph * cth * sps) * dyo
		+  sth * sps * dzo;
	*dy = (-cph * sps - sph * cth * cps) * dxo
		+ (-sph * sps + cph * cth * cps) * dyo
		+  sth * cps * dzo;
	*dz = sth * sph * dxo - sth * cph * dyo + cth * dzo;
	return 0;
}

/* Rotate Move in 3-Dimensional Coordinate ( Cylindorical Coordinate ) */

int crot3(double cro, double cto, double czo, double phai, double thet,
	double psai, double *cr, double *ct, double *cz)
{
	double dx, dy, dz, rdx, rdy, rdz;

	ctod3(cro, cto, czo, &dx, &dy, &dz);
	drot3(dx, dy, dz, phai, thet, psai, &rdx, &rdy, &rdz);
	dtoc3(rdx, rdy, rdz, cr, ct, cz);
	return 0;
}

/* Rotate Move in 3-Dimensional Coordinate ( Polar Coordinate ) */

int prot3(double pro, double peo, double ppo, double phai, double thet,
	double psai, double *pr, double *pe, double *pp)
{
	double dx, dy, dz, rdx, rdy, rdz;

	ptod3(pro, peo, ppo, &dx, &dy, &dz);
	drot3(dx, dy, dz, phai, thet, psai, &rdx, &rdy, &rdz);
	dtop3(rdx, rdy, rdz, pr, pe, pp);
	return 0;
}
