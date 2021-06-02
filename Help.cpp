#include "Help.h"

void polar_perenos(const double& x1, const double& y1, const double& x2, const double& y2, double& u, double& v)
{
	double phi1 = polar_angle(x1, y1);
	double phi2 = polar_angle(x2, y2);
	double fr = u * cos(phi1) + v * sin(phi1);
	double ff = -u * sin(phi1) + v * cos(phi1);
	u = fr * cos(phi2) - ff * sin(phi2);
	v = fr * sin(phi2) + ff * cos(phi2);
	return;
}


double polar_angle(const double& x, const double& y)
{
	if (x < 0)
	{
		return atan(y / x) + 1.0 * pi_;
	}
	else if (x > 0 && y >= 0)
	{
		return atan(y / x);
	}
	else if (x > 0 && y < 0)
	{
		return atan(y / x) + 2.0 * pi_;
	}
	else if (y > 0 && x >= 0 && x <= 0)
	{
		return pi_ / 2.0;
	}
	else if (y < 0 && x >= 0 && x <= 0)
	{
		return  3.0 * pi_ / 2.0;
	}
	return 0.0;
}

double minmod(const double& x, const double& y)
{
	if (sign(x) + sign(y) == 0)
	{
		return 0.0;
	}
	else
	{
		return   ((sign(x) + sign(y)) / 2.0) * min(fabs(x), fabs(y));  ///minmod
		//return (2*x*y)/(x + y);   /// vanleer
	}
}

int sign(const double& x)
{
	if (x > 0)
	{
		return 1;
	}
	else if (x < 0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& x3, const double& t3, const double& y)
// Главное значение с параметрами 2
{
	double d = minmod((t1 - t2) / (x1 - x2), (t2 - t3) / (x2 - x3));
	return  (d * (y - x2) + t2);
}

double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& y)
// Главное значение с параметрами 2
{
	double d = (t1 - t2) / (x1 - x2);
	return  (d * (y - x2) + t2);
}