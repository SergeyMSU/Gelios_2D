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

void polar_provorot(const double& phi, double& u, double& v)
{
	double uu = u, vv = v;
	u = uu * cos(phi) - vv * sin(phi);
	v = uu * sin(phi) + vv * cos(phi);
	return;
}

//double max(const double& x, const double& y)
//{
//	if (x >= y)
//	{
//		return x;
//	}
//	else
//	{
//		return y;
//	}
//}


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

double sign(const double& x)
{
	if (x > 0)
	{
		return 1.0;
	}
	else if (x < 0)
	{
		return -1.0;
	}
	else
	{
		return 0.0;
	}
}

double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& x3, const double& t3, const double& y)
// Главное значение с параметрами 2
// Строим линии между 1 и 2,  2 и 3, потом находим минмодом значение в y
{
	if (true)
	{
		double d = minmod((t1 - t2) / (x1 - x2), (t2 - t3) / (x2 - x3));
		return  (d * (y - x2) + t2);
	}
	else
	{
		// Новая процедура с сжатием
		double dUl = (t2 - t1) / (x2 - x1);
		double dUr = (t3 - t2) / (x3 - x2);
		return t2 + 0.5 * ((1.0 - eta_) * minmod(dUl, betta_ * dUr) + (1.0 + eta_) * minmod(betta_ * dUl, dUr)) * (y - x2);
	}

}

double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& y)
// Главное значение с параметрами 2
{
	double d = (t1 - t2) / (x1 - x2);
	return  (d * (y - x2) + t2);
}

void dekard_skorost(double x, double y, double z, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz)
{
	double r_2 = sqrt(x * x + y * y + z * z);
	double the_2 = acos(z / r_2);
	double phi_2 = polar_angle(x, y);

	Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
	Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
	Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
}

void dekard_skorost2(double r2, double the_2, double phi_2, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz)
{
	Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
	Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
	Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
}

void spherical_skorost(const double& x, const double& y, const double& z, const double& Vx,//
	const double& Vy, const double& Vz, double& Vr, double& Vphi, double& Vtheta)
{
	double r_1 = sqrt(x * x + y * y + z * z);
	double the_1 = acos(z / r_1);
	double phi_1 = polar_angle(x, y);

	Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
	Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
	Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
}


double Godunov_squere_rad(const double& x1, const double& r1, const double& x2, const double& r2,//
	const double& x3, const double& r3, const double& x4, const double& r4)
{
	return (pi_ / 3.0) * ((r4 - r2) * (x3 * (r2 + r4 + r3) - x1 * (r2 + r4 + r1)) - //
		(r3 - r1) * (x4 * (r1 + r3 + r4) - x2 * (r1 + r3 + r2)));
}