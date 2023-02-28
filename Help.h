#pragma once


// ����� �������:
//  1) ���������� ����������� ���� ��� �����
//  2) ���������� ���� ����������
//  3) ���������� ���������� �������� Kn
//  4) ������� �������� ������� (�����-����� \ ������� ��������
//  5) ���������� ����� �����\����������
//
//

#define USEMPI false            // ����� �� ������������ MPI ?

#define Kn_inf true            // ������, ����� ������� �������
 
#include <iostream>
#include <iomanip>

#define pop_atom 6         // ����� ��������� ������

#define low_type double
#define high_type double

// ��� S+  � S-
#define n_S 100
#define max_S 100.0           //  (i + 1) * (max_S/(n_S + 1))   - �� ����� �������� i-�� ������

#define I_ 8 //5 // ���������� ��� �� �������
#define J_ 11 //5 // ���������� ��� �� ����

#define RR_ 390.222 // 389.8557              // ����������� ������ ������

#define ga (5.0/3.0)          // ���������� ��������
#define ggg (5.0/3.0)
#define g1 (ga - 1.0)
#define kv(x) ( (x)*(x) )
#define pow3(x) ( (x)*(x)*(x) )
#define pow4(x) ( (x)*(x)*(x)*(x) )
#define pow5(x) ( (x)*(x)*(x)*(x)*(x) )
#define pow6(x) ( (x)*(x)*(x)*(x)*(x)*(x) )
#define pow7(x) ( (x)*(x)*(x)*(x)*(x)*(x)*(x) )
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))

#define pi_ 3.14159265358979323846
#define sqrtpi_ 1.77245385

#define a_2 0.1307345665  // 0.102578  // 0.10263


#include <initializer_list>
#include <random>
#include <ctime> // ���������� clock
#include <stdint.h>
#include <xmmintrin.h>

#include "Setka.h"
#include "Rail.h"
#include "Point.h"
#include "Cell.h"
#include "Gran.h"
#include "Solvers.h"
#include "sensor.h"
#include "sensor2.h"
#include "MKmethod.h"
#include "Dist_func.h"
#include <mutex>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <algorithm>

#if USEMPI 
#include "mpi.h"
#endif


inline double sign(const double& x);
inline int sgn(const double& val);
inline double minmod(const double& x, const double& y);
inline double polar_angle(const double& x, const double& y);
double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& x3, const double& t3, const double& y);
double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& y);
void polar_perenos(const double& x1, const double& y1, const double& x2, const double& y2, double& u, double& v);
inline void dekard_skorost(const double& x, const double& y, const double& z,
	const double& Vr, const double& Vphi, const double& Vtheta, double& Vx,
	double& Vy, double& Vz);
inline double MyRandom(int& s1, int& s2, int& s3);

void Change(const double& a, const double& b);


void dekard_skorost2(double r2, double the_2, double phi_2, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz);
inline void spherical_skorost(const double& x, const double& y, const double& z, const double& Vx,//
	const double& Vy, const double& Vz, double& Vr, double& Vphi, double& Vtheta);
double Godunov_squere_rad(const double& x1, const double& r1, const double& x2, const double& r2,//
	const double& x3, const double& r3, const double& x4, const double& r4);
void polar_provorot(const double& phi, double& u, double& v);
void Vector_product(const double& a1, const double& a2, const double& a3,//
	const double& b1, const double& b2, const double& b3,//
	double& x, double& y, double& z);
//double max(const double& x, const double& y);

#define Q_gran -90.0              // ���� Q
#define Velosity_inf -2.54338 //-4.54127 //-2.54127
#define M_inf 1.97009
#define chi_ 1.0 //36.1059 // 36.1059
#define chi_real 36.1275
#define kurant  0.1    // 0.1
//#define Kn_  0.622171
//#define Kn_  0.4326569808  // 0.622171 // 0.5  // 242.785									// ����� ��������
#define Kn_  12.0 // 0.4326569808 // 0.4326569808 // 6.0			
//#define Kn_  0.2	                                            // ����� ��������
//#define a_2 0.102578  // 0.10263
#define n_p_LISM_ (3.0) 
#define n_H_LISM_ (1.0)
#define sigma(x) (kv(1.0 - a_2 * log(x)))               // ���������������� ������� �����������
//#define sigma2(x, y) (kv(1.0 - (a_2/(1.0 - a_2 * log(y))) * log(x)))  // ��� ������� ���������������� �������� �� cp
#define sigma2(x, y) (kv(1.0 - a_2 * log((x) * (y))))  // ��� ������� ���������������� �������� �� cp


#define n_inner 20  // ������� ����� �� ���������� ����
#define n_outer_shock 13  // �� ����� ����� �������� ������� ����� (������� �� ���������� �����)
#define n_inner_shock 5  // �� ����� ����� �������� ������� ����� (������ �� 90 ��������)
//#define m_inner 16  // ������� ����� �� ���������� ����  (�� ����� ������ ���������� ���� ���������)
#define zazor 3.0   // ����� ������ ����� ��������� � ����������� ������� 
#define zazor2 7.0   // ����� ������ ����� ������� ������ � ����������� ������� 
#define s_k 0.025   // ������ ������� 4 � ������������
#define alpha_rot 0.1  // ���� �������� (������ ������), �� ���� � ���� ������� �� �������� ���� � � ������
#define R1_ (1.0/RR_)  // 1.0
#define R11_ (38.0/RR_)  // 38.0    // ���������! �� ����� ���������� ������� �� ������� ������
#define R111_ (50.0/RR_)  // 50.0 
#define R2_ (100.0/RR_)  // 100.0
#define R3_ (150.0/RR_)  // 150.0
#define R4_ (600.0/RR_)  // 600.0
#define R5_ (2200.0/RR_)  // 2200.0
//#define Rmax_ (1300.0/RR_)  // 1300.0  // (1300.0)   // ������������ ������ ��� �����-�����
#define Rmax_ (1450.0/RR_)  // 1800.0  // (1300.0)   // ������������ ������ ��� �����-�����
#define Left_ (-1500.0/RR_)  // -1500
#define H_pow 0.5 //0.45  // ���������� �������� ��������� ��� ������ ���������� ��������
#define weight_ 0.25     // ���� ��� ������
#define gam(i) ( 1.0/(kv(R5_ - 2.0/RR_)/kv(Ri[(i)]) - 1.0) )
#define R_m 4.0   // ��� ��������� ���������� �� ��� ���������� ������ � ����� ����� �������
#define Ur_m 0.0   // �������� �����������
#define geo_step 100.0   // ������� ��� ���� �������� � ���� ������
#define geo_accur 100.0    // �������� ��������� ����� ����������
#define Weyght 1000000.0  // ��� ����� ����� �� ���� ������ ��������
#define eta_ 0.5    // �������� ��������
#define betta_ 1.0  // ������

#define polusum false  // ���� �� ��������� ���������� ��� ������ ����
#define mu_statistic false  // ������� �� ���������� ����� �� �����
#define func_stat true     // ������� �� ������� ������������� � ������� �� ����� �������?
#define R_stat 80.0        // ����� �� ������� ������� ����������
#define Al_stat 54         // ������ �������� �� ����, �� ������� ������ ����


#define L_Igor 200.0         // � ������������� ����� �� ������ ������ �������� ������� ���������
#define k_Igor 1200         // � ������������� ����� �� ������ ������ �������� ������� ���������



// �������� �������� ���������� ���� �����

#define pred(i, j, al) (Mu[i][j] * mu_start * 0.3 * sin(al))
#define pred2(i, j, area, k, l) (min(Mu_[area][i][j],Mu_[area][k][l]))

//#define pred(i, j, al) (Mu[i][j] * sin(al))

inline double polar_angle(const double& x, const double& y)
{
	if (fabs(x) + fabs(y) < 0.000001 / RR_)
	{
		return 0.0;
	}

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

inline void spherical_skorost(const double& x, const double& y, const double& z, const double& Vx,//
	const double& Vy, const double& Vz, double& Vr, double& Vphi, double& Vtheta)
{
	double r_1 = sqrt(x * x + y * y + z * z);
	double the_1 = acos(z / r_1);
	double phi_1 = polar_angle(x, y);

	Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
	Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
	Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
}

inline void dekard_skorost(const double& x, const double& y, const double& z,
	const double& Vr, const double& Vphi, const double& Vtheta, double& Vx,
	double& Vy, double& Vz)
{
	double r_2 = sqrt(x * x + y * y + z * z);
	double the_2 = acos(z / r_2);
	double phi_2 = polar_angle(x, y);

	Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
	Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
	Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
}

inline void fff_velocity(const double& a, const double& b) // inline 
{
	double c = a + b + sin(a) + sin(b) + log(a * b);
}

inline int sgn(const double& val) {
	return (double(0) < val) - (val < double(0));          // ���������������� ����� ������� ����������
}

inline double sign(const double& x)
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

inline double minmod(const double& x, const double& y)
{
	if (sgn(x) + sgn(y) == 0)
	{
		return 0.0;
	}
	else
	{
		return   ((sgn(x) + sgn(y)) / 2.0) * min(fabs(x), fabs(y));  ///minmod
		//return (2*x*y)/(x + y);   /// vanleer
	}
}


inline float time_to_peregel_axis(const double& y0, const double& z0, const double& Vy, const double& Vz)
// ������ ��� ���� �����������, �.�. ��� ��������� ��������
{
	float a = -Vz * Vz * y0 * y0 + 2.0 * Vy * Vz * y0 * z0 - Vy * Vy * z0 * z0;
	if (a > 0.0)
	{
		// float t1 = (-Vy * y0 - Vz * z0 - sqrt(a)) / (Vy * Vy + Vz * Vz);
		float t2 = (-Vy * y0 - Vz * z0 + sqrt(a)) / (Vy * Vy + Vz * Vz);
		if (t2 > 0.0)
		{
			return t2;
		}
	}
	else
	{
		return 0.0;
	}

	return 0.0;
}


inline double MyRandom(int& s1, int& s2, int& s3)
{
	int ic15 = 32768, ic10 = 1024;
	int mz = 710, my = 17784, mx = 11973;
	double xi = 9.0949470177292824E-13, c = 1.073741824E9;
	double b;
	int i13, i12, i11, ii;
	i13 = mz * s1 + my * s2 + mx * s3;
	i12 = my * s1 + mx * s2;
	i11 = mx * s1;
	ii = i11 / ic15;
	i12 = i12 + ii;
	s1 = i11 - ic15 * ii;
	ii = i12 / ic15;
	i13 = i13 + ii;
	s2 = i12 - ic15 * ii;
	s3 = i13 % ic10;
	b = xi * (c * s3 + ic15 * s2 + s1);
	return b;
}
