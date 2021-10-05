#pragma once

#include "Setka.h"
#include "Rail.h"
#include "Point.h"
#include "Cell.h"
#include "Gran.h"
#include "Solvers.h"
#include "sensor.h"
#include "sensor2.h"
#include <mutex>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <algorithm>


double sign(const double& x);
double minmod(const double& x, const double& y);
double polar_angle(const double& x, const double& y);
double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& x3, const double& t3, const double& y);
double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& y);
void polar_perenos(const double& x1, const double& y1, const double& x2, const double& y2, double& u, double& v);
void dekard_skorost(double x, double y, double z, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz);
void dekard_skorost2(double r2, double the_2, double phi_2, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz);
void spherical_skorost(const double& x, const double& y, const double& z, const double& Vx,//
	const double& Vy, const double& Vz, double& Vr, double& Vphi, double& Vtheta);
//double max(const double& x, const double& y);

#define ga (5.0/3.0)          // ���������� ��������
#define ggg (5.0/3.0)
#define g1 (ga - 1.0)
#define kv(x) ( (x)*(x) )
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))

#define pi_ 3.14159265358979323846
#define sqrtpi_ 1.77245385

#define Velosity_inf -2.54127 //-4.54127 //-2.54127
#define M_inf 1.968
#define chi_ 2.0 //36.1059 // 36.1059
#define chi_real 36.1059
#define kurant  0.3
#define Kn_  242.944										// ����� ��������
#define a_2 0.10263
#define n_p_LISM_ (3.0) 
#define n_H_LISM_ (1.0)
#define sigma(x) (kv(1.0 - a_2 * log(x)))               // ���������������� ������� �����������
//#define sigma2(x, y) (kv(1.0 - (a_2/(1.0 - a_2 * log(y))) * log(x)))  // ��� ������� ���������������� �������� �� cp
#define sigma2(x, y) (kv(1.0 - a_2 * log(x * y)))  // ��� ������� ���������������� �������� �� cp


#define n_inner 10  // ������� ����� �� ���������� ����
#define m_inner 16  // ������� ����� �� ���������� ����  (�� ����� ������ ���������� ���� ���������)
#define zazor 1.3   // ����� ������ ����� ��������� � ����������� ������� 
#define zazor2 2.0   // ����� ������ ����� ������� ������ � ����������� ������� 
#define s_k 0.03   // ������ ������� 4 � ������������
#define alpha_rot 1.0  // ���� �������� (������ ������), �� ���� � ���� ������� �� �������� ���� � � ������
#define R1_ 1.0
#define R11_ 27.0
#define R111_ 37.0
#define R2_ 100.0
#define R3_ 150.0
#define R4_ 600.0
#define R5_ 2200.0
#define Left_ -1500
#define H_pow 0.5 //0.45  // ���������� �������� ��������� ��� ������ ���������� ��������
#define I_ 7 //5 // ���������� ������ ������
#define geo_accur 0.1    // �������� ��������� ����� ����������
#define weight_ 0.25     // ���� ��� ������
#define gam(i) ( 1.0/(kv(R5_ - 2.0)/kv(Ri[(i)]) - 1.0) )
#define R_m 4.0   // ��� ��������� ���������� �� ��� ���������� ������ � ����� ����� �������
#define Ur_m 0.5   // �������� �����������