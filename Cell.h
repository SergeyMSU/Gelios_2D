#pragma once
#include <vector>
#include "Help.h"
#include <mutex>

using namespace std;

enum Cell_type  // ��� ����� ����� ��� ��������� �������
{
	C_1,
	C_2,
	C_3,
	C_4,
	C_5,
	C_centr,
	C_no,
};

struct Parametr
{                    // ��������� ���������� ������ � �����!!!!
	double ro = 0.0;
	double p = 0.0;
	double u = 0.0;
	double v = 0.0;
	double Q = 0.0;
	double ro_H1 = 0.0;
	double p_H1 = 0.0;
	double u_H1 = 0.0;
	double v_H1 = 0.0;
	double ro_H2 = 0.0;
	double p_H2 = 0.0;
	double u_H2 = 0.0;
	double v_H2 = 0.0;
	double ro_H3 = 0.0;
	double p_H3 = 0.0;
	double u_H3 = 0.0;
	double v_H3 = 0.0;
	double ro_H4 = 0.0;
	double p_H4 = 0.0;
	double u_H4 = 0.0;
	double v_H4 = 0.0;
	double F_n = 0.0;
	double F_u = 0.0;
	double F_v = 0.0;
	double F_T = 0.0;
	double I_u = 0.0;
	double I_v = 0.0;
	double I_T = 0.0;
	double II_u = 0.0;
	double II_v = 0.0;
	double II_T = 0.0;
	double M_u = 0.0;   // ����������� ���������
	double M_v = 0.0;
	double M_T = 0.0;
	double H_n[4];
	double H_u[4];
	double H_v[4];
	double H_T[4];
	double H_n2[4];
	double H_u2[4];
	double H_v2[4];
	double H_T2[4];
	double H_n3[4];
	double H_u3[4];
	double H_v3[4];
	double H_T3[4];
	double H_uu[4];
	double H_vv[4];
	double H_uv[4];
	double H_uuu[4];
	double k_u = 0.0;
	double k_v = 0.0;
	double k_T = 0.0;
	double k_u2 = 0.0;
	double k_v2 = 0.0;
	double k_T2 = 0.0;
	double k_u3 = 0.0;
	double k_v3 = 0.0;
	double k_T3 = 0.0;
	int num_atoms[4];   // ����� ����������� � ������
	double w_m[7];       // ������� ���� �� ������
	double I1_mf[4];
	double I2_mf[4];
	double I3_mf[4];
	double I1_mc[4];
	double I2_mc[4];
	double I3_mc[4];

	double npui = 0.0;
	double Tpui = 0.0;
	double pp = 0.0;
	double ppui = 0.0;
	double divV = 0.0;
};

class Point;
class Gran;
class Sensor;

class Cell
{
public:
	double Potok[5];
	double Potok_H1[4];
	double Potok_H2[4];
	double Potok_H3[4];
	double Potok_H4[4];
	double L;                   // ����������� ������ ������
	Parametr par[2];
	vector <Point*> contour;    // ������������� ������������ ����� �� �����
	vector <Gran*> Grans;

	vector <double> fpui;
	double fpui_max = 1.0;      // ������������ �������� fpui  
	double Wmin = 0.0;
	double Wmax = 10.0;         // ������������ �������� ��� ������� fpui �� �������
	vector <double> nu_pui;     // ������� ����������� �� L, ��� L �� 0 �� 20 � ����� 100
	vector <double> nu2_pui;     // �������� �������� ����������� �� L, ��� L �� 0 �� 20 � ����� 100
	vector <double> nu3_pui;     // �������� �������� ����������� �� L, ��� L �� 0 �� 20 � ����� 100
	bool pui_ = false;                  // ����� �� ������ ������� ��� PUI?

	int number;
	int zona;
	int zona_alpha;
	Cell_type type;
	mutex mut;                 // ������� ��� ������ � ������
	double x_min;
	double x_max;
	double y_min;
	double y_max;

	double r_istoch;           // ������ ���������� (� ����� ����� ��� �����������, ����� ��� ����������������)
	Cell* Back;                // ������ �� ����������������
	Cell* Next;

    // ���� ���������� ����� ��� �����-�����
	bool axis_; // �������� �� ������ ��������� � ���� (��� ����� �������� � �����-�����
	double y_ax; // ������ ������ ��������� ������ ��� ����

	Cell(Point* A, Point* B, Point* C, Point* D);
	Cell(void);

	void Initial(void);

	void Get_Center(double& x, double& y);
	void Get_Center_posle(double& x, double& y);
	void Get_Center2(double& x, double& y);
	void Get_Center_posle2(double& x, double& y);
	double Get_Volume(void);
	double Get_Volume_rotate(const double& angle);
	double Get_Volume_posle(void);
	double Get_Volume_posle_rotate(const double& angle);

	double get_nu_pui(double L);
	double get_nu_pui2(double L);
	double get_nu_pui3(double L);
	double get_fpui(const double& W, const double& Wmin, const double& Wmax);
	bool Change_Velosity_PUI(Sensor* sens, const double& Vh1, const double& Vh2, const double& Vh3, //
		const double& Vp1, const double& Vp2, const double& Vp3, double& W1, double& W2, double& W3, int Nw, //
		const double& wmin, const double& wmax);
	bool Change_Velosity_PUI2(Sensor* sens, const double& Vh1, const double& Vh2, const double& Vh3, //
		const double& Vp1, const double& Vp2, const double& Vp3, double& W1, double& W2, double& W3, int Nw, //
		const double& wmin, const double& wmax);

	// ��������� ��������� �����-�����
	void Get_Sourse_MK1(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p);
	void Get_Sourse_MK2(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p);
	void Get_Sourse_MK3(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p);

	bool belong(const double& x, const double& y);
	void renew(void); // �������� �������� L

	void Calc_Sourse(void);  // ������� ����������� ���������
};

