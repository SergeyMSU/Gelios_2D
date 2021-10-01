#pragma once

#include "Help.h"
#include <vector>
#include <string>
#include <mutex>
#include "sensor2.h"

class Rail;
class Point;
class Cell;
class Gran;
class Sensor;
using namespace std;

class Setka
{
public:
	vector <Rail*> A_Rails;
	vector <Rail*> B_Rails;
	vector <Rail*> C_Rails;
	vector <Rail*> D_Rails;

	vector <Point*> All_Points;
	vector <Gran*> All_Gran;           // ����� ��������
	vector <Gran*> All_Gran_copy;      // ����� ��������� (������� � ������ �������)

	vector <Cell*> All_Cells;          // ��� ������
	vector <Cell*> All_Cells_Inner;    // ������ ������ ���������� ������� (��� �������� �������)

	vector <Gran*> Line_Contact;     // �������
	vector <Gran*> Line_Inner;       // ���������� �����
	vector <Gran*> Line_Outer;		 // ������� �����

	vector <Point*> Contact;	     // �������
	vector <Point*> Inner;		     // ���������� �����
	vector <Point*> Outer;		     // ������� �����

	vector <Cell*> Cell_sphere;      // ������ ������� 
	vector <Cell*> Cell_side;      // ������ �������
	vector <Cell*> Cell_back;      // ������ ������
	Cell* Cell_m;                  // ��� 4-�� ���� ������ ������

	vector<Sensor*> Sensors;
	vector<sensor2*> Sensors2;
	vector <double> Ri;            // �������������� ���� ��� ����������� �����-�����
	vector <double> Mu;            // �������������� ���� ��� ����������� �����-�����


	// ���� ���������� �����
	double mmu1;
	double mmu2;
	double mmu3;
	double mmu4;
	double mmu5;
	int mn1;
	int mn2;
	int mn3;
	int mn4;
	int mn5;



	int N1;
	int N2;
	int N3;
	int N4;
	int M1;
	int M2;
	int M3;
	int M4;


	int Number1;
	int Number2;
	int Number3;
	int Number4;
	int AllNumber;
	double sqv_1;
	double sqv_2;
	double sqv_3;
	double sqv_4;
	double sum_s;

	mutex Smut;


	Setka(int N1, int N2, int N3, int N4, int M1, int M2, int M3, int M4);
	Setka();
	void Inizialization(void);

	void Save_Setka_ALL_ALPHA(string name); // ������� � ������� ������� ���������� ������ �����
	void Download_Setka_ALL_ALPHA(string name);
	void Download_Setka_ALL_ALPHA_2_0(string name);


	void Print_point();      // �������� ����� � ����� (�� ������, � ����)
	void Print_Gran();
	void Print_Gran_type();
	void Print_cell(void);   // �������� ������ (�� ����� ������ ����� ����������, �.�. ������ ������ ������ ������ � ����� �������������
	void Print_cell_type(void);
	void Print_cell2(void);  // � ������� �� ���������� �� ������� ������ �����, ����� ��� �������������� �������� ������������ ���������
	void Print_connect(void); // ����������, ��� ������� ������
	void Print_point_connect(void); // ����������, ��� ������� ���� � ��������
	void Print_Tecplot(void);
	void Print_Tecplot_MK(void);



	void Proverka(void);   // �������� ������������ ���������� ����� (��������� � ������). ���� ����� ��������� ��� ����� ����� � ��������

	// �������� �����
	void Move_Setka_Calculate(const double& dt);
	void Move_surface(int ii);  // ���������� ��������� ������������   ii - ����� ��������� �������� par[ii]
	void Move_surface_hand(void);  // ������ �������� �����

	void TVD_prepare(void);
	void M_K_prepare(void);
	void Print_TVD(void);

	// ������� ��������

	void Save_G_D(void);
	void Download_G_D(void);
	void Save_G_D_5_komponent(void);
	void Download_G_D_5_komponent(void);
	void Go_stationary_5_komponent_inner(int step);
	void Go_stationary_5_komponent_inner_MK(int step);  // ��������� ������� �� �����-�����
	void Go_stationary(int step);
	void Go_stationary_TVD(int step);
	void Go_stationary_5_komponent(int step);
	void Go(int step); // ��������� ��������  ������� ��������
	void Go_5_komponent(int step);
	void Go_5_komponent_MK(int step);  // ������ � ����������� �� �����-�����
	double HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
		const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
		double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp, bool nul_potok = false);
	void Init_conditions(void);


	// �����-�����

	void Init_Velosity(sensor2* sens, const double& A2, vector <double>& mu, vector <double>& Wt, vector <double>& Wp, vector <double>& Wr, const double& the);
	void Velosity_initial2(sensor2* s, double& Vx, double& Vy, double& Vz); //   ��� ������ �����
	void Velosity_initial(sensor2* s, double& Vx, double& Vy, double& Vz);  //  ��� ������ � ��������� ������� (4 ���)
	void Init_Pozision(sensor2* sens, const double& A1, double& phi, double& the);    // ������������� ��������� ��������� �� ���������
	double F_mk(const double& gamma, const double& Yr);
	void MK_start(void);
	Cell* Belong_point(int b, const double& x, const double& y);   // ������� ��������� ������, ������� ����������� �����
	void Fly_exchenge(sensor2* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, //
		double mu, const double& mu_0, bool ExCh);
	void Fly_exchenge_Split(sensor2* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
		const double& mu_0, bool ExCh, int zone);
	void Change_Velosity(sensor2* s, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, double& X, double& Y, double& Z, const double& cp);
	void Change_Velosity_Split(sensor2* s, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr, vector <double>& Wthe,//
		vector <double>& Wphi, vector <double>& mu_, const double& cp, const double& r, int I);
	double w_c_v_s(const double& r, const double& v, int i);
	void Fly_exchenge_Imit(sensor2* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
		 double KSI, double I_do, int area, const double& mu_start); // ������������ �����
	void Fly_exchenge_Imit_Korol(sensor2* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
		 int area, bool ExCh, const double& mu_start);  // ������ �������� ������� � ���� �������
	// ����������� � ������������ �� ����������
	double Velosity_1(const double& u, const double& cp);
	double Velosity_2(const double& u, const double& cp);
	double Velosity_3(const double& u, const double& cp);


private:

};

