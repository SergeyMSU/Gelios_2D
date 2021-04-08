#pragma once

#include "Help.h"
#include <vector>
#include <string>

class Rail;
class Point;
class Cell;
class Gran;
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

	vector <Cell*> All_Cells;

	vector <Gran*> Line_Contact;
	vector <Gran*> Line_Inner;       // ���������� �����
	vector <Gran*> Line_Outer;		 // ������� �����

	vector <Point*> Contact;
	vector <Point*> Inner;       // ���������� �����
	vector <Point*> Outer;		 // ������� �����

	int N1;
	int N2;
	int N3;
	int N4;
	int M1;
	int M2;
	int M3;
	int M4;


	Setka(int N1, int N2, int N3, int N4, int M1, int M2, int M3, int M4);
	Setka();

	void Save_Setka_ALL_ALPHA(string name); // ������� � ������� ������� ���������� ������ �����
	void Download_Setka_ALL_ALPHA(string name);


	void Print_point();      // �������� ����� � ����� (�� ������, � ����)
	void Print_Gran();
	void Print_Gran_type();
	void Print_cell(void);   // �������� ������ (�� ����� ������ ����� ����������, �.�. ������ ������ ������ ������ � ����� �������������
	void Print_cell_type(void);
	void Print_cell2(void);  // � ������� �� ���������� �� ������� ������ �����, ����� ��� �������������� �������� ������������ ���������
	void Print_connect(void); // ����������, ��� ������� ������
	void Print_point_connect(void); // ����������, ��� ������� ���� � ��������
	void Print_Tecplot(void);

	void Save_G_D(void);
	void Download_G_D(void);
	void Save_G_D_5_komponent(void);
	void Download_G_D_5_komponent(void);

	void Proverka(void);   // �������� ������������ ���������� ����� (��������� � ������). ���� ����� ��������� ��� ����� ����� � ��������

	// �������� �����
	void Move_Setka_Calculate(const double& dt);
	void Move_surface(int ii);  // ���������� ��������� ������������   ii - ����� ��������� �������� par[ii]
	void Move_surface_hand(void);  // ������ �������� �����


	// ������� ��������

	void Go_stationary(int step);
	void Go_stationary_5_komponent(int step);
	void Go(int step); // ��������� ��������  ������� ��������
	void Go_5_komponent(int step);
	double HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
		const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
		double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp);
	void Init_conditions(void);

private:
	double polar_angle(const double& x, const double& y);


};

