#pragma once
#include "Help.h"

class MKmethod
{
public:
	// ����������
	double R_[8];       // ������� ��������� ��� (� �� ����������)
	double alpha_[8];       // ������������ ��� ��� �����������(� �� ����������)
	double gam_[8];       // ��� ����� ��� ������� � �������
	double A0_;         // ������� ��������� ����� ����������� 
	double A1_;         // ������� ��������� ����� ����������� 
	double A2_;         // ������� ��������� ����� ����������� 
	int num_area;      // ���������� ��� (������ ���������� �������������� ���)
	double Int_[501];
	double Int_002[51][50];
	double Int_00625[51][50];
	double Int_02[51][50];
	double Int_055[51][50];

	MKmethod(void);


    // ����� ���� ��������� ���������

	// �������� ��������� ���������� �������
	bool Init_Parametrs(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_);
	bool Init_Parametrs_2s(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_);
	// ���������� false, ���� �� ����� ��������� �������� ����
	bool Init_Parametrs_mini(Sensor* sens, double& Wt_, double& Wp_, double& Wr_, double& X_);
	bool Init_Parametrs2(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_);
	int Init(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, double& X_);
	double FF(const double& gam, const double& Yr);

	// �������� �������� ��� �����������
	bool Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // ������ �������� ������
	bool Change_Velosity2(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // ����� ������ �������� ������
	bool Change_Velosity3(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // ��� �������� (�� ������� ������) � ��������� ����������� �����
	bool Change_Velosity4(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // ��� �������� (�� ������� ������) � ��������� ����������� �����
	// ���������� false, ���� �� ����� ��������� �������� ����

	double for_Wr_1(const double& Z, const double& gam, const double& ur);
	double H_Wr_1(const double& gam1, const double& gam2, const double& V, const double& ur, const double& p, const double& ksi);
	double for_Wr_2(const double& Z, const double& gam, const double& ur, const double& ut);
	double H_Wr_2(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& p, const double& ksi);


	void TEST(void);   // ������� - ����������� ������ ��������������� �������, ����� ���� ��� �������
	double play_mho(Sensor* sens, const double& c);
	double play_mho2(Sensor* sens, const double& c, const double& p, const double& p1);
	// ����������� \mho 
	double norm_mho(const double& c);  // ��� ������� �������
	double norm_mho2(const double& c); // ��� ������� Struvel
	double mho_H1(const double& c, const double& k, const double& ksi);
	double mho_H2(const double& c, const double& k, const double& ksi);
	double h_mho(const double& x, const double& c);


	// ������������ �����������

	double Int_cp_1(const double& x);
	double Int_cp_2(const double& x);

	// ��������� ���������� �����������

	double Get_Int(const double& uu);
	double Get_Int002(const double& Ur, const double& uu);
	double Get_Int00625(const double& Ur, const double& uu);
	double Get_Int02(const double& Ur, const double& uu);
	double Get_Int055(const double& Ur, const double& uu);

	// ��������� ���������� ����������� (�����)
	double int_1(const double& x, const double& cp);
	double int_2(const double& x, const double& cp);
	double int_3(const double& x, const double& cp);
	double int_1_f1(const double& x);
	double int_1_f2(const double& x);
	double int_1_f3(const double& x);
	double int_2_f1(const double& x);
	double int_2_f2(const double& x);
	double int_2_f3(const double& x);
	double int_3_f1(const double& x);
	double int_3_f2(const double& x);
	double int_3_f3(const double& x);

	double Lin_Interpolate(const double& x1, const double& y1, const double& x2, const double& y2, const double& x);

	// ��� ��������� �� �����
	double A0(const double& Y);  // ����� ������� ������������� (�������� ��������� � ������)    ���������!
	double G(const double& gam, const double& Y);  // ����� ���������     ���������!
	double F(const double& X, const double& gam, const double& Y);  // ������� ������������� ��� X   // ���������!
	double F0(const double& X, const double& Y);  // ������� ������������� ��� X   // ���������!
	double FI(const double& Z, const double& X, const double& gam, const double& Y);   // ���������!
	double R(const double& X, const double& Y);              // ���������!

	//��� �������� �����
	double B0(const double& k);
	double R_sphere(const double& k, const double& Z);
	double H_sphere(const double& k, const double& Z, const double& B, const double& ksi);

	// ��� ��������� ����������� ������ �������
	double f2(const double& V, const double& gam, const double& ur, const double& ut);
	double f2k(const double& V, const double& gam, const double& ur, const double& ut);  // ������ ��� ����� ����

	// ��������������� �������
	double Hx(const double& gam1, const double& gam2, const double& X, const double& Y, const double& ksi);
	double Hx2(const double& gam1, const double& gam2, const double& X, const double& Y, const double& ksi);
	double Hwr(const double& gam1, const double& gam2, const double& Z, const double& X, const double& Y, const double& ksi);

	double Hvr(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi);
	double Hvrk(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi);


private:
	

};

