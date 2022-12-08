#pragma once
#include "Help.h"
#include <string>
class Dist_func
{
public:
	// ����� ��� ���������� ���-������ ������� ������������� 
	unsigned short n1; // ������ ������� �� 2!
	unsigned short n2;
	unsigned short n3;
	double c1;          // ����� �������������
	double d1;          // ������ �������������
	double a1;
	double c2;
	double d2;
	double a2;
	double c3;
	double d3;
	double a3;
	string name;

	double*** V;

	Dist_func(const unsigned short& a1, const unsigned short& a2, const unsigned short& a3, //
		const double& c1, const double& d1, //
		const double& c2, const double& d2, //
		const double& c3, const double& d3);
	bool call_name(string name);         // ������ ��� (����� ����� ������)
	bool Add_point(const double& V1, const double& V2, const double& V3, const double& mu);
	bool normir(const double& ccc);

	bool print_1d(int koord); // �� ����� ���������� ����

	bool print_3d(void); // �� ����� ���������� ����

private:
};

