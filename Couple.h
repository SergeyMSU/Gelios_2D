#pragma once
#include <vector>
#include <map>
#include "Help.h"
using namespace std;
class Point;
class Gran;


class Couple
{
public:
	Cell* A1;
	Cell* A2;
	double dist;       // ���������� ����� ������� ��������
	double n1;
	double n2;
	double n3;
	double dvig1 = 0.0;                   // ������ ������ ��� ������ �����
	double dvig2 = 0.0;                   // ������ ������ ��� ������ �����
	double dvig3 = 0.0;                   // ������ ������ ��� ������ �����
	double mo1 = 0.0;                   // ������ ������ ��� ������ �����
	double mo2 = 0.0;                   // ������ ������ ��� ������ �����
	double mo3 = 0.0;                   // ������ ������ ��� ������ �����
	double tension_x = 0.0;
	double tension_y = 0.0;
	double tension_z = 0.0;
	double d_sosed;
	int number = -1;
	bool extern_boundary = false;            // �������� �� ������ � �������� (������� �� ����������, �� �����.

	Couple(Cell* A1, Cell* A2, const double& dist);

	void Resolve(void); // ������������� ���������� ����� ��������
	void get_centr(double& x, double& y, double& z);
	void get_normal(double& x, double& y, double& z);  // ��������� ������� ���� �� ������� A1 � A2, � �� ������ ���� �.
	void orient(void); // ������������� ���� �� � ������� (������ ���������� ����)
	void move(const double& m1, const double& m2, const double& m3);


};

