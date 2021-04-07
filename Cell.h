#pragma once
#include <vector>
#include "Help.h"

using namespace std;

enum Cell_type  // Тип грани нужен для граничных условий
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
{
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
};

class Point;
class Gran;

class Cell
{
public:
	double Potok[5];
	double Potok_H1[4];
	double Potok_H2[4];
	double Potok_H3[4];
	double Potok_H4[4];
	Parametr par[2];
	vector <Point*> contour;
	vector <Gran*> Grans;
	int number;
	Cell_type type;

	Cell(Point* A, Point* B, Point* C, Point* D);
	Cell(void);

	void Get_Center(double& x, double& y);
	double Get_Volume(void);
	double Get_Volume_posle(const double& time);
};

