#pragma once
#include "Help.h"
#include <vector>

class Cell;

enum Point_type  // Тип грани нужен для граничных условий
{
	P_Usualy,  
	P_U_1,
	P_U_2,
	P_U_3,
	P_U_4,
	P_U_5,
	P_Contact,
	P_Inner_shock,
	P_Outer_shock,
	P_Inner_Boandary,
	P_Outer_Boandary,
	P_No,
	P_Null,
};

class Point
{
public:
	double x;
	double y;
	double Vx;
	double Vy;
	int number;
	Point_type type;


	vector < Cell * > my_cell;

	Point(const double& x, const double& y);
};

