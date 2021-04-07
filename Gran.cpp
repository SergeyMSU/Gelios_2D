#include "Gran.h"
#include <math.h>

Gran::Gran(Point* A, Point* B, Gran_type type)
{
	this->A = A;
	this->B = B;
	this->type = type;
	this->Gran_copy = nullptr;
	this->Master = nullptr;
	this->Sosed = nullptr;
	this->number = -1;
	this->main_gran = true;
}


void Gran::Get_Center(double& x, double& y)
{
	x = (this->A->x + this->B->x) / 2.0;
	y = (this->A->y + this->B->y) / 2.0;
}

void Gran::Get_normal(double& n1, double& n2)
{
	double t1 = this->B->x - this->A->x;
	double t2 = this->B->y - this->A->y;
	double n = sqrt(kv(t1) + kv(t2));
	t1 = t1 / n;
	t2 = t2 / n;
	n1 = t2;
	n2 = -t1;
	return;
}

double Gran::Get_square(void)
{
	return sqrt(kv(this->B->x - this->A->x) + kv(this->B->y - this->A->y));
}

void Gran::Get_par(Parametr& par, int i)  // Здесь задаются граничные условия
{

	if (this->type == Usualy)  // Не надо менять
	{
		par = this->Sosed->par[i];
	}
	else if (this->type == Extern)  // Не надо менять
	{
		par = this->Master->par[i];
		if (par.u >= -2.0)
		{
			par.u = -2.0;
		}
	}
	else if (this->type == Axis)  // Не надо менять
	{
		par = this->Master->par[i];
		par.v = -par.v;
		par.v_H1 = -par.v_H1;
		par.v_H2 = -par.v_H2;
		par.v_H3 = -par.v_H3;
		par.v_H4 = -par.v_H4;
	}
	else if (this->type == Inner_sphere)
	{
		auto par2 = this->Sosed->par[i];
		double x, y;
		this->Get_Center(x, y);
		double dist = sqrt(x * x + y * y);
		double ro = (389.988 * 389.988) / (chi_real * chi_real);
		double P_E = ro * chi_real * chi_real / (ggg * 5.0 * 5.0);
		par = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_real * x / dist, chi_real * y / dist, ro / (dist * dist),//
		0.000001, 0.000001, chi_real* x / dist, chi_real* y / dist, par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
			par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
			par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4};
	}
	else if (this->type == Input || this->type == Upper_wall)
	{
		auto par2 = this->Master->par[i];
		par = {1.0, 1.0, Velosity_inf, 0.0, 100.0, par2.ro_H1, par2.p_H1, par2.u_H1, par2.v_H1,//
		par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		1.0, 0.5, Velosity_inf, 0.0};
	}

	return;
}