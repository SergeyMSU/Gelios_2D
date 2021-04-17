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
	this->Sosed_down = nullptr;
	this->Sosed_up = nullptr;
}


void Gran::Get_Center(double& x, double& y)
{
	x = (this->A->x + this->B->x) / 2.0;
	y = (this->A->y + this->B->y) / 2.0;
}

void Gran::Get_Center_posle(double& x, double& y)
{
	x = (this->A->x2 + this->B->x2) / 2.0;
	y = (this->A->y2 + this->B->y2) / 2.0;
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
		if (par.u >= 0.0)
		{
			par.u = -0.01;
		}
	}
	else if (this->type == Upper_wall)  // Не надо менять
	{
		par = this->Master->par[i];
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
		double ro = (389.988 * 389.988) / (chi_ * chi_);
		double P_E = ro * chi_ * chi_ / (ggg * 5.0 * 5.0);
		par = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x / dist, chi_ * y / dist, ro / (dist * dist),//
		0.000001, 0.0000001, chi_real* x / dist, chi_real* y / dist, par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
			par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
			par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4};


		//cout << par.ro << " " << dist << endl;
	}
	else if (this->type == Input)
	{
		auto par2 = this->Master->par[i];
		par = {1.0, 1.0, Velosity_inf, 0.0, 100.0, par2.ro_H1, par2.p_H1, par2.u_H1, par2.v_H1,//
		par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		1.0, 0.5, Velosity_inf, 0.0};
	}

	return;
}


void Gran::Get_par_TVD(Parametr& par, int i)  // Здесь задаются граничные условия
// Сносит значение в текущей ячейке на текущую гарнь
{
	if (this->type == Usualy)  // Не надо менять
	{
		if (this->Sosed_down == nullptr)
		{
			par = this->Master->par[i];
			return;
		}
		double x, y, x2, y2, dist1, dist2, dist3;
		this->Get_Center(x, y);
		auto par2 = this->Master->par[i];
		auto par1 = this->Sosed_down->par[i];
		auto par3 = this->Sosed->par[i];
		this->Master->Get_Center(x2, y2);
		dist2 = sqrt(kv(x - x2) + kv(y - y2));
		this->Sosed->Get_Center(x2, y2);
		dist3 = sqrt(kv(x - x2) + kv(y - y2));
		this->Sosed_down->Get_Center(x2, y2);
		dist1 = sqrt(kv(x - x2) + kv(y - y2));

		par.ro = linear(-dist1, par1.ro, -dist2, par2.ro, dist3, par3.ro, 0.0);
		par.p = linear(-dist1, par1.p, -dist2, par2.p, dist3, par3.p, 0.0);
		par.u = linear(-dist1, par1.u, -dist2, par2.u, dist3, par3.u, 0.0);
		par.v = linear(-dist1, par1.v, -dist2, par2.v, dist3, par3.v, 0.0);
		par.Q = linear(-dist1, par1.Q, -dist2, par2.Q, dist3, par3.Q, 0.0);
		if (par.ro <= 0.0)
		{
			cout << x << " " << y << " " << x2 << " " << y2 << " " << par1.ro << " " << par2.ro << " " << par3.ro << endl;
			exit(-1);
		}

		if (par.p <= 0.0)
		{
			cout << x << " " << y << " " << x2 << " " << y2 << " " << par1.p << " " << par2.p << " " << par3.p << endl;
			exit(-1);
		}
		par.ro_H1 = linear(-dist1, par1.ro_H1, -dist2, par2.ro_H1, dist3, par3.ro_H1, 0.0);
		par.p_H1 = linear(-dist1, par1.p_H1, -dist2, par2.p_H1, dist3, par3.p_H1, 0.0);
		par.u_H1 = linear(-dist1, par1.u_H1, -dist2, par2.u_H1, dist3, par3.u_H1, 0.0);
		par.v_H1 = linear(-dist1, par1.v_H1, -dist2, par2.v_H1, dist3, par3.v_H1, 0.0);
		par.ro_H2 = linear(-dist1, par1.ro_H2, -dist2, par2.ro_H2, dist3, par3.ro_H2, 0.0);
		par.p_H2 = linear(-dist1, par1.p_H2, -dist2, par2.p_H2, dist3, par3.p_H2, 0.0);
		par.u_H2 = linear(-dist1, par1.u_H2, -dist2, par2.u_H2, dist3, par3.u_H2, 0.0);
		par.v_H2 = linear(-dist1, par1.v_H2, -dist2, par2.v_H2, dist3, par3.v_H2, 0.0);
		par.ro_H3 = linear(-dist1, par1.ro_H3, -dist2, par2.ro_H3, dist3, par3.ro_H3, 0.0);
		par.p_H3 = linear(-dist1, par1.p_H3, -dist2, par2.p_H3, dist3, par3.p_H3, 0.0);
		par.u_H3 = linear(-dist1, par1.u_H3, -dist2, par2.u_H3, dist3, par3.u_H3, 0.0);
		par.v_H3 = linear(-dist1, par1.v_H3, -dist2, par2.v_H3, dist3, par3.v_H3, 0.0);
		par.ro_H4 = linear(-dist1, par1.ro_H4, -dist2, par2.ro_H4, dist3, par3.ro_H4, 0.0);
		par.p_H4 = linear(-dist1, par1.p_H4, -dist2, par2.p_H4, dist3, par3.p_H4, 0.0);
		par.u_H4 = linear(-dist1, par1.u_H4, -dist2, par2.u_H4, dist3, par3.u_H4, 0.0);
		par.v_H4 = linear(-dist1, par1.v_H4, -dist2, par2.v_H4, dist3, par3.v_H4, 0.0);
	}
	else if (this->type == Extern)  // Не надо менять
	{
		par = this->Master->par[i];
		if (par.u >= -2.5)
		{
			par.u = -2.5;
		}
	}
	else if (this->type == Upper_wall)  // Не надо менять
	{
		par = this->Master->par[i];
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
		double ro = (389.988 * 389.988) / (chi_ * chi_);
		double P_E = ro * chi_ * chi_ / (ggg * 5.0 * 5.0);
		par = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x / dist, chi_ * y / dist, ro / (dist * dist),//
		0.000001, 0.000001, chi_real * x / dist, chi_real * y / dist, par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
			par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
			par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4 };


		//cout << par.ro << " " << dist << endl;
	}
	else if (this->type == Input)
	{
		auto par2 = this->Master->par[i];
		par = { 1.0, 1.0, Velosity_inf, 0.0, 100.0, par2.ro_H1, par2.p_H1, par2.u_H1, par2.v_H1,//
		par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		1.0, 0.5, Velosity_inf, 0.0 };
	}

	return;
}


void Gran::Get_par_TVD_radial(Parametr& par, int i)  // Здесь задаются граничные условия
// Сносит значение в текущей ячейке на текущую гарнь
{
	if (this->type == Usualy)  // Не надо менять
	{
		if (this->Sosed_down == nullptr)
		{
			par = this->Master->par[i];
			return;
		}
		double x, y, x2, y2, dist1, dist2, dist3;
		this->Get_Center(x, y);
		double phi = polar_angle(x, y);
		auto par2 = this->Master->par[i];
		auto par1 = this->Sosed_down->par[i];
		auto par3 = this->Sosed->par[i];
		this->Master->Get_Center(x2, y2);
		double phi2 = polar_angle(x2, y2);
		dist2 = sqrt(kv(x - x2) + kv(y - y2));
		this->Sosed->Get_Center(x2, y2);
		dist3 = sqrt(kv(x - x2) + kv(y - y2));
		double phi3 = polar_angle(x2, y2);
		this->Sosed_down->Get_Center(x2, y2);
		double phi1 = polar_angle(x2, y2);
		dist1 = sqrt(kv(x - x2) + kv(y - y2));

		double fr1 = par1.u_H1 * cos(phi1) + par1.v_H1 * sin(phi1);
		double ff1 = -par1.u_H1 * sin(phi1) + par1.v_H1 * cos(phi1);

		double fr2 = par2.u_H1 * cos(phi2) + par2.v_H1 * sin(phi2);
		double ff2 = -par2.u_H1 * sin(phi2) + par2.v_H1 * cos(phi2);

		double fr3 = par3.u_H1 * cos(phi3) + par3.v_H1 * sin(phi3);
		double ff3 = -par3.u_H1 * sin(phi3) + par3.v_H1 * cos(phi3);

		
		double fr = linear(-dist1, fr1, -dist2, fr2, dist3, fr3, 0.0);
		par.u_H1 = fr * cos(phi) - ff2 * sin(phi);
		par.v_H1 = fr * sin(phi) + ff2 * cos(phi);
	}

	return;
}