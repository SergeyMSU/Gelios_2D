#include "Cell.h"
#include <math.h>


Cell::Cell(Point* A, Point* B, Point* C, Point* D)
{
	this->contour.push_back(A);
	this->contour.push_back(B);
	this->contour.push_back(C);
	this->contour.push_back(D);
	this->number = -1;
	this->type = C_no;
}

Cell::Cell(void)
{
	this->number = -1;
	this->type = C_no;
}

void Cell::Get_Center(double& x, double& y)
{
	if (this->contour.size() == 4)
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x;
			y += i->y;
		}
		x = x / (1.0 * this->contour.size());
		y = y / (1.0 * this->contour.size());
	}
	else
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x;
			y += i->y;
		}
		x = x / (1.35 * this->contour.size());
		y = y / (1.35 * this->contour.size());
	}
}

double Cell::Get_Volume(void)
{
	if (this->contour.size() == 4)
	{
		double a1 = this->contour[0]->x - this->contour[2]->x;
		double a2 = this->contour[0]->y - this->contour[2]->y;
		double a = sqrt(kv(a1) + kv(a2));

		double b1 = this->contour[1]->x - this->contour[3]->x;
		double b2 = this->contour[1]->y - this->contour[3]->y;
		double b = sqrt(kv(b1) + kv(b2));

		a1 = a1 / a;
		a2 = a2 / a;
		b1 = b1 / b;
		b2 = b2 / b;
		double c = a1 * b1 + a2 * b2;
		double s = sqrt(1 - kv(c));
		return 0.5 * a * b * s;
	}
	else
	{
		double x = 0.0;
		double y = 0.0;
		
		this->Get_Center(x, y);

		double a1, a2, a3, p;
		double S = 0.0;
		for (auto& i : this->Grans)
		{
			a1 = sqrt(kv(i->A->x - x) + kv(i->A->y - y));
			a2 = sqrt(kv(i->B->x - x) + kv(i->B->y - y));
			a3 = sqrt(kv(i->A->x - i->B->x) + kv(i->A->y - i->B->y));
			p = (a1 + a2 + a3) / 2.0;
			S = S +  sqrt(p * (p - a1) * (p - a2) * (p - a3));
		}
		return S;
	}
	return 0.0;
}

double Cell::Get_Volume_posle(const double& time)
{
	double x1 = this->contour[0]->x + this->contour[0]->Vx * time;
	double y1 = this->contour[0]->y + this->contour[0]->Vy * time;

	double x2 = this->contour[1]->x + this->contour[1]->Vx * time;
	double y2 = this->contour[1]->y + this->contour[1]->Vy * time;

	double x3 = this->contour[2]->x + this->contour[2]->Vx * time;
	double y3 = this->contour[2]->y + this->contour[2]->Vy * time;

	double x4 = this->contour[3]->x + this->contour[3]->Vx * time;
	double y4 = this->contour[3]->y + this->contour[3]->Vy * time;

	double a1 = x1 - x3;
	double a2 = y1 - y3;
	double a = sqrt(kv(a1) + kv(a2));

	double b1 = x2 - x4;
	double b2 = y2 - y4;
	double b = sqrt(kv(b1) + kv(b2));

	a1 = a1 / a;
	a2 = a2 / a;
	b1 = b1 / b;
	b2 = b2 / b;
	double c = a1 * b1 + a2 * b2;
	double s = sqrt(1 - kv(c));
	return 0.5 * a * b * s;
}
