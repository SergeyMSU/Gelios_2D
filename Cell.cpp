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
		/*if (this->contour[0]->y < 0.0001)
		{
			y = y * (4.0 / 3.0);
		}*/
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
		x = x / (1.1 * this->contour.size());
		y = y / (1.1 * this->contour.size());
	}
}

void Cell::Get_Center_posle(double& x, double& y)
{
	if (this->contour.size() == 4)
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x2;
			y += i->y2;
		}
		x = x / (1.0 * this->contour.size());
		y = y / (1.0 * this->contour.size());

		/*if (this->contour[0]->y2 < 0.0001)
		{
			y = y * (4.0 / 3.0);
		}*/
	}
	else
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x2;
			y += i->y2;
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

double Cell::Get_Volume_rotate(const double& angle)
{
	double A = 0.0;
	double G = 0.0;
	int k = 0;
	for(int i = 0; i < this->contour.size(); i++)
	{
		auto a = this->contour[i];
		k = i + 1;
		if (k == this->contour.size()) k = 0;
		auto b = this->contour[k];
		A = A + (a->x * b->y - a->y * b->x);
		G = G + (a->y + b->y) * (a->x * b->y - a->y * b->x);
	}
	A = A * 0.5;

	return fabs(A) * (1.0 / (6.0 * A)) * G * pi_ * angle / 180.0;
	
}

double Cell::Get_Volume_posle(void)
{
	if (this->contour.size() == 4)
	{
		double a1 = this->contour[0]->x2 - this->contour[2]->x2;
		double a2 = this->contour[0]->y2 - this->contour[2]->y2;
		double a = sqrt(kv(a1) + kv(a2));

		double b1 = this->contour[1]->x2 - this->contour[3]->x2;
		double b2 = this->contour[1]->y2 - this->contour[3]->y2;
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

		for (auto& i : this->contour)
		{
			x += i->x2;
			y += i->y2;
		}
		x = x / (1.0 * this->contour.size());
		y = y / (1.0 * this->contour.size());

		double a1, a2, a3, p;
		double S = 0.0;
		for (auto& i : this->Grans)
		{
			a1 = sqrt(kv(i->A->x2 - x) + kv(i->A->y2 - y));
			a2 = sqrt(kv(i->B->x2 - x) + kv(i->B->y2 - y));
			a3 = sqrt(kv(i->A->x2 - i->B->x2) + kv(i->A->y2 - i->B->y2));
			p = (a1 + a2 + a3) / 2.0;
			S = S + sqrt(p * (p - a1) * (p - a2) * (p - a3));
		}
		return S;
	}
	return 0.0;
}

double Cell::Get_Volume_posle_rotate(const double& angle)
{
	double A = 0.0;
	double G = 0.0;
	int k = 0;
	for (int i = 0; i < this->contour.size(); i++)
	{
		auto a = this->contour[i];
		k = i + 1;
		if (k == this->contour.size()) k = 0;
		auto b = this->contour[k];
		A = A + (a->x2 * b->y2 - a->y2 * b->x2);
		G = G + (a->y2 + b->y2) * (a->x2 * b->y2 - a->y2 * b->x2);
	}

	A = A * 0.5;

	return fabs(A) * (1.0 / (6.0 * A)) * G * pi_ * angle / 180.0;

}


bool Cell::belong(const double& x, const double& y)
{
	for (auto& i : this->Grans)
	{
		if (i->belong(x, y) == false)
		{
			return false;
		}
	}
	return true;
}

void Cell::renew(void)
{
	double d = 100000000000.0;
	for (auto& i : this->Grans)
	{
		d = min(d, i->Get_lenght());
	}
	d = min(d, sqrt(kv(this->contour[0]->x - this->contour[2]->x) + kv(this->contour[0]->y - this->contour[2]->y)));
	this->L = d;
}
