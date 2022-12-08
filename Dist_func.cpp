#include "Dist_func.h"

Dist_func::Dist_func(const unsigned short& a1, const unsigned short& a2, const unsigned short& a3, //
	const double& c1, const double& d1, //
	const double& c2, const double& d2, //
	const double& c3, const double& d3)
{
	this->n1 = a1;
	this->n2 = a2;
	this->n3 = a3;
	this->name = "null_name";

	this->V = new double** [this->n1];
	for (int i = 0; i < this->n1; i++)
	{
		V[i] = new double* [this->n2];
	}

	for (int i = 0; i < this->n1; i++)
	{
		for (int j = 0; j < this->n2; j++)
		{
			V[i][j] = new double [this->n3];
		}
	}

	for (int i = 0; i < this->n1; i++)
	{
		for (int j = 0; j < this->n2; j++)
		{
			for (int k = 0; k < this->n3; k++)
			{
				V[i][j][k] = 0.0;
			}
		}
	}

	this->c1 = c1;
	this->d1 = d1;

	this->c2 = c2;
	this->d2 = d2;

	this->c3 = c3;
	this->d3 = d3;

	if ((a1 % 2 != 0) || (a2 % 2 != 0) || (a3 % 2 != 0))
	{
		cout << "ERROR nreufvyv34753247  27637424  " << endl;
	}

	this->a1 = pow(this->d1 + 1.0, 2.0 / this->n1);
	this->a2 = pow(this->d2 + 1.0, 2.0 / this->n2);
	this->a3 = pow(this->d3 + 1.0, 2.0 / this->n3);
}

bool Dist_func::call_name(string name)
{
	this->name = name;
	return true;
}

bool Dist_func::Add_point(const double& V1, const double& V2, const double& V3, const double& mu)
{
	// Функция добавляет точку в массив

	//сначала нужно определить в какую ячейку массива записать число.
	unsigned short i = 0, j = 0, k = 0;

	// double a = pow(this->d1 + 1.0 - this->c1, 2.0 / this->n1);

	i = (int)(log(fabs(V1 - this->c1) + 1) / log(this->a1));
	if (i >= this->n1 / 2)
	{
		i = this->n1 / 2 - 1;
	}
	if (V1 - this->c1 < 0.0)
	{
		i = i + this->n1 / 2;
	}

	//cout << this->a1 << " " << i << endl;

	j = (int)(log(fabs(V2 - this->c2) + 1) / log(this->a2));
	if (j >= this->n2 / 2)
	{
		j = this->n2 / 2 - 1;
	}
	if (V2 - this->c2 < 0.0)
	{
		j = j + this->n2 / 2;
	}

	k = (int)(log(fabs(V3 - this->c3) + 1) / log(this->a3));
	if (k >= this->n1 / 2)
	{
		k = this->n1 / 2 - 1;
	}
	if (V3 - this->c3 < 0.0)
	{
		k = k + this->n3 / 2;
	}

	V[i][j][k] += mu;
	return true;
}

bool Dist_func::normir(const double& ccc)
{
	for (int i = 0; i < this->n1; i++)
	{
		for (int j = 0; j < this->n2; j++)
		{
			for (int k = 0; k < this->n3; k++)
			{
				this->V[i][j][k] = this->V[i][j][k] * ccc;
			}
		}
	}

	// Нормировка из-за неравномерности функции распределения

	for (int i = 0; i < this->n1; i++)
	{
		double l1 = this->c1 - 1.0 + pow(this->a1, i);
		double r1 = this->c1 - 1.0 + pow(this->a1, i + 1);
		if (i >= this->n1 / 2)
		{
			l1 = -(this->c1 - 1.0 + pow(this->a1, i - this->n1 / 2));
			r1 = -(this->c1 - 1.0 + pow(this->a1, i + 1 - this->n1 / 2));
		}

		for (int j = 0; j < this->n2; j++)
		{
			double l2 = this->c2 - 1.0 + pow(this->a2, j);
			double r2 = this->c2 - 1.0 + pow(this->a2, j + 1);
			if (j >= this->n2 / 2)
			{
				l2 = -(this->c2 - 1.0 + pow(this->a2, j - this->n2 / 2));
				r2 = -(this->c2 - 1.0 + pow(this->a2, j + 1 - this->n2 / 2));
			}

			for (int k = 0; k < this->n3; k++)
			{
				double l3 = this->c3 - 1.0 + pow(this->a3, k);
				double r3 = this->c3 - 1.0 + pow(this->a3, k + 1);
				if (k >= this->n3 / 2)
				{
					l3 = -(this->c3 - 1.0 + pow(this->a3, k - this->n3 / 2));
					r3 = -(this->c3 - 1.0 + pow(this->a3, k + 1 - this->n3 / 2));
				}

				double S = fabs((r1 - l1) * (r2 - l2) * (r3 - l3));
				this->V[i][j][k] = this->V[i][j][k] / S;
			}
		}
	}


	return true;
}

bool Dist_func::print_1d(int koord)
{
	ofstream fout;
	string name_f = to_string(koord) + "_Dist_func_" + this->name + ".txt";
	fout.open(name_f);

	fout << "TITLE = \"HP\"  VARIABLES = \"v\", \"f\"," << "ZONE T = \"HP\"" << endl;

	if (koord == 1)
	{
		for (int i = 0; i < this->n1 / 2; i++)
		{
			double l = this->c1 - 1.0 + pow(this->a1, i);
			double r = this->c1 - 1.0 + pow(this->a1, i + 1);
			double S = 0.0;

			for (int j = 0; j < this->n2; j++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << 0.5 * (l + r) << " " << S << endl;
		}

		for (int i = 0 + this->n1 / 2; i < this->n1; i++)
		{
			double l = - 1.0 + pow(this->a1, i - this->n1 / 2);
			double r = - 1.0 + pow(this->a1, i + 1 - this->n1 / 2);
			double S = 0.0;

			for (int j = 0; j < this->n2; j++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << this->c1 - 0.5 * (l + r) << " " << S << endl;
		}
	}
	else if (koord == 2)
	{
		for (int j = 0; j < this->n2 / 2; j++)
		{
			double l = this->c2 - 1.0 + pow(this->a2, j);
			double r = this->c2 - 1.0 + pow(this->a2, j + 1);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << 0.5 * (l + r) << " " << S << endl;
		}

		for (int j = 0 + this->n2 / 2; j < this->n2; j++)
		{
			double l = - 1.0 + pow(this->a2, j - this->n2 / 2);
			double r = - 1.0 + pow(this->a2, j + 1 - this->n2 / 2);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << this->c2 -0.5 * (l + r) << " " << S << endl;
		}
	}
	else if (koord == 3)
	{
		for (int k = 0; k < this->n3 / 2; k++)
		{
			double l = this->c3 - 1.0 + pow(this->a3, k);
			double r = this->c3 - 1.0 + pow(this->a3, k + 1);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int j = 0; j < this->n2; j++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << 0.5 * (l + r) << " " << S << endl;
		}

		for (int k = 0 + this->n3 / 2; k < this->n3; k++)
		{
			double l = - 1.0 + pow(this->a3, k - this->n3 / 2);
			double r = - 1.0 + pow(this->a3, k + 1 - this->n3 / 2);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int j = 0; j < this->n2; j++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << this->c3 -0.5 * (l + r) << " " << S << endl;
		}
	}

	return true;
}

bool Dist_func::print_3d(void)
{
	ofstream fout;
	string name_f = "3D_Dist_func_" + this->name + ".txt";
	fout.open(name_f);

	fout << "TITLE = \"HP\"  VARIABLES = \"v1\", \"v2\", \"v3\", \"f\"," << "ZONE T = \"HP\"" << endl;


	for (int i = 0; i < this->n1; i++)
	{
		double l1 = this->c1 - 1.0 + pow(this->a1, i);
		double r1 = this->c1 - 1.0 + pow(this->a1, i + 1);
		if (i >= this->n1 / 2)
		{
			l1 = -(this->c1 - 1.0 + pow(this->a1, i - this->n1 / 2));
			r1 = -(this->c1 - 1.0 + pow(this->a1, i + 1 - this->n1 / 2));
		}
		
		for (int j = 0; j < this->n2; j++)
		{
			double l2 = this->c2 - 1.0 + pow(this->a2, j);
			double r2 = this->c2 - 1.0 + pow(this->a2, j + 1);
			if (j >= this->n2 / 2)
			{
				l2 = -(this->c2 - 1.0 + pow(this->a2, j - this->n2 / 2));
				r2 = -(this->c2 - 1.0 + pow(this->a2, j + 1 - this->n2 / 2));
			}

			for (int k = 0; k < this->n3; k++)
			{
				double l3 = this->c3 - 1.0 + pow(this->a3, k);
				double r3 = this->c3 - 1.0 + pow(this->a3, k + 1);
				if (k >= this->n3 / 2)
				{
					l3 = -(this->c3 - 1.0 + pow(this->a3, k - this->n3 / 2));
					r3 = -(this->c3 - 1.0 + pow(this->a3, k + 1 - this->n3 / 2));
				}

				fout << 0.5 * (l1 + r1) << " " << 0.5 * (l2 + r2) << " " << 0.5 * (l3 + r3) << " " << this->V[i][j][k] << endl;
			}
		}
	}


	return true;
}
