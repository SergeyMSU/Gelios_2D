#include "MKmethod.h"

using namespace std;

MKmethod::MKmethod(void)
{
	double Yr = fabs(Velosity_inf);
	this->A0_ = this->A0(Yr);     // Сразу посчитаем константу, чтобы постоянно не считать. 
	this->A1_ = 1.0 + (1.0 + 1.0 / (2.0 * kv(Yr))) * erf(Yr) + exp(-kv(Yr)) / (sqrtpi_ * Yr);
	this->num_area = I_; // I_;
	this->R_[0] = 1.0/ RR_;                 // Здесь задаются радиусы
	this->R_[1] = 2.61 / RR_;
	this->R_[2] = 6.84 / RR_;
	this->R_[3] = 17.92 / RR_;
	this->R_[4] = 46.9 / RR_;
	this->R_[5] = 122.7 / RR_;
	this->R_[6] = 321.23 / RR_;                 // Здесь задаются радиусы
	this->R_[7] = 840.65 / RR_;


	for (int i = 0; i < this->num_area; i++)
	{
		this->gam_[i] = 1.0 / (kv(Rmax_ / this->R_[i]) - 1.0);
	}

	ifstream fout;
	fout.open("cp_1_osn.txt");
	ifstream fout2;
	fout2.open("py_cp_1_al_0.02.txt");
	ifstream fout3;
	fout3.open("py_cp_1_al_0.0625.txt");
	ifstream fout4;
	fout4.open("py_cp_1_al_0.2.txt");
	ifstream fout5;
	fout5.open("py_cp_1_al_0.55.txt");
	double a, b, c;

	for (int i = 0; i < 501; i++)
	{
		fout >> a >> b;
		this->Int_[i] = b;
	}

	/*for (int i = 0; i < 51; i++)
	{
		for (int j = 0; j < 50; j++)
		{
			fout2 >> a >> b >> c;
			this->Int_002[i][j] = c;
		}
	}

	for (int i = 0; i < 51; i++)
	{
		for (int j = 0; j < 50; j++)
		{
			fout3 >> a >> b >> c;
			this->Int_00625[i][j] = c;
		}
	}

	for (int i = 0; i < 51; i++)
	{
		for (int j = 0; j < 50; j++)
		{
			fout4 >> a >> b >> c;
			this->Int_02[i][j] = c;
		}
	}

	for (int i = 0; i < 51; i++)
	{
		for (int j = 0; j < 50; j++)
		{
			fout5 >> a >> b >> c;
			this->Int_055[i][j] = c;
		}
	}*/

}

double MKmethod::Get_Int(const double& uu)
{
	if (uu < 5.0)
	{
		int k = (int)(uu / 0.01);
		double a = this->Int_[k];
		double b = this->Int_[k + 1];
		return Lin_Interpolate(k * 0.01, a, (k + 1.0) * 0.01, b, uu);
	}
	else
	{
		return Lin_Interpolate(495 * 0.01, this->Int_[495], 500 * 0.01, this->Int_[500], uu);
	}
}

double MKmethod::Get_Int002(const double& Ur, const double& uu)
{
	if (Ur >= -5.0 && uu <= 4.9)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = (int)(uu / 0.1);
		//cout << "k1 , k2 = " << k1 << " " << k2 << endl;
		double a = this->Int_002[k1][k2 + 1];
		double b = this->Int_002[k1 + 1][k2 + 1];
		double c = this->Int_002[k1][k2];
		double d = this->Int_002[k1 + 1][k2];
		//cout << a << " " << b << " " << c << " " << d << endl;
		//cout << (k1) * 0.1 << " " << (k1 + 1.0) * 0.1 << " " << (Ur + 5.0) << endl;
		double ab = Lin_Interpolate((k1) * 0.1, a, (k1 + 1.0) * 0.1, b, (Ur + 5.0));
		double cd = Lin_Interpolate((k1) * 0.1, c, (k1 + 1.0) * 0.1, d, (Ur + 5.0));
		//cout << ab << " " << cd << endl;
		return Lin_Interpolate((k2 + 1.0) * 0.1, ab, (k2) * 0.1, cd, uu);
	}
	else if (Ur >= -5.0)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = 49;
		return Lin_Interpolate((k2) * 0.1, this->Int_002[k1][k2], (k2 - 3.0) * 0.1, this->Int_002[k1][k2 - 3], uu);
	}
	else if (uu <= 4.9)
	{
		int k1 = 0;
		int k2 = (int)(uu / 0.1);
		return Lin_Interpolate((k1) * 0.0, this->Int_002[k1][k2], (k1 + 1.0) * 0.1, this->Int_002[k1 + 1][k2], Ur + 5.0);
	}
	else
	{
		cout << "ERROR  88  efrfecer" << endl;
		exit(-1);
	}
}

double MKmethod::Get_Int00625(const double& Ur, const double& uu)
{
	if (Ur >= -5.0 && uu <= 4.9)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = (int)(uu / 0.1);
		//cout << "k1 , k2 = " << k1 << " " << k2 << endl;
		double a = this->Int_00625[k1][k2 + 1];
		double b = this->Int_00625[k1 + 1][k2 + 1];
		double c = this->Int_00625[k1][k2];
		double d = this->Int_00625[k1 + 1][k2];
		///cout << a << " " << b << " " << c << " " << d << endl;
		double ab = Lin_Interpolate((k1) * 0.1, a, (k1 + 1.0) * 0.1, b, (Ur + 5.0));
		double cd = Lin_Interpolate((k1) * 0.1, c, (k1 + 1.0) * 0.1, d, (Ur + 5.0));
		//cout << ab << " " << cd << endl;
		return Lin_Interpolate((k2 + 1.0) * 0.1, ab, (k2) * 0.1, cd, uu);
	}
	else if (Ur >= -5.0)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = 49;
		return Lin_Interpolate((k2) * 0.1, this->Int_00625[k1][k2], (k2 - 3.0) * 0.1, this->Int_00625[k1][k2 - 3], uu);
	}
	else if (uu <= 4.9)
	{
		int k1 = 0;
		int k2 = (int)(uu / 0.1);
		return Lin_Interpolate((k1) * 0.0, this->Int_00625[k1][k2], (k1 + 1.0) * 0.1, this->Int_00625[k1 + 1][k2], Ur + 5.0);
	}
	else
	{
		cout << "ERROR  157  efrfecer" << endl;
		exit(-1);
	}
}

double MKmethod::Get_Int02(const double& Ur, const double& uu)
{
	if (Ur >= -5.0 && uu <= 4.9)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = (int)(uu / 0.1);
		//cout << "k1 , k2 = " << k1 << " " << k2 << endl;
		double a = this->Int_02[k1][k2 + 1];
		double b = this->Int_02[k1 + 1][k2 + 1];
		double c = this->Int_02[k1][k2];
		double d = this->Int_02[k1 + 1][k2];
		//cout << a << " " << b << " " << c << " " << d << endl;
		double ab = Lin_Interpolate((k1) * 0.1, a, (k1 + 1.0) * 0.1, b, (Ur + 5.0));
		double cd = Lin_Interpolate((k1) * 0.1, c, (k1 + 1.0) * 0.1, d, (Ur + 5.0));
		//cout << ab << " " << cd << endl;
		return Lin_Interpolate((k2 + 1.0) * 0.1, ab, (k2) * 0.1, cd, uu);
	}
	else if (Ur >= -5.0)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = 49;
		return Lin_Interpolate((k2) * 0.1, this->Int_02[k1][k2], (k2 - 3.0) * 0.1, this->Int_02[k1][k2 - 3], uu);
	}
	else if (uu <= 4.9)
	{
		int k1 = 0;
		int k2 = (int)(uu / 0.1);
		return Lin_Interpolate((k1) * 0.0, this->Int_02[k1][k2], (k1 + 1.0) * 0.1, this->Int_02[k1 + 1][k2], Ur + 5.0);
	}
	else
	{
		cout << "ERROR  88  efrfecer" << endl;
		exit(-1);
	}
}

double MKmethod::Get_Int055(const double& Ur, const double& uu)
{
	if (Ur >= -5.0 && uu <= 4.9)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = (int)(uu / 0.1);
		//cout << "k1 , k2 = " << k1 << " " << k2 << endl;
		double a = this->Int_055[k1][k2 + 1];
		double b = this->Int_055[k1 + 1][k2 + 1];
		double c = this->Int_055[k1][k2];
		double d = this->Int_055[k1 + 1][k2];
		//cout << a << " " << b << " " << c << " " << d << endl;
		double ab = Lin_Interpolate((k1) * 0.1, a, (k1 + 1.0) * 0.1, b, (Ur + 5.0));
		double cd = Lin_Interpolate((k1) * 0.1, c, (k1 + 1.0) * 0.1, d, (Ur + 5.0));
		//cout << ab << " " << cd << endl;
		return Lin_Interpolate((k2 + 1.0) * 0.1, ab, (k2) * 0.1, cd, uu);
	}
	else if (Ur >= -5.0)
	{
		int k1 = (int)((Ur + 5.0) / 0.1);
		int k2 = 49;
		return Lin_Interpolate((k2) * 0.1, this->Int_055[k1][k2], (k2 - 3.0) * 0.1, this->Int_055[k1][k2 - 3], uu);
	}
	else if (uu <= 4.9)
	{
		int k1 = 0;
		int k2 = (int)(uu / 0.1);
		return Lin_Interpolate((k1) * 0.0, this->Int_055[k1][k2], (k1 + 1.0) * 0.1, this->Int_055[k1 + 1][k2], Ur + 5.0);
	}
	else
	{
		cout << "ERROR  88  efrfecer" << endl;
		exit(-1);
	}
}

double MKmethod::Lin_Interpolate(const double& x1, const double& y1, const double& x2, const double& y2, const double& x)
{
	double k = (y1 - y2) / (x1 - x2);
	double b = (y1 * x2 - y2 * x1) / (x2 - x1);
	return (k * x + b);
}

double MKmethod::A0(const double& Y)
{
	return (Y + 1.0 / (2.0 * Y)) * erf(Y) + exp(-kv(Y)) / sqrtpi_ + Y;
}

double MKmethod::G(const double& gam, const double& Y)
{
	double b = sqrt(1.0 + gam);
	return (exp(-kv(Y)) / (2.0 * sqrtpi_ * b * kv(b) * Y)) * //
		(-sqrtpi_ * (kv(b) * (gam + 4.0) + (3.0 * gam + 2.0) * kv(Y)) - 2.0 * b * gam * Y - 4.0 * b * Y + //
			sqrtpi_ * exp(kv(Y) / kv(b)) * (gam * (gam + 2.0 * kv(Y) + 5.0) + 4.0) * (erf(Y / b) + 1.0));
}

double MKmethod::F(const double& X, const double& gam, const double& Y)
{
	double b = sqrt(1.0 + gam);
	double gg = 1.0 + gam;
	double Y2 = kv(Y);
	double Y4 = pow4(Y);
	double X2 = kv(X);
	double A = exp(-Y2) / (2.0 * sqrtpi_ * Y * pow7(b));
	double B = 2 * X * Y * pow3(b) * (2.0 + gam * (3.0 + (X2 - 1.0) * Y2 + gam));
	double C = sqrtpi_ * gg * (gg * (gg * (4.0 + gam) + Y2 * (2.0 + 3.0 * gam)) + //
		exp(X2 * Y2 / gg) * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - kv(gg) * (4.0 + gam) + //
			Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))));
	double D = sqrtpi_ * gg * exp(X2 * Y2 / gg) * //
		(-2.0 * X2 * (X2 - 1.0) * Y4 * gam + kv(gg) * (4.0 + gam) - Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))) * //
		erf(X * Y / b);
	//cout << A << endl;
	//cout << B << endl;
	//cout << C << endl;
	//cout << D << endl;
	return A * (B + C - D);
}

double MKmethod::FF(const double& gam, const double& Yr)
// Из препринта Маламы
{
	return gam * ((0.5 + kv(Yr)) * (1.0 + erf(-Yr)) - exp(-kv(Yr)) * Yr / sqrtpi_);
}

double MKmethod::F0(const double& X, const double& Y)
{

	return -(2.0 * X * Y * exp(-kv(X * Y)) * (5.0 - 2.0 * (-2.0 + kv(X)) * kv(Y)) + //
		sqrtpi_ * (-4.0 * pow4(X * Y) + 8.0 * kv(X * Y) * (1.0 + kv(Y)) + (3.0 + 4.0 * kv(Y) - 4.0 * pow4(X * Y) //
			+ 8.0 * kv(X * Y) * (1.0 + kv(Y))) * erf(X * Y)))//
		/ (8.0 * sqrtpi_ * Y);
}

double MKmethod::FI(const double& Z, const double& X, const double& gam, const double& Y)
{
	double b = sqrt(1.0 + gam);
	double gg = 1.0 + gam;
	double X2 = kv(X);
	double Y2 = kv(Y);
	double Y4 = pow4(Y);
	double Z2 = kv(Z);
	double A = exp(-Y2 - 2.0 * X * Y * Z - Z2 * gg) / (sqrtpi_ * pow7(b));
	double ee = exp(kv(X * Y + Z + Z * gam) / gg);
	double B = ee * sqrtpi_ * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * kv(gg) + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam));
	double C = 2.0 * b * (X2 * (X2 - 1.0) * Y4 * gam - X * (X2 - 1.0) * pow3(Y) * Z * gam * gg - kv(gg) + //
		(X2 - 1.0) * Y2 * gg * (1 + gam * (2.0 + Z2 * gg)));
	double D = ee * sqrtpi_ * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * kv(gg) + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam)) * //
		erf((X * Y + Z + Z * gam) / b);
	return A * (B + C + D);
}

void MKmethod::TEST(void)
{
	cout << "test F:" << endl;
	cout << this->F(0.9, 0.1, 2.54) << "  -   " << 0.260062 << endl;
	cout << this->F(0.5, 0.01, 2.3) << "  -   " << 0.0327696 << endl;
	cout << this->F(0.6, 0.01, 0.3) << "  -   " << 0.346545 << endl << endl;

	cout << "test FI:" << endl;
	cout << this->FI(-1.1, 0.5, 0.01, 2.54) << "  -   " << 10.5833 << endl;
	cout << this->FI(-2.1, 0.9, 0.1, 2.0) << "  -   " << 1.74939 << endl << endl;

	cout << "test R:" << endl;
	cout << this->R(0.9, 2.2) << "  -   " << 0.819301  << endl;
	cout << this->R(0.6, 2.5) << "  -   " << 0.384307  << endl << endl;

	cout << "test f2:" << endl;
	cout << this->f2(0.0, 0.12, -1.5, 0.7) << "  -   " << 1.17911 << endl;
	cout << this->f2(-2.0, 0.12, -1.5, 0.7) << "  -   " << 0.218135 << endl;
	cout << this->f2(-1.1, 0.02, -1.2, 0.1) << "  -   " << 0.5287 << endl;
	cout << this->f2(0.0, 0.02, -1.2, 0.1) << "  -   " << 0.927412 << endl;
	cout << this->f2(-1.2, 0.12, -1.5, 0.7) << "  -   " << 0.730397 << endl << endl;


	cout << "test h_mho:" << endl;
	cout << this->h_mho(0.2, 2.1) << "  -   " << 1.05157 << endl;
	cout << this->h_mho(2.2, 1.7) << "  -   " << 0.705761 << endl;
	cout << this->h_mho(1.2, 1.7) << "  -   " << 0.437421 << endl << endl;


}

double MKmethod::R(const double& X, const double& Y)
{
	return ((X * X * Y + 0.5 / (Y)) * erf(X * Y) + X * (exp(-X * X * Y * Y) / sqrtpi_ + X * Y)) / this->A0(Y);
}

double MKmethod::f2(const double& V, const double& gam, const double& ur, const double& ut)
{
	double b = sqrt(1.0 + gam);
	double b4 = pow4(b);
	double b2 = kv(b);

	double A = (exp(2.0 * V * ur - kv(ur) - kv(V) * b2) / (2.0 * sqrtpi_ * b4));
	if (A <= 0.0000000000001)
	{
		return 0.0;
	}
	return A * (-gam * kv(ut) * (ur + gam * V + V) + //
		(sqrtpi_ / (2.0 * b)) * exp(kv(V + gam * V - ur) / b2) * (2.0 * b4 + kv(ut) * (b2 * (3.0 * gam + 2.0) + 2.0 * gam * kv(ur))) *//
		(1.0 + erf((-ur + gam * V + V) / b))); 
}

double MKmethod::f2k(const double& V, const double& gam, const double& ur, const double& ut)
{
	double b = sqrt(1.0 + gam);
	double b4 = pow4(b);
	double b2 = kv(b);
	double V2 = kv(V);
	double g2 = kv(gam);
	double gg = (1.0 + gam);

	double ee = exp(kv(-ur + V + V * gam) / (1.0 + gam));
	double A = exp(-kv(ur - V) - V2 * gam) / (32.0 * pow(b, 7));
	double B = -4.0 * ee * sqrtpi_ * pow4(ur * ut) * b * g2 + 2.0 * pow3(ur * ut) * ut * g2 * gg;
	double C = V * gam * (1.0 + gam) - 2.0 * ee * sqrtpi_ * b * (2.0 + 5.0 * gam);
	double D = V * gam * (1.0 + gam) - ee * sqrtpi_ * b * (2.0 + 3.0 * gam);
	double E = 2.0 * pow3(V) * g2 * kv(1.0 + gam) + V * gam * b2 * (4.0 + 7.0 * gam) - //
		ee * sqrtpi_ * b * (8.0 + 5.0 * gam * (4.0 + 3.0 * gam));
	double G = 8.0 * (1.0 + gam) + kv(ut) * (4.0 + gam * (9.0 + 2.0 * V2 * (1.0 + gam)));
	double H = 8.0 + gam * (36.0 + 4.0 * gam * pow4(ur) + 4.0 * kv(ur) * (1.0 + gam) * (2.0 + 5.0 * gam) + //
		gam * (63.0 + 5.0 * gam * (10.0 + 3.0 * gam)));
	double K = -16.0 * ee * sqrtpi_ * pow5(b) + 8.0 * kv(ut) * (1.0 + gam) * D + pow4(ut) * E;
	double L = -8.0 * ee * sqrtpi_ * pow3(b) + kv(ut) * C;
	double M = B + 2.0 * kv(ur * ut) * gam * (1.0 + gam) * L + ur * ut * ut * gam * kv(1.0 + gam) * G + kv(1.0 + gam) * K;
	double N = 16.0 * pow4(1.0 + gam) + 8.0 * kv(ut * (1.0 + gam)) * (2.0 + gam * (5.0 + 2.0 * kv(ur) + 3.0 * gam)) + pow4(ut) * H;

	//cout << G << endl;

	return A * (1.0 / (sqrtpi_ * pow3(b)) * 2.0 * M + 1.0 / (1.0 + gam) * ee * N * erfc((-ur + V + V * gam) / b));
}

double MKmethod::h_mho(const double& x, const double& c)
{
	if (x <= pi_ / 2.0)
	{
		return exp(2.0 * c * cos(x))/(exp(2.0 * c) + (2.0 * x * (1.0 - exp(2.0 * c))) / (pi_));
	}
	else
	{
		return exp(2.0 * c * cos(x)) / (pow(pi_ / (2.0 * x), 2.0 * c / (log(2.0))));
	}
}

double MKmethod::play_mho(Sensor* sens, const double& c)
{
	//cout << c << endl;
	double x;
	do
	{
		x = sens->MakeRandom() * pi_;
	} while (exp(2.0 * c * cos(x)) < sens->MakeRandom() * exp(2.0 * fabs(c)) );

	if (sens->MakeRandom() < 0.5)
	{
		return x;
	}
	else
	{
		return (2.0 * pi_ - x);
	}
}

double MKmethod::play_mho2(Sensor* sens, const double& c, const double& p, const double& p1)
{
	//cout << c << endl;
	double cc = c;
	bool minus = false;
	if (c < 0.0)
	{
		minus = true;
		cc = -c;
	}

	double x = 0.0;
	double ksi, k1, k2, k0;

	if (sens->MakeRandom() <= p1 / p)
	{
		do
		{
			ksi = sens->MakeRandom();
			k1 = 0.0;
			k2 = pi_ / 2.0;
			k0 = (k1 + k2) / 2.0;
			while (fabs(k2 - k1) > 0.0001)     // Деление пополам, иначе разваливается
			{
				k0 = (k1 + k2) / 2.0;
				if (this->mho_H1(c, k0, ksi) < 0)
				{
					k1 = k0;
				}
				else
				{
					k2 = k0;
				}
			}
			x = k0;
		} while (sens->MakeRandom() > exp(2.0 * c * cos(x)) / exp(-8.0 * c * kv(x / pi_) + 2.0 * c));
	}
	else
	{
		do
		{
			ksi = sens->MakeRandom();
			k1 = pi_ / 2.0;
			k2 = pi_;
			k0 = (k1 + k2) / 2.0;
			while (fabs(k2 - k1) > 0.0001)     // Деление пополам, иначе разваливается
			{
				k0 = (k1 + k2) / 2.0;
				if (this->mho_H2(c, k0, ksi) < 0)
				{
					k1 = k0;
				}
				else
				{
					k2 = k0;
				}
			}
			x = k0;
		} while (sens->MakeRandom() > exp(2.0 * c * cos(x)) / (exp(-4.0 * c / pi_ * (x - pi_ / 2.0))));
	}
	

	if (minus)
	{
		x = pi_ - x;
	}

	//return x;

	if (sens->MakeRandom() < 0.5)
	{
		return x;
	}
	else
	{
		return (2.0 * pi_ - x);
	}
}

double MKmethod::mho_H1(const double& c, const double& k, const double& ksi)
{
	return pi_ * sqrt(pi_ / (2.0 * c)) * exp(2.0 * c) / (4.0) * (erf(2.0 * sqrt(2.0 * c) * k / pi_) - ksi * erf(sqrt(2.0 * c)));
}

double MKmethod::mho_H2(const double& c, const double& k, const double& ksi)
{
	return -pi_ / (4.0 * c) * ((exp(c * (2.0 - 4.0 * k / pi_)) - 1.0) - ksi * (exp(-2.0 * c) - 1.0));
}


double MKmethod::norm_mho(const double& c)
{
	double s = 1.0;
	double d = 0.0;
	double cc = 1.0;
	double nn = 1.0;
	for (int i = 1; i < 100; i++)
	{
		cc = cc * c;
		nn = nn * i;
		d = kv(cc / nn);
		s = s + d;
		if (d < 0.00001)
		{
			break;
		}
	}

	return s;
}

double MKmethod::norm_mho2(const double& c)
{
	double s = 4.0 * c;
	double d = 4.0 * c;
	double cc = 1.0;
	double nn = 1.0;
	double fo = 1.0;
	for (int i = 0; i < 100; i++)
	{
		d = d * 4.0 * c * c/ kv(2.0 * i + 3.0);
		s = s + d;
		if (d < 0.00001)
		{
			break;
		}
	}

	return s;
}

bool MKmethod::Init_Parametrs(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_)
{
	double Y = fabs(Velosity_inf);
	double ksi, ksi1, ksi2;
	double X1;
	double X2;
	double X0 = 1.0;                // для метода хорд
	double split, gam1, gam2;
	vector <double> Wa_(this->num_area);
	vector <double> Mho_(this->num_area);

	// Разыгрываем  X
	for (int i = 0; i < this->num_area; i++)
	{
		ksi = sens->MakeRandom();
		//cout << "ksi  " << ksi << endl;
		X1 = 0.0;
		X2 = 1.0;
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = this->gam_[i];
		}
		else
		{
			gam1 = this->gam_[i - 1];
			gam2 = this->gam_[i];
		}

		while (fabs(X2 - X1) > 0.0000001)     // Деление пополам, иначе разваливается
		{
			X0 = (X1 + X2) / 2.0;
			if (this->Hx(gam1, gam2, X0, Y, ksi) < 0.0)
			{
				X1 = X0;
			}
			else
			{
				X2 = X0;
			}
			//k++;
		}
		/*if (std::fpclassify(X0) != FP_NORMAL && std::fpclassify(X0) != FP_ZERO)
		{
			cout << X0 << "    ERROR  229" << endl;
			exit(-1);
		}*/
		X_[i] = X0;
		//cout << X0 << " " << ksi << " " << gam1 << " " << gam2 << endl;
		//cout << "X2 = " << X2 << "    " << k << endl;
	}

	double Wr1, Wr2, Wr0;
	// Разыгрываем  Wr
	for (int i = 0; i < this->num_area; i++)
	{
		ksi = sens->MakeRandom();
		//cout << "ksi  " << ksi << endl;
		Wr1 = -5.0;
		Wr2 = 0.0;
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = this->gam_[i];
		}
		else
		{
			gam1 = this->gam_[i - 1];
			gam2 = this->gam_[i];
		}
		while (this->Hwr(gam1, gam2, Wr1, X_[i], Y, ksi) >= 0.0)
		{
			Wr1 = Wr1 - 1.0;
		}
		//int k = 0;
		while (fabs(Wr2 - Wr1) > 0.000001)     // Деление пополам, иначе разваливается
		{
			Wr0 = (Wr1 + Wr2) / 2.0;
			if (this->Hwr(gam1, gam2, Wr0, X_[i], Y, ksi) < 0)
			{
				Wr1 = Wr0;
			}
			else
			{
				Wr2 = Wr0;
			}
			//k++;
		}

		/*if (std::fpclassify(Wr2) != FP_NORMAL && std::fpclassify(Wr2) != FP_ZERO)
		{
			cout << Wr2 << "    ERROR  Wr2  276" << endl;
			exit(-1);
		}*/

		Wr_[i] = Wr0;
		//cout << Wr0 << " " << ksi << " " << X_[i] << endl;
		//cout << "Wr2 = " << Wr2 << "    " << k << endl;
	}

	double W1, W2, Wa, Yt;
	// Разыгрываем  Wa
	for (int i = 0; i < this->num_area; i++)
	{
		Yt = Y * sqrt(1.0 - kv(X_[i]));
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = this->gam_[i];
		}
		else
		{
			gam1 = this->gam_[i - 1];
			gam2 = this->gam_[i];
		}

		do
		{
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			W1 = sqrt(gam1 * kv(Wr_[i]));
			W2 = sqrt(gam2 * kv(Wr_[i]));
			Wa = sqrt(-log(exp(-kv(W1)) - ksi1 * (exp(-kv(W1)) - exp(-kv(W2)))));
		} while ((1.0 + kv(Yt * Wa)) / (1.0 + kv(Yt * W2)) < ksi2);

		if (std::fpclassify(Wa) != FP_NORMAL && std::fpclassify(Wa) != FP_ZERO)
		{
			cout << Wa << "    ERROR  Wa  311" << endl;
			cout << ksi1 << endl;
			exit(-1);
		}

		Wa_[i] = Wa;
		//cout << Wa << " " << ksi1 << " " << Wr_[i] << " " << X_[i] << endl;
	}

	//exit(-1);

	// Разыгрываем  Mho
	//for (int i = 0; i < this->num_area; i++)
	//{
	//	Mho_[i] = play_mho(sens, Y * sqrt(1.0 - kv(X_[i])) * Wa_[i]);
	//	/*if (std::fpclassify(Mho_[i]) != FP_NORMAL && std::fpclassify(Mho_[i]) != FP_ZERO)
	//	{
	//		cout << Mho_[i] << "    ERROR  Mho_[i]  324   " << Y * sqrt(1.0 - kv(X_[i])) * Wa_[i] << endl;
	//		exit(-1);
	//	}*/
	//}

	// Считаем веса
	for (int i = 0; i < this->num_area; i++)
	{
		double c = Y * sqrt(1.0 - kv(X_[i])) * Wa_[i];
		double p = this->norm_mho(c);
		if (false)//(fabs(c) > 0.1)
		{
			double p1 = (pi_ / 2.0) * (p + this->norm_mho2(c) / pi_);
			Mho_[i] = play_mho2(sens, c, p * pi_, p1);
		}
		else
		{
			Mho_[i] = play_mho(sens, c);
		}

		Wt_[i] = Wa_[i] * cos(Mho_[i]);
		Wp_[i] = Wa_[i] * sin(Mho_[i]);
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = this->gam_[i];
		}
		else
		{
			gam1 = this->gam_[i - 1];
			gam2 = this->gam_[i];
		}
		mu_[i] = (this->F(1.0, gam2, Y) - this->F(1.0, gam1, Y)) * (p) / (this->A0_ * (1.0 + kv(c)));
		//cout << mu_[i] << "     -  " << i <<  endl;
		//cout << (this->F(1.0, gam2, Y) - this->F(1.0, gam1, Y)) << "  " << this->A0_ << endl;
		/*if (std::fpclassify(mu_[i]) != FP_NORMAL && std::fpclassify(mu_[i]) != FP_ZERO)
		{
			cout << mu_[i] << "    ERROR  mu_[i]  348" << endl;
			cout << gam1 << " " << gam2 << endl;
			cout << this->F(1.0, gam2, Y) << endl;
			cout << this->F(1.0, gam1, Y) << endl;
			cout << c << endl;
			cout << this->A0_ << endl;
			cout << Wa_[i] << endl;
			exit(-1);
		}*/
	}
	//exit(-1);

	// Разыгрываем основной атом

	// Розыгрыш X целиком из препринта
	double p1 = erf(Y) / (this->A1_ * kv(Y));
	double ksi3, ksi4;
	double z, h;
	ksi1 = sens->MakeRandom();
	if (p1 < ksi1)
	{
		do
		{
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			X2 = sqrt(ksi2);
			h = (1.0 + erf(X2 * Y)) / (1.0 + erf(Y));
		} while (h < ksi3);
	}
	else
	{
		do
		{
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			X2 = (1.0 / Y) * sqrt(-log(ksi2)) * cos(pi_ * ksi3 / 2.0);
		} while (X2 >= 1.0);
	}

	X_[this->num_area] = X2;

	double gg = 0.0;
	if (this->num_area > 0)
	{
		gg = this->gam_[this->num_area - 1];
	}
	
	double p4_ = sqrtpi_ * (Y * X2) / (1.0 + sqrtpi_ * (Y * X2));                  // Для розыгрыша основного атома на границе по препринту Маламы
	//do
	//{
		do
		{
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			if (p4_ > ksi1)
			{
				z = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			}
			else
			{
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}

			Wr_[this->num_area] = z - (Y * X2);
			h = fabs(-(Y * X2) + z) / ((Y * X2) + fabs(z));
			ksi4 = sens->MakeRandom();
		} while (h <= ksi4 || z >= (Y * X2));

		double ksi5, ksi6;
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		//ksi7 = sens->MakeRandom();

		Wt_[this->num_area] = Y * sqrt(1.0 - kv(X_[this->num_area])) + sqrt(-log(ksi5)) * cos(2.0 * pi_ * ksi6);
		Wp_[this->num_area] = sqrt(-log(ksi5)) * sin(2.0 * pi_ * ksi6);
	//} while(Wr_[this->num_area] < 0.0 && kv(Wt_[this->num_area]) + kv(Wp_[this->num_area]) < gg * kv(Wr_[this->num_area]));

	/*mu_[this->num_area] = 1.0;
	for (int i = 0; i < this->num_area; i++)
	{
		mu_[this->num_area] = mu_[this->num_area] - mu_[i];
	}*/

	/*for (int i = 0; i <= this->num_area; i++)
	{
		cout << mu_[i] << endl;
	}
	cout << "----------------------------------" << endl;*/

	if (this->num_area == 0)
	{
		mu_[this->num_area] = 1.0;
		//cout << "A " << endl;
		return true;
	}
	else
	{
		if (Wr_[this->num_area] >= 0.0 || kv(Wt_[this->num_area]) + kv(Wp_[this->num_area]) > this->gam_[this->num_area - 1] * kv(Wr_[this->num_area]))
		{
			mu_[this->num_area] = 1.0;
			//cout << "B" << endl;
		}
		else
		{
			mu_[this->num_area] = 0.0;  // Чтобы не запускать этот атом
			//cout << "C" << endl;
			return false;
		}
	}

	return true;
}

bool MKmethod::Init_Parametrs2(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_)
// Алгоритм Маламы из препринта
{
	double Y = fabs(Velosity_inf);

	double p1 = erf(Y) / (this->A1_ * kv(Y));
	double ksi1, ksi2, ksi3;
	double X, h;

	// Разыграли X
	ksi1 = sens->MakeRandom();
	if (p1 >= ksi1)
	{
		do
		{
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			X = (1.0 / Y) * sqrt(-log(ksi2)) * cos(ksi3 * pi_ / 2.0);
		} while (X >= 1.0);
	}
	else
	{
		do
		{
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			X = sqrt(ksi2);
			h = (1.0 + erf(X * Y)) / (1.0 + erf(Y));
		} while (h <= ksi3);
	}

	for (int i = 0; i <= this->num_area; i++)
	{
		X_[i] = X;
	}

	double Yr = -Y * X;
	double Yt = Y * sqrt(1.0 - kv(X));
	double e1, e2, e3;

	this->A2_ = exp(-kv(Y * X)) / sqrtpi_ + Y * X * (1.0 + erf(X * Y));

	e1 = sqrtpi_ * kv(Yr);
	e2 = 2.0 * fabs(Yr);
	e3 = 0.5 * sqrtpi_;

	double p2, p4, p5;

	p1 = e1 / (e1 + e2 + e3);
	p2 = e2 / (e1 + e2 + e3);

	double ksi4, ksi5, ksi6, ksi7, ksi8;
	double z = 0.0;
	double Wa = 0.0;
	double V = 0.0;
	double gam1, gam2;

	for (int i = 0; i < this->num_area; i++)
	{
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = this->gam_[i];
		}
		else
		{
			gam1 = this->gam_[i - 1];
			gam2 = this->gam_[i];
		}

		// Разыграли Wr
		do
		{
			ksi1 = sens->MakeRandom();
			if (p1 > ksi1)
			{
				ksi2 = sens->MakeRandom();
				ksi3 = sens->MakeRandom();
				z = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			}
			else if (p1 + p2 > ksi1)
			{
				ksi2 = sens->MakeRandom();
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}
			else
			{
				ksi2 = sens->MakeRandom();
				ksi3 = sens->MakeRandom();
				ksi4 = sens->MakeRandom();
				ksi5 = sens->MakeRandom();
				z = sign(ksi5 - 0.5) * sqrt(-log(ksi2) - log(ksi3) * kv(cos(pi_ * ksi4)));
			}

			Wr_[i] = z + Yr;
			h = kv(Yr + z) / kv(fabs(Yr) + fabs(z));
			ksi6 = sens->MakeRandom();
		} while (h <= ksi6 || z > -Yr);

		ksi7 = sens->MakeRandom();
		ksi8 = sens->MakeRandom();
		V = 2.0 * pi_ * ksi7;
		Wa = sqrt(gam1 * kv(Wr_[i]) + ksi8 * (gam2 - gam1) * kv(Wr_[i]));

		Wt_[i] = Wa * cos(V);
		Wp_[i] = Wa * sin(V);

		mu_[i] = ((this->FF(gam2, Yr) - this->FF(gam1, Yr)) / this->A2_) * fabs(Wr_[i]) * exp(-kv(Wt_[i] - Yt) - kv(Wp_[i]));
	}

	// Разыгрываем основной атом

	double gg = 0.0;
	if (this->num_area > 0)
	{
		gg = this->gam_[this->num_area - 1];
	}

	double p4_ = sqrtpi_ * fabs(Yr) / (1.0 + sqrtpi_ * fabs(Yr));                  // Для розыгрыша основного атома на границе по препринту Маламы
	do
	{
		do
		{
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			if (p4_ > ksi1)
			{
				z = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			}
			else
			{
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}

			Wr_[this->num_area] = z + Yr;
			h = fabs(Yr + z) / (fabs(Yr) + fabs(z));
			ksi4 = sens->MakeRandom();
		} while (h <= ksi4 || z >= -Yr);

		double ksi5, ksi6;
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		//ksi7 = sens->MakeRandom();

		Wt_[this->num_area] = Yt + sqrt(-log(ksi5)) * cos(2.0 * pi_ * ksi6);
		Wp_[this->num_area] = sqrt(-log(ksi5)) * sin(2.0 * pi_ * ksi6);
	} while (Wr_[this->num_area] < 0.0 && kv(Wt_[this->num_area]) + kv(Wp_[this->num_area]) < gg * kv(Wr_[this->num_area]));


	mu_[this->num_area] = 1.0;

	for (int i = 0; i < this->num_area; i++)
	{
		mu_[this->num_area] = mu_[this->num_area] - mu_[i];
		//cout << i << " " << mu_[i] << endl;
	}

	return true;
}

int MKmethod::Init(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, double& X_)
{
	double Y = fabs(Velosity_inf);
	double ksi, ksi1, ksi2;
	double X1;
	double X2;
	double X0 = 1.0;                // для метода хорд
	double split;

	// Разыгрываем основной атом
	ksi = sens->MakeRandom();
	X1 = 0.0;
	X2 = 1.0;
	//int k = 0;
	while (fabs(X2 - X1) > 0.00001)     // Алгоритм метода секущих (быстрее, но опаснее)
	{
		split = X2;
		X2 = X1 - (X2 - X1) * (this->R(X1, Y) - ksi) / ((this->R(X2, Y) - ksi) - (this->R(X1, Y) - ksi));
		if (X2 > 1.0)  // Важная проверка, иногда X2 выходит за пределы
		{
			X2 = 1.0;
		}
		X1 = split;
		//k++;
	}
	X_ = X2;

	int num = 1;
	if (X_ > 0.989)
	{
		//cout << "Drobim" << endl;
		num = 10;
	}
	else if (X_ > 0.98)
	{
		//cout << "Drobim" << endl;
		num = 5;
	}
	else if (X_ > 0.97)
	{
		//cout << "Drobim" << endl;
		num = 3;
	}
	mu_.resize(num);
	Wt_.resize(num);
	Wp_.resize(num);
	Wr_.resize(num);

	//cout << X2 << endl;
	double p4_ = sqrtpi_ * (Y * X2) / (1.0 + sqrtpi_ * (Y * X2));                  // Для розыгрыша основного атома на границе по препринту Маламы
	double ksi3, ksi4, z, h;
	double ksi5, ksi6;

	for (int i = 0; i < num; i++)
	{
		do
		{
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			if (p4_ > ksi1)
			{
				z = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			}
			else
			{
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}

			Wr_[i] = z - (Y * X2);
			h = fabs(-(Y * X2) + z) / ((Y * X2) + fabs(z));
			ksi4 = sens->MakeRandom();
		} while (h <= ksi4 || z >= (Y * X2));

		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();

		Wt_[i] = Y * sqrt(1.0 - kv(X_)) + sqrt(-log(ksi5)) * cos(2.0 * pi_ * ksi6);
		Wp_[i] = sqrt(-log(ksi5)) * sin(2.0 * pi_ * ksi6);

		mu_[i] = 1.0 / num;
	}


	return num;
}

double MKmethod::Hx(const double& gam1, const double& gam2, const double& X, const double& Y, const double& ksi)
{
	return F(X, gam2, Y) - F(X, gam1, Y) - ksi * (F(1.0, gam2, Y) - F(1.0, gam1, Y));
}

double MKmethod::Hwr(const double& gam1, const double& gam2, const double& Z, const double& X, const double& Y, const double& ksi)
{
	return FI(Z, X, gam2, Y) - FI(Z, X, gam1, Y) - ksi * (FI(0.0, X, gam2, Y) - FI(0.0, X, gam1, Y));
}

bool MKmethod::Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
	vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex,//
	const double& y_ex, const double& z_ex)
{
	double X = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	//double uu = exp(-kv(X)) / sqrtpi_ + (X + 1.0 / (2.0 * X)) * erf(X);

	//double IS = this->Int_cp_1(X);
	vector <double> gamma_(I);
	vector <double> Wa_(I);
	vector <double> Mho_(I);

	for (int i = 0; i < I; i++)
	{
		gamma_[i] = 1.0 / (kv(r / this->R_[i]) - 1.0);
		/*if (gamma_[i] < 0.0)
		{
			cout << "EROR  882   gam < 0" << endl;
		}*/
	}

	/*for (int i = 0; i < I; i++)
	{
		gamma_[i] = 1.0 / (kv(1.0 / this->alpha_[i]) - 1.0);
	}*/

	double ksi, gam1, gam2, Wr1, Wr2, Wr0 = -1.0;
	// Разыграем Wr
	for (int i = 0; i < I; i++)
	{
		ksi = sens->MakeRandom();
		Wr1 = -3.0;
		Wr2 = 0.0;
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		while (this->Hvr(gam1, gam2, Wr1, Ur, Uthe, ksi) >= 0.0)
		{
			Wr1 = Wr1 - 1.0;
		}
		int k = 0;
		while (fabs(Wr2 - Wr1) > 0.0001)     // Деление пополам, иначе разваливается
		{
			Wr0 = (Wr1 + Wr2) / 2.0;
			if (this->Hvr(gam1, gam2, Wr0, Ur, Uthe, ksi) < 0)
			{
				Wr1 = Wr0;
			}
			else
			{
				Wr2 = Wr0;
			}
			k++;
		}
		Wr_[i] = Wr1;
		/*if (Wr_[i] >= 0)
		{
			cout << "ERROR 541 hfgfh   Wr >= 0" << endl;
			exit(-1);
		}*/
		/*if (i == 1)
		{
			cout << Wr0 << endl;
		}*/
		/*if (std::fpclassify(Wr1) != FP_NORMAL && std::fpclassify(Wr1) != FP_ZERO)
		{
			cout << Wr1 << "    ERROR  Wr0 537" << endl;
			exit(-1);
		}*/
	}

	double W1, W2, Wa, ksi1, ksi2;
	// Разыгрываем  Wa
	for (int i = 0; i < I; i++)
	{
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		do
		{
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			W1 = sqrt(gam1 * kv(Wr_[i]));
			W2 = sqrt(gam2 * kv(Wr_[i]));
			Wa = sqrt(-log(exp(-kv(W1)) - ksi1 * (exp(-kv(W1)) - exp(-kv(W2)))));
		} while ((1.0 + kv(Uthe * Wa)) / (1.0 + kv(Uthe * W2)) < ksi2);
		Wa_[i] = Wa;
		/*if (std::fpclassify(Wa) != FP_NORMAL && std::fpclassify(Wa) != FP_ZERO)
		{
			cout << Wa << "    ERROR  Wa 568" << endl;
			exit(-1);
		}*/
	}

	// Разыгрываем  Mho
	//for (int i = 0; i < I; i++)
	//{
	//	Mho_[i] = play_mho(sens, Uthe * Wa_[i]);
	//	/*if (std::fpclassify(Mho_[i]) != FP_NORMAL && std::fpclassify(Mho_[i]) != FP_ZERO)
	//	{
	//		cout << Mho_[i] << "    ERROR  Mho_[i] 579" << endl;
	//		cout << "Uthe  " << Uthe << endl;
	//		cout << "Wa_[i]  " << Wa_[i] << endl;
	//		cout << Uthe * Wa_[i] << endl;
	//		exit(-1);
	//	}*/
	//}

	// Считаем веса
	for (int i = 0; i < I; i++) 
	{
		double c = Uthe * Wa_[i];
		double p = this->norm_mho(c);
		if (false)//(fabs(c) > 0.1)
		{
			double p1 = (pi_ / 2.0) * (p + this->norm_mho2(c) / pi_);
			Mho_[i] = play_mho2(sens, c, p * pi_, p1);
		}
		else
		{
			Mho_[i] = play_mho(sens, c);
		}

		Wthe_[i] = Wa_[i] * cos(Mho_[i]);
		Wphi_[i] = Wa_[i] * sin(Mho_[i]);
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		double u = sqrt(kvv(Vr - Wr_[i], Vthe - Wthe_[i], Vphi - Wphi_[i]));

		if (X > 7)
		{
			double uu = exp(-kv(X)) / sqrtpi_ + (X + 1.0 / (2.0 * X)) * erf(X);
			mu_[i] = (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (p) / (1.0 + kv(c));
		}
		else
		{
			mu_[i] = (u * sigma2(u, cp) / (this->int_1(X * cp, cp))) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (p) / (1.0 + kv(c));
		}
		//mu_[i] = (u * sigma2(u, cp) / (IS)) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c));
		//cout << IS << " " << uu * sigma2(uu, cp) << endl;
		
		//cout << "mu = " << mu_[i] << endl;
		//cout << u << " " << uu << " " << gam1 << " " << gam2 << " " << Uthe << " " << Ur << " " << Wa_[i] << "  " << cp << endl;
		//cout << f2(0.0, gam1, Ur, Uthe) -  f2(0.0, gam2, Ur, Uthe) << endl;
		//cout << (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) << endl;
		//exit(-4);

		/*if (std::fpclassify(mu_[i]) != FP_NORMAL && std::fpclassify(mu_[i]) != FP_ZERO)
		{
			cout << mu_[i] << "    ERROR  mu_[i] 604" << endl;
			cout << "f2(0.0, gam1, Ur, Uthe) = " << f2(0.0, gam1, Ur, Uthe) << endl;
			cout << "f2(0.0, gam2, Ur, Uthe) = " << f2(0.0, gam2, Ur, Uthe) << endl;
			cout << "gam1 = " << gam1 << endl;
			cout << "gam2 = " << gam2 << endl;
			cout << "Ur = " << Ur << endl;
			cout << "Uthe = " << Uthe << endl;

			exit(-1);
		}*/

		//cout << mu_[i] << "  " << i << endl;
	}
	

	// Розыгрыш основного атома
	double p4 = 0.5 * sqrtpi_ * X / (1.0 + 0.5 * sqrtpi_ * X);
	double om1, om2, om3, lo;
	double y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h;
	double ksi3, ksi4, ksi5, ksi6;
	double D, ko;

	double gg;
	if (I > 0)
	{
		gg = gamma_[I - 1];
	}
	else
	{
		gg = 0.0;
	}

	do
	{
		ksi1 = sens->MakeRandom();
		ksi2 = sens->MakeRandom();
		ksi3 = sens->MakeRandom();
		ksi4 = sens->MakeRandom();
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм   --  выйгрыша нет вроде от него
			/*do
			{
				om2 = 1.0 - 2.0 * sens->MakeRandom();
				om3 = 1.0 - 2.0 * sens->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1.0);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * pi_ * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * pi_ * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uuu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
		h = ((uuu * sigma2(uuu, cp)) / (sigma2(X, cp) * (X + yy)));
	} while (h < ksi6); // || (v1 < 0.0 && kv(v2) + kv(v3) < gg * kv(v1)));


	Wr_[I] = v1;
	Wthe_[I] = v2;
	Wphi_[I] = v3;

	/*mu_[I] = 1.0;
	for (int i = 0; i < I; i++)   
	{
		mu_[I] = mu_[I] - mu_[i];
	}*/

	//if (mu_[I] < 0.0)
	//{
	//	for (int i = 0; i < I; i++)
	//	{
	//		cout << "gamma_[i]  " << gamma_[i] << endl;
	//		cout << "mu_[i]  " << mu_[i] << endl;
	//	}
	//	
	//	for (int i = 0; i < I; i++)
	//	{

	//		double u = sqrt(kvv(Vr - Wr_[i], Vthe - Wthe_[i], Vphi - Wphi_[i]));
	//		cout << i << " " << (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) << endl;
	//	}
	//	exit(-1);
	//	//cout << "ERROR  710" << endl;
	//	//cout << r << endl;
	//	//cout << mu_[I] << endl;
	//	for (int i = 0; i <= I; i++)
	//	{
	//		double aa, bb, cc;
	//		dekard_skorost(y_ex, z_ex, x_ex, Wr_[i], Wphi_[i], Wthe_[i], bb, cc, aa);
	//		for (int ij = 0; ij < 27000; ij++)
	//		{
	//			cout << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << i << endl;
	//		}
	//	}
	//	exit(-1);
	//}

	if (Wr_[I] >= 0.0 || kv(Wthe_[I]) + kv(Wphi_[I]) > gg * kv(Wr_[I]))
	{
		mu_[I] = 1.0;
		return true;
		//cout << 1 << endl;
	}
	else
	{
		mu_[I] = 0.0;  // Чтобы не запускать этот атом
		//cout << 0 << endl;
		return false;
	}

	return true;
}

bool MKmethod::Change_Velosity2(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
	vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex,//
	const double& y_ex, const double& z_ex)
// Алгоритм Маламы из препринта
{
	vector <double> gamma_(I);

	for (int i = 0; i < I; i++)
	{
		gamma_[i] = 1.0 / (kv(r / this->R_[i]) - 1.0);
	}

	double e1 = sqrtpi_ * kv(Ur);
	double e2 = 2.0 * fabs(Ur);
	double e3 = 0.5 * sqrtpi_;
	double p1 = e1 / (e1 + e2 + e3);
	double p2 = e2 / (e1 + e2 + e3);
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7, ksi8;
	double z = 0.0, h = 0.0;

	for (int i = 0; i < I; i++)
	{
		do
		{
			ksi1 = sens->MakeRandom();
			if (p1 > ksi1)
			{
				ksi2 = sens->MakeRandom();
				ksi3 = sens->MakeRandom();
				z = sqrt(-log(ksi2)) * cos(ksi3 * pi_);
			}
			else if (p1 + p2 > ksi1)
			{
				ksi2 = sens->MakeRandom();
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}
			else
			{
				ksi2 = sens->MakeRandom();
				ksi3 = sens->MakeRandom();
				ksi4 = sens->MakeRandom();
				ksi5 = sens->MakeRandom();
				z = sign(ksi5 - 0.5) * sqrt(-log(ksi2) - log(ksi3) * //
					kv(cos(pi_ * ksi4)));
			}

			Wr_[i] = Ur + z;
			h = kv(Ur + z) / kv(fabs(Ur) + fabs(z));
			ksi6 = sens->MakeRandom();
		} while (h <= ksi6 || z >= -Ur);
	}

	//vector <double> w_(I + 1);
	vector <double> F_(I);

	for (int i = 0; i < I; i++)
	{
		F_[i] = 0.5 * gamma_[i] * ((0.5 + kv(Ur)) * (1.0 - sign(Ur) * erf(fabs(Ur))) - Ur * exp(-kv(Ur)) / sqrtpi_);
	}

	double Val;
	double The;
	double gam1, gam2;
	for (int i = 0; i < I; i++)
	{
		ksi7 = sens->MakeRandom();
		ksi8 = sens->MakeRandom();
		The = 2.0 * pi_ * ksi7;
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		Val = sqrt(gam1 * kv(Wr_[i]) + ksi8 * (gam2 * kv(Wr_[i]) - gam1 * kv(Wr_[i])));

		Wthe_[i] = Val * cos(The);
		Wphi_[i] = Val * sin(The);
	}

	double x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double uu = exp(-kv(x)) / sqrtpi_ + (x + 1.0 / (2.0 * x)) * erf(x);
	double u, F1, F2;

	// Считаем веса
	for (int i = 0; i < I; i++)
	{
		u = sqrt(kvv(Vr - Wr_[i], Vthe - Wthe_[i], Vphi - Wphi_[i]));
		if (i == 0)
		{
			F1 = 0.0;
			F2 = F_[i];
		}
		else
		{
			F1 = F_[i - 1];
			F2 = F_[i];
		}
		mu_[i] = (F2 - F1) * ( (u * sigma2(u, cp)) / (uu * sigma2(uu, cp))) * exp(-kv(Wthe_[i] - Uthe) - kv(Wphi_[i]));
	}


	// Расчёт основного атома
	//x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi_ * x / (1.0 + 0.5 * sqrtpi_ * x);
	double om1, om2, om3, lo;
	double y1, y2, y3;
	double v1, v2, v3, u1, u2, u3;
	double yy, vv, D, ko, uuu;

	double gg;
	if (I > 0)
	{
		gg = gamma_[I - 1];
	}
	else
	{
		gg = 0.0;
	}

	do
	{
		ksi1 = sens->MakeRandom();
		ksi2 = sens->MakeRandom();
		ksi3 = sens->MakeRandom();
		ksi4 = sens->MakeRandom();
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		//cout << "sd " << endl;
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм
			/*do
			{
				om2 = 1.0 - 2.0 * sens->MakeRandom();
				om3 = 1.0 - 2.0 * sens->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * pi_ * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * pi_ * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uuu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
		h = ((uuu * sigma2(uuu, cp)) / (sigma2(x, cp) * (x + yy)));
	} while (h < ksi6 || (v1 <= 0 && kvv(v2, v3, 0.0) < gg * kv(v1)));
	//cout << v2 << endl;
	Wr_[I] = v1;
	Wthe_[I] = v2;
	Wphi_[I] = v3;

	mu_[I] = 1.0;
	for (int i = 0; i < I; i++)
	{
		mu_[I] = mu_[I] - mu_[i];
	}

	//if (Wr_[I] >= 0.0 || kv(Wthe_[I]) + kv(Wphi_[I]) > gg * kv(Wr_[I]))
	//{
	//	mu_[I] = 1.0;
	//	return true;
	//	//cout << 1 << endl;
	//}
	//else
	//{
	//	mu_[I] = 0.0;  // Чтобы не запускать этот атом
	//	//cout << 0 << endl;
	//	return false;
	//}
	
	return true;
}

bool MKmethod::Change_Velosity3(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
	vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex,//
	const double& y_ex, const double& z_ex)
{
	double X = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double uu = exp(-kv(X)) / sqrtpi_ + (X + 1.0 / (2.0 * X)) * erf(X);

	//double IS = this->Int_cp_1(X);
	vector <double> gamma_(I);
	vector <double> Wa_(I);
	vector <double> Mho_(I);

	for (int i = 0; i < I; i++)
	{
		gamma_[i] = 1.0 / (kv(r / this->R_[i]) - 1.0);
		/*if (gamma_[i] < 0.0)
		{
			cout << "EROR  882   gam < 0" << endl;
		}*/
	}

	/*for (int i = 0; i < I; i++)
	{
		gamma_[i] = 1.0 / (kv(1.0 / this->alpha_[i]) - 1.0);
	}*/

	double ksi, gam1, gam2, Wr1, Wr2, Wr0 = -1.0;
	// Разыграем Wr
	for (int i = 0; i < I; i++)
	{
		ksi = sens->MakeRandom();
		//cout << "ksi  " << ksi << endl;
		Wr1 = -3.0;
		Wr2 = 0.0;
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		while (this->Hvrk(gam1, gam2, Wr1, Ur, Uthe, ksi) >= 0.0)
		{
			Wr1 = Wr1 - 1.0;
		}
		int k = 0;
		while (fabs(Wr2 - Wr1) > 0.00001)     // Деление пополам, иначе разваливается
		{
			Wr0 = (Wr1 + Wr2) / 2.0;
			if (this->Hvrk(gam1, gam2, Wr0, Ur, Uthe, ksi) < 0)
			{
				Wr1 = Wr0;
			}
			else
			{
				Wr2 = Wr0;
			}
			k++;
		}
		Wr_[i] = Wr1;
	}

	double W1, W2, Wa, ksi1, ksi2;
	// Разыгрываем  Wa
	for (int i = 0; i < I; i++)
	{
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		do
		{
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			W1 = sqrt(gam1 * kv(Wr_[i]));
			W2 = sqrt(gam2 * kv(Wr_[i]));
			Wa = sqrt(-log(exp(-kv(W1)) - ksi1 * (exp(-kv(W1)) - exp(-kv(W2)))));
		} while ((1.0 + kv(Uthe * Wa) + pow4(Uthe * Wa)/4.0) / (1.0 + kv(Uthe * W2) + pow4(Uthe * W2) / 4.0) < ksi2);
		Wa_[i] = Wa;
	}

	// Разыгрываем  Mho
	for (int i = 0; i < I; i++)
	{
		Mho_[i] = play_mho(sens, Uthe * Wa_[i]);
	}

	// Считаем веса
	for (int i = 0; i < I; i++)
	{
		Wthe_[i] = Wa_[i] * cos(Mho_[i]);
		Wphi_[i] = Wa_[i] * sin(Mho_[i]);
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}
		double c = Uthe * Wa_[i];
		double u = sqrt(kvv(Vr - Wr_[i], Vthe - Wthe_[i], Vphi - Wphi_[i]));

		mu_[i] = (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * (f2k(0.0, gam2, Ur, Uthe) - f2k(0.0, gam1, Ur, Uthe)) * exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c) + pow4(c)/4.0);
		
		//if (exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c)) > 1.0)
		//{
		//	cout << "A - " << exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c)) << " " << exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c) + pow4(c) / 4.0) << endl;
		//	cout << Uthe << endl;
		//	cout << Wa_[i] << endl;
		//	//exit(-1);
		//}

		//mu_[i] = (u * sigma2(u, cp) / (IS)) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c));
		//cout << IS << " " << uu * sigma2(uu, cp) << endl;
	}


	// Розыгрыш основного атома
	double p4 = 0.5 * sqrtpi_ * X / (1.0 + 0.5 * sqrtpi_ * X);
	double om1, om2, om3, lo;
	double y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h;
	double ksi3, ksi4, ksi5, ksi6;
	double D, ko;

	double gg;
	if (I > 0)
	{
		gg = gamma_[I - 1];
	}
	else
	{
		gg = 0.0;
	}

	do
	{
		ksi1 = sens->MakeRandom();
		ksi2 = sens->MakeRandom();
		ksi3 = sens->MakeRandom();
		ksi4 = sens->MakeRandom();
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			//om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			//om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм
			do
			{
				om2 = 1.0 - 2.0 * sens->MakeRandom();
				om3 = 1.0 - 2.0 * sens->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * pi_ * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * pi_ * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uuu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
		h = ((uuu * sigma2(uuu, cp)) / (sigma2(X, cp) * (X + yy)));
	} while (h < ksi6);// || (v1 < 0.0 && kv(v2) + kv(v3) < gg * kv(v1)));


	Wr_[I] = v1;
	Wthe_[I] = v2;
	Wphi_[I] = v3;


	if (Wr_[I] >= 0.0 || kv(Wthe_[I]) + kv(Wphi_[I]) > gg * kv(Wr_[I]))
	{
		mu_[I] = 1.0;
		return true;
		//cout << 1 << endl;
	}
	else
	{
		mu_[I] = 0.0;  // Чтобы не запускать этот атом
		//cout << 0 << endl;
		return false;
	}
	return true;
}

bool MKmethod::Change_Velosity4(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
	vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex,//
	const double& y_ex, const double& z_ex)
	// Как первая часть, но розыгрышь идёт по-частям
{
	double X = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	//double uu = exp(-kv(X)) / sqrtpi_ + (X + 1.0 / (2.0 * X)) * erf(X);

	//double IS = this->Int_cp_1(X);
	vector <double> gamma_(I);
	vector <double> Wa_(I);
	vector <double> Mho_(I);

	for (int i = 0; i < I; i++)
	{
		gamma_[i] = 1.0 / (kv(r / this->R_[i]) - 1.0);
	}

	double ksi, gam1, gam2, Wr1, Wr2, Wr0 = -1.0, ksi1, ksi2, W1, W2, Wa;
	int met = 0;
	// Разыграем Wr  
	for (int i = 0; i < I; i++)
	{
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}
		double pp1 = this->for_Wr_1(0.0, gam2, Ur) - this->for_Wr_1(0.0, gam1, Ur);
		double pp2 = this->for_Wr_2(0.0, gam2, Ur, Uthe) - this->for_Wr_2(0.0, gam1, Ur, Uthe);

		if (sens->MakeRandom() < pp1 / (pp1 + pp2))
		{
			ksi = sens->MakeRandom();
			Wr1 = -3.0;
			Wr2 = 0.0;

			while (this->H_Wr_1(gam1, gam2, Wr1, Ur, pp1, ksi) >= 0.0)
			{
				Wr1 = Wr1 - 1.0;
			}
			int k = 0;
			while (fabs(Wr2 - Wr1) > 0.0001)     // Деление пополам, иначе разваливается
			{
				Wr0 = (Wr1 + Wr2) / 2.0;
				if (this->H_Wr_1(gam1, gam2, Wr0, Ur, pp1, ksi) < 0)
				{
					Wr1 = Wr0;
				}
				else
				{
					Wr2 = Wr0;
				}
				k++;
			}
			Wr_[i] = Wr1;
		}
		else
		{
			ksi = sens->MakeRandom();
			Wr1 = -3.0;
			Wr2 = 0.0;

			while (this->H_Wr_2(gam1, gam2, Wr1, Ur, Uthe, pp2, ksi) >= 0.0)
			{
				Wr1 = Wr1 - 1.0;
			}
			int k = 0;
			while (fabs(Wr2 - Wr1) > 0.0001)     // Деление пополам, иначе разваливается
			{
				Wr0 = (Wr1 + Wr2) / 2.0;
				if (this->H_Wr_2(gam1, gam2, Wr0, Ur, Uthe, pp2, ksi) < 0)
				{
					Wr1 = Wr0;
				}
				else
				{
					Wr2 = Wr0;
				}
				k++;
			}
			Wr_[i] = Wr1;
		}
	}

	// Разыгрываем  Wa
	for (int i = 0; i < I; i++)
	{
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		do
		{
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			W1 = sqrt(gam1 * kv(Wr_[i]));
			W2 = sqrt(gam2 * kv(Wr_[i]));
			Wa = sqrt(-log(exp(-kv(W1)) - ksi1 * (exp(-kv(W1)) - exp(-kv(W2)))));
		} while ((1.0 + kv(Uthe * Wa)) / (1.0 + kv(Uthe * W2)) < ksi2);
		Wa_[i] = Wa;
	}

	// Считаем веса и Mho
	for (int i = 0; i < I; i++)
	{
		double c = Uthe * Wa_[i];
		double p = this->norm_mho(c);
		if (false)//(fabs(c) > 0.1)
		{
			double p1 = (pi_ / 2.0) * (p + this->norm_mho2(c) / pi_);
			Mho_[i] = play_mho2(sens, c, p * pi_, p1);
		}
		else
		{
			Mho_[i] = play_mho(sens, c);
		}

		Wthe_[i] = Wa_[i] * cos(Mho_[i]);
		Wphi_[i] = Wa_[i] * sin(Mho_[i]);
		if (i == 0)
		{
			gam1 = 0.0;
			gam2 = gamma_[i];
		}
		else
		{
			gam1 = gamma_[i - 1];
			gam2 = gamma_[i];
		}

		double u = sqrt(kvv(Vr - Wr_[i], Vthe - Wthe_[i], Vphi - Wphi_[i]));

		if (X > 7)
		{
			double uu = exp(-kv(X)) / sqrtpi_ + (X + 1.0 / (2.0 * X)) * erf(X);
			mu_[i] = (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (p) / (1.0 + kv(c));
		}
		else
		{
			mu_[i] = (u * sigma2(u, cp) / (this->int_1(X * cp, cp))) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (p) / (1.0 + kv(c));
		}
	}


	// Розыгрыш основного атома
	double p4 = 0.5 * sqrtpi_ * X / (1.0 + 0.5 * sqrtpi_ * X);
	double om1, om2, om3, lo;
	double y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h;
	double ksi3, ksi4, ksi5, ksi6;
	double D, ko;

	double gg;
	if (I > 0)
	{
		gg = gamma_[I - 1];
	}
	else
	{
		gg = 0.0;
	}

	do
	{
		ksi1 = sens->MakeRandom();
		ksi2 = sens->MakeRandom();
		ksi3 = sens->MakeRandom();
		ksi4 = sens->MakeRandom();
		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм   --  выйгрыша нет вроде от него
			/*do
			{
				om2 = 1.0 - 2.0 * sens->MakeRandom();
				om3 = 1.0 - 2.0 * sens->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1.0);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * pi_ * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * pi_ * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uuu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
		h = ((uuu * sigma2(uuu, cp)) / (sigma2(X, cp) * (X + yy)));
	} while (h < ksi6); // || (v1 < 0.0 && kv(v2) + kv(v3) < gg * kv(v1)));


	Wr_[I] = v1;
	Wthe_[I] = v2;
	Wphi_[I] = v3;

	if (Wr_[I] >= 0.0 || kv(Wthe_[I]) + kv(Wphi_[I]) > gg * kv(Wr_[I]))
	{
		mu_[I] = 1.0;
		return true;
		//cout << 1 << endl;
	}
	else
	{
		mu_[I] = 0.0;  // Чтобы не запускать этот атом
		//cout << 0 << endl;
		return false;
	}

	return true;
}


double MKmethod::for_Wr_1(const double& Z, const double& gam, const double& ur)
// Для розыгрыша Wr в перезарядке по-частям (первая часть)
{
	return -exp(-gam * kv(ur) / (gam + 1.0)) * (1.0 + erf((-ur + gam * Z + Z) / sqrt(gam + 1.0))) / (2.0 * sqrt(gam + 1.0));
}

double MKmethod::H_Wr_1(const double& gam1, const double& gam2, const double& V, const double& ur, const double& p, const double& ksi)
// Для розыгрыша Wr в перезарядке по-частям (первая часть)
{
	return for_Wr_1(V, gam2, ur) - for_Wr_1(V, gam1, ur) - ksi * p;
}

double MKmethod::for_Wr_2(const double& Z, const double& gam, const double& ur, const double& ut)
// Для розыгрыша Wr в перезарядке по-частям (первая часть)
{
	double gam1 = gam + 1.0;  // gam + 1
	double ur2 = kv(ur);
	return  1.0 / (4.0 * gam1 * gam1 * gam1 * sqrtpi_) * kv(ut) * exp((-1.0 - gam / gam1) * ur2 - gam1 * Z * Z) * //
		(2.0 * exp(ur * (gam * ur / gam1 + 2.0 * Z)) * gam * gam1 * (ur + Z + gam * Z) + //
			exp(ur2 + gam1 * Z * Z) * sqrt(gam1) * sqrtpi_ * (2.0 + gam * (5.0 + 3.0 * gam + 2.0 * ur2)) * //
			(-1.0 - erf((-ur + Z + gam * Z) / sqrt(gam1))));
}

double MKmethod::H_Wr_2(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& p, const double& ksi)
// Для розыгрыша Wr в перезарядке по-частям (первая часть)
{
	return for_Wr_2(V, gam2, ur, ut) - for_Wr_2(V, gam1, ur, ut) - ksi * p;
}

double MKmethod::Hvr(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi)
{
	return f2(V, gam1, ur, ut) - f2(V, gam2, ur, ut) - ksi * (f2(0.0, gam1, ur, ut) - f2(0.0, gam2, ur, ut));
}

double MKmethod::Hvrk(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi)
{
	return f2k(V, gam2, ur, ut) - f2k(V, gam1, ur, ut) - ksi * (f2k(0.0, gam2, ur, ut) - f2k(0.0, gam1, ur, ut));
}

double MKmethod::Int_cp_1(const double& x)
{
	if (x < 0.0)
	{
		cout << "Error  1789 fuihferfrefdewed   " << x << endl;
		exit(-1);
	}
	else if (x <= 1.0)
	{
		return 1.0818657860101877 + 0.285218386210128 * kv(x) - //
			0.0013036140991678873 * pow(x, 3) - 0.027948426972744425 * pow(x, 4) - //
			0.02488624181111631 * pow(x, 5) + 0.05903086910650146 * pow(x, 6) - //
			0.07375065435261574 * pow(x, 7) + 0.061880698880672345 * pow(x, 8) - //
			0.02992298377445923 * pow(x, 9) + 0.006218848236075065 * pow(x, 10);
	}
	else if (x <= 10.0)
	{
		return 1.163413866981175 - 0.2948788733781214 * x + // 
			0.7222333466034818 * kv(x) - 0.33744822137483277 * pow(x, 3) + //
			0.1011120159809275 * pow(x, 4) - 0.020694315456416357 * pow(x, 5) + //
			0.0029240572198959364 * pow(x, 6) - 0.00028083329445766233 * pow(x, 7) + //
			0.000017508447325531657 * pow(x, 8) - 638.9792015974883 * pow((0.1 * x), 9) + //
			103.58193963558761 * pow((0.1 * x), 10);
	}
	else
	{
		cout << "Error  1805 fuihferfrefdewed   " << x << endl;
		exit(-1);
	}
}

double MKmethod::Int_cp_2(const double& x)
{
	if (x < 0.0)
	{
		cout << "Error  1771 dwcewrc   " << x << endl;
		exit(-1);
	}
	else if (x <= 1.0)
	{
		return 1.08187 + 0.0000148503 * x + 0.284859 * kv(x) + 0.00139356 * pow(x, 3) - //
			0.0392891 * pow(x, 4) + 0.00585294 * pow(x, 5) + 0.00170631 * pow(x, 6);
	}
	else if (x <= 10.0)
	{
		return 1.1427 - 0.227066 * x + 0.631423 * kv(x) - 0.271463 * pow(x, 3) + 0.0720498 * pow(x, 4) - //
			0.0125149 * pow(x, 5) + 0.00142189 * pow(x, 6) - 0.000101751 * pow(x, 7) + //
			415.925 * pow(0.1 * x, 8) - 74.0093 * pow(0.1 * x, 9);
	}
	else
	{
		cout << "Error  1893 vertrwx   " << x << endl;
		exit(-1);
	}
}


double MKmethod::int_1_f1(const double& x)
{
	if (x <= 1.0)
	{
		return 6.283185155644284 + 0.000024846677279866114 * x +//
			2.0934329078277405 * x * x + 0.008055998193903208 * x * x * x -//
			0.2355169235647438 * x * x * x * x + 0.03820480582423355 * x * x * x * x * x +//
			0.006992274370591744 * x * x * x * x * x * x;
	}
	else if(x <= 3.0)
	{
		return 6.437524091973454 - 0.6331520099380095 * x +//
			3.1348881317268997 * x * x - 0.8454201478027856 * x * x * x +//
			0.1004702004260311 * x * x * x * x + 0.0009895488638964746 * x * x * x * x * x -//
			0.000920750276197054 * x * x * x * x * x * x;
	}
	else if (x <= 5)
	{
		return 4.4920780630505135 + 2.5133093267020654 * x +//
			1.1327223176567935 * x * x - 0.24648691152318875 * x * x * x +//
			0.031326738629523766 * x * x * x * x - 0.0021366031960331384 * x * x * x * x * x +//
			0.00005954097505746697 * x * x * x * x * x * x;
	}
	else if (x <= 7)
	{
		return 1.9138683588136232 + 5.350374732905213 * x -//
			0.16380205801427633 * x * x + 0.06765898334856263 * x * x * x -//
			0.011071118267864083 * x * x * x * x + 0.0008673476933852199 * x * x * x * x * x -//
			0.00002691859374483661 * x * x * x * x * x * x;
	}
	else if (x <= 50.0)
	{
		return 1.3138472469154294 + 5.336877156136497 * x +//
			0.020286308991329983 * x * x - 0.9780973533544137 * pow(x / 10.0, 3) +//
			0.26354051936651874 * pow(x / 10.0, 4) - 0.03711733070841712 * pow(x / 10.0, 5) +//
			0.002120935433043921 * pow(x / 10.0, 6);
	}
	else
	{
		cout << "Error  int_f1 > 7  =  " << x << endl;
		return 1.3138472469154294 + 5.336877156136497 * x +//
			0.020286308991329983 * x * x - 0.9780973533544137 * pow(x / 10.0, 3) +//
			0.26354051936651874 * pow(x / 10.0, 4) - 0.03711733070841712 * pow(x / 10.0, 5) +//
			0.002120935433043921 * pow(x / 10.0, 6);
	}
	return 0.0;
}

double MKmethod::int_1_f2(const double& x)
{
	if (x <= 1.0)
	{
		return 1.328216167939543 - 0.000004545681954848391 * x +//
			2.537368073155103 * x * x - 0.0020584991728545624 * x * x * x -//
			0.03742568018912792 * x * x * x * x - 0.010312136385277346 * x * x * x * x * x +//
			0.002767736179209713 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return 1.2959616295629246 + 0.1533684067037866 * x +//
			2.2354849981206106 * x * x + 0.3113395567715921 * x * x * x -//
			0.21656309882941488 * x * x * x * x + 0.041957500887605075 * x * x * x * x * x -//
			0.0029978773724628604 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 1.903643456971281 - 1.4801836911099535 * x + 3.973958664572268 * x * x -//
			0.6482729779428982 * x * x * x + 0.07665007314658864 * x * x * x * x -//
			0.005369758193338703 * x * x * x * x * x + 0.00016605531871992049 * x * x * x * x * x * x;
	}
	else if (x <= 7.0)
	{
		return -4.484415105552316 + 5.3747429756690766 * x +//
			0.8892806582308143 * x * x + 0.09767316152573671 * x * x * x -//
			0.025704749778475783 * x * x * x * x + 0.0021937998296249206 * x * x * x * x * x -//
			0.00006928845984076111 * x * x * x * x * x * x;
	}
	return 0.0;
}

double MKmethod::int_1_f3(const double& x)
{
	if (x <= 1.0)
	{
		return 1.2938345594193854 - 0.000031719847351174835 * x +//
			1.3183710041280094 * x * x - 0.014150512069488197 * x * x * x +//
			0.4226114681928129 * x * x * x * x - 0.06985750969880078 * x * x * x * x * x -//
			0.015347864048406958 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return 0.9667460440956788 + 1.336271810704016 * x -//
			0.8687355257991665 * x * x + 1.7676868273627229 * x * x * x -//
			0.2731222764016417 * x * x * x * x + 0.004801770033831665 * x * x * x * x * x +//
			0.001780776080720323 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 4.760566734123174 - 5.048204299463048 * x + 3.332342585228025 * x * x +//
			0.47584339615235993 * x * x * x - 0.12072786272726124 * x * x * x * x +//
			0.011870955604980658 * x * x * x * x * x - 0.0004580199652402304 * x * x * x * x * x * x;//
	}
	else if (x <= 7.0)
	{
		return 9.370493362469261 - 10.848615619383615 * x + 6.423326878282571 * x * x -//
			0.4148977656870439 * x * x * x + 0.025300923044176957 * x * x * x * x -//
			0.0010108688120876522 * x * x * x * x * x + 0.00001864423130429156 * x * x * x * x * x * x;
	}
	cout << " ERRRROROROROROROO MKMETHOD 2251" << endl;
	return 0.0;
}

double MKmethod::int_1(const double& x, const double& cp)
{
	double b = 1.0 - a_2 * log(cp);
	return (cp / (sqrtpi_ * sqrtpi_ * sqrtpi_)) * (b * b * int_1_f1(x/cp) - 2.0 * a_2 * b * int_1_f2(x/cp) + kv(a_2) * int_1_f3(x/cp));
}

double MKmethod::int_2_f1(const double& x)
{
	if (x <= 1.0)
	{
		return 8.377571213788123 * x + 0.00047608508679086725 * x * x +//
			1.6710478320575737 * x * x * x + 0.016857530811432916 * x * x * x * x -//
			0.15132474125724812 * x * x * x * x * x + 0.030723378194358945 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return -0.11788367995598747 + 8.937936705157014 * x -//
			1.119886471634001 * x * x + 2.8831031948885917 * x * x * x -//
			0.735250146386514 * x * x * x * x + 0.10356311378423572 * x * x * x * x * x -//
			0.006231417172309398 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 2.9044497739429858 + 2.757712415967557 * x + 4.239161941189675 * x * x +//
			0.36198838294786784 * x * x * x - 0.05737777787138304 * x * x * x * x +//
			0.004956250079677106 * x * x * x * x * x - 0.0001809238236975877 * x * x * x * x * x * x;
	}
	else if (x <= 7.0)
	{
		return 41.6323689028892 - 38.118317864344135 * x + 22.211189528076645 * x * x -//
			3.8547348524931246 * x * x * x + 0.5000517174807501 * x * x * x * x -//
			0.03446294709493891 * x * x * x * x * x + 0.0009860204070962582 * x * x * x * x * x * x;
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

double MKmethod::int_2_f2(const double& x)
{
	if (x <= 1.0)
	{
		return 3.8653461103376667 * x + 0.0001975300512691014 * x * x +//
			2.4468141895384012 * x * x * x + 0.005987984681429616 * x * x * x * x -//
			0.06453987836713967 * x * x * x * x * x + 0.0066920981111229004 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return -0.10983889480661446 + 4.321087890898017 * x -//
			0.7707850845797619 * x * x + 3.1237901158486583 * x * x * x -//
			0.31485222316123385 * x * x * x * x + 0.010270249760261791 * x * x * x * x * x +//
			0.0008259803934338584 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 1.8468011847509729 + 1.1986396743254275 * x +//
			1.1421489029509448 * x * x + 2.606316149569781 * x * x * x -//
			0.2788783468089509 * x * x * x * x + 0.019815317035281846 * x * x * x * x * x -//
			0.0006379970557448899 * x * x * x * x * x * x;
	}
	else if (x <= 7.0)
	{
		return 9.480707804348185 - 8.022988228952784 * x + 5.823555900242488 * x * x +//
			1.3277220473440972 * x * x * x - 0.08074921612981413 * x * x * x * x +//
			0.0033058587723954185 * x * x * x * x * x - 0.00006041810279926061 * x * x * x * x * x * x;
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

double MKmethod::int_2_f3(const double& x)
{
	if (x <= 1.0)
	{
		return 2.6106039258326 * x - 0.0008357997793049243 * x * x +//
			2.0764571907368174 * x * x * x - 0.03182275644841273 * x * x * x * x +//
			0.26310521962808975 * x * x * x * x * x - 0.06034325471643871 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return 0.20784760901369737 + 1.5920325291316857 * x +//
			2.0985329535259014 * x * x - 0.26286255221171206 * x * x * x +//
			1.4610329096120886 * x * x * x * x - 0.25626862029131897 * x * x * x * x * x +//
			0.01684969647300594 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return -6.284115352064703 + 15.665162343948523 * x //
			- 10.766105772158252 * pow(x, 2) + 6.0821074614870465 * pow(x, 3) -//
			0.3181501403196319 * pow(x, 4) + 0.01232319194701587 * pow(x, 5)//
			- 0.00017890661550876597 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return -4.355962170454177 + 13.835332665069274 * x //
			- 10.12766071646978 * pow(x, 2) + 5.999392227482686 * pow(x, 3) -//
			0.32171318647989955 * pow(x, 4) + 0.014181987261856027 * pow(x, 5)//
			- 0.00030579035551447497 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

double MKmethod::int_2(const double& x, const double& cp)
{
	double b = 1.0 - a_2 * log(cp);
	return -(cp * cp / (sqrtpi_ * sqrtpi_ * sqrtpi_)) * (b * b * int_2_f1(x / cp) - 2.0 * a_2 * b * int_2_f2(x / cp) + kv(a_2) * int_2_f3(x / cp));
}

double MKmethod::int_3(const double& x, const double& cp)
{
	double b = 1.0 - a_2 * log(cp);
	//cout << "1 = " << int_3_f1(x / cp) << endl;
	//cout << "2 = " << int_3_f2(x / cp) << endl;
	//cout << "3 = " << int_3_f3(x / cp) << endl;
	return ( pow(cp, 3.0) / (sqrtpi_ * sqrtpi_ * sqrtpi_)) * (b * b * int_3_f1(x / cp) - 2.0 * a_2 * b * int_3_f2(x / cp) + kv(a_2) * int_3_f3(x / cp));
}


double MKmethod::int_3_f1(const double& x)
{
	if (x <= 1.0)
	{
		return 12.566370586001975 - 0.00001944816384202852 * x + 12.567558607381049 * pow(x, 2) -//
			0.010507444068349692 * pow(x, 3) + 1.2911398125420694 * pow(x, 4) - 0.05048482363937502 * pow(x, 5) -//
			0.029947937607835744 * pow(x, 6);
	}
	else if (x <= 3.0)
	{
		return 12.451555448799724 + 0.40252674353016715 * x + 12.081033298182223 * pow(x, 2) + 0.12939193776331415 * pow(x, 3) +//
			1.478876561367302 * pow(x, 4) - 0.2237491583356496 * pow(x, 5) + 0.014474521138620033 * pow(x, 6);
	}
	else if (x <= 5.0)
	{
		return 8.281962852844913 + 9.783032429527339 * x + 3.165848344887614 * pow(x, 2) + 4.711462666119614 * pow(x, 3) +//
			0.13739453130827392 * pow(x, 4) - 0.012096015315889195 * pow(x, 5) + 0.0004514225555943018 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return 17.29966669035025 + 2.1457227895928668 * x + 5.572082327920818 * pow(x, 2) + 4.4083748449004645 * pow(x, 3) +//
			0.13640200155890422 * pow(x, 4) - 0.00854302147917508 * pow(x, 5) + 0.00022205921430504255 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

double MKmethod::int_3_f2(const double& x)
{
	if (x <= 1.0)
	{
		return 5.798024979296493 - 0.00001772671478406096 * x + 9.98769025405073 * pow(x, 2) - 0.0073593516014156535 * pow(x, 3) +//
			2.27820901023418 * pow(x, 4) - 0.03135086958655956 * pow(x, 5) - 0.030403716978821237 * pow(x, 6);
	}
	else if (x <= 3.0)
	{
		return 5.864728705834779 - 0.34799875480550213 * x + 10.760249318358127 * pow(x, 2) - 0.9493738978205943 * pow(x, 3) +//
			2.948810798239708 * pow(x, 4) - 0.29690823082625284 * pow(x, 5) + 0.015284639719532296 * pow(x, 6);
	}
	else if (x <= 5.0)
	{
		return -3.21810405152587 + 17.08054504204813 * x - 3.43427364926184 * pow(x, 2) + 5.342946243936558 * pow(x, 3) +//
			1.3459970589832742 * pow(x, 4) - 0.07448103623442687 * pow(x, 5) + 0.0021616888853670958 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return -59.42860988375287 + 79.24709914283142 * x - 32.304295129821256 * pow(x, 2) + 12.558114389999929 * pow(x, 3) +//
			0.32138208180756495 * pow(x, 4) + 0.003976327853989966 * pow(x, 5) - 0.0003703205838879336 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

double MKmethod::int_3_f3(const double& x)
{
	if (x <= 1.0)
	{
		return 3.915885322797866 + 0.000044284982651632276 * x + 7.778946280540595 * pow(x, 2) + 0.01990833982643192 * pow(x, 3) +//
			2.710531013102172 * pow(x, 4) + 0.09480498076917274 * pow(x, 5) + 0.026454483470905545 * pow(x, 6);
	}
	else if (x <= 3.0)
	{
		return 4.371710139648201 - 1.787430730965795 * x + 10.530386754306065 * pow(x, 2) - 1.9945949141813966 * pow(x, 3) +//
			3.310710912254515 * pow(x, 4) + 0.1356460090430464 * pow(x, 5) - 0.019853464614841342 * pow(x, 6);
	}
	else if (x <= 5.0)
	{
		return 4.6940208896644435 - 6.423081982183987 * x + 18.137713177954524 * pow(x, 2) - 7.298291873534797 * pow(x, 3) +//
			5.202095367169497 * pow(x, 4) - 0.20620743501245353 * pow(x, 5) + 0.005094092519200813 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return 194.2826492092263 - 195.77327124311523 * x + 95.83253215419927 * pow(x, 2) - 23.99281397751645 * pow(x, 3) +//
			7.170549820159753 * pow(x, 4) - 0.32568903981020303 * pow(x, 5) + 0.007955090180048264 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}
