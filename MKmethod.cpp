#include "MKmethod.h"

MKmethod::MKmethod(void)
{
	double Yr = fabs(Velosity_inf);
	this->A0_ = this->A0(Yr);     // ����� ��������� ���������, ����� ��������� �� �������. 
	this->A1_ = 1.0 + (1.0 + 1.0 / (2.0 * kv(Yr))) * erf(Yr) + exp(-kv(Yr)) / (sqrtpi_ * Yr);
	this->num_area = I_; // I_;
	this->R_[0] = 1.0;                 // ����� �������� �������
	this->R_[1] = 3.0;
	this->R_[2] = 10.0;
	this->R_[3] = 30.0;
	this->R_[4] = 100.0;
	this->R_[5] = 250.0;
	this->R_[6] = 500.0;                 // ����� �������� �������
	this->R_[7] = 1000.0;

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
	double cc = c;
	bool minus = false;
	if (c < 0.0)
	{
		minus = true;
		cc = -c;
	}

	double ec = exp(2.0 * cc);
	double ln2 = log(2.0);
	double p1 = 0.25 * pi_ * ec + pi_/4.0;
	double p2;
	if (fabs(cc - ln2 / 2.0) < 0.001)
	{
		p2 = 0.5 * pi_ * ln2;
	}
	else
	{
		p2 = pi_ * ln2 * (1.0 - 2.0 * exp(-2.0 * cc)) / (2.0 * (2.0 * cc - ln2));
	}
	double p = p1 + p2;
	//cout << p1 << " " << p2 << endl;
	//cout << p1 << " " << p2 << " " << p << endl;
	//cout << c - ln2 / 2.0 << endl;
	double x = 0.0;
	double k = 0.0;
	double ksi1, ksi2;
	//int tt = 0;
	if (cc < 0.0001)
	{
		x = sens->MakeRandom() * pi_ / 2.0;
	}
	else
	{
		do
		{
			ksi1 = sens->MakeRandom();
			//ksi1 = 0.8;
			if (ksi1 < p1 / p)
			{
				//cout << "A" << endl;
				//tt = 1;
				//kappa = 1.13;
				ksi2 = sens->MakeRandom();
				//ksi2 = 0.4;
				x = (pi_ * ec - sqrtpi_ * sqrt(-4.0 * ec * ksi2 * p1 + pi_ * kv(ec) + 4.0 * ksi2 * p1)) / (2.0 * (exp(2.0 * cc) - 1.0));
				//cout << x << endl;
			}
			else
			{
				//cout << "B" << endl;
				//tt = 2;
				//kappa = 1.0;
				ksi2 = sens->MakeRandom(); // -p1 / p;
				//ksi2 = 0.4;
				if (fabs(cc - ln2 / 2.0) < 0.001)
				{
					//cout << "B1" << endl;
					x = 0.5 * pi_ * exp(2.0 * ksi2 * p2 / pi_);
					//cout << x << endl;
				}
				else
				{
					//cout << "B2" << endl;
					x = pow(2.0 * ln2, 1.0 / (2.0 * cc / ln2 - 1.0)) * pow(ec * pow(pi_, -2.0 * cc / ln2) * ((4.0 * cc - 2.0 * ln2) * (-ksi2 * p2 + pi_ * ln2/(4.0 * cc - 2.0 * ln2))),//
						2.0 * cc / (ln2 - 2.0 * cc) + 1.0);
					//cout << x << endl;
				}

				/*if (std::fpclassify(x) != FP_NORMAL && std::fpclassify(x) != FP_ZERO)
				{
					cout << x << "    ERROR  144     " << ksi1 << "  " << c <<  "  " << ec << "  " << p1 <<  endl;
					cout << pow(2.0 * ln2, 1.0 / (2.0 * c / ln2 - 1.0)) << endl;
					exit(-1);
				}*/
				//cout << x << endl;
			}

			if (this->h_mho(x, cc) / 1.13 > 1.0)
			{
				cout << "ERROR  179" << endl;
				cout << this->h_mho(x, cc) / 1.13 << endl;
				cout << "cc = " << cc << endl;
				exit(-1);
			}

		} while (this->h_mho(x, cc) / 1.13 < sens->MakeRandom());
	}
	//cout << tt << endl;

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

bool MKmethod::Init_Parametrs(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_)
{
	double Y = fabs(Velosity_inf);
	double ksi, ksi1, ksi2;
	double X1;
	double X2;
	double X0 = 1.0;                // ��� ������ ����
	double split, gam1, gam2;
	vector <double> Wa_(this->num_area);
	vector <double> Mho_(this->num_area);

	// �����������  X
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

		while (fabs(X2 - X1) > 0.0000001)     // ������� �������, ����� �������������
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
	// �����������  Wr
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
		while (fabs(Wr2 - Wr1) > 0.000001)     // ������� �������, ����� �������������
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
	// �����������  Wa
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

	// �����������  Mho
	for (int i = 0; i < this->num_area; i++)
	{
		Mho_[i] = play_mho(sens, Y * sqrt(1.0 - kv(X_[i])) * Wa_[i]);
		/*if (std::fpclassify(Mho_[i]) != FP_NORMAL && std::fpclassify(Mho_[i]) != FP_ZERO)
		{
			cout << Mho_[i] << "    ERROR  Mho_[i]  324   " << Y * sqrt(1.0 - kv(X_[i])) * Wa_[i] << endl;
			exit(-1);
		}*/
	}

	// ������� ����
	for (int i = 0; i < this->num_area; i++)
	{
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
		double c = Y * sqrt(1.0 - kv(X_[i])) * Wa_[i];
		mu_[i] = (this->F(1.0, gam2, Y) - this->F(1.0, gam1, Y)) * (this->norm_mho(c)) / (this->A0_ * (1.0 + kv(c)));
		//cout << mu_[i] << "     -  " << i <<  endl;
		//cout << (this->F(1.0, gam2, Y) - this->F(1.0, gam1, Y)) << "  " << this->A0_ << endl;
		if (std::fpclassify(mu_[i]) != FP_NORMAL && std::fpclassify(mu_[i]) != FP_ZERO)
		{
			cout << mu_[i] << "    ERROR  mu_[i]  348" << endl;
			cout << gam1 << " " << gam2 << endl;
			cout << this->F(1.0, gam2, Y) << endl;
			cout << this->F(1.0, gam1, Y) << endl;
			cout << c << endl;
			cout << this->A0_ << endl;
			cout << Wa_[i] << endl;
			exit(-1);
		}
	}
	//exit(-1);

	// ����������� �������� ����

	// �������� X ������� �� ���������
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
	
	double p4_ = sqrtpi_ * (Y * X2) / (1.0 + sqrtpi_ * (Y * X2));                  // ��� ��������� ��������� ����� �� ������� �� ��������� ������
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
	//} while (kv(Wt_[this->num_area]) + kv(Wp_[this->num_area]) < gg * kv(Wr_[this->num_area]));

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
		return true;
	}
	else
	{
		if (kv(Wt_[this->num_area]) + kv(Wp_[this->num_area]) > this->gam_[this->num_area - 1] * kv(Wr_[this->num_area]))
		{
			mu_[this->num_area] = 1.0;
			//cout << 1 << endl;
		}
		else
		{
			mu_[this->num_area] = 0.0;  // ����� �� ��������� ���� ����
			//cout << 0 << endl;
			return false;
		}
	}

	return true;
}

int MKmethod::Init(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, double& X_)
{
	double Y = fabs(Velosity_inf);
	double ksi, ksi1, ksi2;
	double X1;
	double X2;
	double X0 = 1.0;                // ��� ������ ����
	double split;

	// ����������� �������� ����
	ksi = sens->MakeRandom();
	X1 = 0.0;
	X2 = 1.0;
	//int k = 0;
	while (fabs(X2 - X1) > 0.00001)     // �������� ������ ������� (�������, �� �������)
	{
		split = X2;
		X2 = X1 - (X2 - X1) * (this->R(X1, Y) - ksi) / ((this->R(X2, Y) - ksi) - (this->R(X1, Y) - ksi));
		if (X2 > 1.0)  // ������ ��������, ������ X2 ������� �� �������
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
	double p4_ = sqrtpi_ * (Y * X2) / (1.0 + sqrtpi_ * (Y * X2));                  // ��� ��������� ��������� ����� �� ������� �� ��������� ������
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
	// ��������� Wr
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

		while (this->Hvr(gam1, gam2, Wr1, Ur, Uthe, ksi) >= 0.0)
		{
			Wr1 = Wr1 - 1.0;
		}
		int k = 0;
		while (fabs(Wr2 - Wr1) > 0.00001)     // ������� �������, ����� �������������
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
		if (std::fpclassify(Wr1) != FP_NORMAL && std::fpclassify(Wr1) != FP_ZERO)
		{
			cout << Wr1 << "    ERROR  Wr0 537" << endl;
			exit(-1);
		}
	}

	double W1, W2, Wa, ksi1, ksi2;
	// �����������  Wa
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

	// �����������  Mho
	for (int i = 0; i < I; i++)
	{
		Mho_[i] = play_mho(sens, Uthe * Wa_[i]);
		if (std::fpclassify(Mho_[i]) != FP_NORMAL && std::fpclassify(Mho_[i]) != FP_ZERO)
		{
			cout << Mho_[i] << "    ERROR  Mho_[i] 579" << endl;
			cout << "Uthe  " << Uthe << endl;
			cout << "Wa_[i]  " << Wa_[i] << endl;
			cout << Uthe * Wa_[i] << endl;
			exit(-1);
		}
	}

	// ������� ����
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

		mu_[i] = (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c));
		//mu_[i] = (u * sigma2(u, cp) / (IS)) * (f2(0.0, gam1, Ur, Uthe) - f2(0.0, gam2, Ur, Uthe)) * exp(-kv(Uthe)) * (this->norm_mho(c)) / (1.0 + kv(c));
		//cout << IS << " " << uu * sigma2(uu, cp) << endl;
		
		//cout << "mu = " << mu_[i] << endl;
		//cout << u << " " << uu << " " << gam1 << " " << gam2 << " " << Uthe << " " << Ur << " " << Wa_[i] << "  " << cp << endl;
		//cout << f2(0.0, gam1, Ur, Uthe) -  f2(0.0, gam2, Ur, Uthe) << endl;
		//cout << (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) << endl;
		//exit(-4);

		if (std::fpclassify(mu_[i]) != FP_NORMAL && std::fpclassify(mu_[i]) != FP_ZERO)
		{
			cout << mu_[i] << "    ERROR  mu_[i] 604" << endl;
			cout << "f2(0.0, gam1, Ur, Uthe) = " << f2(0.0, gam1, Ur, Uthe) << endl;
			cout << "f2(0.0, gam2, Ur, Uthe) = " << f2(0.0, gam2, Ur, Uthe) << endl;
			cout << "gam1 = " << gam1 << endl;
			cout << "gam2 = " << gam2 << endl;
			cout << "Ur = " << Ur << endl;
			cout << "Uthe = " << Uthe << endl;

			exit(-1);
		}

		//cout << mu_[i] << "  " << i << endl;
	}
	

	// �������� ��������� �����
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
			// ����� ����������� ��������
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
		mu_[I] = 0.0;  // ����� �� ��������� ���� ����
		//cout << 0 << endl;
		return false;
	}
	return true;
}

bool MKmethod::Change_Velosity2(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
	vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex,//
	const double& y_ex, const double& z_ex)
// �������� ������ �� ���������
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

	// ������� ����
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


	// ������ ��������� �����
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
			// ����� ����������� ��������
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
	} while (h <= ksi6 || (v1 <= 0 && kvv(v2, v3, 0.0) < gg * kv(v1)));
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
	//	mu_[I] = 0.0;  // ����� �� ��������� ���� ����
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
	// ��������� Wr
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
		while (fabs(Wr2 - Wr1) > 0.00001)     // ������� �������, ����� �������������
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
	// �����������  Wa
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

	// �����������  Mho
	for (int i = 0; i < I; i++)
	{
		Mho_[i] = play_mho(sens, Uthe * Wa_[i]);
	}

	// ������� ����
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


	// �������� ��������� �����
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
			// ����� ����������� ��������
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
		mu_[I] = 0.0;  // ����� �� ��������� ���� ����
		//cout << 0 << endl;
		return false;
	}
	return true;
}

bool MKmethod::Change_Velosity4(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
	vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex,//
	const double& y_ex, const double& z_ex)
{
	double X = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double UUr = Ur - Vr;
	double UUa = sqrt(kvv(0.0, Vthe - Uthe, Vphi - Uphi));
	double uu = exp(-kv(X)) / sqrtpi_ + (X + 1.0 / (2.0 * X)) * erf(X);
	double vh = sqrt(kvv(Vr, Vthe, Vphi));
	double up = sqrt(kvv(Ur, Uthe, Uphi));

	vector <double> gamma_(I);
	vector <double> Fint_(I);
	double IS;

	for (int i = 0; i < I; i++)
	{
		gamma_[i] = 1.0 / (kv(1.0 / this->alpha_[i]) - 1.0);
	}

	if (I == 1)
	{
		IS = this->Get_Int(X);
		Fint_[0] = this->Get_Int002(UUr, UUa);
	}
	else if (I == 2)
	{
		IS = this->Get_Int(X);
		Fint_[0] = this->Get_Int002(UUr, UUa);
		Fint_[1] = this->Get_Int00625(UUr, UUa) - Fint_[0];
	}
	else if (I == 3)
	{
		IS = this->Get_Int(X);
		Fint_[0] = this->Get_Int002(UUr, UUa);
		Fint_[1] = this->Get_Int00625(UUr, UUa) - Fint_[0];
		Fint_[2] = this->Get_Int02(UUr, UUa) - Fint_[1] - Fint_[0];
	}
	else if (I == 4)
	{
		IS = this->Get_Int(X);
		Fint_[0] = this->Get_Int002(UUr, UUa);
		Fint_[1] = this->Get_Int00625(UUr, UUa) - Fint_[0];
		Fint_[2] = this->Get_Int02(UUr, UUa) - Fint_[1] - Fint_[0];
		Fint_[3] = this->Get_Int055(UUr, UUa) - Fint_[2] - Fint_[1] - Fint_[0];
	}

	double e1 = sqrtpi_ * kv(Ur);
	double e2 = 2.0 * fabs(Ur);
	double e3 = 0.5 * sqrtpi_;
	double p1 = e1 / (e1 + e2 + e3);
	double p2 = e2 / (e1 + e2 + e3);
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7, ksi8;
	double z = 0.0, h = 0.0;
	double gam1, gam2;
	double The, Val;
	// ��������� �������������� �����

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
			} while (h < ksi6 || z >= -Ur);

			ksi7 = sens->MakeRandom();
			ksi8 = sens->MakeRandom();
			The = 2.0 * pi_ * ksi7;

			Val = sqrt(gam1 * kv(Wr_[i]) + ksi8 * (gam2 * kv(Wr_[i]) - gam1 * kv(Wr_[i])));

			Wthe_[i] = Val * cos(The);
			Wphi_[i] = Val * sin(The);

			//double c = Uthe * Wa_[i];
			double u = sqrt(kvv(Vr - Wr_[i], Vthe - Wthe_[i], Vphi - Wphi_[i]));
			//h = (u * sigma2(u, cp) / ((2.0 * vh + up) * sigma2(X, cp))) * (this->norm_mho(c)) / (1.0 + kv(c));
			//h = 0.6 * (u * sigma2(u, cp) / ((vh + up) * sigma2(X, cp))) * exp(-kv(Wthe_[i] - Uthe) - kv(Wphi_[i] - Uphi));
			h = (1.0/97.0) * (u * sigma2(u, cp)) * exp(-kv(Wthe_[i] - Uthe) - kv(Wphi_[i] - Uphi));
			if (h > 1.0)
			{
				cout << "h > 1  1555   ERRORR" << endl;
				cout << Wr_[i] << " " << Wthe_[i] << " " << Wphi_[i] << endl;
				cout << sqrt(kvv(0.0, Wthe_[i], Wphi_[i])) << endl;
				cout << gam1 << " " << gam2 << endl;
				cout << Val << endl;
				cout << h << endl;
				cout << u << endl;
				cout << (2.0 * vh + up) << endl;
				cout << sigma2(u, cp) << endl;
				cout << sigma2(X, cp) << endl;
				cout << Uthe << endl;
				cout << (Uthe * X + 0.1) << endl;
				exit(-2);
			}
		} while (h < sens->MakeRandom());
		mu_[i] = Fint_[i] / IS;
		if (mu_[i] < 0.0)
		{
			cout << " < 0  Error   1573  frvbjkutrex" << endl;
			cout << UUr << " " << UUa << endl;
			cout << IS << endl;
			cout << this->Get_Int002(UUr, UUa) << endl;
			cout << this->Get_Int00625(UUr, UUa) << endl;
			cout << this->Get_Int02(UUr, UUa) << endl;
			cout << this->Get_Int055(UUr, UUa) << endl;
			for (int k = 0; k < I; k++)
			{
				cout << Fint_[k] << endl;
			}
			cout << mu_[i] << endl;
			exit(-5);
		}
		if (std::fpclassify(mu_[i]) != FP_NORMAL && std::fpclassify(mu_[i]) != FP_ZERO)
		{
			cout << mu_[i] << "    ERROR  mu_[i] 1331" << endl;
			exit(-1);
		}
	}

	// �������� ��������� �����
	double p4 = 0.5 * sqrtpi_ * X / (1.0 + 0.5 * sqrtpi_ * X);
	double om1, om2, om3, lo;
	double y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy;

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
	} while (h < ksi6 || (v1 < 0.0 && kv(v2) + kv(v3) < gg * kv(v1)));


	Wr_[I] = v1;
	Wthe_[I] = v2;
	Wphi_[I] = v3;

	mu_[I] = 1.0;
	for (int i = 0; i < I; i++)
	{
		mu_[I] = mu_[I] - mu_[i];
	}

	if (mu_[I] < 0.0)
	{
		cout << " < 0  Error   1383  dgvewjkl;kjh" << endl;
		cout << mu_[I] << endl;
		exit(-5);
	}

	/*if (true)
	{
		for (int i = 0; i <= I; i++)
		{
			double aa, bb, cc;
			dekard_skorost(y_ex, z_ex, x_ex, Wr_[i], Wphi_[i], Wthe_[i], bb, cc, aa);
			for (int ij = 0; ij < 27000; ij++)
			{
				cout << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << i << endl;
			}
		}
		exit(-1);
	}*/

	return true;
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