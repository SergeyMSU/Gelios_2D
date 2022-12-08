#include "voro++.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "Help.h"
#include "Setka.h"
#include "Cell.h"


// ��� �� ��������

// ������ 1: �������� � ������ ����� ����� �� ���������, ����� ����� � ������ �������� ������ � ������� ����������� ��� ������ �����
// ����� ��� �� ���������� � ������� ����
// ������! ��� �� �� ������ ��������� �������
// ����� ��� ��������� �������� ������ � ������ ���������, � ����� ��������� � �����������, ������� ����� ����� ��� ���� ��������!
// ������ �� ����� �� ���������� ���������, � ������� �����������, ������ ��� ��� ���� ����� ���������! ��� ���������, �� ������ �� ������
// ����� ���������� ���������! ��. ������ 3

// ������ 2: ��� � ���������� �������, �� ������������ ���������, ������ �� ������� � ������. ����� ������� �� �����
// ����� ���������� ��������� � ����
// ����� ���������� ������� - �������

// ������ 3: �� �� �����, ��� ������ 1, �� ����� ��������� ����������� ������

// �����: 
// ������ �� ��������� ������ ���������� �����, ������� ��� ��������� � ����������  (������ ����� ��� �� ����� ��������� � ������)

using namespace voro;
using namespace std;

double rnd() { return double(rand()) / RAND_MAX; }

int sign(double x)
{
	if (x > 0)
	{
		return 1;
	}
	else if (x < 0)
	{
		return -1;
	}
	return 0;
}


int main()
{
	Setka A = Setka();
	A.Construct_start();
	A.save_soseds();
	//A.HLLD_Test();
	//return 0;
	A.Add_couple_start();
	//A.Print_points();
	A.Construct_start();
	A.Disable_cells();
	//A.Construct_start();
	A.Reconstruct_couple(); 
	A.Reconstruct_couple();
	A.Reconstruct_couple(true);
	A.set_normal();
	A.move_par();
	A.Reconstruct_couple();
	A.Reconstruct_couple();
	A.Reconstruct_couple();
	A.Calc_normal();
	A.Reconstruct_couple();


	A.Cut_Plane_z(1.0);
	A.Cut_Plane_y(1.0);
	A.Cut_Surface();
	//A.Cut_Plane_y_Tecplot(1.0);
	A.Init();
	A.Download_MHD("Vers_1.txt");
	//A.Go_MHD(30);
	A.Calc_normal();                           // ��������� � ����� ������, ��� ������� ��� ���-�� �� ���������
	//A.Start_MHD(100);
	A.Zapusk();
	//A.Reconstruct_couple();
	A.Cut_Plane_z_Tecplot(1.0);
	A.Cut_Plane_y_Tecplot(1.0);
	A.Save_MHD("Vers_2.txt");
	A.Cut_Plane_z(1.0);
	A.Cut_Plane_y(1.0);
	A.Cut_Surface();

	/*for (auto& i : A.All_Cell[80000]->Grans)
	{
		cout << i->S << " " << i->Sosed->number << " " << i->n1 << " " << i->n2 << " " << i->n3 << " " <<  endl;
	}*/

	return 0;
	double x, y, z, rsq, r;
	vector<int> neigh, f_vert, ney, edge, vershin;
	vector<double> X, Y, Ver, Versh;

	X.push_back(-1.0);
	Y.push_back(0.0);

	X.push_back(0.0);
	Y.push_back(0.0);

	X.push_back(1.0);
	Y.push_back(0.0);

	X.push_back(-1.0);
	Y.push_back(1.0);

	X.push_back(0.0);
	Y.push_back(1.0);

	X.push_back(1.0);
	Y.push_back(1.0);

	X.push_back(0.5);
	Y.push_back(0.5);

	X.push_back(-0.5);
	Y.push_back(0.5);

	X.push_back(-1.0);
	Y.push_back(2.0);

	X.push_back(0.0);
	Y.push_back(0.5);

	X.push_back(0.0);
	Y.push_back(2.0);

	double nx, ny, nz, nn;

	cout << "START " << endl;

	ofstream fout;
	string name_f = "setka.txt";
	fout.open(name_f);

	ofstream fout2;
	name_f = "points.txt";
	fout2.open(name_f);

	int N = 0;
	int E = 0;

	for (int i = 0; i < X.size(); i++) // X.size()
	{
		fout2 << X[i] << " " << Y[i] << endl;
	}

	for (int i = 0; i < X.size(); i++) // X.size()
	{
		cout << "D = " << i << endl;
		voronoicell_neighbor c;
		c.init(-20.0, 20.0, -20.0, 20.0, -20.0, 20.0);

		for (int j = 0; j < X.size(); j++)
		{
			if (i == j)
			{
				continue;
			}

			// ������� ����������� ������, �� ��� ����. ����� ��������� ������� �������� � ��������� ���������� �� ���� (� �������������, �� ��� ����)
			double A = (X[j] - X[i]);
			double B = (Y[j] - Y[i]);
			double C = 0.0;
			c.nplane(A, B, C, (A * A + B * B + C * C), j);
			// ������ ����������� ����� �������� ����� ���� ��������� ������ ��� �� ������ (����� ���������� ������)
		}

		// ���� ��� ����� ��� ������ ����� � ���� ��������. ��������� � ���, ��� ����� �� ������� ����������� ��� ���� � ����� ����� (������� ��� ��������)
		c.nplane(0.0, 0.0, 1.0, 0.01, -7);

		c.neighbors(ney); // �������� ������ ������� ������� �� �������
		int k = -1;  // ����� ������� ������
		// ���� ����� ������ "-7" - ������ ����� ��� ��������� ����� (���� ���������)  ����� - ���!!!
		for (int j = 0; j < ney.size(); j++)
		{
			k++;
			if (ney[j] == -7)
			{
				break;
			}
		}

		if (k == -1)
		{
			continue;
		}
		cout << "K = " << k << endl;
		//���� ������ �������:
		c.face_vertices(f_vert);
		int m = 0;
		int n = 0;
		int r = 0;
		while (r < k)
		{
			n = f_vert[m];
			m = m + n + 1;
			r++;
		}

		cout << "m = " << m << endl;

		for (int j = m + 1; j < m + 1 + f_vert[m] - 1; j++)
		{
			vershin.push_back(N + f_vert[j] + 1);  // ������ ������� ������ ������
			vershin.push_back(N + f_vert[j + 1] + 1);  // ������ ������� ������ ������
		}

		vershin.push_back(N + f_vert[m + 1 + f_vert[m] - 1] + 1);
		vershin.push_back(N + f_vert[m + 1] + 1);

		c.vertices(Ver);

		for (int j = 0; j < Ver.size(); j++)
		{
			if (j % 3 == 0)
			{
				Versh.push_back(Ver[j] + X[i]);  // ������ ���� ������
			}
			else if (j % 3 == 1)
			{
				Versh.push_back(Ver[j] + Y[i]);  // ������ ���� ������
			}
			else if(j % 3 == 2)
			{
				Versh.push_back(Ver[j]);  // ������ ���� ������
			}
		}

		N += Ver.size() / 3;
		E += f_vert[m];

		if (true)
		{
			printf("Total vertices      : %d\n", c.p);
			printf("Vertex positions    : "); c.output_vertices(); puts("");
			printf("Vertex orders       : "); c.output_vertex_orders(); puts("");
			printf("Max rad. sq. vertex : %g\n\n", 0.25 * c.max_radius_squared());
			cout << endl;
			printf("Total faces         : %d\n", c.number_of_faces());
			printf("Surface area        : %g\n", c.surface_area());
			printf("Face freq. table    : "); c.output_face_freq_table(); puts("");
			printf("Face orders         : "); c.output_face_orders(); puts("");
			printf("Face areas          : "); c.output_face_areas(); puts("");
			printf("Face normals        : "); c.output_normals(); puts("");
			printf("Face vertices       : "); c.output_face_vertices(); puts("\n");
			c.neighbors(ney);
			c.vertices(Ver);
			c.face_vertices(f_vert);
			for (int j = 0; j < Ver.size(); j++)
			{
				cout << Ver[j] << " ";
			}
			cout << endl;
			for (int j = 0; j < ney.size(); j++)
			{
				cout << ney[j] << " ";
			}
			cout << endl;
			for (int j = 0; j < f_vert.size(); j++)
			{
				cout << f_vert[j] << " ";
			}
			cout << endl;
		}

	}

	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=" << N << ", E=" << E << ", F=FEPOINT, ET=LINESEG " << endl;

	int jkl = 0;
	for (int j = 0; j < Versh.size(); j++)
	{
		jkl++;
		fout << Versh[j] << " ";
		if (jkl % 3 == 0)
		{
			fout << endl;
		}
	}

	jkl = 0;
	for (int j = 0; j < vershin.size(); j++)
	{
		jkl++;
		fout << vershin[j] << " ";
		if (jkl % 2 == 0)
		{
			fout << endl;
		}
	}




}