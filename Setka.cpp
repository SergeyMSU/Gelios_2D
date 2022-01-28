#include "Setka.h"
#include <algorithm>
#include <omp.h>
#include <mutex>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

Setka::Setka()
{
	this->Inizialization();
}

Setka::Setka(int N1, int N2, int N3, int N4, int M1, int M2, int M3, int M4)
// Начальная инициализация сетки,  эту функцию лучше не трогать, она писалась поэтапно и результаты постоянно проверялись
// Изменение чего-либо может привести к ошибке
{
	this->Inizialization();
	this->N1 = N1;
	this->N2 = N2;
	this->N3 = N3;
	this->N4 = N4;
	this->M1 = M1;
	this->M2 = M2;
	this->M3 = M3;
	this->M4 = M4;
	N1++;
	N3++;
	M1++;
	double s = 0.0;
	int kk, kkk, kkkk;
	int n = M1 + M2 + M3 + M4 - 1; // Число граней по лучу. 
	kk = n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1);
	kkk = kk + (N4 - 1) * (M2 + M3 + M4);

	auto Cell_Centr = new Cell();; // Центральная особая ячейка (создаётся здесь отдельно, чтобы положить в общий массив только в конце связки всего)

	for (int i = 0; i < N1; i++)
	{
		//double x = 1.0 * i / (N1 - 1);
		//s = (0.35 * x * (1.0 - x) + kv(x)) * (pi_ / 2.0);  //i * (pi_ / 2.0) / (N1 - 1.0);
		s = i * (pi_ / 2.0) / (N1 - 1.0);
		auto R = new Rail(s);
		R->type = A;
		this->A_Rails.push_back(R);
		R->M1 = M1;
		R->M2 = M2;
		R->M3 = M3;
		R->M4 = M4;

		R->Init_start(R);
		for (auto& j : R->All_point)
		{
			this->All_Points.push_back(j);
		}
	}

	for (int i = 0; i < N2; i++)
	{
		s = (pi_ / 2.0) + (i + 1) * (pi_ / 4.0) / (N2);
		auto R = new Rail(s);
		R->type = B;
		this->B_Rails.push_back(R);
		R->M1 = M1;
		R->M2 = M2;
		R->M3 = M3;
		R->M4 = M4;
		R->Init_start(R);
		for (auto& j : R->All_point)
		{
			this->All_Points.push_back(j);
		}
	}

	for (int i = 0; i < N3; i++)
	{
		s = (3.0 * pi_ / 4.0) + (i) * (pi_ / 4.0) / (N3 - 1);
		auto R = new Rail(s);
		R->type = C;
		this->C_Rails.push_back(R);
		R->M1 = M1;
		R->M2 = N4;
		R->M3 = 0.0;
		R->M4 = 0.0;
		if (i > 0)
		{
			R->Init_start(R);
		}
	}

	for (int i = 0; i < N4; i++)
	{
		s = -R1_ - (i + 1) * (-R1_ - Left_) / (N4);
		auto R = new Rail(s);
		R->type = D;
		this->D_Rails.push_back(R);
		R->M1 = N3;
		R->M2 = M2;
		R->M3 = M3;
		R->M4 = M4;
	}

	// Отдельно нужно обработать С-type линию от тройной точки
	auto C = this->C_Rails[0];
	for (int i = 0; i < M1; i++)
	{
		C->All_point.push_back(this->B_Rails[N2 - 1]->All_point[i]);
		if (i == M1 - 1)
		{
			C->Key_point.push_back(this->B_Rails[N2 - 1]->All_point[i]);
		}
	}
	double l;
	for (int i = 0; i < N4; i++)
	{
		double x = pow(Left_ / (R2_ * cos(3.0 * pi_ / 4.0)), 1.0 / (N4));
		l = R2_ * cos(3.0 * pi_ / 4.0) * pow(x, i + 1); ; // R2_* cos(3.0 * pi_ / 4.0) - (i + 1) * (R2_ * cos(3.0 * pi_ / 4.0) - Left_) / N4;
		auto K = new Point(l, R2_ * sin(3.0 * pi_ / 4.0));
		C->All_point.push_back(K);
		this->All_Points.push_back(K);
		K->type = P_U_2;
		if (i == N4 - 1)
		{
			C->Key_point.push_back(K);
			K->type = P_Outer_Boandary;
		}
	}
	// Обработка завершена

	for (int i = 1; i < N3; i++)
	{
		for (auto& j : this->C_Rails[i]->All_point)
		{
			this->All_Points.push_back(j);
		}
	}

	// Отдельно обрабатываем четвёртый тип линий
	for (int i = 0; i < N4; i++)
	{
		auto P = C->All_point[M1 + i];
		this->D_Rails[i]->All_point.push_back(P);
		this->D_Rails[i]->Key_point.push_back(P);
		for (int j = 0; j < M2; j++)
		{
			l = P->y + (j + 1) * (R3_ - P->y) / (M2);
			auto K = new Point(P->x, l);
			this->D_Rails[i]->All_point.push_back(K);
			this->All_Points.push_back(K);
			K->type = P_U_3;
			if (j == M2 - 1)
			{
				this->D_Rails[i]->Key_point.push_back(K);
				K->type = P_Contact;
			}
		}
		for (int j = 0; j < M3; j++)
		{
			l = R3_ + (j + 1) * (R4_ - R3_) / M3;
			auto K = new Point(P->x, l);
			this->D_Rails[i]->All_point.push_back(K);
			this->All_Points.push_back(K);
			K->type = P_U_4;
			if (j == M3 - 1)
			{
				this->D_Rails[i]->Key_point.push_back(K);
				K->type = P_Outer_shock;
			}
		}
		double x = pow(R5_ / R4_, 1.0 / (M4));
		for (int j = 0; j < M4; j++)
		{
			l = R4_ * pow(x, j + 1);  //R4_ + (j + 1) * (R5_ - R4_) / M4;
			auto K = new Point(P->x, l);
			this->D_Rails[i]->All_point.push_back(K);
			this->All_Points.push_back(K);
			K->type = P_U_5;
			if (j == M4 - 1)
			{
				this->D_Rails[i]->Key_point.push_back(K);
				K->type = P_Outer_Boandary;
			}
		}
	}

	// Теперь начинается создание самих ячеек
	for (int i = 0; i < this->A_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->A_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->A_Rails[i]->All_point[j], this->A_Rails[i]->All_point[j + 1], //
				this->A_Rails[i + 1]->All_point[j + 1], this->A_Rails[i + 1]->All_point[j]);
			this->A_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->A_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->A_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->A_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}

	for (int i = 0; i < this->B_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->B_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->B_Rails[i]->All_point[j], this->B_Rails[i]->All_point[j + 1], //
				this->B_Rails[i + 1]->All_point[j + 1], this->B_Rails[i + 1]->All_point[j]);
			this->B_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->B_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->B_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->B_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}

	for (int i = 0; i < this->C_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->C_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->C_Rails[i]->All_point[j], this->C_Rails[i]->All_point[j + 1], //
				this->C_Rails[i + 1]->All_point[j + 1], this->C_Rails[i + 1]->All_point[j]);
			this->C_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->C_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->C_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->C_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}

	for (int i = 0; i < this->D_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->D_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->D_Rails[i]->All_point[j], this->D_Rails[i]->All_point[j + 1], //
				this->D_Rails[i + 1]->All_point[j + 1], this->D_Rails[i + 1]->All_point[j]);
			this->D_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->D_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->D_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->D_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}


	// Связка по-середине
	auto P1 = this->A_Rails[this->A_Rails.size() - 1];
	auto P2 = this->B_Rails[0];
	for (int j = 0; j < P1->All_point.size() - 1; j++)
	{
		auto C = new Cell(P1->All_point[j], P1->All_point[j + 1], //
			P2->All_point[j + 1], P2->All_point[j]);
		P1->All_point[j]->my_cell.push_back(C);
		P1->All_point[j + 1]->my_cell.push_back(C);
		P2->All_point[j + 1]->my_cell.push_back(C);
		P2->All_point[j]->my_cell.push_back(C);
		this->All_Cells.push_back(C);
	}

	// Связка от тройной точки
	P1 = this->B_Rails[this->B_Rails.size() - 1];
	P2 = this->D_Rails[0];
	for (int j = 0; j < P2->All_point.size() - 1; j++)
	{
		auto C = new Cell(P1->All_point[j + M1 - 1], P1->All_point[j + M1], //
			P2->All_point[j + 1], P2->All_point[j]);
		P1->All_point[j + M1 - 1]->my_cell.push_back(C);
		P1->All_point[j + M1]->my_cell.push_back(C);
		P2->All_point[j + 1]->my_cell.push_back(C);
		P2->All_point[j]->my_cell.push_back(C);
		this->All_Cells.push_back(C);
	}

	// На этом этапе все узлы связаны, ячейки созданы. Теперь нужно связать ячейки. Создать грани. ---------------------------------------------

	// Пронумеруем узлы и ячейки.
	int num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;  
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}

	n = M1 + M2 + M3 + M4 - 1; // Число граней по лучу. 

	if (true)
	{

		for (int i = 0; i < n * (N1 - 1); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1) % (n) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Input);
				K->Grans.push_back(G);
				G->main_gran = true;
				G->Master = K;
				this->All_Gran.push_back(G);
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = 0; i < n * (N1 - 2); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + n];
			G2->Master = this->All_Cells[i + n];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + n]->Grans.push_back(G2);
		}


		for (int i = n * (N1 - 1); i < n * (N2 + N1 - 2); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1) % (n) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
				G->main_gran = true;
				this->All_Gran.push_back(G);
				K->Grans.push_back(G);
				G->Master = K;
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = n * (N1 - 1); i < n * (N2 + N1 - 3); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + n];
			G2->Master = this->All_Cells[i + n];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + n]->Grans.push_back(G2);
		}


		for (int i = n * (N2 + N1 - 2); i < n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1 - n * (N2 + N1 - 2)) % (N4 + M1 - 1) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Extern);
				G->main_gran = true;
				this->All_Gran.push_back(G);
				K->Grans.push_back(G);
				G->Master = K;
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = n * (N2 + N1 - 2); i < n * (N2 + N1 - 2) + (N3 - 2) * (N4 + M1 - 1); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + (N4 + M1 - 1)];
			G2->Master = this->All_Cells[i + (N4 + M1 - 1)];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + (N4 + M1 - 1)]->Grans.push_back(G2);
		}


		kk = n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1);

		for (int i = kk; i < kk + (N4 - 1) * (M2 + M3 + M4); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1 - kk) % (M2 + M3 + M4) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
				G->main_gran = true;
				K->Grans.push_back(G);
				G->Master = K;
				this->All_Gran.push_back(G);
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = kk; i < kk + (N4 - 2) * (M2 + M3 + M4); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + (M2 + M3 + M4)];
			G2->Master = this->All_Cells[i + (M2 + M3 + M4)];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + (M2 + M3 + M4)]->Grans.push_back(G2);
		}


		// Связываем шов поцентру
		kkk = kk + (N4 - 1) * (M2 + M3 + M4);

		for (int i = kkk; i < kkk + M1 + M2 + M3 + M4 - 2; i++)
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = kkk; i < kkk + M1 + M2 + M3 + M4 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[n * (N1 - 2) + i - kkk];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		for (int i = kkk; i < kkk + M1 + M2 + M3 + M4 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[n * (N1 - 1) + i - kkk];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}


		kkkk = n * (N2 + N1 - 2);
		for (int i = kkkk; i < kkkk + M1 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[i - n];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}


		for (int i = this->All_Cells.size() - M2 - M3 - M4; i < this->All_Cells.size() - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = this->All_Cells.size() - M2 - M3 - M4; i < this->All_Cells.size(); i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[kk + i - this->All_Cells.size() + M2 + M3 + M4];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		for (int i = this->All_Cells.size() - M2 - M3 - M4; i < this->All_Cells.size(); i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[kkkk + M1 - n + i - (this->All_Cells.size() - M2 - M3 - M4) - 1];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		if (true)  // Ручное связывание одного из узлов
		{
			auto K = this->All_Cells[this->All_Cells.size() - M2 - M3 - M4];
			auto F = this->All_Cells[kkkk + M1 - 1];
			auto G = new Gran(K->contour[3], K->contour[0], Usualy);
			auto G2 = new Gran(K->contour[0], K->contour[3], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		for (int i = kkkk + M1; i < kkkk + M1 + N4 - 1; i++)  // Последний шов от тройной точки вбок
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[kk + (i - kkkk - M1) * (M2 + M3 + M4)];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		// Добавляем грани - границы

		// Начнём с оси симметрии

		for (int i = 0; i < M1 + M2 + M3 + M4 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[0], K->contour[1], Axis);
			auto G2 = new Gran(K->contour[1], K->contour[0], Axis);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Master = nullptr;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G2->Sosed = K;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

		for (int i = n * (N2 + N1 - 2) + (N3 - 2) * (N4 + M1 - 1); i < n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Axis);
			auto G2 = new Gran(K->contour[3], K->contour[2], Axis);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Master = nullptr;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}


		for (int i = 0; i < N1 + N2 - 2; i++)
		{
			auto K = this->All_Cells[i * (M1 + M2 + M3 + M4 - 1)];
			auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
			auto G2 = new Gran(K->contour[0], K->contour[3], Inner_sphere);  // Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			Cell_Centr->Grans.push_back(G2);
			G->Master = K;
			G2->Master = Cell_Centr;
			G->Sosed = Cell_Centr;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}


		// Шов по центру с верхней стенкой
		if (true)
		{
			auto K = this->All_Cells[kkk + M1 + M2 + M3 + M4 - 2];
			auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
			auto G2 = new Gran(K->contour[2], K->contour[1], Upper_wall);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

		if (true)
		{
			auto K = this->All_Cells[this->All_Cells.size() - 1];
			auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
			auto G2 = new Gran(K->contour[2], K->contour[1], Upper_wall);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

	}

	int h = (N1 + N2 - 2) * (M1 + M2 + M3 + M4 - 1);
	for (int i = 0; i < N3 - 1; i++)
	{
		auto K = this->All_Cells[i * (M1 + N4 - 1) + h];
		auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
		auto G2 = new Gran(K->contour[0], K->contour[3], Inner_sphere);  //Usualy);
		G->main_gran = true;
		G2->main_gran = false;
		K->Grans.push_back(G);
		Cell_Centr->Grans.push_back(G2);
		G->Master = K;
		G2->Master = Cell_Centr;
		G->Sosed = Cell_Centr;
		G2->Sosed = K;
		G->Gran_copy = G2;
		G2->Gran_copy = G;
		this->All_Gran.push_back(G);
		this->All_Gran_copy.push_back(G2);
	}

	// Шов по центру связываем с внутренней границей
	auto K = this->All_Cells[kkk];
	auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
	auto G2 = new Gran(K->contour[0], K->contour[3], Inner_sphere); // Usualy);
	G->main_gran = true;
	G2->main_gran = false;
	K->Grans.push_back(G);
	Cell_Centr->Grans.push_back(G2);
	G->Master = K;
	G2->Master = Cell_Centr;
	G->Sosed = Cell_Centr;
	G2->Sosed = K;
	G->Gran_copy = G2;
	G2->Gran_copy = G;
	this->All_Gran.push_back(G);
	this->All_Gran_copy.push_back(G2);

	for (int i = kk + (N4 - 2) * (M2 + M3 + M4); i < kk + (N4 - 1) * (M2 + M3 + M4); i++) // Связываем ячейки по лучу (продольная связка)
	{
		auto K = this->All_Cells[i];
		auto G = new Gran(K->contour[2], K->contour[3], Extern);
		auto G2 = new Gran(K->contour[3], K->contour[2], Extern);
		G->main_gran = true;
		G2->main_gran = true;
		K->Grans.push_back(G);
		G->Master = K;
		G2->Sosed = K;
		G->Gran_copy = G2;
		G2->Gran_copy = G;
		this->All_Gran.push_back(G);
		this->All_Gran_copy.push_back(G2);
	}

	for (auto& i : this->All_Cells)  // Нужно для выделения контактной поверхности
	{
		if (i->contour[0]->type == P_Inner_Boandary || i->contour[0]->type == P_U_1)
		{
			i->type = C_1;
		}
		else if (i->contour[3]->type == P_U_2 || i->contour[2]->type == P_U_2)
		{
			i->type = C_2;
		}
		
		if (i->contour[0]->type == P_U_3 || i->contour[1]->type == P_U_3 || i->contour[2]->type == P_U_3)
		{
			i->type = C_3;
		}
		else if (i->contour[0]->type == P_U_4 || i->contour[1]->type == P_U_4 || i->contour[2]->type == P_U_4)
		{
			i->type = C_4;
		}
		else if(i->contour[0]->type == P_U_5 || i->contour[1]->type == P_U_5 || i->contour[2]->type == P_U_5)
		{
			i->type = C_5;
		}
	}


	// Хотим добавить особые ячейки.
	vector<Cell*> centr_cell;
	if (true)
	{
		sort(Cell_Centr->Grans.begin(), Cell_Centr->Grans.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});

		int kl = 8;
		while (Cell_Centr->Grans.size() % kl != 0)
		{
			kl++;
		}
		int how = Cell_Centr->Grans.size() / kl;

		for (int i = 0; i < kl; i++)
		{
			auto A = new Cell();
			A->type = C_centr;
			for (int j = i * how; j < (i + 1) * how; j++)
			{
				A->Grans.push_back(Cell_Centr->Grans[j]);
				Cell_Centr->Grans[j]->Master = A;
				Cell_Centr->Grans[j]->Gran_copy->Sosed = A;
			}
			centr_cell.push_back(A);
		}

		auto O = new Point(0.0, 0.0);
		O->type = P_Null;
		this->All_Points.push_back(O);

		for (auto& i : centr_cell)
		{
			i->contour.push_back(O);
			O->my_cell.push_back(i);
			for (auto& j : i->Grans)
			{
				i->contour.push_back(j->A);
				j->A->my_cell.push_back(i);
			}
			i->contour.push_back(i->Grans[i->Grans.size() - 1]->B);
			i->Grans[i->Grans.size() - 1]->B->my_cell.push_back(i);
		}

		for (int i = 0; i < centr_cell.size() - 1; i++)
		{
			auto K = centr_cell[i];
			int vm = K->contour.size();
			auto G = new Gran(K->contour[vm - 1], K->contour[0], Usualy);
			auto G2 = new Gran(K->contour[0], K->contour[vm - 1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			centr_cell[i+1]->Grans.push_back(G2);
			G->Master = K;
			G2->Master = centr_cell[i + 1];
			G->Sosed = centr_cell[i + 1];
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

		auto PO = centr_cell[0];
		auto G2 = new Gran(PO->contour[0], PO->contour[1], Axis);
		auto G3 = new Gran(PO->contour[1], PO->contour[0], Axis);
		G2->main_gran = true;
		G3->main_gran = false;
		PO->Grans.push_back(G2);
		G2->Master = PO;
		G3->Sosed = PO;
		G3->Gran_copy = G2;
		G2->Gran_copy = G3;
		this->All_Gran.push_back(G2);
		this->All_Gran_copy.push_back(G3);

		PO = centr_cell[centr_cell.size() - 1];
		G2 = new Gran(PO->contour[PO->contour.size() - 1], PO->contour[0], Axis);
		G3 = new Gran(PO->contour[0], PO->contour[PO->contour.size() - 1], Axis);
		G2->main_gran = true;
		G3->main_gran = false;
		PO->Grans.push_back(G2);
		G2->Master = PO;
		G3->Sosed = PO;
		G3->Gran_copy = G2;
		G2->Gran_copy = G3;
		this->All_Gran.push_back(G2);
		this->All_Gran_copy.push_back(G3);


		/*this->All_Cells.push_back(Cell_Centr);
		Cell_Centr->type = C_centr;
		for (int i = 0; i < this->A_Rails.size(); i++)
		{
			Cell_Centr->contour.push_back(this->A_Rails[i]->All_point[0]);
		}
		for (int i = 0; i < this->B_Rails.size(); i++)
		{
			Cell_Centr->contour.push_back(this->B_Rails[i]->All_point[0]);
		}
		for (int i = 1; i < this->C_Rails.size(); i++)
		{
			Cell_Centr->contour.push_back(this->C_Rails[i]->All_point[0]);
		}

		auto G2 = new Gran(Cell_Centr->contour[Cell_Centr->contour.size() - 1], Cell_Centr->contour[0], Axis);
		G2->main_gran = true;
		Cell_Centr->Grans.push_back(G2);
		G2->Master = Cell_Centr;
		this->All_Gran.push_back(G2);*/
	}

	for (auto& i : centr_cell)
	{
		this->All_Cells.push_back(i);
	}


	delete Cell_Centr;


	// Нумерация всего
	num = 0;
	for (auto& i : this->All_Gran)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
		if (i->A->type == P_Contact && i->B->type == P_Contact)
		{
			this->Line_Contact.push_back(i);
		}
		sort(this->Line_Contact.begin(), this->Line_Contact.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});

		if (i->A->type == P_Inner_shock && i->B->type == P_Inner_shock)
		{
			this->Line_Inner.push_back(i);
		}

		sort(this->Line_Inner.begin(), this->Line_Inner.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});

		if (i->A->type == P_Outer_shock && i->B->type == P_Outer_shock)
		{
			this->Line_Outer.push_back(i);
		}

		sort(this->Line_Outer.begin(), this->Line_Outer.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});
	}

	num = 0;
	for (auto& i : this->All_Gran_copy)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}


	cout << "Vsego tochek = " << this->All_Points.size() << endl;
	cout << "Vsego Yacheek = " << this->All_Cells.size() << endl;
}

void Setka::Inizialization(void)
{
	mmu1 = 0.0;
	mmu2 = 0.0;
	mmu3 = 0.0;
	mmu4 = 0.0;
	mmu5 = 0.0;
	mn1 = 0;
	mn2 = 0;
	mn3 = 0;
	mn4 = 0;
	mn5 = 0;
	f_num = 0;
	this->Cell_m = nullptr;
	f_way.open("file_work.txt");

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
			for (int k = 0; k < J_ + 1; k++)
			{
				Mu_stat[i][j][k] = 0.0;
				Num_stat[i][j][k] = 0;
				Sinus[k] = sin(k * pi_ / 12.0);
			}
		}
	}
	Sinus[0] = sin(pi_ / 12.0);
	Sinus[11] = sin(pi_ / 12.0);

	// Считаем веса
	int i1, i2, i3;
	ifstream fout;
	fout.open("weyght_1.txt");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
			for (int k = 0; k < J_ + 1; k++)
			{
				fout >> i1 >> i2 >> i3 >> Mu_stat[i][j][k];
			}
		}
	}

	Ri.resize(I_ + 1);
	//Ri[0] = R5_ - 2.0;
	Ri[0] = 1.0;                 // Здесь задаются радиусы
	Ri[1] = 3.0;
	Ri[2] = 10.0;
	Ri[3] = 30.0;
	Ri[4] = 100.0;
	Ri[5] = 250.0;
	Ri[6] = 500.0;                 // Здесь задаются радиусы
	Ri[7] = 1000.0;
	Ri[8] = Rmax_;
}

//void Setka::TVD_prepare(void)
//{
//	double n1, n2, n;
//	double x, y;
//	double x1, y1;
//	double x3, y3;
//	double a, b, N;
//	Cell* A = nullptr;
//	cout << "TVD  prepare" << endl;
//	for (auto& i : this->All_Gran)
//	{
//		N = 1.0;
//		A = nullptr;
//		if (i->type == Axis)
//		{
//			i->Get_Center(x, y);
//			i->Get_normal(n1, n2);
//			n1 = -n1;
//			n2 = -n2;
//			for (auto& j : i->Master->Grans)
//			{
//				if (j->type == Usualy)
//				{
//					j->Sosed->Get_Center(x1, y1);
//
//					a = x1 - x;
//					b = y1 - y;
//					n = sqrt(kv(a) + kv(b));
//					a = a / n;
//					b = b / n;
//					//cout << fabs(a * n1 + b * n2) << endl;
//					if (fabs(1.0 - fabs(a * n1 + b * n2)) < N)
//					{
//						N = fabs(1.0 - fabs(a * n1 + b * n2));
//						A = j->Sosed;
//					}
//				}
//
//			}
//			i->Sosed_down = A;
//			continue;
//		}
//		if (i->type != Usualy)
//		{
//			continue;
//		}
//		i->Sosed->Get_Center(x3, y3);
//		i->Get_Center(x, y);
//		n1 = x3 - x;                            // Направление от центра грани к центру соседа для грани
//		n2 = y3 - y;
//		n = sqrt(kv(n1) + kv(n2));
//		n1 = n1 / n;
//		n2 = n2 / n;
//
//		for (auto& j : i->Sosed->Grans)
//		{
//			if (j->type != Usualy || j->Sosed->number != i->Master->number)
//			{
//				if (j->Sosed != nullptr)
//				{
//					j->Sosed->Get_Center(x1, y1);
//				}
//				else
//				{
//					j->Get_Center(x1, y1);
//				}
//
//				a = x1 - x;
//				b = y1 - y;
//				n = sqrt(kv(a) + kv(b));
//				a = a / n;
//				b = b / n;
//				//cout << fabs(a * n1 + b * n2) << endl;
//				if (fabs(1.0 - fabs(a * n1 + b * n2)) < N)
//				{
//					N = fabs(1.0 - fabs(a * n1 + b * n2));
//					A = j->Sosed;
//				}
//			}
//			
//		}
//		i->Sosed_up = A;
//	}
//
//	for (auto& i : this->All_Gran_copy)
//	{
//		A = nullptr;
//		if (i->type != Usualy)
//		{
//			continue;
//		}
//		N = 1.0;
//		i->Sosed->Get_Center(x3, y3);
//		i->Get_Center(x, y);
//		n1 = x3 - x;
//		n2 = y3 - y;
//		n = sqrt(kv(n1) + kv(n2));
//		n1 = n1 / n;
//		n2 = n2 / n;
//		for (auto& j : i->Sosed->Grans)
//		{
//			if (j->type != Usualy || j->Sosed->number != i->Master->number)
//			{
//				if (j->Sosed != nullptr)
//				{
//					j->Sosed->Get_Center(x1, y1);
//				}
//				else
//				{
//					j->Get_Center(x1, y1);
//				}
//
//				a = x1 - x;
//				b = y1 - y;
//				n = sqrt(kv(a) + kv(b));
//				a = a / n;
//				b = b / n;
//				//cout << fabs(a * n1 + b * n2) << endl;
//				if (fabs(1.0 - fabs(a * n1 + b * n2)) < N)
//				{
//					N = fabs(1.0 - fabs(a * n1 + b * n2));
//					A = j->Sosed;
//				}
//			}
//
//		}
//		i->Sosed_up = A;
//	}
//
//	for (auto& i : this->All_Gran)
//	{
//		if (i->Gran_copy != nullptr)
//		{
//			i->Sosed_down = i->Gran_copy->Sosed_up;
//		}
//	}
//
//	for (auto& i : this->All_Gran_copy)
//	{
//		if (i->Gran_copy != nullptr)
//		{
//			i->Sosed_down = i->Gran_copy->Sosed_up;
//		}
//	}
//
//	cout << "TVD  prepare  end" << endl;
//}

void Setka::TVD_prepare(void)
{
	double n1, n2, n;
	double x, y;
	double x1, y1;
	double x3, y3;
	double a, b, N;
	Cell* A = nullptr;
	Gran* B = nullptr;
	cout << "TVD  prepare" << endl;

	// Находим Sosed_up для этой грани
	for (auto& i : this->All_Gran)
	{
		N = -1.0;
		B = nullptr;
		if (i->Sosed == nullptr)
		{
			i->Sosed_up = nullptr;
			continue;
		}

		i->Sosed->Get_Center(x3, y3);
		i->Get_Center(x, y);
		n1 = x3 - x;                            // Направление от центра грани к центру соседа для грани
		n2 = y3 - y;
		n = sqrt(kv(n1) + kv(n2));
		n1 = n1 / n;
		n2 = n2 / n;

		for (auto& j : i->Sosed->Grans)
		{
			if (j->type == Usualy)
			{
				if (j->Gran_copy->number == i->number)
				{
					continue;
				}
			}

			j->Get_Center(x1, y1);
			a = x1 - x3;
			b = y1 - y3;
			n = sqrt(kv(a) + kv(b));
			a = a / n;
			b = b / n;
			if (a * n1 + b * n2 > N)
			{
				N = a * n1 + b * n2;
				B = j;
			}
		}

		if (B->Sosed != nullptr)
		{
			i->Sosed_up = B->Sosed;
		}
		else
		{
			i->Sosed_up = nullptr;
		}
	}

	for (auto& i : this->All_Gran_copy)
	{
		N = -1.0;
		B = nullptr;
		if (i->Sosed == nullptr)
		{
			if (i->type == Axis)
			{
				cout << "Axis " << endl;
			}
			i->Sosed_up = nullptr;
			continue;
		}

		i->Sosed->Get_Center(x3, y3);
		i->Get_Center(x, y);
		n1 = x3 - x;                            // Направление от центра грани к центру соседа для грани
		n2 = y3 - y;
		n = sqrt(kv(n1) + kv(n2));
		n1 = n1 / n;
		n2 = n2 / n;

		for (auto& j : i->Sosed->Grans)
		{
			if (j->type == Usualy)
			{
				if (j->Gran_copy->number == i->number)
				{
					continue;
				}
			}

			j->Get_Center(x1, y1);
			a = x1 - x3;
			b = y1 - y3;
			n = sqrt(kv(a) + kv(b));
			a = a / n;
			b = b / n;
			if (a * n1 + b * n2 > N)
			{
				N = a * n1 + b * n2;
				B = j;
			}
		}

		if (B->Sosed != nullptr)
		{
			i->Sosed_up = B->Sosed;
		}
		else
		{
			i->Sosed_up = nullptr;
		}
	}

	for (auto& i : this->All_Gran)
	{
		if (i->Gran_copy != nullptr)
		{
			i->Sosed_down = i->Gran_copy->Sosed_up;
		}
		else
		{
			i->Sosed_down = nullptr;
		}
	}

	for (auto& i : this->All_Gran_copy)
	{
		if (i->Gran_copy != nullptr)
		{
			i->Sosed_down = i->Gran_copy->Sosed_up;
		}
		else
		{
			i->Sosed_down = nullptr;
		}
	}

	cout << "TVD  prepare  end" << endl;
}

void Setka::Print_TVD(void)
{
	int ll = this->All_Gran.size();
	ofstream fout;
	fout.open("TVD.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"C\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	double x1, y1;
	double x2, y2;
	int nn = 0;
	for (auto& i : this->All_Gran)
	{
		nn++;
		i->Get_Center(x1, y1);
		if (i->Sosed_down != nullptr)
		{
			i->Sosed_down->Get_Center(x2, y2);
		}
		else
		{
			x2 = x1;
			y2 = y1;
		}
		fout << x1 << " " << y1 << " " << nn << endl;
		fout << x2 << " " << y2 << " " << nn << endl;
	}

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_point()
{
	ofstream fout;
	fout.open("point.txt");

	for (auto& i : this->All_Points)
	{
		fout << i->x << " " << i->y << " " << i->type << endl;
	}
}

void Setka::Print_Gran()
{
	int ll = this->Line_Contact.size() + this->Line_Inner.size() + this->Line_Outer.size();
	ofstream fout;
	fout.open("surfaces.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	double x1, y1;
	for (auto& i : this->Line_Contact)
	{
		fout << i->A->x << " " << i->A->y << endl;
		fout << i->B->x << " " << i->B->y << endl;
	}
	for (auto& i : this->Line_Inner)
	{
		fout << i->A->x << " " << i->A->y << endl;
		fout << i->B->x << " " << i->B->y << endl;
	}
	for (auto& i : this->Line_Outer)
	{
		fout << i->A->x << " " << i->A->y << endl;
		fout << i->B->x << " " << i->B->y << endl;
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_cell_type(void)
{
	ofstream fout;
	fout.open("Cell_type.txt");
	double x, y;
	for (auto& i : this->All_Cells)
	{
		i->Get_Center(x, y);
		fout << x << " " << y << " " << i->type << endl;
	}
}

void Setka::Print_Gran_type(void)
{
	ofstream fout;
	fout.open("Gran_type.txt");
	double x, y;
	for (auto& i : this->All_Gran)
	{
		i->Get_Center(x, y);
		fout << x << " " << y << " " << i->type << endl;
	}
}

void Setka::Print_cell(void)
{
	int ll = this->All_Cells.size();
	ofstream fout;
	fout.open("Setka_all_cell.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 4 * ll;
	fout << " , E= " << 4 * ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Cells)
	{
		for (auto& j : i->contour)
		{
			fout << j->x << " " << j->y << endl;
		}
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 4 * i + 1 << " " << 4 * i + 2 << endl;
		fout << 4 * i + 2 << " " << 4 * i + 3 << endl;
		fout << 4 * i + 3 << " " << 4 * i + 4 << endl;
		fout << 4 * i + 4 << " " << 4 * i + 1 << endl;
	}

	fout.close();
}

void Setka::Print_cell2(void)
{
	int ll = 0;
	for (auto& i : this->All_Cells)
	{
		ll += i->contour.size();
	}

	ofstream fout;
	fout.open("Setka_all_cell.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << this->All_Points.size();
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Points)
	{
		fout << i->x << " " << i->y << endl;
	}

	for (auto& i : this->All_Cells)
	{
		for (int j = 0; j < i->contour.size() - 1; j++)
		{
			fout << i->contour[j]->number + 1 << " " << i->contour[j + 1]->number + 1 << endl;
		}
		fout << i->contour[i->contour.size() - 1]->number + 1 << " " << i->contour[0]->number + 1 << endl;
	}

	fout.close();
}

void Setka::Print_connect(void)
{
	int ll = 0;
	for (auto& i : this->All_Cells)
	{
		ll += i->Grans.size();
	}
	double x, y;
	double x2, y2;
	ofstream fout;
	fout.open("Setka_connect.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << ll * 2;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Cells)
	{
		for (auto& j : i->Grans)
		{
			if (false) // (j->type == Usualy)
			{
				i->Get_Center(x, y);
				j->Sosed->Get_Center(x2, y2);
				fout << x << " " << y << endl;
				fout << (x + x2) / 2.0 << " " << (y + y2) / 2.0 << endl;
			}
			else
			{
				i->Get_Center(x, y);
				j->Get_Center(x2, y2);
				fout << x << " " << y << endl;
				fout << x2 << " " << y2 << endl;
			}
		}
	}

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_point_connect(void)
{
	int ll = 0;
	for (auto& i : this->All_Points)
	{
		ll += i->my_cell.size();
	}
	double x, y;
	double x2, y2;
	ofstream fout;
	fout.open("Setka_connect_point.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << ll * 2;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Points)
	{
		x = i->x;
		y = i->y;
		for (auto& j : i->my_cell)
		{
			j->Get_Center(x2, y2);
			fout << x << " " << y << endl;
			fout << x2 << " " << y2 << endl;
		}
	}

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_Tecplot(void)
{
	double r_o = 1.0;    // Размер расстояния
	double ro_o = 0.06;   // Размер плотности
	double ro_o_H = 0.18;   // Размер плотности
	double p_o = 1.0;   // Размер давления
	double u_o = 10.38;   // Размер давления

	ofstream fout;
	string name_f = "2D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\",\"Ro_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", \"RO_H\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->All_Cells)
	{
		double kk;
		if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		{
			kk = (chi_real / chi_);
		}
		else
		{
			kk = 1.0;
		}
		//double kk = 1.0;
		double Max = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			QQ = i->par[0].Q / i->par[0].ro;
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}

		fout << x * r_o << " " << y * r_o << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o/kv(kk) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " " //
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << //
			" " << i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
			(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << endl;

	}
	fout.close();

	name_f = "1D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"Ro_r_r\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\",\"Ro_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"M_H1\", \"T_H1\", \"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", \"RO_H\", ZONE T = \"HP\"" << endl;
	int num = 0;
	//double ro = (389.988 * 389.988) / (chi_ * chi_);


	fout << 1 * r_o << " " << 0.0 * r_o << " " << 1.0 << //
		" " << 0.06 * (389.988 * 389.988) / (chi_real * chi_real) << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0  << //
		" " << 0.001666666 * ro_o_H << " " << (0.001666666 * chi_real * chi_real / (ggg * 50.0 * 50.0)) * pow(1.0 / 1.0, 2.0 * ggg) << " " //
		<< chi_real * u_o << " " << 0.0 << //
		" " << 50.0 << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << //
		" " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << //
		" " << 0.0 << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << " " <<//
		0.0 << endl;


	for (auto& i : this->All_Cells)
	{
		double kk;
		if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		{
			kk = (chi_real / chi_);
		}
		else
		{
			kk = 1.0;
		}
		//kk = 1.0;
		num++;
		if (num > this->M1 + this->M2 + this->M3 + this->M4)
		{
			break;
		}
		double Max = 0.0;
		double Max_H1 = 0.0;
		double T_H1 = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			QQ = i->par[0].Q / i->par[0].ro;
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}
		if (i->par[0].ro_H1 > 0.000000000001)
		{
			Max_H1 = sqrt(kvv(i->par[0].u_H1, i->par[0].v_H1, 0.0) / (ggg * i->par[0].p_H1 / i->par[0].ro_H1));
			T_H1 = 2.0 * (i->par[0].p_H1 / i->par[0].ro_H1) * 6527.0;
		}
		double dist = sqrt(kv(x) + kv(y));
		fout << x * r_o << " " << y * r_o << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o / kv(kk) << " " << 0.06 * (389.988 * 389.988) / (chi_real * chi_real)/ (dist * dist) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " "//
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << " " << Max_H1 << " " << T_H1 << " "//
			<< i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
			(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << endl;

	}
	fout.close();
}

void Setka::Print_Sourse(void)
{
	string name_f = "1D_Sourse.txt";
	ofstream fout;
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"I1_mc_1\",  \"I1_mc_2\", \"I1_mc_3\", \"I1_mc_4\", \"I2_mc_1\",  \"I2_mc_2\", \"I2_mc_3\", \"I2_mc_4\", \"I3_mc_1\",  \"I3_mc_2\", \"I3_mc_3\", \"I3_mc_4\", \"I1_mf_1\",  \"I1_mf_2\", \"I1_mf_3\", \"I1_mf_4\", \"I2_mf_1\",  \"I2_mf_2\", \"I2_mf_3\", \"I2_mf_4\", \"I3_mf_1\",  \"I3_mf_2\", \"I3_mf_3\", \"I3_mf_4\",  " << //
		"ZONE T = \"HP\"" << endl;
	int num = 0;
	for (auto& i : this->All_Cells)
	{
		num++;
		if (num > this->M1 + this->M2 + this->M3 + this->M4)
		{
			break;
		}
		
		double x, y;
		i->Get_Center(x, y);
		i->Calc_Sourse();
		double dist = sqrt(kv(x) + kv(y));
		fout << sqrt(kvv(0.0, x, y)) << //
			" " << i->par[0].I1_mc[0] << " " << i->par[0].I1_mc[1] << " " << i->par[0].I1_mc[2] << " " << i->par[0].I1_mc[3] << //
			" " << i->par[0].I2_mc[0] << " " << i->par[0].I2_mc[1] << " " << i->par[0].I2_mc[2] << " " << i->par[0].I2_mc[3] << //
			" " << i->par[0].I3_mc[0] << " " << i->par[0].I3_mc[1] << " " << i->par[0].I3_mc[2] << " " << i->par[0].I3_mc[3] << //
			" " << i->par[0].I1_mf[0] << " " << i->par[0].I1_mf[1] << " " << i->par[0].I1_mf[2] << " " << i->par[0].I1_mf[3] << //
			" " << i->par[0].I2_mf[0] << " " << i->par[0].I2_mf[1] << " " << i->par[0].I2_mf[2] << " " << i->par[0].I2_mf[3] << //
			" " << i->par[0].I3_mf[0] << " " << i->par[0].I3_mf[1] << " " << i->par[0].I3_mf[2] << " " << i->par[0].I3_mf[3] << endl;

	}
	fout.close();
}

void Setka::Print_Tecplot_MK(void)
{
	double r_o = 1.0;    // Размер расстояния
	double ro_o = 0.06;   // Размер плотности
	double ro_o_H = 0.18;    //0.18;   // Размер плотности
	double p_o = 1.0;    //1.0;   // Размер давления
	double u_o = 10.38;   // Размер скорости

	ofstream fout;
	string name_f = "2D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Rho\", \"P\", \"Vx\", \"Vy\", \"M\",\"Q\",\"Rho_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"Rho_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Rho_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Rho_H4\", \"P_H4\"," << //
		" \"Vx_H4\", \"Vy_H4\", \"RO_H\", \"F_n\", \"F_u\", \"F_v\", \"F_T\", \"I_u\", \"I_v\", \"I_T\", \"II_u\", \"II_v\", \"II_T\", \"M_u\", \"M_v\", \"M_T\",\"H1_n\", \"H1_u\", \"H1_v\", \"H1_T\"," << //
		"\"H2_n\", \"H2_u\", \"H2_v\", \"H2_T\", \"H3_n\", \"H3_u\", \"H3_v\", \"H3_T\", \"H4_n\", \"H4_u\", \"H4_v\", \"H4_T\", \"k_u\", \"k_v\", \"k_T\", \"m_1\",  \"m_2\",  \"m_3\",  \"m_4\",  \"m_5\",  \"m_6\",  \"m_7\",  \"num1\", \"num2\", \"num3\", \"num4\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->All_Cells)
	{
		double kk = 1.0;
		/*if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		{
			kk = (chi_real / chi_);
		}
		else
		{
			kk = 1.0;
		}*/
		double Max = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			QQ = i->par[0].Q / i->par[0].ro;
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}

		fout << x * r_o << " " << y * r_o << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o / kv(kk) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " " //
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << //
			" " << i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
			(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << //
			" " << i->par[0].F_n << " " << i->par[0].F_u << " " << i->par[0].F_v << " " << i->par[0].F_T << " " //
			<< i->par[0].I_u << " " << i->par[0].I_v << " " << i->par[0].I_T << " " << //
			i->par[0].II_u << " " << i->par[0].II_v << " " << i->par[0].II_T << " " << //
			i->par[0].M_u << " " << i->par[0].M_v << " " << i->par[0].M_T << " " << //
			i->par[0].H_n[0] << " " << i->par[0].H_u[0] << " " << i->par[0].H_v[0] << " " << i->par[0].H_T[0] << " " //
			<< i->par[0].H_n[1] << " " << i->par[0].H_u[1] << " " << i->par[0].H_v[1] << " " << i->par[0].H_T[1] << " " //
			<< i->par[0].H_n[2] << " " << i->par[0].H_u[2] << " " << i->par[0].H_v[2] << " " << i->par[0].H_T[2] << " " //
			<< i->par[0].H_n[3] << " " << i->par[0].H_u[3] << " " << i->par[0].H_v[3] << " " << i->par[0].H_T[3] << //
			" " << i->par[0].k_u << " " << i->par[0].k_v << " " << i->par[0].k_T << //
			" " << i->par[0].w_m[0] << //
			" " << i->par[0].w_m[1] << //
			" " << i->par[0].w_m[2] << //
			" " << i->par[0].w_m[3] << //
			" " << i->par[0].w_m[4] << //
			" " << i->par[0].w_m[5] << //
			" " << i->par[0].w_m[6] << //
			" " << i->par[0].num_atoms[0] << //
			" " << i->par[0].num_atoms[1] << //
			" " << i->par[0].num_atoms[2] << //
			" " << i->par[0].num_atoms[3] << endl;

	}
	fout.close();

	name_f = "1D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"Ro_r_r\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\",\"Ro_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"M_H1\", \"T_H1\", \"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", \"RO_H\", \"I_u\", \"I_v\", \"I_T\"," << //
		"ZONE T = \"HP\"" << endl;
	int num = 0;
	//double ro = (389.988 * 389.988) / (chi_ * chi_);


	fout << 1 * r_o << " " << 0.0 * r_o << " " << 1.0 << //
		" " << 0.06 * (389.988 * 389.988) / (chi_real * chi_real) << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << //
		" " << 0.001666666 * ro_o_H << " " << (0.001666666 * chi_real * chi_real / (ggg * 50.0 * 50.0)) * pow(1.0 / 1.0, 2.0 * ggg) << " " //
		<< chi_real * u_o << " " << 0.0 << //
		" " << 50.0 << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << //
		" " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << //
		" " << 0.0 << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << " " <<//
		0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 <<  endl;


	for (auto& i : this->All_Cells)
	{
		double kk;
		/*if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		{
			kk = (chi_real / chi_);
		}
		else
		{
			kk = 1.0;
		}*/
		kk = 1.0;
		num++;
		if (num > this->M1 + this->M2 + this->M3 + this->M4)
		{
			break;
		}
		double Max = 0.0;
		double Max_H1 = 0.0;
		double T_H1 = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			QQ = i->par[0].Q / i->par[0].ro;
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}
		if (i->par[0].ro_H1 > 0.000000000001)
		{
			Max_H1 = sqrt(kvv(i->par[0].u_H1, i->par[0].v_H1, 0.0) / (ggg * i->par[0].p_H1 / i->par[0].ro_H1));
			T_H1 = 2.0 * (i->par[0].p_H1 / i->par[0].ro_H1) * 6527.0;
		}
		double dist = sqrt(kv(x) + kv(y));
		fout << x * r_o << " " << y * r_o << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o / kv(kk) << " " << 0.06 * (389.988 * 389.988) / (chi_real * chi_real) / (dist * dist) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " "//
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << " " << Max_H1 << " " << T_H1 << " "//
			<< i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
			(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << //
			" " << i->par[0].I_u << " " << i->par[0].I_v << " " << i->par[0].I_T <<  endl;

	}
	fout.close();

	ofstream f_weyght;
	f_weyght.open("weyght.txt");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
			for (int k = 0; k < J_ + 1; k++)
			{
				f_weyght << i << " " << j << " " << k << " " << kv(Rmax_/ Ri[j]) * Mu_stat[i][j][k] * 12.0/(411.0 * 2) << endl;
			}
		}
	}
	f_weyght.close();
}

void Setka::Proverka(void)
{
	cout << "Proverka  Start" << endl;
	// Выписываем ячейки у которых не 4 соседа - особые
	for (auto& i : this->All_Cells)  
	{
		if (i->Grans.size() != 4)
		{
			double x, y;
			i->Get_Center(x, y);
			cout << "Gran  82753   in  point:  " << x << " " <<  y << "   Have  sosed: " << i->Grans.size() << endl;
		}
	}


	// Проверяем правильное определение нормалей к граням
	for (auto& i : this->All_Gran)
	{
		if (i->type == Usualy)
		{
			double xc, yc;
			double n1, n2;
			double x2, y2;
			double x, y;
			i->Master->Get_Center(x, y);
			i->Get_Center(xc, yc);
			i->Sosed->Get_Center(x2, y2);
			i->Get_normal(n1, n2);
			if (n1 * (x2 - xc) + n2 * (y2 - yc) < 0)
			{
				cout << "ERROR  gran  2345656435434543234:  " << xc << " " << yc << endl;
			}
			if (n1 * (x - xc) + n2 * (y - yc) > 0)
			{
				cout << "ERROR  gran  wecwcefwcfwcfwcwcfwf:  " << xc << " " << yc << endl;
			}
		}

	}

	cout << "Proverka TVD" << endl;

	for (auto& i : this->All_Gran)
	{
		if (i->type == Usualy)
		{
			Parametr par;
			i->Get_par_TVD(par, 0);
			if ((par.ro < i->Master->par[0].ro && par.ro < i->Sosed->par[0].ro) || (par.ro > i->Master->par[0].ro && par.ro > i->Sosed->par[0].ro))
			{
				cout << "PROBLEM BIG" << endl;
			}
			if ((par.p < i->Master->par[0].p && par.p < i->Sosed->par[0].p) || (par.p > i->Master->par[0].p && par.p > i->Sosed->par[0].p))
			{
				cout << "PROBLEM BIG" << endl;
			}
			if ((par.u < i->Master->par[0].u && par.u < i->Sosed->par[0].u) || (par.u > i->Master->par[0].u && par.u > i->Sosed->par[0].u))
			{
				cout << "PROBLEM BIG" << endl;
			}
			if ((par.v < i->Master->par[0].v && par.v < i->Sosed->par[0].v) || (par.v > i->Master->par[0].v && par.v > i->Sosed->par[0].v))
			{
				cout << "PROBLEM BIG" << endl;
			}
			if ((par.Q < i->Master->par[0].Q && par.Q < i->Sosed->par[0].Q) || (par.Q > i->Master->par[0].Q && par.Q > i->Sosed->par[0].Q))
			{
				cout << "PROBLEM BIG" << endl;
			}
		}

	}

	cout << "Proverka TVD 2" << endl;

	for (auto& i : this->All_Gran_copy)
	{
		if (i->type == Usualy)
		{
			Parametr par;
			i->Get_par_TVD(par, 0);
			if ((par.ro < i->Master->par[0].ro && par.ro < i->Sosed->par[0].ro) || (par.ro > i->Master->par[0].ro && par.ro > i->Sosed->par[0].ro))
			{
				cout << "PROBLEM BIG" << endl;
			}
		}

	}
	
	cout << "Proverka TVD  3" << endl;

	// Проверяем расположение точек в ячейке по кругу!  (как грани расположены тоже надо будет проверить)
	for (auto& i : this->All_Cells)
	{
		bool t = false;
		for (int j = 0; j < i->contour.size() - 1; j++)
		{
			auto A = i->contour[j];
			auto B = i->contour[j + 1];
			for (auto& k : i->Grans)
			{
				if ((k->A == A && k->B == B) || (k->A == B && k->B == A))
				{
					t = true;
				}
			}
		}

		if (t == false)
		{
			cout << "EROR hygwyfwuhlduwhguydcw" << endl;
			for (auto& j: i->contour)
			{
				cout << j->x << " " << j->y << endl;
			}
			cout << "EROR hygwyfwuhlduwhguydcw" << endl;
		}
	}


	cout << "Proverka  End" << endl;
}

void Setka::Copy(Setka* S)
{
	for (auto& i : this->A_Rails)
	{
		for (auto& j : S->A_Rails)
		{
			if (i->s <= j->s)
			{
				//cout << i->s << " " << j->s << endl;
				i->Key_point[0]->Vx = j->Key_point[0]->x - i->Key_point[0]->x;
				i->Key_point[0]->Vy = j->Key_point[0]->y - i->Key_point[0]->y;
				i->Key_point[1]->Vx = j->Key_point[1]->x - i->Key_point[1]->x;
				i->Key_point[1]->Vy = j->Key_point[1]->y - i->Key_point[1]->y;
				i->Key_point[2]->Vx = j->Key_point[2]->x - i->Key_point[2]->x;
				i->Key_point[2]->Vy = j->Key_point[2]->y - i->Key_point[2]->y;
				break;
			}
		}
	}

	for (auto& i : this->B_Rails)
	{
		//cout << i->s << endl;
		for (auto& j : S->B_Rails)
		{
			if (i->s - 0.1 <= j->s)
			{
				//cout << i->s << " " << j->s << endl;
				i->Key_point[0]->Vx = j->Key_point[0]->x - i->Key_point[0]->x;
				i->Key_point[0]->Vy = j->Key_point[0]->y - i->Key_point[0]->y;
				i->Key_point[1]->Vx = j->Key_point[1]->x - i->Key_point[1]->x;
				i->Key_point[1]->Vy = j->Key_point[1]->y - i->Key_point[1]->y;
				i->Key_point[2]->Vx = j->Key_point[2]->x - i->Key_point[2]->x;
				i->Key_point[2]->Vy = j->Key_point[2]->y - i->Key_point[2]->y;
				break;
			}
		}
	}

	for (auto& i : this->C_Rails)
	{
		for (auto& j : S->C_Rails)
		{
			if (i->s - 0.001 <= j->s)
			{
				//cout << i->s << " " << j->s << endl;
				i->Key_point[0]->Vx = j->Key_point[0]->x - i->Key_point[0]->x;
				i->Key_point[0]->Vy = j->Key_point[0]->y - i->Key_point[0]->y;
				break;
			}
		}
	}

	for (auto& i : this->D_Rails)
	{
		//cout << i->s << endl;
		for (auto& j : S->D_Rails)
		{
			if (i->s + 1.0 >= j->s)
			{
				//cout << i->s << " " << j->s << endl;
				i->Key_point[1]->Vx = j->Key_point[1]->x - i->Key_point[1]->x;
				i->Key_point[1]->Vy = j->Key_point[1]->y - i->Key_point[1]->y;
				i->Key_point[2]->Vx = j->Key_point[2]->x - i->Key_point[2]->x;
				i->Key_point[2]->Vy = j->Key_point[2]->y - i->Key_point[2]->y;
				break;
			}
		}
	}

	this->Move_Setka_Calculate_2(1.0);

	for (auto& i : this->All_Points)
	{
		i->x = i->x2;
		i->y = i->y2;
		i->Vx = 0.0;
		i->Vy = 0.0;
		i->count = 0;
	}

	this->Line_Outer[0]->A->Vx = this->Line_Outer[0]->B->Vx + (this->Line_Outer[0]->B->x - this->Line_Outer[0]->A->x);
	this->Line_Outer[this->Line_Outer.size() - 1]->B->Vy = this->Line_Outer[this->Line_Outer.size() - 1]->A->y - this->Line_Outer[this->Line_Outer.size() - 1]->B->y;
	this->Line_Contact[0]->A->Vx = this->Line_Contact[0]->B->Vx + (this->Line_Contact[0]->B->x - this->Line_Contact[0]->A->x);

	this->Line_Contact[this->Line_Contact.size() - 1]->B->Vy = this->Line_Contact[this->Line_Contact.size() - 1]->A->Vy + //
		(this->Line_Contact[this->Line_Contact.size() - 1]->A->y - this->Line_Contact[this->Line_Contact.size() - 1]->B->y);
	this->Move_Setka_Calculate_2(1.0);
	
	for (auto& i : this->All_Points)
	{
		i->x = i->x2;
		i->y = i->y2;
		i->Vx = 0.0;
		i->Vy = 0.0;
		i->count = 0;
	}

	cout << "Filling data" << endl;
	double x, y;
	for (auto& i : this->All_Cells)
	{
		i->Get_Center(x, y);
		for (auto& j : S->All_Cells)
		{
			if (j->belong(x, y) == true)
			{
				i->par[0] = j->par[0];
				i->par[1] = j->par[1];
				break;
			}
		}
	}
}

void Setka::Save_G_D(void)
{
	ofstream fout;
	fout.open("Setka_all.txt");

	for (auto& i : this->All_Cells)
	{
		fout << i->par[0].ro << " " << i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].Q << endl;
	}
}

void Setka::Save_G_D_5_komponent(void)
{
	ofstream fout;
	fout.open("Setka_all_5_component.txt");

	for (auto& i : this->All_Cells)
	{
		fout << i->par[0].ro << " " << i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].Q << " " << //
			i->par[0].ro_H1 << " " << i->par[0].p_H1 << " " << i->par[0].u_H1 << " " << i->par[0].v_H1 << " " <<//
			i->par[0].ro_H2 << " " << i->par[0].p_H2 << " " << i->par[0].u_H2 << " " << i->par[0].v_H2 << " " <<//
			i->par[0].ro_H3 << " " << i->par[0].p_H3 << " " << i->par[0].u_H3 << " " << i->par[0].v_H3 << " " <<//
			i->par[0].ro_H4 << " " << i->par[0].p_H4 << " " << i->par[0].u_H4 << " " << i->par[0].v_H4 << " " <<//
			endl;
	}
}

void Setka::Save_Setka_ALL_ALPHA(string name)
{
	int num = 0;
	for (auto& i : this->All_Gran)  
	{
		i->number = num;
		num++;
	}
	num = 0;
	for (auto& i : this->All_Gran_copy)  
	{
		i->number = num;
		num++;
	}


	ofstream fout;
	fout.open(name);

	fout << this->N1 << " " << this->N2 << " " << this->N3 << " " << this->N4 << " " << this->M1 << " " << this->M2 << " " << this->M3 << " " << this->M4 << " " << endl;
	
	fout << this->All_Points.size() << endl;
	for (auto& i : this->All_Points)
	{
		fout << i->x << " " << i->y << " " << i->type << endl;
	}

	fout << this->All_Gran.size() << endl;
	for (auto& i : this->All_Gran)
	{
		fout << i->A->number << " " << i->B->number << " " << i->type << endl;
	}

	fout << this->All_Gran_copy.size() << endl;
	for (auto& i : this->All_Gran_copy)
	{
		fout << i->A->number << " " << i->B->number << " " << i->type << endl;
	}

	fout << this->A_Rails.size() << endl;
	for (auto& i : this->A_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->B_Rails.size() << endl;
	for (auto& i : this->B_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->C_Rails.size() << endl;
	for (auto& i : this->C_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->D_Rails.size() << endl;
	for (auto& i : this->D_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->All_Cells.size() << endl;
	for (auto& i : this->All_Cells)
	{
		fout << i->type << endl;
		fout << i->contour.size() << endl;
		for (auto& j : i->contour)
		{
			fout << j->number << endl;
		}

		fout << i->Grans.size() << endl;
		for (auto& j : i->Grans)
		{
			fout << j->number << " " << j->main_gran << endl;
		}
	}

	fout << this->All_Gran.size() << endl;
	for (auto& i : this->All_Gran)
	{
		int i1, i2, i3;
		if (i->Sosed != nullptr)
		{
			i1 = i->Sosed->number;
		}
		else
		{
			i1 = -1;
		}

		if (i->Gran_copy != nullptr)
		{
			i2 = i->Gran_copy->number;
		}
		else
		{
			i2 = -1;
		}

		fout << i->Master->number << " " << i1 << " " << i2 << endl;
	}

	fout << this->All_Gran_copy.size() << endl;
	for (auto& i : this->All_Gran_copy)
	{
		int i_1, i_2, i_3;
		if (i->Master != nullptr)
		{
			i_1 = i->Master->number;
		}
		else
		{
			i_1 = -1;
		}

		if (i->Sosed != nullptr)
		{
			i_2 = i->Sosed->number;
		}
		else
		{
			i_2 = -1;
		}

		fout << i_1 << " " << i_2 << " " << i->Gran_copy->number << endl;
	}

	fout << this->All_Points.size() << endl;
	for (auto& i : this->All_Points)
	{
		fout << i->my_cell.size() << endl;
		for (auto& j : i->my_cell)
		{
			fout << j->number << endl;
		}
	}

	fout << this->Line_Contact.size() << endl;
	for (auto& i : this->Line_Contact)
	{
		fout << i->number << endl;
	}

	fout << this->Line_Inner.size() << endl;
	for (auto& i : this->Line_Inner)
	{
		fout << i->number << endl;
	}

	fout << this->Line_Outer.size() << endl;
	for (auto& i : this->Line_Outer)
	{
		fout << i->number << endl;
	}

	fout << this->All_Cells.size() << endl;
	for (auto& i : this->All_Cells)
	{
		fout << i->par[0].ro << " " << i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].Q << endl;
		fout << i->par[0].ro_H1 << " " << i->par[0].p_H1 << " " << i->par[0].u_H1 << " " << i->par[0].v_H1 << endl;
		fout << i->par[0].ro_H2 << " " << i->par[0].p_H2 << " " << i->par[0].u_H2 << " " << i->par[0].v_H2 << endl;
		fout << i->par[0].ro_H3 << " " << i->par[0].p_H3 << " " << i->par[0].u_H3 << " " << i->par[0].v_H3 << endl;
		fout << i->par[0].ro_H4 << " " << i->par[0].p_H4 << " " << i->par[0].u_H4 << " " << i->par[0].v_H4 << endl;
		fout << i->par[0].I_u << " " << i->par[0].I_v << " " << i->par[0].I_T << endl;
		for (int j = 0; j < 4; j++)
		{
			fout << i->par[0].H_n[j] << " " << i->par[0].H_u[j] << " " << i->par[0].H_v[j] << " " << i->par[0].H_T[j] << endl;
		}
		fout << i->par[0].k_u << " " << i->par[0].k_v << " " << i->par[0].k_T << endl;
	}

	fout.close();
}

void Setka::Download_Setka_ALL_ALPHA(string name)
{
	int n, m, k, type;
	double x, y;
	int num;
	int m1, m2, m3, m4;
	double s;
	bool b;
	ifstream fout;
	fout.open(name);

	fout >> this->N1 >> this->N2 >> this->N3 >> this->N4 >> this->M1 >> this->M2 >> this->M3 >> this->M4;

	fout >> n;
	this->All_Points.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> x >> y; 
		auto P = new Point(x, y);
		fout >> type;
		P->type = static_cast<Point_type>(type);
		this->All_Points.push_back(P);
	}

	fout >> n;
	this->All_Gran.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = true;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran.push_back(G);
	}

	fout >> n;
	this->All_Gran_copy.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = false;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran_copy.push_back(G);
	}

	fout >> n;
	this->A_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->A_Rails.push_back(R);
	}

	fout >> n;
	this->B_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->B_Rails.push_back(R);
	}

	fout >> n;
	this->C_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->C_Rails.push_back(R);
	}

	fout >> n;
	this->D_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->D_Rails.push_back(R);
	}

	fout >> n;
	this->All_Cells.reserve(n);
	for (int i = 0; i < n; i++)
	{
		auto C = new Cell();
		fout >> type;
		C->type = static_cast<Cell_type>(type);

		fout >> m;
		C->contour.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			C->contour.push_back(this->All_Points[k]);
		}

		fout >> m;
		C->Grans.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k >> b;
			if (b == true)
			{
				C->Grans.push_back(this->All_Gran[k]);
			}
			else
			{
				C->Grans.push_back(this->All_Gran_copy[k]);
			}
		}

		this->All_Cells.push_back(C);
	}

	fout >> n;
	for (auto& i : this->All_Gran)
	{
		fout >> m1 >> m2 >> m3;
		i->Master = this->All_Cells[m1];
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
			i->Gran_copy = this->All_Gran_copy[m3];
		}
		else
		{
			i->Sosed = nullptr;
			i->Gran_copy = nullptr;
		}
	}

	fout >> n;
	for (auto& i : this->All_Gran_copy)
	{
		fout >> m1 >> m2 >> m3;
		i->Master = this->All_Cells[m1];
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
			i->Gran_copy = this->All_Gran[m3];
		}
		else
		{
			cout << "Syda ne dolgny popadat  hrgrfwgydwfy2e2443" << endl;
			i->Sosed = nullptr;
			i->Gran_copy = nullptr;
		}
	}


	fout >> n;
	for (auto& i : this->All_Points)
	{
		fout >> m;
		i->my_cell.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			i->my_cell.push_back(this->All_Cells[k]);
		}
	}

	fout >> n;
	this->Line_Contact.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Contact.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Inner.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Inner.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Outer.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Outer.push_back(this->All_Gran[k]);
	}

	fout >> n;
	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >> i->par[0].p >> i->par[0].u >> i->par[0].v >> i->par[0].Q;
		fout >> i->par[0].ro_H1 >> i->par[0].p_H1 >> i->par[0].u_H1 >> i->par[0].v_H1;
		fout >> i->par[0].ro_H2 >> i->par[0].p_H2 >> i->par[0].u_H2 >> i->par[0].v_H2;
		fout >> i->par[0].ro_H3 >> i->par[0].p_H3 >> i->par[0].u_H3 >> i->par[0].v_H3;
		fout >> i->par[0].ro_H4 >> i->par[0].p_H4 >> i->par[0].u_H4 >> i->par[0].v_H4;

		i->par[1] = i->par[0];
	}

	// Нумерация всего
	num = 0;
	for (auto& i : this->All_Gran)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Gran_copy)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}


	cout << "Vsego tochek = " << this->All_Points.size() << endl;
	cout << "Vsego Yacheek = " << this->All_Cells.size() << endl;

	fout.close();
}

void Setka::Download_Setka_ALL_ALPHA_2_0(string name)
{
	int n, m, k, type;
	double x, y;
	int num;
	int m1, m2, m3, m4;
	double s;
	bool b;
	ifstream fout;
	fout.open(name);
	cout << "Downloading start" << endl;
	fout >> this->N1 >> this->N2 >> this->N3 >> this->N4 >> this->M1 >> this->M2 >> this->M3 >> this->M4;

	fout >> n;
	this->All_Points.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> x >> y;
		auto P = new Point(x, y);
		fout >> type;
		P->type = static_cast<Point_type>(type);
		this->All_Points.push_back(P);
	}

	fout >> n;
	this->All_Gran.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = true;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran.push_back(G);
	}

	fout >> n;
	this->All_Gran_copy.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = false;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran_copy.push_back(G);
	}

	fout >> n;
	this->A_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->A_Rails.push_back(R);
	}

	fout >> n;
	this->B_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->B_Rails.push_back(R);
	}

	fout >> n;
	this->C_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->C_Rails.push_back(R);
	}

	fout >> n;
	this->D_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->D_Rails.push_back(R);
	}
	cout << "AAA1" << endl;
	fout >> n;
	this->All_Cells.reserve(n);
	for (int i = 0; i < n; i++)
	{
		auto C = new Cell();
		fout >> type;
		C->type = static_cast<Cell_type>(type);

		fout >> m;
		C->contour.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			C->contour.push_back(this->All_Points[k]);
		}

		fout >> m;
		C->Grans.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k >> b;
			if (b == true)
			{
				C->Grans.push_back(this->All_Gran[k]);
			}
			else
			{
				if (this->All_Gran_copy.size() - 1 < k)
				{
					cout << "ERROR  1741 " << endl;
					cout << k << endl;
					exit(-2);
				}
				C->Grans.push_back(this->All_Gran_copy[k]);
			}
		}

		this->All_Cells.push_back(C);
	}
	cout << "AAA2" << endl;
	fout >> n;
	for (auto& i : this->All_Gran)
	{
		fout >> m1 >> m2 >> m3;
		i->Master = this->All_Cells[m1];
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
		}
		else
		{
			i->Sosed = nullptr;
		}

		if (m3 != -1)
		{
			i->Gran_copy = this->All_Gran_copy[m3];
		}
		else
		{
			i->Gran_copy = nullptr;
		}
	}

	fout >> n;
	cout << "AAA3" << endl;
	for (auto& i : this->All_Gran_copy)
	{
		fout >> m1 >> m2 >> m3;
		if (m1 != -1)
		{
			i->Master = this->All_Cells[m1];
		}
		else
		{
			i->Master = nullptr;
		}
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
			i->Gran_copy = this->All_Gran[m3];
		}
		else
		{
			cout << "Syda ne dolgny popadat  hrgrfwgydwfy2e2443" << endl;
			i->Sosed = nullptr;
			i->Gran_copy = nullptr;
		}
	}


	fout >> n;
	for (auto& i : this->All_Points)
	{
		fout >> m;
		i->my_cell.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			i->my_cell.push_back(this->All_Cells[k]);
		}
	}

	fout >> n;
	this->Line_Contact.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Contact.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Inner.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Inner.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Outer.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Outer.push_back(this->All_Gran[k]);
	}

	fout >> n;
	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >> i->par[0].p >> i->par[0].u >> i->par[0].v >> i->par[0].Q;
		fout >> i->par[0].ro_H1 >> i->par[0].p_H1 >> i->par[0].u_H1 >> i->par[0].v_H1;
		fout >> i->par[0].ro_H2 >> i->par[0].p_H2 >> i->par[0].u_H2 >> i->par[0].v_H2;
		fout >> i->par[0].ro_H3 >> i->par[0].p_H3 >> i->par[0].u_H3 >> i->par[0].v_H3;
		fout >> i->par[0].ro_H4 >> i->par[0].p_H4 >> i->par[0].u_H4 >> i->par[0].v_H4;

		fout >> i->par[0].I_u >> i->par[0].I_v >> i->par[0].I_T;

		for (int j = 0; j < 4; j++)
		{
			fout >> i->par[0].H_n[j] >> i->par[0].H_u[j] >> i->par[0].H_v[j] >> i->par[0].H_T[j];
		}

		fout >> i->par[0].k_u >> i->par[0].k_v >> i->par[0].k_T;

		i->par[1] = i->par[0];
	}




	// Нумерация всего
	num = 0;
	for (auto& i : this->All_Gran)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	for (auto& i : this->All_Gran_copy)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}


	cout << "Vsego tochek = " << this->All_Points.size() << endl;
	cout << "Vsego Yacheek = " << this->All_Cells.size() << endl;
	cout << "Downloading end" << endl;
	fout.close();

	for (auto& i : this->All_Cells)
	{
		i->renew();
	}
}

void Setka::Download_G_D_5_komponent(void)
{
	ifstream fout;
	fout.open("Setka_all_5_component.txt");

	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >>i->par[0].p >>i->par[0].u >>i->par[0].v >>i->par[0].Q >>//
			i->par[0].ro_H1 >>i->par[0].p_H1 >>i->par[0].u_H1 >>i->par[0].v_H1 >>//
			i->par[0].ro_H2 >>i->par[0].p_H2 >>i->par[0].u_H2 >>i->par[0].v_H2 >>//
			i->par[0].ro_H3 >>i->par[0].p_H3 >>i->par[0].u_H3 >>i->par[0].v_H3 >>//
			i->par[0].ro_H4 >>i->par[0].p_H4 >>i->par[0].u_H4 >>i->par[0].v_H4;
	}
}

void Setka::Download_G_D(void)
{
	ifstream fout;
	fout.open("Setka_all.txt");

	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >> i->par[0].p >> i->par[0].u >> i->par[0].v >> i->par[0].Q;
		i->par[1] = i->par[0];
	}
}

void Setka::Move_surface_hand(void)
{
	/*for (auto& i : this->Line_Contact)
	{
		i->A->Vx = 0.0;
		i->A->Vy = 0.0;
		i->B->Vx = 0.0;
		i->B->Vy = 0.0;
	}*/
	
	double x, y, d, dist;

	for (auto& i : this->Line_Contact)
	{
		d = 1000;
		x = i->A->x;
		y = i->A->y;
		for (auto& j : this->Contact)
		{
			dist = sqrt(kv(j->x - x) + kv(j->y - y));
			if (dist < d)
			{
				d = dist;
				i->A->Vx = (j->x - x) * 0.1;
				i->A->Vy = (j->y - y) * 0.1;
			}
		}
	}

	this->Line_Contact[0]->A->Vx = this->Line_Contact[0]->B->Vx + (this->Line_Contact[0]->B->x - this->Line_Contact[0]->A->x);

	this->Line_Contact[this->Line_Contact.size() - 1]->B->Vy = this->Line_Contact[this->Line_Contact.size() - 1]->A->Vy + //
		(this->Line_Contact[this->Line_Contact.size() - 1]->A->y - this->Line_Contact[this->Line_Contact.size() - 1]->B->y);

	double VY;

	for (auto& i : this->Line_Outer)
	{
		d = 1000;
		x = i->A->x;
		y = i->A->y;
		if (x > 0)
		{
			for (auto& j : this->Outer)
			{
				dist = sqrt(kv(j->x - x) + kv(j->y - y));
				if (dist < d)
				{
					d = dist;
					i->A->Vx = (j->x - x) * 0.1;
					i->A->Vy = (j->y - y) * 0.1;
					VY = i->A->Vy;
				}
			}
		}
		else
		{
			i->A->Vx = 0.0;
			i->A->Vy = VY;
		}
	}

	this->Line_Outer[0]->A->Vx = this->Line_Outer[0]->B->Vx + (this->Line_Outer[0]->B->x - this->Line_Outer[0]->A->x);

	this->Line_Outer[this->Line_Outer.size() - 1]->B->Vy = this->Line_Outer[this->Line_Outer.size() - 1]->A->Vy + //
		(this->Line_Outer[this->Line_Outer.size() - 1]->A->y - this->Line_Outer[this->Line_Outer.size() - 1]->B->y);

	for (auto& i : this->Line_Inner)
	{
		d = 1000;
		x = i->A->x;
		y = i->A->y;
		for (auto& j : this->Inner)
		{
			dist = sqrt(kv(j->x - x) + kv(j->y - y));
			if (dist < d)
			{
				d = dist;
				i->A->Vx = (j->x - x) * 0.1;
				i->A->Vy = (j->y - y) * 0.1;
			}
		}
	}

	this->Line_Inner[0]->A->Vx = this->Line_Inner[0]->B->Vx + (this->Line_Inner[0]->B->x - this->Line_Inner[0]->A->x);

	this->Line_Inner[this->Line_Inner.size() - 1]->B->Vx = this->Line_Inner[this->Line_Inner.size() - 1]->A->Vx + //
		(this->Line_Inner[this->Line_Inner.size() - 1]->A->x - this->Line_Inner[this->Line_Inner.size() - 1]->B->x);


}

void Setka::Move_surface(int ii, const double& dt = 1.0)
{
	double koef = 1.0; // 1.0;
	// Разбираемся с контактом

	//for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
	//{
	//	auto i = this->Line_Contact[j];
	//	i->A->Vx = 0.3;
	//	i->A->Vy = 0.0;
	//	i->B->Vx = 0.3;
	//	i->B->Vy = 0.0;

	//}
	
	// Контакт
	if (true)
	{
		//int bb = -1;
		for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
		{
			auto S = Solvers();
			auto i = this->Line_Contact[j];
			Parametr par1;// = i->Master->par[ii];
			Parametr par2;// = i->Sosed->par[ii];
			i->Get_par_TVD(par1, ii);
			i->Gran_copy->Get_par_TVD(par2, ii);
			double n1, n2;
			vector<double> qqq1(5);
			vector<double> qqq2(5);
			vector<double> qqq(5);
			vector<double> n(3);
			i->Get_normal(n1, n2);
			double VV, Vl, Vp;
			double P[4];
			double PQ;

			n[0] = n1;
			n[1] = n2;
			n[2] = 0.0;

			qqq1[0] = par1.ro;
			qqq1[1] = par1.p;
			qqq1[2] = par1.u;
			qqq1[3] = par1.v;
			qqq1[4] = 0.0;

			qqq2[0] = par2.ro;
			qqq2[1] = par2.p;
			qqq2[2] = par2.u;
			qqq2[3] = par2.v;
			qqq2[4] = 0.0;

			//this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
			//	par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, 1.0, 1, Vl, VV, Vp);

			S.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV);

			double Max = sqrt((kv(par1.u) + kv(par1.v)) / (ggg * par1.p / par1.ro));
			VV = VV * koef;// *0.05;

			/*if (i->A->x < -400)
			{
				VV = 0.0;
			}*/

			/*if (i->A->x < -200)
			{
				VV = par1.u * n1 + par1.v * n2;
			}*/



			/*if ( (i->A->x < -200 && i->A->y < 200 && VV <= 0)|| (i->B->x < -200 && i->B->y < 200 && VV <= 0) )
			{
				VV = 0.0;
			}*/

			double t1 = -n2;
			double t2 = n1;

			if (false)//(Max < 1)
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->A->Vx += VV * n1;
				i->A->Vy += VV * n2;
				i->B->count++;
				i->A->count++;
			}
			else
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->B->count++;
			}

			/*if (par1.u * t1 + par1.v * t2 > 0)
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->B->count++;
			}
			else if (par1.u * t1 + par1.v * t2 < 0)
			{
				i->A->Vx += VV * n1;
				i->A->Vy += VV * n2;
				i->A->count++;
			}
			else
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->A->Vx += VV * n1;
				i->A->Vy += VV * n2;
				i->B->count++;
				i->A->count++;
			}*/

			/*if (i->A->x < -200 && i->A->y < 400 && i->B->x < -200 && i->B->y < 400)
			{
				VV = 0.05;
			}*/
			//cout << i->A->x << " " << i->A->y << " " << VV << " " << n1 << " " << n2 <<  endl;
			/*if (Max > 1)
			{
				if (bb == -1)
				{
					bb = j;
				}
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
			}
			else
			{
				i->A->Vx += VV * n1 * 0.5;
				i->A->Vy += VV * n2 * 0.5;
				i->B->Vx += VV * n1 * 0.5;
				i->B->Vy += VV * n2 * 0.5;
			}*/
		}

		for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
		{
			auto i = this->Line_Contact[j];
			if (i->A->count > 0)
			{
				i->A->Vx /= (1.0 * i->A->count);
				i->A->Vy /= (1.0 * i->A->count);
			}
		}

		this->Line_Contact[0]->A->Vx = this->Line_Contact[0]->B->Vx + (this->Line_Contact[0]->B->x - this->Line_Contact[0]->A->x)/dt;

		this->Line_Contact[this->Line_Contact.size() - 1]->B->Vy = this->Line_Contact[this->Line_Contact.size() - 1]->A->Vy + //
			(this->Line_Contact[this->Line_Contact.size() - 1]->A->y - this->Line_Contact[this->Line_Contact.size() - 1]->B->y)/dt;

		/*if (bb >= 0)
		{
			this->Line_Contact[bb]->A->Vx *= 2.0;
			this->Line_Contact[bb]->A->Vy *= 2.0;
		}*/
	}

	// Движение внутренней ударной волны
	if (true)
	{
		// Движение внутренней ударной волны
		for (int j = 0; j < this->Line_Inner.size(); j++)
		{
			auto i = this->Line_Inner[j];
			Parametr par1 = i->Master->par[ii];
			double x, y, x2, y2;
			i->Master->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			i->Get_Center(x2, y2);
			double dis = sqrt(kv(x2) + kv(y2));
			Parametr par2;// = i->Sosed->par[ii];
			//par2 = i->Gran_copy->Master->par[ii];
			i->Gran_copy->Get_par_TVD(par2, ii);
			double n1, n2;
			i->Get_normal(n1, n2);
			double VV, Vl, Vp;
			double P[4];
			double PQ;
			par1.ro = par1.ro * kv(radius) / kv(dis);
			par1.p = par1.p * pow(radius / dis, 2.0 * ggg);
			polar_perenos(x, y, x2, y2, par1.u, par1.v);

			auto S = Solvers();
			vector<double> qqq1(5);
			vector<double> qqq2(5);
			vector<double> qqq(5);
			vector<double> n(3);

			n[0] = n1;
			n[1] = n2;
			n[2] = 0.0;

			qqq1[0] = par1.ro;
			qqq1[1] = par1.p;
			qqq1[2] = par1.u;
			qqq1[3] = par1.v;
			qqq1[4] = 0.0;

			qqq2[0] = par2.ro;
			qqq2[1] = par2.p;
			qqq2[2] = par2.u;
			qqq2[3] = par2.v;
			qqq2[4] = 0.0;

			S.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV);

			//this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
			//	par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, 1.0, 1, Vl, VV, Vp);
			Vl = Vl * koef;
			double t1 = -n2;
			double t2 = n1;

			if (fabs(par2.u * t1 + par2.v * t2) < 0.1)
			{
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->A->count++;
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}

			/*if ( fabs(par2.u * t1 + par2.v * t2) < 0.01)
			{
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->A->count++;
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else if (par2.u * t1 + par2.v * t2 > 0)
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else
			{
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->A->count++;
			}*/

			/*if (i->A->x < 0)
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->B->count++;
				i->A->count++;
			}*/

			//if (true)// (j > 12)
			//{
			//	if (par1.u * t1 + par1.v * t2 > 0)
			//	{
			//		i->B->Vx += Vl * n1;
			//		i->B->Vy += Vl * n2;
			//		i->B->count++;
			//	}
			//	else if (par1.u * t1 + par1.v * t2 < 0)
			//	{
			//		i->A->Vx += Vl * n1;
			//		i->A->Vy += Vl * n2;
			//		i->A->count++;
			//	}
			//	else
			//	{
			//		i->B->Vx += Vl * n1;
			//		i->B->Vy += Vl * n2;
			//		i->A->Vx += Vl * n1;
			//		i->A->Vy += Vl * n2;
			//		i->B->count++;
			//		i->A->count++;
			//	}
			//}
			/*else if(j < 12)
			{
				if (par1.u * t1 + par1.v * t2 < 0)
				{
					i->B->Vx += Vl * n1;
					i->B->Vy += Vl * n2;
					i->B->count++;
				}
				else if (par1.u * t1 + par1.v * t2 > 0)
				{
					i->A->Vx += Vl * n1;
					i->A->Vy += Vl * n2;
					i->A->count++;
				}
				else
				{
					i->B->Vx += Vl * n1;
					i->B->Vy += Vl * n2;
					i->A->Vx += Vl * n1;
					i->A->Vy += Vl * n2;
					i->B->count++;
					i->A->count++;
				}
			}
			else
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->B->count++;
				i->A->count++;
			}*/


		}

		for (int j = 0; j < this->Line_Inner.size(); j++)  // Вычисляем скорость контакта
		{
			auto i = this->Line_Contact[j];
			// Можно что-то ещё придумать с узлом где count = 0
			if (i->A->count > 0)
			{
				i->A->Vx /= (1.0 * i->A->count);
				i->A->Vy /= (1.0 * i->A->count);
			}
			/*else if(j > 0 && j < this->Line_Inner.size() - 1)
			{
				i->A->Vx = (i->B->count + this->Line_Contact[j - 1]->A->Vx)/2.0;
			}*/
		}

		this->Line_Inner[0]->A->Vx = this->Line_Inner[0]->B->Vx + (this->Line_Inner[0]->B->x - this->Line_Inner[0]->A->x)/dt;

		this->Line_Inner[this->Line_Inner.size() - 1]->B->Vx = this->Line_Inner[this->Line_Inner.size() - 1]->A->Vx + //
			(this->Line_Inner[this->Line_Inner.size() - 1]->A->x - this->Line_Inner[this->Line_Inner.size() - 1]->B->x)/dt;

	}

	// Двигаем ли внешнюю волну
	if (true)
	{
		//double Vconst = 0.0;
		for (int j = 0; j < this->Line_Outer.size(); j++)
		{
			auto i = this->Line_Outer[j];
			Parametr par1;// = i->Master->par[ii];
			Parametr par2;// = i->Sosed->par[ii];
			i->Get_par_TVD(par1, ii);
			i->Gran_copy->Get_par_TVD(par2, ii);
			double n1, n2;
			i->Get_normal(n1, n2);
			double VV, Vl, Vp;
			double P[4];
			double PQ;

			vector<double> qqq1(5);
			vector<double> qqq2(5);
			vector<double> qqq(5);
			vector<double> n(3);
			n[0] = n1;
			n[1] = n2;
			n[2] = 0.0;

			qqq1[0] = par1.ro;
			qqq1[1] = par1.p;
			qqq1[2] = par1.u;
			qqq1[3] = par1.v;
			qqq1[4] = 0.0;

			qqq2[0] = par2.ro;
			qqq2[1] = par2.p;
			qqq2[2] = par2.u;
			qqq2[3] = par2.v;
			qqq2[4] = 0.0;

			auto S = Solvers();
			S.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV);

			//this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
			//	par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, 1.0, 1, Vl, VV, Vp);
			Vp = Vp * koef;
			//cout << Vp << endl;

			double t1 = -n2;
			double t2 = n1;

			if (kvv(par1.u, par1.v, 0.0) / (ggg * par1.p / par1.ro) > 1.0)
			{
				i->B->Vx += Vp * n1;
				i->B->Vy += Vp * n2;
				i->B->count++;
				/*i->A->Vx += Vp * n1;
				i->A->Vy += Vp * n2;
				i->A->count++;*/
			}
			else
			{
				if (fabs(par1.u * t1 + par1.v * t2) < 0.01)
				{
					i->A->Vx += Vp * n1;
					i->A->Vy += Vp * n2;
					i->A->count++;
					i->B->Vx += Vp * n1;
					i->B->Vy += Vp * n2;
					i->B->count++;
				}
				else if (par1.u * t1 + par1.v * t2 > 0)
				{
					i->B->Vx += Vp * n1;
					i->B->Vy += Vp * n2;
					i->B->count++;
				}
				else
				{
					i->A->Vx += Vp * n1;
					i->A->Vy += Vp * n2;
					i->A->count++;
				}
			}
			
			
		}

		for (int j = 0; j < this->Line_Outer.size(); j++)  // Вычисляем скорость контакта
		{
			auto i = this->Line_Contact[j];
			if (i->A->count > 0)
			{
				i->A->Vx /= (1.0 * i->A->count);
				i->A->Vy /= (1.0 * i->A->count);
			}
		}

		this->Line_Outer[0]->A->Vx = this->Line_Outer[0]->B->Vx + (this->Line_Outer[0]->B->x - this->Line_Outer[0]->A->x)/dt;

		this->Line_Outer[this->Line_Outer.size() - 1]->B->Vy = this->Line_Outer[this->Line_Outer.size() - 1]->A->Vy + //
			(this->Line_Outer[this->Line_Outer.size() - 1]->A->y - this->Line_Outer[this->Line_Outer.size() - 1]->B->y)/dt;
	}
}

void Setka::Move_Setka_Calculate_2(const double& dt)
{
	double Vx, Vy, V;
	double R2, r, R3, R4;

	//double y_Outer = this->D_Rails[16]->Key_point[1]->y;

	for (int jj = 0; jj < this->A_Rails.size(); jj++)
	{
		auto i = this->A_Rails[jj];
		// Подвинем ключевые точки
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

		V = i->Key_point[1]->Vx * cos(i->s) + i->Key_point[1]->Vy * sin(i->s);
		i->Key_point[1]->x2 = i->Key_point[1]->x + dt * V * cos(i->s);
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V * sin(i->s);

		V = i->Key_point[2]->Vx * cos(i->s) + i->Key_point[2]->Vy * sin(i->s);
		i->Key_point[2]->x2 = i->Key_point[2]->x + dt * V * cos(i->s);
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V * sin(i->s);

		i->Key_point[3]->x2 = i->Key_point[3]->x;
		i->Key_point[3]->y2 = i->Key_point[3]->y;


		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;

		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++) 
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++) 
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		//x = pow(R2 / R1_, 1.0 / (i->M1 - 1));
		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ +  pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		//R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2)) - zazor;
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		for (int j = 0; j < i->M2 - 1; j++)  // Передвинули точки до контакта
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);         // Равномерное распределение точек
			i->All_point[i->M1 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + j]->y2 = r * sin(i->s);
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = R3 * cos(i->s);   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 * sin(i->s);
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		//i->All_point[i->M1 + i->M2]->x2 = R3 * cos(i->s);   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 * sin(i->s);

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));


		for (int j = 0; j < i->M3 - 1; j++)// Передвинули точки до внешней волны
		{
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);         // Равномерное распределение точек
			i->All_point[i->M1 + i->M2 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + j]->y2 = r * sin(i->s);
		}

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));

		//x = log(R5_ / R4) / (log(1.0 * i->M4 + 1) * (i->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.03 * x) * (R5_ - R4);
			//r = R4 *  pow( 1.0 * (j + 2), x * (j + 1));  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			//r = R4 + (R5_ - R4) * (j) / (i->M4 - 1);   // равномерно
			x = (j + 1) * 1.0 / (i->M4);
			double dd = 1.5708/sqrt(1.0 - 1.0/1.76);
			r = R4 + pow(x, 1.76 * (1.0 - kv(i->s/dd))) * (R5_ - R4);
			//r = r + (j + 1) * dd;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r * sin(i->s);
		}
	}


	auto i_ = this->A_Rails[this->A_Rails.size() - 1];
	double R_distant_shok_1 = sqrt(kv(i_->Key_point[0]->x2) + kv(i_->Key_point[0]->y2));

	for (int jj = 0; jj < this->B_Rails.size(); jj++)
	{
		auto i = this->B_Rails[jj];
		// Подвинем ключевые точки   0 - внутренняя ударная волна
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = R_distant_shok_1 * cos(i->s);
		i->Key_point[0]->y2 = R_distant_shok_1 * sin(i->s);

		V = i->Key_point[1]->Vx * 0.0 + i->Key_point[1]->Vy * 1.0;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;

		V = i->Key_point[2]->Vx * 0.0 + i->Key_point[2]->Vy * 1.0;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		
		
		R2 = i->Key_point[0]->y2;
		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < i->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[i->M1 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + j]->y2 = r;
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 ;
		R3 = i->Key_point[1]->y2;

		//i->All_point[i->M1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 ;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < i->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[i->M1 + i->M2 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + j]->y2 = r;
		}

		//x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;

		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * (0.03 + 0.19 * jj / (this->B_Rails.size() - 1.0)) * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r;
		}

	}

	for (auto& i : this->C_Rails)
	{
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = R_distant_shok_1 * cos(i->s);
		i->Key_point[0]->y2 = R_distant_shok_1 * sin(i->s);

		i->Key_point[1]->x2 = Left_; // i->Key_point[1]->x;
		i->Key_point[1]->y2 = i->Key_point[1]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		R2 = i->Key_point[0]->x2;

		R3 = i->Key_point[1]->x2;

		x = pow(R3 / (R2), 1.0 / (i->M2));

		for (int j = i->M1; j < i->All_point.size(); j++)
		{
			r = R2  * pow(x, j + 1 - i->M1);
			i->All_point[j]->x2 = r;
			i->All_point[j]->y2 = i->Key_point[0]->y2;
		}
	}

	for (int jj = 0; jj < this->D_Rails.size(); jj++)
	{
		double x;
		// Подвинем ключевые точки
		auto i = this->D_Rails[jj];
		V = i->Key_point[1]->Vy;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		/*if (jj > 16)
		{
			i->Key_point[1]->y2 = y_Outer + (247.0 - y_Outer) * (jj - 16.0) / (this->D_Rails.size() - 1.0 - 16.0);
		}
		else
		{
			i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		}*/
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		/*if (jj > 22)
		{
			i->Key_point[1]->y2 = this->D_Rails[22]->Key_point[1]->y2;
		}*/

		V = i->Key_point[2]->Vy;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;

		if (jj > n_outer_shock)
		{
			i->Key_point[2]->y2 = this->D_Rails[jj - 1]->Key_point[2]->y2;
		}

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = i->Key_point[0]->y2;

		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < this->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[j + 1]->y2 = r;
		}

		//i->All_point[1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[1 + i->M2 - 2]->y2 = R3;
		R3 = i->Key_point[1]->y2;

		//i->All_point[1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[1 + i->M2]->y2 = R3;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < this->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[this->M2 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + j + 1]->y2 = r;
		}

		//double x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < this->M4; j++)
		{
			//double x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.22 * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[this->M2 + this->M3 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + this->M3 + j + 1]->y2 = r;
		}


	}
}

void Setka::Move_Setka_Calculate(const double& dt)
{
	double Vx, Vy, V;
	double R2, r, R3, R4;

	//double y_Outer = this->D_Rails[16]->Key_point[1]->y;

	for (int jj = 0; jj < this->A_Rails.size(); jj++)
	{
		auto i = this->A_Rails[jj];
		// Подвинем ключевые точки
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

		V = i->Key_point[1]->Vx * cos(i->s) + i->Key_point[1]->Vy * sin(i->s);
		i->Key_point[1]->x2 = i->Key_point[1]->x + dt * V * cos(i->s);
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V * sin(i->s);

		V = i->Key_point[2]->Vx * cos(i->s) + i->Key_point[2]->Vy * sin(i->s);
		i->Key_point[2]->x2 = i->Key_point[2]->x + dt * V * cos(i->s);
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V * sin(i->s);

		i->Key_point[3]->x2 = i->Key_point[3]->x;
		i->Key_point[3]->y2 = i->Key_point[3]->y;


		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;

		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++) 
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++) 
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		//x = pow(R2 / R1_, 1.0 / (i->M1 - 1));
		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		//R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2)) - zazor;
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		for (int j = 0; j < i->M2 - 1; j++)  // Передвинули точки до контакта
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);         // Равномерное распределение точек
			i->All_point[i->M1 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + j]->y2 = r * sin(i->s);
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = R3 * cos(i->s);   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 * sin(i->s);
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		//i->All_point[i->M1 + i->M2]->x2 = R3 * cos(i->s);   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 * sin(i->s);

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));


		for (int j = 0; j < i->M3 - 1; j++)// Передвинули точки до внешней волны
		{
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);         // Равномерное распределение точек
			i->All_point[i->M1 + i->M2 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + j]->y2 = r * sin(i->s);
		}

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));

		//x = log(R5_ / R4) / (log(1.0 * i->M4 + 1) * (i->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.03 * x) * (R5_ - R4);
			//r = R4 *  pow( 1.0 * (j + 2), x * (j + 1));  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			//r = R4 + (R5_ - R4) * (j) / (i->M4 - 1);   // равномерно
			x = (j + 1) * 1.0 / (i->M4);
			double dd = 1.5708 / sqrt(1.0 - 1.0 / 1.76);
			r = R4 + pow(x, 1.76 * (1.0 - kv(i->s / dd))) * (R5_ - R4);
			//r = r + (j + 1) * dd;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r * sin(i->s);
		}
	}

	for (int jj = 0; jj < this->B_Rails.size(); jj++)
	{
		auto i = this->B_Rails[jj];
		// Подвинем ключевые точки
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

		V = i->Key_point[1]->Vx * 0.0 + i->Key_point[1]->Vy * 1.0;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;

		V = i->Key_point[2]->Vx * 0.0 + i->Key_point[2]->Vy * 1.0;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		R2 = i->Key_point[0]->y2;
		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < i->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[i->M1 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + j]->y2 = r;
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 ;
		R3 = i->Key_point[1]->y2;

		//i->All_point[i->M1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 ;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < i->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[i->M1 + i->M2 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + j]->y2 = r;
		}

		//x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;

		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * (0.03 + 0.19 * jj / (this->B_Rails.size() - 1.0)) * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r;
		}

	}

	for (auto& i : this->C_Rails)
	{
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

		i->Key_point[1]->x2 = Left_; // i->Key_point[1]->x;
		i->Key_point[1]->y2 = i->Key_point[1]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		R2 = i->Key_point[0]->x2;

		R3 = i->Key_point[1]->x2;

		x = pow(R3 / (R2), 1.0 / (i->M2));

		for (int j = i->M1; j < i->All_point.size(); j++)
		{
			r = R2 * pow(x, j + 1 - i->M1);
			i->All_point[j]->x2 = r;
			i->All_point[j]->y2 = i->Key_point[0]->y2;
		}
	}

	for (int jj = 0; jj < this->D_Rails.size(); jj++)
	{
		double x;
		// Подвинем ключевые точки
		auto i = this->D_Rails[jj];
		V = i->Key_point[1]->Vy;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		/*if (jj > 16)
		{
			i->Key_point[1]->y2 = y_Outer + (247.0 - y_Outer) * (jj - 16.0) / (this->D_Rails.size() - 1.0 - 16.0);
		}
		else
		{
			i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		}*/
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		/*if (jj > 22)
		{
			i->Key_point[1]->y2 = this->D_Rails[22]->Key_point[1]->y2;
		}*/

		V = i->Key_point[2]->Vy;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = i->Key_point[0]->y2;

		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < this->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[j + 1]->y2 = r;
		}

		//i->All_point[1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[1 + i->M2 - 2]->y2 = R3;
		R3 = i->Key_point[1]->y2;

		//i->All_point[1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[1 + i->M2]->y2 = R3;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < this->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[this->M2 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + j + 1]->y2 = r;
		}

		//double x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < this->M4; j++)
		{
			//double x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.22 * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[this->M2 + this->M3 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + this->M3 + j + 1]->y2 = r;
		}


	}
}

void Setka::Init_conditions(void)
{
	double x, y, dist;
	double ro = (389.988 * 389.988) / (chi_ * chi_);
	double P_E = ro * chi_ * chi_ / (ggg * 10.0 * 10.0);

	for (auto& i : this->All_Cells)
	{
		i->Get_Center(x, y);
		dist = sqrt(kv(x) + kv(y));
		double T_p = (P_E * pow(1.0 / dist, 2.0 * ggg)) / (2.0 * ro / (dist * dist));
		//if (dist < 1)//(i->type == C_centr)
		//{
		//	i->par[1] = i->par[0] = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x / dist, chi_ * y / dist, ro / (dist * dist),//
		//	0.001666666, (0.001666666 * chi_real * chi_real / (ggg * 50.0 * 50.0))* pow(1.0 / dist, 2.0 * ggg), chi_real* x / dist, chi_real* y / dist, 0.000001, 0.000001, 0.0, 0.0, 0.000001, 0.000001, 0.0, 0.0,//
		//	0.000001, 0.000001, 0.0, 0.0};
		//}
		//if (dist <= 110)
		//{
		//	i->par[0] = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x / dist, chi_ * y / dist, ro / (dist * dist),//
		//(2.0 * (sigma(chi_real) * 389.988) / (2.0 * Kn_ * chi_real))* (1.0 / dist - 0.1 / kv(dist)), 0.000001, chi_real * x / dist, chi_real * y / dist, 0.000001, 0.000001, 0.0, 0.0, 0.000001, 0.000001, 0.0, 0.0,//
		//	0.000001, 0.000001, 0.0, 0.0};
		//	i->par[1] = i->par[0];
		//}
		//else
		//{
		i->par[0] = { 1.0, 1.0, Velosity_inf, 0.0, 100.0,//
		0.000001, 0.000001, 0.0, 0.0, 0.000001, 0.000001, 0.0, 0.0, 0.000001, 0.000001, 0.0, 0.0,//
		0.000001, 0.000001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		i->par[1] = i->par[0];
		//}
	}
}

void Setka::Go_stationary(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st ++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double rad = sqrt(kv(x) + kv(y));

			if (K->type == C_centr)
			{
				continue;
			}

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11; 
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				int met = 1;
				if (rad < 95 && i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(rad) / (kv(x2) + kv(y2));
					par2.ro = par2.ro * kv(dd) / (kv(x2) + kv(y2));
					par11.p = par11.p * pow(rad / sqrt(kv(x2) + kv(y2)), 2 * ggg);
					par2.p = par2.p * pow(dd / sqrt(kv(x2) + kv(y2)), 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(rad) / (kv(x2) + kv(y2));
					par11.p = par11.p * pow(rad / sqrt(kv(x2) + kv(y2)), 2 * ggg);
				}
				else if (rad < 95 && i->type == Axis)
				{
					//par11.ro = par11.ro * kv(rad) / (kv(x2) + kv(y2));
					//par11.p = par11.p * pow(rad / sqrt(kv(x2) + kv(y2)), 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par11.v = 0.0;
					par2 = par11;
					//par2.v = -par2.v;
				}
				else if (i->type == Axis)
				{
					par11.v = 0.0;
					par2 = par11;
				}
				/*if (num_cell == 1)
				{
					cout << "-----" << endl;
					cout << n1 << " " << n2 << endl;
					cout << par11.ro << " " << par11.u << " " << par11.v << endl;
					cout << par2.ro << " " << par2.u << " " << par2.v << endl;
				}*/
				/*if (step > 20000)
				{
					met = 1;
				}*/
				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
							par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok[k] = K->Potok[k] + P[k] * S;
				}
				//cout << "B5" << endl;
				K->Potok[4] = K->Potok[4] + PQ * S;
				//cout << "B6" << endl;
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			double ro3, p3, u3, v3, Q33;
			double Volume = K->Get_Volume();

			ro3 = par1.ro - T[now1] * (K->Potok[0] / Volume + par1.ro * par1.v / y);
			Q33 = par1.Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * par1.Q * par1.v / y;
			if (ro3 <= 0)
			{
				printf("Problemsssss  par1.ro < 0!   %lf,   %lf,   %lf\n", x, y, par1.ro);
				ro3 = 0.00001;
			}
			u3 = (par1.ro * par1.u - T[now1] * (K->Potok[1] / Volume + par1.ro * par1.v * par1.u / y)) / ro3;
			v3 = (par1.ro * par1.v - T[now1] * (K->Potok[2] / Volume + par1.ro * par1.v * par1.v / y)) / ro3;
			p3 = (((par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) - T[now1] * (K->Potok[3] / Volume +  //
				+ par1.v * (ggg * par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) / y)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;
		}
	}
}

void Setka::Go_stationary_TVD(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			if (K->type == C_centr)
			{
				continue;
			}

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2;
				Parametr par3;
				if (i->type == Usualy)
				{
					i->Get_par_TVD(par3, now1);
					i->Gran_copy->Get_par_TVD(par2, now1);
				}
				else
				{
					par3 = par1;
					i->Get_par(par2, now1);
				}
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				int met = 0;
				/*if (step > 20000)
				{
					met = 1;
				}*/
				/*if (num_cell == 1)
				{
					cout << par3.ro << " " << par2.ro << " " << x2 << " " << y2 << endl;
				}*/
				time = min(time, this->HLLC_2d_Korolkov_b_s(par3.ro, par3.Q, par3.p, par3.u, par3.v, par2.ro, par2.Q, //
					par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok[k] = K->Potok[k] + P[k] * S;
				}
				//cout << "B5" << endl;
				K->Potok[4] = K->Potok[4] + PQ * S;
				//cout << "B6" << endl;
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			double ro3, p3, u3, v3, Q33;
			double Volume = K->Get_Volume();

			ro3 = par1.ro - T[now1] * (K->Potok[0] / Volume + par1.ro * par1.v / y);
			Q33 = par1.Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * par1.Q * par1.v / y;
			if (ro3 <= 0)
			{
				printf("Problemsssss  par1.ro < 0!   %lf,   %lf,   %lf\n", x, y, par1.ro);
				ro3 = 0.00001;
			}
			u3 = (par1.ro * par1.u - T[now1] * (K->Potok[1] / Volume + par1.ro * par1.v * par1.u / y)) / ro3;
			v3 = (par1.ro * par1.v - T[now1] * (K->Potok[2] / Volume + par1.ro * par1.v * par1.v / y)) / ro3;
			p3 = (((par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) - T[now1] * (K->Potok[3] / Volume +  //
				+par1.v * (ggg * par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) / y)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;
		}
	}
}

void Setka::Go_stationary_5_komponent(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min, y_min;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));
				if (radius < 85 && i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (radius < 85 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}
				

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;
			if (par1.Q / par1.ro < 90)
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}
			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
			q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}
			
			if (par1.Q / par1.ro < 90)
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}

			int k1 = 0, k2 = 0, k3 = 0, k4 = 0;

			if (kappa < 90)
			{
				k3 = 0;
				k4 = 0;
				if (((kv(u) + kv(v)) / (ggg * p / ro) > 1.0) && ((radius <= 400.0)))
				{
					k1 = 1;
					k2 = 0;
				}
				else
				{
					k1 = 0;
					k2 = 1;
				}
			}
			else
			{
				k1 = 0;
				k2 = 0;
				if (p / ro > 1.8)   // ( ((y < 8.0)&&(p / ro > (y * (-0.0238) + 0.36)))||( (y >= 8.0)&&(p / ro > 0.17) ) )
				{
					k4 = 0;
					k3 = 1;
				}
				else
				{
					k3 = 0;
					k4 = 1;
				}
			}

			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));


			if (K->type != C_centr)
			{
				ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume + ro_H1 * v_H1 / y - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!\n");
					ro3 = 0.0000001;
				}
				u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 - T[now1] * (K->Potok_H1[2] / Volume + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume + //
					+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume + ro_H2 * v_H2 / y - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 - T[now1] * (K->Potok_H2[2] / Volume + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume + //
				+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;

			
			ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume + ro_H3 * v_H3 / y - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 - T[now1] * (K->Potok_H3[2] / Volume + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume + //
				+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;

			
			ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume + ro_H4 * v_H4 / y - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 - T[now1] * (K->Potok_H4[2] / Volume + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume + //
				+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;
			

		}
	}
}

void Setka::Go_stationary_5_komponent_inner(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 10000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			if (K->type == C_centr)
			{
				continue;
			}
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			//cout << "Do = " << y << endl;
			//y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * K->Get_Volume());
			//y = K->Get_Volume_rotate(5.0) / (2.0 * pi_ * K->Get_Volume() / 72.0);
			//cout << "Posle = " << y << endl;

			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
					//par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					//par2.p_H1 = par2.p_H1 * pow(dd / dis, 2 * ggg);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par11.p_H1 = pp1.p_H1;
						par2.p_H1 = pp2.p_H1;
						par11.ro_H1 = pp1.ro_H1;
						par2.ro_H1 = pp2.ro_H1;

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					/*par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);*/

					//cout << "inner " << dis << endl;
					/*par2.ro_H3 = par1.ro_H3;
					par2.p_H3 = par1.p_H3;
					par2.u_H3 = par1.u_H3;
					par2.v_H3 = par1.v_H3;
					par2.ro_H4 = par1.ro_H4;
					par2.p_H4 = par1.p_H4;
					par2.u_H4 = par1.u_H4;
					par2.v_H4 = par1.v_H4;*/
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
					par11.ro_H1 = par2.ro_H1;
					par11.p_H1 = par2.p_H1;
					par11.u_H1 = par2.u_H1;
					par11.v_H1 = par2.v_H1;
				}
				else if (i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}


				if (K->type != C_centr)
				{
					/*double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, 0.0, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;*/


					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);

			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
			q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));// / (chi_real / chi_);


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);
			

			int k1 = 1, k2 = 0, k3 = 0, k4 = 0;


			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));

			if (K->type == C_centr)
			{
				q1_H1 = 0.0;
				q1_H2 = 0.0;
				q1_H3 = 0.0;
				q1_H4 = 0.0;

				q21_H1 = 0.0;
				q21_H2 = 0.0;
				q21_H3 = 0.0;
				q21_H4 = 0.0;

				q22_H1 = 0.0;
				q22_H2 = 0.0;
				q22_H3 = 0.0;
				q22_H4 = 0.0;

				q3_H1 = 0.0;
				q3_H2 = 0.0;
				q3_H3 = 0.0;
				q3_H4 = 0.0;
			}

			//cout << x << " " << y << " " << q21_H1 << " " << q22_H1 << " " << q3_H1 << " " << S2 << " " << nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)) << endl;
			//cout << x << " " << y << " " << par1.ro_H3 << " " << par1.ro_H4 << endl;
			//cout << nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro)) << " " << //
			//	nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro)) << " " << //
			//	nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro)) << " " << //
			//	nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro)) << endl;
			//cout << nu_H1 << " " << nu_H2 << " " << nu_H3 << " " << nu_H4 << endl;
			//cout << U_H1 / U_M_H1 << " " << U_H2 / U_M_H2 << " " << U_H3 / U_M_H3 << " " << U_H4 / U_M_H4 << endl;
			/*if (radius < 100)
			{
				 q3_H1 = q22_H1 = q21_H1 = q1_H1 = 0.0;
			}*/


			if (K->type != C_centr)
			{
				ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume + ro_H1 * v_H1 / y - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!  %lf,  %lf\n", x, y);
					ro3 = 0.0003;
				}
				u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 - T[now1] * (K->Potok_H1[2] / Volume + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume + //
					+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				//if (num_cell == 0)
				//{
				//	cout << ro3 << " " << u3 << " " << v3 << " " << q1_H1 << endl;
				//	//ro3 = 0.0003;
				//}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume + ro_H2 * v_H2 / y - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 - T[now1] * (K->Potok_H2[2] / Volume + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume + //
				+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume + ro_H3 * v_H3 / y - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 - T[now1] * (K->Potok_H3[2] / Volume + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume + //
				+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume + ro_H4 * v_H4 / y - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 - T[now1] * (K->Potok_H4[2] / Volume + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume + //
				+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;


		}
	}
}

void Setka::Go_stationary_5_komponent_inner_2(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 10000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			if (K->type == C_centr)
			{
				continue;
			}
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			//cout << "Do = " << y << endl;
			//y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * K->Get_Volume());
			//y = K->Get_Volume_rotate(5.0) / (2.0 * pi_ * K->Get_Volume() / 72.0);
			//cout << "Posle = " << y << endl;

			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume_rotate(alpha_rot);

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square_rotate(alpha_rot);
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
					//par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					//par2.p_H1 = par2.p_H1 * pow(dd / dis, 2 * ggg);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par11.p_H1 = pp1.p_H1;
						par2.p_H1 = pp2.p_H1;
						par11.ro_H1 = pp1.ro_H1;
						par2.ro_H1 = pp2.ro_H1;

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					/*par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);*/

					//cout << "inner " << dis << endl;
					/*par2.ro_H3 = par1.ro_H3;
					par2.p_H3 = par1.p_H3;
					par2.u_H3 = par1.u_H3;
					par2.v_H3 = par1.v_H3;
					par2.ro_H4 = par1.ro_H4;
					par2.p_H4 = par1.p_H4;
					par2.u_H4 = par1.u_H4;
					par2.v_H4 = par1.v_H4;*/
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
					par11.ro_H1 = par2.ro_H1;
					par11.p_H1 = par2.p_H1;
					par11.u_H1 = par2.u_H1;
					par11.v_H1 = par2.v_H1;
				}
				else if (i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}


				if (K->type != C_centr)
				{
					/*double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, 0.0, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;*/


					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);

			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
			q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));// / (chi_real / chi_);


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				/*ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;*/
				ro3 = ro - T[now1] * (K->Potok[0] / Volume);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1];
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				//u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				//v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
				//	+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * ((K->Potok[2] - 2.0 * p * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume  - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);


			int k1 = 1, k2 = 0, k3 = 0, k4 = 0;


			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));

			if (K->type == C_centr)
			{
				q1_H1 = 0.0;
				q1_H2 = 0.0;
				q1_H3 = 0.0;
				q1_H4 = 0.0;

				q21_H1 = 0.0;
				q21_H2 = 0.0;
				q21_H3 = 0.0;
				q21_H4 = 0.0;

				q22_H1 = 0.0;
				q22_H2 = 0.0;
				q22_H3 = 0.0;
				q22_H4 = 0.0;

				q3_H1 = 0.0;
				q3_H2 = 0.0;
				q3_H3 = 0.0;
				q3_H4 = 0.0;
			}

			//cout << x << " " << y << " " << q21_H1 << " " << q22_H1 << " " << q3_H1 << " " << S2 << " " << nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)) << endl;
			//cout << x << " " << y << " " << par1.ro_H3 << " " << par1.ro_H4 << endl;
			//cout << nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro)) << " " << //
			//	nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro)) << " " << //
			//	nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro)) << " " << //
			//	nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro)) << endl;
			//cout << nu_H1 << " " << nu_H2 << " " << nu_H3 << " " << nu_H4 << endl;
			//cout << U_H1 / U_M_H1 << " " << U_H2 / U_M_H2 << " " << U_H3 / U_M_H3 << " " << U_H4 / U_M_H4 << endl;
			/*if (radius < 100)
			{
				 q3_H1 = q22_H1 = q21_H1 = q1_H1 = 0.0;
			}*/


			if (K->type != C_centr)
			{
				//ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume + ro_H1 * v_H1 / y - q1_H1);
				ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!  %lf,  %lf\n", x, y);
					ro3 = 0.0003;
				}
				//u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 - T[now1] * (K->Potok_H1[2] / Volume + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				//p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume + //
				//	+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 - T[now1] * ((K->Potok_H1[2] - 2.0 * p_H1 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				//if (num_cell == 0)
				//{
				//	cout << ro3 << " " << u3 << " " << v3 << " " << q1_H1 << endl;
				//	//ro3 = 0.0003;
				//}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			//ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume + ro_H2 * v_H2 / y - q1_H2);
			ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 - T[now1] * (K->Potok_H2[2] / Volume + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			//p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume + //
			//	+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 - T[now1] * ((K->Potok_H2[2] - 2.0 * p_H2 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			//ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume + ro_H3 * v_H3 / y - q1_H3);
			ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 - T[now1] * (K->Potok_H3[2] / Volume + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			//p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume + //
			//	+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 - T[now1] * ((K->Potok_H3[2] - 2.0 * p_H3 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			//ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume + ro_H4 * v_H4 / y - q1_H4);
			ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 - T[now1] * (K->Potok_H4[2] / Volume + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			//p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume + //
			//	+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 - T[now1] * ((K->Potok_H4[2] - 2.0 * p_H4 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;


		}
	}
}


void Setka::Go_stationary_5_komponent_inner_MK(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
				}
				else if (i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par2 = par11;
					par2.v = -par2.v;
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}

			/////////////////////////////////////////////////////////////////////////
			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;


			double u_H1 = K->par[0].H_u[0], v_H1 = K->par[0].H_v[0], ro_H1 = K->par[0].H_n[0], p_H1 = 0.5 * K->par[0].H_T[0] * K->par[0].H_n[0];
			double u_H2 = K->par[0].H_u[1], v_H2 = K->par[0].H_v[1], ro_H2 = K->par[0].H_n[1], p_H2 = 0.5 * K->par[0].H_T[1] * K->par[0].H_n[1];
			double u_H3 = K->par[0].H_u[2], v_H3 = K->par[0].H_v[2], ro_H3 = K->par[0].H_n[2], p_H3 = 0.5 * K->par[0].H_T[2] * K->par[0].H_n[2];
			double u_H4 = K->par[0].H_u[3], v_H4 = K->par[0].H_v[3], ro_H4 = K->par[0].H_n[3], p_H4 = 0.5 * K->par[0].H_T[3] * K->par[0].H_n[3];
			/*u = K->par[0].u;
			v = K->par[0].v;
			ro = K->par[0].ro;
			p = K->par[0].p;*/

			/*u = par1.u * (chi_real / chi_);
			v = par1.v * (chi_real / chi_);
			ro = par1.ro / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q / kv(chi_real / chi_);*/
			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;


			if (ro <= 0.0000001)
			{
				ro = 0.0000001;
				p = 0.0;
			}
			if (ro_H1 <= 0.0000001)
			{
				ro_H1 = 0.0000001;
				p_H1 = 0.0;
			}
			if (ro_H2 <= 0.0000001)
			{
				ro_H2 = 0.0000001;
				p_H2 = 0.0;
			}
			if (ro_H3 <= 0.0000001)
			{
				ro_H3 = 0.0000001;
				p_H3 = 0.0;
			}
			if (ro_H4 <= 0.0000001)
			{
				ro_H4 = 0.0000001;
				p_H4 = 0.0;
			}

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			K->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			K->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));


			//K->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
			//	(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
			//	nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
			//		(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
			//	nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
			//		(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
			//	nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
			//		(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);
			K->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));

			q2_1 = K->par[0].M_u * K->par[0].k_u;
			q2_2 = K->par[0].M_v * K->par[0].k_v;
			q3 = K->par[0].M_T * K->par[0].k_T;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;
			
			double ro3, Q33, u3, v3, p3;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}


		}
	}
}

void Setka::Go(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double xx, yy;

	vector <double> q1(5);
	vector <double> q2(5);
	vector <double> q(5);
	vector <double> n(3);
	double dsl, dsp;

	Solvers Sol;

	for (int st = 0; st < step; st++)
	{
		if (st % 10000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << " " << xx << " " << yy << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//this->Move_surface(now1);
		//this->Move_Setka_Calculate(T[now1]);

#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			//cout << "A1" << endl;
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			double P[4];
			double PQ;
			//cout << "A2" << endl;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double W;
			//cout << "A3" << endl;

			if (sqrt(kv(x) + kv(y)) <= 70.0)
			{
				continue;
			}

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double x3, y3;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x3, y3);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				W = ((x3 - x2) * n1 + (y3 - y2) * n2) / T[now1];
				double Vc, Vl, Vp;


				if (i->A->type == P_Contact && i->B->type == P_Contact)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, 1, Vl, Vc, Vp));
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, 0, Vl, Vc, Vp));
				}

				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok[k] = K->Potok[k] + P[k] * S;
				}
				//cout << "B5" << endl;
				K->Potok[4] = K->Potok[4] + PQ * S;
				//cout << "B6" << endl;
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				xx = x;
				yy = y;
				mut.unlock();
			}

			double ro3, p3, u3, v3, Q33;
			double Volume = K->Get_Volume();
			double Volume2 = K->Get_Volume_posle();

			ro3 = par1.ro * Volume/Volume2 - T[now1] * (K->Potok[0] / Volume2 + par1.ro * par1.v / y);
			Q33 = par1.Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * par1.Q * par1.v / y;
			if (ro3 <= 0)
			{
				printf("Problemsssss  par1.ro < 0!   %lf,  %lf\n", x, y);
				ro3 = 0.00001;
			}
			u3 = (par1.ro * par1.u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + par1.ro * par1.v * par1.u / y)) / ro3;
			v3 = (par1.ro * par1.v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + par1.ro * par1.v * par1.v / y)) / ro3;
			p3 = (((par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 +  //
				+par1.v * (ggg * par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) / y)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;
		}

		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
		}
	}
}

void Setka::Go_5_komponent(int step)
{
	//cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate(T[now1]);
		
		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			//y = K->Get_Volume_rotate(5.0) / (2.0 * pi_ * K->Get_Volume() / 72.0);

			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume(); // K->Get_Volume_rotate(alpha_rot);
			//cout << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle(); // K->Get_Volume_posle_rotate(alpha_rot);
			//cout << "Vol2 = " << Volume2 << endl;
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				//double S = i->Get_square_rotate(alpha_rot); 
				double S = i->Get_square();
				//cout << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;  // 1
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);
						polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						par11.p_H1 = pp1.p_H1;
						par11.ro_H1 = pp1.ro_H1;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par2.p_H1 = pp2.p_H1;
						par2.ro_H1 = pp2.ro_H1;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}
				

				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = false;
				bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}
				
				
				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}

				//np = false;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}

				if (god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;
			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real/chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}
		
			
			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);
			}
			else
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
			}


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			
			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				//ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
				ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);  // В цилиндрических

				//Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
				Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0.0)
				{
					printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
					ro3 = 0.00001;
				}
				//u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2  - q2_1)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * ( (K->Potok[2] - 2.0 * p * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}


			int k1 = 0, k2 = 0, k3 = 0, k4 = 0;

			if (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				k3 = 0;
				k4 = 0;
				if (K->type == C_1 || K->type == C_centr)//(((kv(u) + kv(v)) / (ggg * p / ro) > 1.0) && ((radius <= 250.0)))
				{
					k1 = 1;
					k2 = 0;
				}
				else
				{
					k1 = 0;
					k2 = 1;
				}
			}
			else
			{
				k1 = 0;
				k2 = 0;
				if (K->type == C_4)//(p / ro > 1.8)   // ( ((y < 8.0)&&(p / ro > (y * (-0.0238) + 0.36)))||( (y >= 8.0)&&(p / ro > 0.17) ) )
				{
					k4 = 0;
					k3 = 1;
				}
				else
				{
					k3 = 0;
					k4 = 1;
				}
			}

			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));


			if (K->type != C_centr)
			{
				ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 + ro_H1 * v_H1 / y - q1_H1);
				//ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.0000001;
				}
				//u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2  - q21_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * ((K->Potok_H1[2] - 2.0 * p_H1 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H1)) / ro3;
				//p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2  - q3_H1)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2 + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[2] / Volume2 + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2 + //
					+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;

				/*if (x > 600.0)
				{
					K->par[now2].ro_H1 = 1.0;
					K->par[now2].p_H1 = 1.0;
					K->par[now2].u_H1 = Velosity_inf;
					K->par[now2].v_H1 = 0.0;
				}*/
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 + ro_H2 * v_H2 / y - q1_H2);
			//ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 - q21_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * ((K->Potok_H2[2] - 2.0 * p_H2 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2 - q22_H2)) / ro3;
			//p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2  - q3_H2)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[2] / Volume2 + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2 + //
				+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2 + ro_H3 * v_H3 / y - q1_H3);
			//ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2  - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2  - q21_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * ((K->Potok_H3[2] - 2.0 * p_H3 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H3)) / ro3;
			//p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2  - q3_H3)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2 + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[2] / Volume2 + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2 + //
				+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2 + ro_H4 * v_H4 / y - q1_H4);
			//ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2  - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				cout << ro_H1 << " " << ro_H2 << " " << ro_H3 << " " << ro_H4 << endl;
				cout << p_H1 << " " << p_H2 << " " << p_H3 << " " << p_H4 << endl;
				cout << ro << " " << p << " " << u << " " << v << endl;
				ro3 = 0.000001;
			}
			//u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 - q21_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * ((K->Potok_H4[2] - 2.0 * p_H4 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H4)) / ro3;
			//p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 - q3_H4)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[2] / Volume2 + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 + //
				+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;

			/*if (x > 600.0)
			{
				K->par[now2].ro_H4 = 1.0;
				K->par[now2].p_H4 = 0.5;
				K->par[now2].u_H4 = Velosity_inf;
				K->par[now2].v_H4 = 0.0;
			}*/


		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}


void Setka::Go_5_komponent_2(int step)
// Идея в том, чтобы использовать расчёты в угле (3д задача, сведённая к 2д)
{
	//cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate(T[now1]);

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PP[8];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			//y = K->Get_Volume_rotate(5.0) / (2.0 * pi_ * K->Get_Volume() / 72.0);

			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume =  K->Get_Volume_rotate(alpha_rot);
			//cout << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle_rotate(alpha_rot);
			//cout << "Vol2 = " << Volume2 << endl;
			double Scw = K->Get_Volume();
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				double S = i->Get_square_rotate(alpha_rot); 
				//double S = i->Get_square();
				//cout << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;  // 1
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);
						polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						par11.p_H1 = pp1.p_H1;
						par11.ro_H1 = pp1.ro_H1;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par2.p_H1 = pp2.p_H1;
						par2.ro_H1 = pp2.ro_H1;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}


				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = false;
				bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}


				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}

				//np = false;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}

				if (god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}

			}


			// Блок учёта боковых граней (которые искусственные из-за трёхмерности)

			if (false)
			{
				//cout << K->Potok_H1[0] << " " << K->Potok_H1[1] << " " << K->Potok_H1[2] - 2.0 * par1.p_H1 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H1[3] << endl;
				Solvers SS = Solvers();
				//-0.0436194 0.999048
				double uu = par1.v_H1, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H1[0] = K->Potok_H1[0] + P[0] * Scw;
				//K->Potok_H1[1] = K->Potok_H1[1] + P[1] * Scw;
				K->Potok_H1[2] = K->Potok_H1[2] + P[2] * Scw;
				//K->Potok_H1[3] = K->Potok_H1[3] + P[4] * Scw;

				uu = par1.v_H1, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H1[0] = K->Potok_H1[0] + P[0] * Scw;
				//K->Potok_H1[1] = K->Potok_H1[1] + P[1] * Scw;
				K->Potok_H1[2] = K->Potok_H1[2] + P[2] * Scw;
				//K->Potok_H1[3] = K->Potok_H1[3] + P[4] * Scw;

				//cout << K->Potok_H1[0] << " " << K->Potok_H1[1] << " " << K->Potok_H1[2] << " " << K->Potok_H1[3] << endl;
				//cout << K->Potok_H2[0] << " " << K->Potok_H2[1] << " " << K->Potok_H2[2] - 2.0 * par1.p_H2 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H2[3] << endl;
				uu = par1.v_H2, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H2[0] = K->Potok_H2[0] + P[0] * Scw;
				//K->Potok_H2[1] = K->Potok_H2[1] + P[1] * Scw;
				K->Potok_H2[2] = K->Potok_H2[2] + P[2] * Scw;
				//K->Potok_H2[3] = K->Potok_H2[3] + P[4] * Scw;

				uu = par1.v_H2, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H2[0] = K->Potok_H2[0] + P[0] * Scw;
				//K->Potok_H2[1] = K->Potok_H2[1] + P[1] * Scw;
				K->Potok_H2[2] = K->Potok_H2[2] + P[2] * Scw;
				//K->Potok_H2[3] = K->Potok_H2[3] + P[4] * Scw;
				//cout << K->Potok_H2[0] << " " << K->Potok_H2[1] << " " << K->Potok_H2[2] << " " << K->Potok_H2[3] << endl;
				//cout << K->Potok_H3[0] << " " << K->Potok_H3[1] << " " << K->Potok_H3[2] - 2.0 * par1.p_H3 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H3[3] << endl;
				uu = par1.v_H3, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H3[0] = K->Potok_H3[0] + P[0] * Scw;
				//K->Potok_H3[1] = K->Potok_H3[1] + P[1] * Scw;
				K->Potok_H3[2] = K->Potok_H3[2] + P[2] * Scw;
				//K->Potok_H3[3] = K->Potok_H3[3] + P[4] * Scw;

				uu = par1.v_H3, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H3[0] = K->Potok_H3[0] + P[0] * Scw;
				//K->Potok_H3[1] = K->Potok_H3[1] + P[1] * Scw;
				K->Potok_H3[2] = K->Potok_H3[2] + P[2] * Scw;
				//K->Potok_H3[3] = K->Potok_H3[3] + P[4] * Scw;
				//cout << K->Potok_H3[0] << " " << K->Potok_H3[1] << " " << K->Potok_H3[2] << " " << K->Potok_H3[3] << endl;
				//cout << K->Potok_H4[0] << " " << K->Potok_H4[1] << " " << K->Potok_H4[2] - 2.0 * par1.p_H4 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H4[3] << endl;
				uu = par1.v_H4, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H4[0] = K->Potok_H4[0] + P[0] * Scw;
				//K->Potok_H4[1] = K->Potok_H4[1] + P[1] * Scw;
				K->Potok_H4[2] = K->Potok_H4[2] + P[2] * Scw;
				//K->Potok_H4[3] = K->Potok_H4[3] + P[4] * Scw;

				uu = par1.v_H4, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H4[0] = K->Potok_H4[0] + P[0] * Scw;
				//K->Potok_H4[1] = K->Potok_H4[1] + P[1] * Scw;
				K->Potok_H4[2] = K->Potok_H4[2] + P[2] * Scw;
				//K->Potok_H4[3] = K->Potok_H4[3] + P[4] * Scw;

				//cout << K->Potok_H4[0] << " " << K->Potok_H4[1] << " " << K->Potok_H4[2] << " " << K->Potok_H4[3] << endl;
				//exit(-1);

				//cout << K->Potok[0] << " " << K->Potok[1] << " " << K->Potok[2] - 2.0 * par1.p * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok[3] << " " << K->Potok[4] << endl;
				uu = par1.v, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro, par1.p, par1.u, uu, vv, //
					par1.ro, par1.p, par1.u, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok[0] = K->Potok[0] + PP[0] * Scw;
				//K->Potok[1] = K->Potok[1] + PP[1] * Scw;
				K->Potok[2] = K->Potok[2] + P[2] * Scw;
				//K->Potok[3] = K->Potok[3] + PP[4] * Scw;
				//K->Potok[4] = K->Potok[4] + PQ * Scw;

				uu = par1.v, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro, par1.p, par1.u, uu, vv, //
					par1.ro, par1.p, par1.u, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok[0] = K->Potok[0] + PP[0] * Scw;
				//K->Potok[1] = K->Potok[1] + PP[1] * Scw;
				K->Potok[2] = K->Potok[2] + P[2] * Scw;
				//K->Potok[3] = K->Potok[3] + PP[4] * Scw;
				//K->Potok[4] = K->Potok[4] + PQ * Scw;

				//cout << K->Potok[0] << " " << K->Potok[1] << " " << K->Potok[2] << " " << K->Potok[3] << " " << K->Potok[4] << endl;
				//cout << "A  " << endl;
				//cout << -2.0 * par1.p * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << endl;
				//exit(-1);
			}


			// Блок законов сохранения в ячейке

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;
			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}


			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);
			}
			else
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
			}


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;


			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
				//ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);  // В цилиндрических

				Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
				//Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0.0)
				{
					printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
					ro3 = 0.00001;
				}

				u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2  - q2_1)) / ro3;
				v3 = (ro * v * Volume / Volume2 - T[now1] * ( (K->Potok[2] - 2.0 * p * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q2_2)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				
				//u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y - q2_1)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
				//	+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}


			int k1 = 0, k2 = 0, k3 = 0, k4 = 0;

			if (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				k3 = 0;
				k4 = 0;
				if (K->type == C_1 || K->type == C_centr)//(((kv(u) + kv(v)) / (ggg * p / ro) > 1.0) && ((radius <= 250.0)))
				{
					k1 = 1;
					k2 = 0;
				}
				else
				{
					k1 = 0;
					k2 = 1;
				}
			}
			else
			{
				k1 = 0;
				k2 = 0;
				if (K->type == C_4)//(p / ro > 1.8)   // ( ((y < 8.0)&&(p / ro > (y * (-0.0238) + 0.36)))||( (y >= 8.0)&&(p / ro > 0.17) ) )
				{
					k4 = 0;
					k3 = 1;
				}
				else
				{
					k3 = 0;
					k4 = 1;
				}
			}

			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));


			if (K->type != C_centr)
			{
				//ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 + ro_H1 * v_H1 / y - q1_H1);
				ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.0000001;
				}
				u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2  - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * ((K->Potok_H1[2] - 2.0 * p_H1 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[2] / Volume2 - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2  - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				//u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2 + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[2] / Volume2 + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				//p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2 + //
				//	+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;

				/*if (x > 600.0)
				{
					K->par[now2].ro_H1 = 1.0;
					K->par[now2].p_H1 = 1.0;
					K->par[now2].u_H1 = Velosity_inf;
					K->par[now2].v_H1 = 0.0;
				}*/
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			//ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 + ro_H2 * v_H2 / y - q1_H2);
			ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * ((K->Potok_H2[2] - 2.0 * p_H2 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2 - q22_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[2] / Volume2 - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2  - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			
			//u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[2] / Volume2 + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			//p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2 + //
			//	+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			//ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2 + ro_H3 * v_H3 / y - q1_H3);
			ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2  - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2  - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * ((K->Potok_H3[2] - 2.0 * p_H3 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[2] / Volume2 - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2  - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			
			//u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2 + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[2] / Volume2 + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			//p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2 + //
			//	+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			//ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2 + ro_H4 * v_H4 / y - q1_H4);
			ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2  - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				cout << ro_H1 << " " << ro_H2 << " " << ro_H3 << " " << ro_H4 << endl;
				cout << p_H1 << " " << p_H2 << " " << p_H3 << " " << p_H4 << endl;
				cout << ro << " " << p << " " << u << " " << v << endl;
				ro3 = 0.000001;
			}
			u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * ((K->Potok_H4[2] - 2.0 * p_H4 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[2] / Volume2 - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

			//u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[2] / Volume2 + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			//p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 + //
			//	+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;

			/*if (x > 600.0)
			{
				K->par[now2].ro_H4 = 1.0;
				K->par[now2].p_H4 = 0.5;
				K->par[now2].u_H4 = Velosity_inf;
				K->par[now2].v_H4 = 0.0;
			}*/


		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}


void Setka::Go_5_komponent_MK(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate(T[now1]);

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume(); // K->Get_Volume_rotate(alpha_rot);
			//cout << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle(); // K->Get_Volume_posle_rotate(alpha_rot);
			//cout << "Vol2 = " << Volume2 << endl;
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				//double S = i->Get_square_rotate(alpha_rot); 
				double S = i->Get_square();
				//cout << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));


				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
				}


				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = true;  // Сглаживание
				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					np = false;
					met = 1;
				}
				//np = true;  // Сглаживание
				
				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;

			double u_H1 = K->par[0].H_u[0], v_H1 = K->par[0].H_v[0], ro_H1 = K->par[0].H_n[0], p_H1 = 0.5 * K->par[0].H_T[0] * K->par[0].H_n[0];
			double u_H2 = K->par[0].H_u[1], v_H2 = K->par[0].H_v[1], ro_H2 = K->par[0].H_n[1], p_H2 = 0.5 * K->par[0].H_T[1] * K->par[0].H_n[1];
			double u_H3 = K->par[0].H_u[2], v_H3 = K->par[0].H_v[2], ro_H3 = K->par[0].H_n[2], p_H3 = 0.5 * K->par[0].H_T[2] * K->par[0].H_n[2];
			double u_H4 = K->par[0].H_u[3], v_H4 = K->par[0].H_v[3], ro_H4 = K->par[0].H_n[3], p_H4 = 0.5 * K->par[0].H_T[3] * K->par[0].H_n[3];


			if (false)//(par1.Q / par1.ro < 90)// (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}

			if (ro <= 0.0)
			{
				ro = 0.0000001;
				p = 0.0;
			}
			if (ro_H1 <= 0.0)
			{
				ro_H1 = 0.0000001;
				p_H1 = 0.0;
			}
			if (ro_H2 <= 0.0)
			{
				ro_H2 = 0.0000001;
				p_H2 = 0.0;
			}
			if (ro_H3 <= 0.0)
			{
				ro_H3 = 0.0000001;
				p_H3 = 0.0;
			}
			if (ro_H4 <= 0.0)
			{
				ro_H4 = 0.0000001;
				p_H4 = 0.0;
			}

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			K->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			K->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

			if (false)//(par1.Q / par1.ro < 90)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//
			{
				// Нужно изменить источник, так как он подставляется в уравнения, где параметры плазмы перенормированы 
				K->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);  
			}
			else
			{
				K->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
			}

			q2_1 = K->par[0].M_u * K->par[0].k_u;
			q2_2 = K->par[0].M_v * K->par[0].k_v;
			q3 = K->par[0].M_T * K->par[0].k_T;

			//cout << q2_1 << " " << q2_2 << " " << q3 << endl;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			double ro3, Q33, u3, v3, p3;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

double Setka::HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
	const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
	double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp, bool nul_potok)
	// BestSeries
	// Лучший работающий 2д распадник
	//
	//  Вывод:
	// P[1]       // Скорости
	// P[2]
	// P[0]       // Масса
	// P[3]       // Энергия
{
	double t1 = n2;
	double t2 = -n1;

	double u1, v1, u2, v2;
	u1 = v1_L * n1 + v2_L * n2;
	v1 = v1_L * t1 + v2_L * t2;
	u2 = v1_R * n1 + v2_R * n2;
	v2 = v1_R * t1 + v2_R * t2;

	double sqrtroL = sqrt(ro_L);
	double sqrtroR = sqrt(ro_R);
	double cL = sqrt(ggg * p_L / ro_L);
	double cR = sqrt(ggg * p_R / ro_R);


	double uu_L = (kv(v1_L) + kv(v2_L)) / 2.0;
	double uu_R = (kv(v1_R) + kv(v2_R)) / 2.0;



	double SL = min(u1, u2) - max(cL, cR);
	double SR = max(u1, u2) + max(cL, cR);

	Vl = SL;
	Vp = SR;

	double suR = (SR - u2);
	double suL = (SL - u1);

	double SM = (suR * ro_R * u2 - suL * ro_L * u1 - p_R + p_L) //
		/ (suR * ro_R - suL * ro_L);

	Vc = SM;

	double PTT = (suR * ro_R * p_L - suL * ro_L * p_R + ro_L * ro_R * suR * suL * (u2 - u1)) / (suR * ro_R - suL * ro_L);

	double UU = max(fabs(SL), fabs(SR));
	double time = kurant * rad / UU;

	double FL[5], FR[5], UL[5], UR[5];

	double e1 = p_L / g1 + ro_L * uu_L;
	double e2 = p_R / g1 + ro_R * uu_R;


	FL[0] = ro_L * u1;
	FL[1] = ro_L * u1 * u1 + p_L;
	FL[2] = ro_L * u1 * v1;
	FL[3] = (e1 + p_L) * u1;
	FL[4] = Q_L * u1;

	FR[0] = ro_R * u2;
	FR[1] = ro_R * u2 * u2 + p_R;
	FR[2] = ro_R * u2 * v2;
	FR[3] = (e2 + p_R) * u2;
	FR[4] = Q_R * u2;


	UL[0] = ro_L;
	UL[1] = ro_L * u1;
	UL[2] = ro_L * v1;
	UL[3] = e1;
	UL[4] = Q_L;

	UR[0] = ro_R;
	UR[1] = ro_R * u2;
	UR[2] = ro_R * v2;
	UR[3] = e2;
	UR[4] = Q_R;

	if (SL >= W)
	{
		P[1] = n1 * (FL[1] - W * UL[1]) + t1 * (FL[2] - W * UL[2]);     // Скорости
		P[2] = n2 * (FL[1] - W * UL[1]) + t2 * (FL[2] - W * UL[2]);
		P[0] = FL[0] - W * UL[0];                       // Масса
		P[3] = FL[3] - W * UL[3];                       // Энергия
		PQ = FL[4] - W * UL[4];
		return time;
	}

	if (SR <= W)
	{
		P[1] = n1 * (FR[1] - W * UR[1]) + t1 * (FR[2] - W * UR[2]);     // Скорости
		P[2] = n2 * (FR[1] - W * UR[1]) + t2 * (FR[2] - W * UR[2]);
		P[0] = FR[0] - W * UR[0];                       // Масса
		P[3] = FR[3] - W * UR[3];                       // Энергия
		PQ = FR[4] - W * UR[4];
		return time;
	}


	double ro_LL = ro_L * (SL - u1) / (SL - SM);
	double ro_RR = ro_R * (SR - u2) / (SR - SM);
	double Q_LL = Q_L * (SL - u1) / (SL - SM);
	double Q_RR = Q_R * (SR - u2) / (SR - SM);


	double UZ0 = (SR * UR[0] - SL * UL[0] + FL[0] - FR[0]) / (SR - SL);
	double UZ1 = (SR * UR[1] - SL * UL[1] + FL[1] - FR[1]) / (SR - SL);
	double UZ2 = (SR * UR[2] - SL * UL[2] + FL[2] - FR[2]) / (SR - SL);
	double UZ3 = (SR * UR[3] - SL * UL[3] + FL[3] - FR[3]) / (SR - SL);
	double UZ4 = (SR * UR[4] - SL * UL[4] + FL[4] - FR[4]) / (SR - SL);
	double vzL, vzR, vLL, vRR, ppLR, ee1, ee2;

	if (metod == 0)
	{
		double  PO[5];
		for (int i = 0; i < 5; i++)
		{
			PO[i] = (SR * FL[i] - SL * FR[i] + SR * SL * (UR[i] - UL[i])) / (SR - SL);
		}

		P[1] = n1 * (PO[1] - W * UZ1) + t1 * (PO[2] - W * UZ2);     // Скорости
		P[2] = n2 * (PO[1] - W * UZ1) + t2 * (PO[2] - W * UZ2);
		P[0] = PO[0] - W * UZ0;                       // Масса
		P[3] = PO[3] - W * UZ3;                       // Энергия
		PQ = PO[4] - W * UZ4;
		return time;
	}


	double suRm = suR / (SR - SM);
	double suLm = suL / (SL - SM);
	double rzR = ro_R * suRm;
	double rzL = ro_L * suLm;

	double ptzR = p_R + ro_R * suR * (SM - u2);
	double ptzL = p_L + ro_L * suL * (SM - u1);
	double ptz = (ptzR + ptzL) / 2.0;


	/*if( fabs(v1 - v2) > 0.1)
	{
		vLL = v1;
		vRR = v2;
	}
	else
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}*/


	if (nul_potok == true)   // Некое сглаживание
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}
	else
	{
		vLL = v1;
		vRR = v2;
	}



	ee2 = e2 * suRm + (ptz * SM - p_R * u2) / (SR - SM);
	ee1 = e1 * suLm + (ptz * SM - p_L * u1) / (SL - SM);


	double  ULL[5], URR[5], PO[5];
	ULL[0] = ro_LL;
	ULL[1] = ro_LL * SM;
	ULL[2] = ro_LL * vLL;
	ULL[3] = ee1;
	ULL[4] = Q_LL;

	URR[0] = ro_RR;
	URR[1] = ro_RR * SM;
	URR[2] = ro_RR * vRR;
	URR[3] = ee2;
	URR[4] = Q_RR;

	if (SL < W && SM >= W)
	{
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FL[i] + SL * ULL[i] - SL * UL[i] - W * ULL[i];
		}
	}
	else if (SR > W && SM < W)
	{
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FR[i] + SR * URR[i] - SR * UR[i] - W * URR[i];
		}
	}

	P[1] = n1 * PO[1] + t1 * PO[2];     // Скорости
	P[2] = n2 * PO[1] + t2 * PO[2];
	P[0] = PO[0];                       // Масса
	P[3] = PO[3];                       // Энергия
	PQ = PO[4];

	return time;
}

double Setka::HLLC_2d_Korolkov_b_s_2(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
	const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
	double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp, bool nul_potok)
	// BestSeries
	// Лучший работающий 2д распадник
	//
	//  Вывод:
	// P[1]       // Скорости
	// P[2]
	// P[0]       // Масса
	// P[3]       // Энергия
{
	double t1 = n2;
	double t2 = -n1;

	double u1, v1, u2, v2;
	u1 = v1_L * n1 + v2_L * n2;
	v1 = v1_L * t1 + v2_L * t2;
	u2 = v1_R * n1 + v2_R * n2;
	v2 = v1_R * t1 + v2_R * t2;

	double sqrtroL = sqrt(ro_L);
	double sqrtroR = sqrt(ro_R);
	double cL = sqrt(ggg * p_L / ro_L);
	double cR = sqrt(ggg * p_R / ro_R);


	double uu_L = (kv(v1_L) + kv(v2_L)) / 2.0;
	double uu_R = (kv(v1_R) + kv(v2_R)) / 2.0;



	double SL = min(u1, u2) - max(cL, cR);
	double SR = max(u1, u2) + max(cL, cR);

	Vl = SL;
	Vp = SR;

	double suR = (SR - u2);
	double suL = (SL - u1);

	double SM = (suR * ro_R * u2 - suL * ro_L * u1 - p_R + p_L) //
		/ (suR * ro_R - suL * ro_L);

	Vc = SM;

	double PTT = (suR * ro_R * p_L - suL * ro_L * p_R + ro_L * ro_R * suR * suL * (u2 - u1)) / (suR * ro_R - suL * ro_L);

	double UU = max(fabs(SL), fabs(SR));
	double time = kurant * rad / UU;

	double FL[5], FR[5], UL[5], UR[5];

	double e1 = p_L / g1 + ro_L * uu_L;
	double e2 = p_R / g1 + ro_R * uu_R;


	FL[0] = ro_L * u1;
	FL[1] = ro_L * u1 * u1 + p_L;
	FL[2] = ro_L * u1 * v1;
	FL[3] = (e1 + p_L) * u1;
	FL[4] = Q_L * u1;

	FR[0] = ro_R * u2;
	FR[1] = ro_R * u2 * u2 + p_R;
	FR[2] = ro_R * u2 * v2;
	FR[3] = (e2 + p_R) * u2;
	FR[4] = Q_R * u2;


	UL[0] = ro_L;
	UL[1] = ro_L * u1;
	UL[2] = ro_L * v1;
	UL[3] = e1;
	UL[4] = Q_L;

	UR[0] = ro_R;
	UR[1] = ro_R * u2;
	UR[2] = ro_R * v2;
	UR[3] = e2;
	UR[4] = Q_R;

	if (SL >= W)
	{
		P[1] = n1 * (FL[1] - W * UL[1]) + t1 * (FL[2] - W * UL[2]);     // Скорости
		P[2] = n2 * (FL[1] - W * UL[1]) + t2 * (FL[2] - W * UL[2]);
		P[0] = FL[0] - W * UL[0];                       // Масса
		P[3] = FL[3] - W * UL[3];                       // Энергия
		PQ = FL[4] - W * UL[4];
		return time;
	}

	if (SR <= W)
	{
		P[1] = n1 * (FR[1] - W * UR[1]) + t1 * (FR[2] - W * UR[2]);     // Скорости
		P[2] = n2 * (FR[1] - W * UR[1]) + t2 * (FR[2] - W * UR[2]);
		P[0] = FR[0] - W * UR[0];                       // Масса
		P[3] = FR[3] - W * UR[3];                       // Энергия
		PQ = FR[4] - W * UR[4];
		return time;
	}


	double ro_LL = ro_L * (SL - u1) / (SL - SM);
	double ro_RR = ro_R * (SR - u2) / (SR - SM);
	double Q_LL = Q_L * (SL - u1) / (SL - SM);
	double Q_RR = Q_R * (SR - u2) / (SR - SM);


	double UZ0 = (SR * UR[0] - SL * UL[0] + FL[0] - FR[0]) / (SR - SL);
	double UZ1 = (SR * UR[1] - SL * UL[1] + FL[1] - FR[1]) / (SR - SL);
	double UZ2 = (SR * UR[2] - SL * UL[2] + FL[2] - FR[2]) / (SR - SL);
	double UZ3 = (SR * UR[3] - SL * UL[3] + FL[3] - FR[3]) / (SR - SL);
	double UZ4 = (SR * UR[4] - SL * UL[4] + FL[4] - FR[4]) / (SR - SL);
	double vzL, vzR, vLL, vRR, ppLR, ee1, ee2;

	if (metod == 0)
	{
		double  PO[5];
		for (int i = 0; i < 5; i++)
		{
			PO[i] = (SR * FL[i] - SL * FR[i] + SR * SL * (UR[i] - UL[i])) / (SR - SL);
		}

		P[1] = n1 * (PO[1] - W * UZ1) + t1 * (PO[2] - W * UZ2);     // Скорости
		P[2] = n2 * (PO[1] - W * UZ1) + t2 * (PO[2] - W * UZ2);
		P[0] = PO[0] - W * UZ0;                       // Масса
		P[3] = PO[3] - W * UZ3;                       // Энергия
		PQ = PO[4] - W * UZ4;
		return time;
	}


	double suRm = suR / (SR - SM);
	double suLm = suL / (SL - SM);
	double rzR = ro_R * suRm;
	double rzL = ro_L * suLm;

	double ptzR = p_R + ro_R * suR * (SM - u2);
	double ptzL = p_L + ro_L * suL * (SM - u1);
	double ptz = (ptzR + ptzL) / 2.0;


	/*if( fabs(v1 - v2) > 0.1)
	{
		vLL = v1;
		vRR = v2;
	}
	else
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}*/


	if (nul_potok == true)   // Некое сглаживание
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}
	else
	{
		vLL = v1;
		vRR = v2;
	}



	ee2 = e2 * suRm + (ptz * SM - p_R * u2) / (SR - SM);
	ee1 = e1 * suLm + (ptz * SM - p_L * u1) / (SL - SM);


	double  ULL[5], URR[5], PO[5];
	ULL[0] = ro_LL;
	ULL[1] = ro_LL * SM;
	ULL[2] = ro_LL * vLL;
	ULL[3] = ee1;
	ULL[4] = Q_LL;

	URR[0] = ro_RR;
	URR[1] = ro_RR * SM;
	URR[2] = ro_RR * vRR;
	URR[3] = ee2;
	URR[4] = Q_RR;

	if (SL < W && SM >= W)
	{
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FL[i] + SL * ULL[i] - SL * UL[i] - W * ULL[i];
		}
	}
	else if (SR > W && SM < W)
	{
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FR[i] + SR * URR[i] - SR * UR[i] - W * URR[i];
		}
	}

	P[1] = n1 * PO[1] + t1 * PO[2];     // Скорости
	P[2] = n2 * PO[1] + t2 * PO[2];
	P[0] = PO[0];                       // Масса
	P[3] = PO[3];                       // Энергия
	PQ = PO[4];

	return time;
}

void Setka::Init_Velosity(Sensor* sens, const double& A2, vector <double>& mu, vector <double>& Wt, vector <double>& Wp, vector <double>& Wr, const double& the)
{
	double Y = fabs(Velosity_inf);
	double Yr = -Y * cos(the); 
	double Yt = Y * sin(the); 
	double e1, e2, e3;

	e1 = sqrtpi_ * kv(Yr);
	e2 = 2.0 * fabs(Yr);
	e3 = 0.5 * sqrtpi_;

	double p1, p2, p4, p5;

	p1 = e1 / (e1 + e2 + e3);
	p2 = e2 / (e1 + e2 + e3);

	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7, ksi8;
	double z = 0.0;
	double h = 0.0;
	double Wa = 0.0;

	for (int i = 1; i < I_; i++)
	{
		//cout << 1 << endl;
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

			Wr[i - 1] = z + Yr;
			h = kv(Yr + z) / kv(fabs(Yr) + fabs(z));
			ksi6 = sens->MakeRandom();
		} while (h <= ksi6 || z > -Yr);

		ksi7 = sens->MakeRandom();
		ksi8 = sens->MakeRandom();
		double V = 2.0 * pi_ * ksi7;
		if (i >= 2)
		{
			Wa = sqrt((gam(i - 2) + ksi8 * ((gam(i - 1) - gam(i - 2)))) * kv(Wr[i - 1]));
		}
		else
		{
			Wa = sqrt((ksi8 * ((gam(i - 1)))) * kv(Wr[i - 1]));
		}
		Wt[i - 1] = Wa * cos(V);
		Wp[i - 1] = Wa * sin(V);
		if (i >= 2)
		{
			mu[i - 1] = ((F_mk(gam(i - 1), Yr) - F_mk(gam(i - 2), Yr)) / A2) * fabs(Wr[i - 1]) * exp(-kv(Wt[i - 1] - Yt) - kv(Wp[i - 1]));
		}
		else
		{
			mu[i - 1] = ((F_mk(gam(i - 1), Yr)) / A2) * fabs(Wr[i - 1]) * exp(-kv(Wt[i - 1] - Yt) - kv(Wp[i - 1]));
		}
		//cout << gam(i - 2) << endl;
	}



	p4 = (sqrtpi_ * fabs(Yr)) / (1.0 + sqrtpi_ * fabs(Yr));
	p5 = 1.0 - p4;
	double gamma_ = 0.0;
	if (I_ >= 2)
	{
		gamma_ = gam(I_ - 2);
	}

	do
	{
		//cout << 2 << endl;
		do
		{
			//cout << 3 << endl;
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			if (p4 > ksi1)
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

			Wr[I_ - 1] = z + Yr;
			h = fabs(Yr + z) / (fabs(Yr) + fabs(z));
			ksi4 = sens->MakeRandom();
		} while (h <= ksi4 || z >= -Yr);

		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		//ksi7 = sens->MakeRandom();

		Wt[I_ - 1] = Yt + sqrt(-log(ksi5)) * cos(2.0 * pi_ * ksi6);
		Wp[I_ - 1] = sqrt(-log(ksi5)) * sin(2.0 * pi_ * ksi6);
	} while (kv(Wt[I_ - 1]) + kv(Wp[I_ - 1]) <= gamma_ * kv(Wr[I_ - 1]));

	mu[I_ - 1] = 1.0;
	for (int k = 0; k < I_ - 1; k++)
	{
		mu[I_ - 1] = mu[I_ - 1] - mu[k];
	}

	if (mu[I_ - 1] <= 0.0)
	{
		cout << "6337 ERROR friuhuyvcyvawecsrygvirserccgr  " << mu[I_ - 1] << endl;
		exit(-1);
	}

}

void Setka::Velosity_initial2(Sensor* s, double& Vx, double& Vy, double& Vz)
{
	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(1.0 - ksi2));
	Vy = a * cos(2.0 * pi_ * ksi1);
	Vz = a * sin(2.0 * pi_ * ksi1);
	//cout << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = 0.5 * fabs(Velosity_inf) * sqrtpi_ / (0.5 + 0.5 * fabs(Velosity_inf) * sqrtpi_);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(pi_ * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			z = sqrt(-log(1.0 - ksi4));

		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) <= ksi6 || z < -Velosity_inf);

	Vx = z + Velosity_inf;
	if (Vx <= 0)
	{
		cout << "dfEEERR 32424442" << endl;
	}
	return;
}

void Setka::Velosity_initial(Sensor* s, double& Vx, double& Vy, double& Vz)
{
	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(ksi2));
	Vy = a * cos(2.0 * pi_ * ksi1);
	Vz = a * sin(2.0 * pi_ * ksi1);
	//cout << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = fabs(Velosity_inf) * sqrtpi_ / (1.0 + fabs(Velosity_inf) * sqrtpi_);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(pi_ * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			if (ksi4 <= 0.5)
			{
				z = -sqrt(-log(2.0 * ksi4));
			}
			else
			{
				z = sqrt(-log(2.0 * (1.0 - ksi4)));
			}
		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) < ksi6 || z > -Velosity_inf);

	Vx = z + Velosity_inf;
	return;
}

double Setka::F_mk(const double& gamma, const double& Yr)   // Нужно задать
{
	return gamma * ((0.5 + kv(Yr)) * (1.0 + erf(-Yr)) - Yr * exp(-kv(Yr)) / sqrtpi_);
}

void Setka::Init_Pozision(Sensor* sens, const double& A1, double& phi, double& the)
{
	double ksi0, ksi1, ksi2, ksi3, ksi4, ksi5, ksi6;
	double X, h;
	double Y = fabs(Velosity_inf);

	ksi0 = sens->MakeRandom();
	ksi1 = sens->MakeRandom();
	phi = 2.0 * pi_ * ksi0;
	double p1 = erf(Y) / (A1 * kv(Y));

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

	the = acos(X);
}
 
Cell* Setka::Belong_point(int b, const double& x, const double& y)
{
	if (b == 1)  // Внешняя сфера
	{
		for (auto& i : this->Cell_sphere)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}
	else if (b == 2)  // Боковая граница
	{
		for (auto& i : this->Cell_side)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}
	else if (b == 3)  // Задняя граница
	{
		for (auto& i : this->Cell_back)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}
	else if (b == 4)  // Передний диск
	{
		for (auto& i : this->Cell_disk)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}

	cout << "ERRORRRORORfireubfvwcefrvjywgkvygdcwkug324h324334" << endl;
	cout << x << " " << y << endl;
	return nullptr;
}

void Setka::M_K_prepare(void)
{
	double xx1, yy1;
	double xx2, yy2;
	bool bb1 = false;
	bool bb2 = false;


	// Заполним граничные ячейки
	//for (auto& i : this->All_Cells)    
	//{
	//	for (auto& j : i->Grans)
	//	{
	//		if (j->type == Input)
	//		{
	//			this->Cell_sphere.push_back(i);
	//			//cout << j->A->x << " " << j->A->y << endl;
	//			break;
	//		}
	//	}
	//}


	for (auto& i : this->All_Cells)    // Заполняем массивы начальными гранями
	{
		bb1 = false;
		bb2 = false;
		i->renew();
		for (auto& j : i->contour)
		{
			if (sqrt(kvv(j->x, j->y, 0.0)) <= Rmax_ && j->x > -0.1)
			{
				bb1 = true;
			}
			else if (sqrt(kvv(j->x, j->y, 0.0)) >= Rmax_ && j->x > -0.1)
			{
				bb2 = true;
			}
		}

		if (bb1 == true && bb2 == true)
		{
			this->Cell_sphere.push_back(i);
		}

		if (bb2 == true)
		{
			this->Cell_other.push_back(i);
		}

		i->par[0].F_n = 0.0;
		i->par[0].F_u = 0.0;
		i->par[0].F_v = 0.0;
		i->par[0].F_T = 0.0;
		i->par[0].I_u = 0.0;
		i->par[0].I_v = 0.0;
		i->par[0].I_T = 0.0;
		i->par[0].II_u = 0.0;
		i->par[0].II_v = 0.0;
		i->par[0].II_T = 0.0;
		i->par[0].H_n[0] = 0.0;
		i->par[0].H_n[1] = 0.0;
		i->par[0].H_n[2] = 0.0;
		i->par[0].H_n[3] = 0.0;
		i->par[0].H_u[0] = 0.0;
		i->par[0].H_u[1] = 0.0;
		i->par[0].H_u[2] = 0.0;
		i->par[0].H_u[3] = 0.0;
		i->par[0].H_v[0] = 0.0;
		i->par[0].H_v[1] = 0.0;
		i->par[0].H_v[2] = 0.0;
		i->par[0].H_v[3] = 0.0;
		i->par[0].H_T[0] = 0.0;
		i->par[0].H_T[1] = 0.0;
		i->par[0].H_T[2] = 0.0;
		i->par[0].H_T[3] = 0.0;
		for (int li = 0; li < 7; li++)
		{
			i->par[0].w_m[li] = 0.0;
		}
	}

	cout << "this->Cell_sphere = " << this->Cell_sphere.size() << endl;

	for (auto& i : this->All_Cells)    // Заполняем массивы начальными гранями
	{
		for (auto& j : i->Grans)
		{
			if (j->type == Upper_wall)
			{
				this->Cell_side.push_back(i);
				break;
			}
		}
	}

	for (auto& i : this->All_Cells)    // Заполняем массивы начальными ячейками
	{
		bb1 = false;
		bb2 = false;

		for (auto& j : i->contour)
		{
			if (j->x <= -1.0 && j->y > Rmax_)
			{
				bb1 = true;
			}
			else if (j->x >= -1.0 && j->y > Rmax_)
			{
				bb2 = true;
			}
		}

		if (bb1 == true && bb2 == true)
		{
			this->Cell_disk.push_back(i);
		}
	}

	cout << "this->Cell_disk = " << this->Cell_disk.size() << endl;

	for (auto& i : this->All_Cells)    // Заполняем массивы начальными ячейками
	{
		for (auto& j : i->Grans)
		{
			if (j->type == Extern)
			{
				this->Cell_back.push_back(i);
				break;
			}
		}
	}

	for (auto& i : this->All_Gran)    // Обновим уравнения граней для правильного нахождения пересечения траекторий с ними
	{
		i->renew();
	}

	for (auto& i : this->All_Gran_copy)
	{
		i->renew();
	}

	for (auto& i : this->All_Cells) 
	{
		//if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		//{
		//	i->par[0].u = i->par[0].u * (chi_real / chi_);       // Перенормировка
		//	i->par[0].v = i->par[0].v * (chi_real / chi_);
		//	i->par[0].ro = i->par[0].ro / kv(chi_real / chi_);
		//	i->par[0].Q = i->par[0].Q / kv(chi_real / chi_);
		//}

		for (auto& j : i->Grans)
		{
			if (j->Sosed != nullptr)
			{
				i->L = min(i->L, j->Sosed->L);
			}
		}
	}


	// Блок загрузки датчиков случайных чисел
	ifstream fin2;
	fin2.open("my_dat_822.txt");
	double d, a1, b1, c;
	for (int i = 0; i < 822; i++)
	{
		fin2 >> a1 >> b1 >> c;
		auto s = new Sensor(a1, b1, c);
		this->Sensors.push_back(s);
	}
	fin2.close();

	fin2.open("Sensors276405.txt");
	unsigned int v1, v2, v3, v4;
	for (int i = 0; i < 276405; i++)
	{
		fin2 >> v1 >> v2 >> v3 >> v4;
		if (i % 1000 == 0)
		{
			auto s = new sensor2(v1, v2, v3, v4);
			this->Sensors2.push_back(s);
		}
	}


	double Y = fabs(Velosity_inf);
	this->sqv_1 = (Rmax_ - 3.0) * (0.5 * (Rmax_ - 3.0) * pi_ * Y * (erf(Y) * (1.0 + 1.0 / (2.0 * kv(Y))) + 1.0 + exp(-kv(Y)) / (Y * sqrtpi_)));
	this->sqv_2 = sqrtpi_ * (R5_) * fabs(Left_);
	this->sqv_3 = pi_ * kv(R5_) * exp(-kv(Velosity_inf)) * (1.0 + exp(kv(Velosity_inf)) * sqrtpi_ * Velosity_inf * (1.0 + erf(Velosity_inf)))/(2.0 * sqrtpi_);
	this->sqv_4 = (sqrtpi_ / 2.0) * (kv((Rmax_ - 3.0)) - kv(R5_ - 0.1)) * exp(-kv(Velosity_inf)) * (sqrtpi_ * Velosity_inf * erfc(Velosity_inf) * exp(kv(Velosity_inf)) - 1.0);
	cout << "this->sqv_1 = " << this->sqv_1 << endl;
	cout << "this->sqv_4 = " << this->sqv_4 << endl;
	this->sum_s = this->sqv_1 + this->sqv_2 + this->sqv_3 + this->sqv_4;
	this->Number1 = 411 * 243 * 10; // 6000;  25
	this->Number2 = 411 * 10 * 2; // 100;
	this->Number3 = 411 * 3; // 100;
	this->Number4 = 411 * 100 * 5; // 411 * 1650; // 135 * 40; // 400;
	this->AllNumber = ((this->Number1) + (this->Number2) + (this->Number3) + (this->Number4));
	cout << "this->AllNumber " << this->AllNumber << endl;

	Ri.resize(I_ + 1);
	//Ri[0] = R5_ - 2.0;
	Ri[0] = 1.0;                 // Здесь задаются радиусы
	Ri[1] = 2.61;
	Ri[2] = 6.84;
	Ri[3] = 17.92;
	Ri[4] = 46.9;
	Ri[5] = 122.7;
	Ri[6] = 321.23;                 // Здесь задаются радиусы
	Ri[7] = 840.65;
	Ri[8] = (Rmax_ - 3.0);

	//Ri[4] = 540.0; // 216.0;
	//Ri[5] = Rmax_;
	//Ri[0] = 1.0;                 // Здесь задаются радиусы
	//Ri[1] = 2.0;
	//Ri[2] = 4.0;


	for (auto& i : this->All_Cells)    // Заполняем зоны ячейки
	{                                    
		double x, y;
		i->Get_Center(x, y);
		double r = sqrt(kvv(x, y, 0.0));
		i->zona = geo_zones(r);
		i->zona_alpha = alpha_zones(x, y);
	}

	double mu1;
	mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);  // Учёт, что для большого количества траекторий исходный вес меньше


	double kas = 0.3; // 0.01 // 0.00025; //0.001;
	cout << "kas = " << kas << endl;
	double ss1 = 0.1; // 0.05;
	double ss2 = 0.01; // 0.003;
	double ss3 = 0.4; // 0.4;

	//double kas = 0.25 * mu1 * 100.0 * 10.0 * 0.25; // 0.00025; //0.001;
	//double ss1 = 1.0;
	//double ss2 = 1.0;
	//double ss3 = 1.0;

	// Основной сорт
	/*Mu[3][0] = 0.0001 * kas;
	Mu[3][1] = 0.001 * kas;
	Mu[3][2] = 0.01 * kas;
	Mu[3][3] = 0.1 * kas;
	Mu[3][4] = 1.0 * kas;
	Mu[3][5] = 10.0 * kas;
	Mu[3][6] = 100.0 * kas;
	Mu[3][7] = 1000.0 * kas;
	Mu[3][8] = 10000.0 * kas;*/

	Mu[3][0] = kv(Ri[0] / Rmax_) * kas;
	Mu[3][1] = kv(Ri[1] / Rmax_) * kas;
	Mu[3][2] = kv(Ri[2] / Rmax_) * kas;
	Mu[3][3] = kv(Ri[3] / Rmax_) * kas;
	Mu[3][4] = kv(Ri[4] / Rmax_) * kas;
	Mu[3][5] = kv(Ri[5] / Rmax_) * kas;
	Mu[3][6] = kv(Ri[6] / Rmax_) * kas;
	Mu[3][7] = kv(Ri[7] / Rmax_) * kas;
	Mu[3][8] = kv(Ri[8] / Rmax_) * kas;

	/*for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
			for (int k = 0; k < J_ + 1; k++)
			{
				this->Mu_stat[i][j][k] = Mu_stat[i][j][k] * kv(Ri[j] / Rmax_);
			}
		}
	}*/

	/*Mu[3][0] = 0.0001 * kas;
	Mu[3][1] = 0.0001 * kas;
	Mu[3][2] = 0.0001 * kas;
	Mu[3][3] = 0.0001 * kas;
	Mu[3][4] = 0.0001 * kas;
	Mu[3][5] = 0.0001 * kas;
	Mu[3][6] = 0.0001 * kas;
	Mu[3][7] = 0.0001 * kas;
	Mu[3][8] = 0.0001 * kas;*/


	Mu[0][0] = Mu[3][0] * ss1;
	Mu[0][1] = Mu[3][1] * ss1;
	Mu[0][2] = Mu[3][2] * ss1;
	Mu[0][3] = Mu[3][3] * ss1;
	Mu[0][4] = Mu[3][4] * ss1;
	Mu[0][5] = Mu[3][5] * ss1;
	Mu[0][6] = Mu[3][6] * ss1;
	Mu[0][7] = Mu[3][7] * ss1;
	Mu[0][8] = Mu[3][8] * ss1;

	Mu[1][0] = Mu[3][0] * ss2;
	Mu[1][1] = Mu[3][1] * ss2;
	Mu[1][2] = Mu[3][2] * ss2;
	Mu[1][3] = Mu[3][3] * ss2;
	Mu[1][4] = Mu[3][4] * ss2;
	Mu[1][5] = Mu[3][5] * ss2;
	Mu[1][6] = Mu[3][6] * ss2;
	Mu[1][7] = Mu[3][7] * ss2;
	Mu[1][8] = Mu[3][8] * ss2;

	Mu[2][0] = Mu[3][0] * ss3;
	Mu[2][1] = Mu[3][1] * ss3;
	Mu[2][2] = Mu[3][2] * ss3;
	Mu[2][3] = Mu[3][3] * ss3;
	Mu[2][4] = Mu[3][4] * ss3;
	Mu[2][5] = Mu[3][5] * ss3;
	Mu[2][6] = Mu[3][6] * ss3;
	Mu[2][7] = Mu[3][7] * ss3;
	Mu[2][8] = Mu[3][8] * ss3;
}

int Setka::geo_zones(const double& r, const double& k)
// В какой геометрической зоне сейчас находится атом (зоны считаются с нуля)
{
	if (r < 0.0)
	{
		cout << "Error  7151 geo_zones" << endl;
	}

	for (int i = 0; i < I_; i++)
	{
		if (r < k * this->Ri[i])
		{
			return i;
		}
	}
	return I_;
}

int Setka::alpha_zones(const double& x, const double& y)
// В какой геометрической зоне сейчас находится атом (зоны считаются с нуля)
{
	double al = polar_angle(x, y);

	for (int i = 0; i < J_; i++)
	{
		if (al < (i + 1) * pi_ / (J_ + 1))
		{
			return i;
		}
	}
	return J_;
}

void Setka::MK_start(void)
{

	mutex mut_1;
	double Y = fabs(Velosity_inf);
	int st = 1;

#pragma omp parallel for
	for (int index = 0; index < 135; index++)
	{
		double A1 = 1.0 + (1.0 + 1.0 / (2.0 * kv(Y))) * erf(Y) + exp(-kv(Y)) / (Y * sqrtpi_);
		double A2 = 0.0; // Нужно задать после нахождения угла theta
		double mu1, mu2, mu3, mu4;
		double Vx, Vy, Vz;
		mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);
		mu2 = ((this->sqv_2) / this->sum_s) * (1.0 * this->AllNumber / this->Number2);
		mu3 = ((this->sqv_3) / this->sum_s) * (1.0 * this->AllNumber / this->Number3);
		mu4 = ((this->sqv_4) / this->sum_s) * (1.0 * this->AllNumber / this->Number4);
		if (index == 0)
		{
			mut_1.lock();
			cout << "Wright weyght   " << mu1 << " " << mu2 << " " << mu3 << " " << mu4 << endl;
			mut_1.unlock();
		}
		Sensor* sens1 = Sensors[2 * index];
		Sensor* sens2 = Sensors[2 * index + 1];
		MKmethod MK = MKmethod();
		mut_1.lock();
		cout << st << " potok  is  135;  index = " << index << endl;
		st++;
		mut_1.unlock();
		double x, y, z;
		double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7;
		double phi, Vphi, Vr;
		double a, b, c;
		double r;

		for (int ii = 0; ii < Number1 / 135; ii++)  //  Вылет с полусферы
		{
			/*if (ii % 100000 == 0 && ii > 1 )
			{
				cout << "ii =  " << ii << endl;
			}*/
			double phi, the;
			// Нахождение начальной траектории
			this->Init_Pozision(sens1, A1, phi, the);
			// перевод позиции в декартовы координаты
			x = (R5_ - 2.0) * cos(the);
			y = (R5_ - 2.0) * sin(the) * cos(phi);
			z = (R5_ - 2.0) * sin(the) * sin(phi);
			//cout << x << " " << sqrt(kv(y)+ kv(z)) << endl;

			Cell* Point = Belong_point(1, x, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			vector <double> mu(I_);
			vector <double> Wt(I_);
			vector <double> Wp(I_);
			vector <double> Wr(I_);

			double X = cos(the);
			A2 = exp(-kv(Y) * kv(X)) / sqrtpi_ + Y * X * (1.0 + erf(Y * X));
			Init_Velosity(sens1, A2, mu, Wt, Wp, Wr, the);

			//Velosity_initial(sens1, Vx, Vy, Vz);
			//Fly_exchenge_Imit(sens2, x, y, z, Vx, Vy, Vz, Point, mu1, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu1, I_ - 1, 5);

			//cout << "Start" << endl;
			if (true)
			{
				for (int i = 0; i < I_; i++)
				{
					dekard_skorost2((R5_ - 2.0), the, phi, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);
					/*for (int k = 0; k < 60000; k++)
					{
						cout << x + 0.3 * k * Vx << " " << sqrt(kv(y + 0.3 * k * Vy) + kv(z + 0.3 * k * Vz)) << endl;
					}*/
					//cout << "START" << endl;
					//Fly_exchenge(sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1, mu[i] * mu1, false);
					//Fly_exchenge_Split(sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1, mu[i] * mu1, false, i);
					//cout << "A" << endl;

					//Fly_exchenge_Imit_Korol(sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1,  3, false, mu1, I_ - 1);
					//cout << ii << endl;
					Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu[i] * mu1, I_ - 1, ii);

				}
			}
		}
		for (int ii = 0; ii < Number2 / 135; ii++)  //  Вылет сверху 
		{
			//cout << "A" << endl;
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();
			ksi3 = sens1->MakeRandom();
			ksi4 = sens1->MakeRandom();
			ksi5 = sens1->MakeRandom();

			x = (Left_ + 0.1) + ksi1 * (- 0.2 - Left_);
			phi = ksi2 * 2.0 * pi_;
			Vphi = cos(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vx = Velosity_inf + sin(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vr = -sqrt(-log(ksi5));
			y = (R5_ - 0.2) * cos(phi);
			z = (R5_ - 0.2) * sin(phi);

			Cell* Point = Belong_point(2, x, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//Fly_exchenge_Imit_Korol(sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
			//	Point, mu2, 3, false, mu2, I_ - 1);
			Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
				Point, mu2, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu2, I_ - 1, ii);
		}
		for (int ii = 0; ii < Number3 / 135; ii++)  //  Вылет сзади
		{
			//cout << "B" << endl;
			Velosity_initial2(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(ksi1 * (R5_ - 2.0) * (R5_ - 2.0));
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(3, Left_ + 2.0, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//cout << "Start" << endl;
			//Fly_exchenge_Imit_Korol(sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1);
			Fly_exchenge_Imit(MK, sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, -log(1.0 - sens1->MakeRandom()), 0.0, //
				3, mu3, I_ - 1, ii);
			//cout << "Stop" << endl;
		}
		for (int ii = 0; ii < Number4 / 135; ii++)  //Number4 / 135
		{
			double a, b, c;
			Velosity_initial(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(-(ksi1 * (kv(Rmax_) - kv(R5_)) - kv(Rmax_)));
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(4, -0.1, sqrt(kv(z) + kv(y)));

			Fly_exchenge_Imit(MK, sens2, -0.1, y, z, a, b, c, Point, mu4, -log(1.0 - sens1->MakeRandom()), 0.0, //
				3, mu4, I_ - 1, ii);
		}
	}


	for (auto& k : this->All_Cells)
	{

		/*for (int il = 0; il < 7; il++)
		{
			if (k->par[0].num_atoms > 0)
			{
				k->par[0].w_m[il] /= k->par[0].num_atoms;
			}
		}*/


		double no = (1.0 * AllNumber * k->Get_Volume_rotate(360.0));

		k->par[0].I_u = sum_s * k->par[0].I_u / no;
		k->par[0].I_v = sum_s * k->par[0].I_v / no;
		k->par[0].I_T = sum_s * k->par[0].I_T / no;

		k->par[0].II_u = sum_s * k->par[0].II_u / no;
		k->par[0].II_v = sum_s * k->par[0].II_v / no;
		k->par[0].II_T = sum_s * k->par[0].II_T / no;

		if (k->par[0].F_n > 0)
		{
			/*k->par[0].I_u = k->par[0].I_u / k->par[0].F_n;
			k->par[0].I_v = k->par[0].I_v / k->par[0].F_n;
			k->par[0].I_T = k->par[0].I_T / k->par[0].F_n;*/

			k->par[0].F_u = k->par[0].F_u / k->par[0].F_n;
			k->par[0].F_v = k->par[0].F_v / k->par[0].F_n;
			k->par[0].F_T = (2.0 / 3.0) * (k->par[0].F_T / k->par[0].F_n - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
		}
		else
		{
			/*k->par[0].I_u = 0.0;
			k->par[0].I_v = 0.0;
			k->par[0].I_T = 0.0;*/

			k->par[0].F_n = 0.0;
			k->par[0].F_u = 0.0;
			k->par[0].F_v = 0.0;
			k->par[0].F_T = 0.0;

		}

		if (k->par[0].F_n > 0)
		{
			for (int i = 0; i < 4; i++)
			{
				if (k->par[0].H_n[i] > 0)
				{
					k->par[0].H_u[i] = k->par[0].H_u[i] / k->par[0].H_n[i];
					k->par[0].H_v[i] = k->par[0].H_v[i] / k->par[0].H_n[i];
					//k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
					k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].H_u[i], k->par[0].H_v[i], 0.0));
				}
				else
				{
					k->par[0].H_u[i] = 0.0;
					k->par[0].H_v[i] = 0.0;
					k->par[0].H_T[i] = 0.0;
				}
			}
		}

		k->par[0].F_n = sum_s * k->par[0].F_n / no;
		for (int i = 0; i < 4; i++)
		{
			k->par[0].H_n[i] = sum_s * k->par[0].H_n[i] / no;
		}

		/*k->par[0].I_u = (n_p_LISM_ / Kn_) * (-k->par[0].ro * k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_ / Kn_) * (-k->par[0].ro * k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_ / Kn_) * (-(1.0 / 2.0) * k->par[0].ro * k->par[0].I_T);*/

		/*k->par[0].I_u = (n_p_LISM_ / Kn_) * (k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_ / Kn_) * (k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_ / Kn_) * (k->par[0].I_T);*/
		k->par[0].I_u = (n_p_LISM_) * (k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_) * (k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_) * (k->par[0].I_T);
		k->par[0].II_u = (n_p_LISM_) * (k->par[0].II_u);
		k->par[0].II_v = (n_p_LISM_) * (k->par[0].II_v);
		k->par[0].II_T = (n_p_LISM_) * (k->par[0].II_T);

		// Считаем мультифдюидные источники

		double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
		double U_H1, U_H2, U_H3, U_H4;
		double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
		double nu_H1, nu_H2, nu_H3, nu_H4;
		double q2_1, q2_2, q3;
		double u, v, ro, p, Q;

		double u_H1 = k->par[0].H_u[0], v_H1 = k->par[0].H_v[0], ro_H1 = k->par[0].H_n[0], p_H1 = 0.5 * k->par[0].H_T[0] * k->par[0].H_n[0];
		double u_H2 = k->par[0].H_u[1], v_H2 = k->par[0].H_v[1], ro_H2 = k->par[0].H_n[1], p_H2 = 0.5 * k->par[0].H_T[1] * k->par[0].H_n[1];
		double u_H3 = k->par[0].H_u[2], v_H3 = k->par[0].H_v[2], ro_H3 = k->par[0].H_n[2], p_H3 = 0.5 * k->par[0].H_T[2] * k->par[0].H_n[2];
		double u_H4 = k->par[0].H_u[3], v_H4 = k->par[0].H_v[3], ro_H4 = k->par[0].H_n[3], p_H4 = 0.5 * k->par[0].H_T[3] * k->par[0].H_n[3];
		u = k->par[0].u;
		v = k->par[0].v;
		ro = k->par[0].ro;
		p = k->par[0].p;

		if (ro <= 0.0)
		{
			ro = 0.0000001;
			p = 0.0;
		}
		if (ro_H1 <= 0.0)
		{
			ro_H1 = 0.0000001;
			p_H1 = 0.0;
		}
		if (ro_H2 <= 0.0)
		{
			ro_H2 = 0.0000001;
			p_H2 = 0.0;
		}
		if (ro_H3 <= 0.0)
		{
			ro_H3 = 0.0000001;
			p_H3 = 0.0;
		}
		if (ro_H4 <= 0.0)
		{
			ro_H4 = 0.0000001;
			p_H4 = 0.0;
		}

		U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H1 / ro_H1));
		U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H2 / ro_H2));
		U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H3 / ro_H3));
		U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H4 / ro_H4));

		U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H1 / ro_H1));
		U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H2 / ro_H2));
		U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H3 / ro_H3));
		U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H4 / ro_H4));

		sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
		sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
		sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
		sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

		nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
		nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
		nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
		nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

		k->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
			+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
		k->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
			+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

		
		k->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
			(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
			nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
				(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
			nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
				(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
			nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
				(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));

		
		k->par[0].k_u = k->par[0].I_u / k->par[0].M_u;
		if (fabs(k->par[0].k_u) > 10.0)
		{
			k->par[0].k_u = 1.0;
		}
		k->par[0].k_v = k->par[0].I_v / k->par[0].M_v;
		if (fabs(k->par[0].k_v) > 10.0)
		{
			k->par[0].k_v = 1.0;
		}
		k->par[0].k_T = k->par[0].I_T / k->par[0].M_T;
		if (fabs(k->par[0].k_T) > 10.0)
		{
			k->par[0].k_T = 1.0;
		}

		
	}
}

void Setka::MK_start_new(void)
{

	mutex mut_1;
	double Y = fabs(Velosity_inf);
	int st = 1;

#pragma omp parallel for
	for (int index = 0; index < 411; index++)
	{
		double A1 = 1.0 + (1.0 + 1.0 / (2.0 * kv(Y))) * erf(Y) + exp(-kv(Y)) / (Y * sqrtpi_);
		double A2 = 0.0; // Нужно задать после нахождения угла theta
		double mu1, mu2, mu3, mu4;
		double Vx, Vy, Vz;
		mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);
		mu2 = ((this->sqv_2) / this->sum_s) * (1.0 * this->AllNumber / this->Number2);
		mu3 = ((this->sqv_3) / this->sum_s) * (1.0 * this->AllNumber / this->Number3);
		mu4 = ((this->sqv_4) / this->sum_s) * (1.0 * this->AllNumber / this->Number4);
		if (index == 0)
		{
			mut_1.lock();
			cout << "Wright weyght   " << mu1 << " " << mu2 << " " << mu3 << " " << mu4 << endl;
			mut_1.unlock();
		}
		Sensor* sens1 = Sensors[2 * index];
		Sensor* sens2 = Sensors[2 * index + 1];
		MKmethod MK = MKmethod();
		mut_1.lock();
		cout << st << " potok  is  411;  index = " << index << endl;
		st++;
		mut_1.unlock();
		double x, y, z;
		double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7;
		double phi, Vphi, Vr;
		double a, b, c;
		double r;

		for (int ii = 0; ii < Number1 / 411; ii++)  //  Вылет с полусферы   Number1 / 135
		{
			if (true)
			{
				vector <double> mu(I_ + 1);
				vector <double> Wt(I_ + 1);
				vector <double> Wp(I_ + 1);
				vector <double> Wr(I_ + 1);
				vector <double> X(I_ + 1);
				double phi, sin_;
				//bool bb = MK.Init_Parametrs(sens1, mu, Wt, Wp, Wr, X);
				//cout << "A " << endl;
				bool bb = MK.Init_Parametrs(sens1, mu, Wt, Wp, Wr, X);
				//cout << "B " << endl;

				for (int i = 0; i <= I_; i++) // I_
				{
					sin_ = sqrt(1.0 - kv(X[i]));
					x = (Rmax_ - 3.0) * X[i];
					phi = 2.0 * pi_ * sens1->MakeRandom();
					y = (Rmax_ - 3.0) * sin_ * cos(phi);
					z = (Rmax_ - 3.0) * sin_ * sin(phi);
					Cell* Point = Belong_point(1, x, sqrt(kv(z) + kv(y)));
					//dekard_skorost2(Rmax_, acos(X[i]), phi, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);
					dekard_skorost(y, z, x, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);

					//Fly_exchenge_Imit_Korol(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1 * mu[i], i);

					if (i != I_ || bb == true)
					{
						double time_do_peregel, peregel;
						int ii_z, ii_alp;
						if (Vx * x + Vy * y + Vz * z < 0.0)
						{
							time_do_peregel = (-Vx * x - Vy * y - Vz * z) / kvv(Vx, Vy, Vz);
							peregel = sqrt(kvv(x + Vx * time_do_peregel, y + Vy * time_do_peregel, z + Vz * time_do_peregel));
							ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
							ii_alp = alpha_zones(x + Vx * time_do_peregel, sqrt(kvv(y + Vy * time_do_peregel, z + Vz * time_do_peregel, 0.0)));
						}
						else
						{
							ii_z = I_;
							ii_alp = alpha_zones(x, sqrt(kvv(y, z, 0.0)));
						}

						Fly_exchenge_Imit_Korol(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1, ii_z, ii_alp, true);
						//Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu1, i, ii); // i
					}
				}
			}

		}
		for (int ii = 0; ii < Number2 / 411; ii++)  //  Вылет сверху 
		{
			//cout << "A" << endl;
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();
			ksi3 = sens1->MakeRandom();
			ksi4 = sens1->MakeRandom();
			ksi5 = sens1->MakeRandom();

			x = (Left_ + 0.1) + ksi1 * (-0.2 - Left_);
			phi = ksi2 * 2.0 * pi_;
			Vphi = cos(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vx = Velosity_inf + sin(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vr = -sqrt(-log(ksi5));
			y = (R5_ - 0.2) * cos(phi);
			z = (R5_ - 0.2) * sin(phi);

			Cell* Point = Belong_point(2, x, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//Fly_exchenge_Imit_Korol(sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
			//	Point, mu2, 3, false, mu2, I_ - 1);
			
			Fly_exchenge_Imit_Korol(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
				Point, mu2, 3, false, mu2, I_ - 1, J_ - 1, false);
			//Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
			//	Point, mu2, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu2, I_ - 1, ii);
		}
		for (int ii = 0; ii < Number3 / 411; ii++)  //  Вылет сзади
		{
			//cout << "B" << endl;
			Velosity_initial2(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(ksi1 * (R5_ - 2.0) * (R5_ - 2.0));
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(3, Left_ + 2.0, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//cout << "Start" << endl;
			//Fly_exchenge_Imit_Korol(sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1);
			
			Fly_exchenge_Imit_Korol(MK, sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1, J_ - 1, false);
			//Fly_exchenge_Imit(MK, sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, -log(1.0 - sens1->MakeRandom()), 0.0, //
			//	3, mu3, I_ - 1, ii);
		}
		for (int ii = 0; ii < Number4 / 411; ii++)  //Number4 / 411
		{
			double a, b, c;
			Velosity_initial(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(ksi1 * (kv(R5_ - 0.1) - kv(Rmax_)) + kv(Rmax_));
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(4, -0.01, sqrt(kv(z) + kv(y)));

			//Fly_exchenge_Imit_Korol(MK, sens2, -0.01, y, z, a, b, c, Point, mu4, 3, false, mu4, I_ - 1, J_ - 1, true);
			Fly_exchenge_Imit(MK, sens2, -0.01, y, z, a, b, c, Point, mu4, -log(1.0 - sens1->MakeRandom()), 0.0, //
				3, mu4, I_ - 1, ii);
		}
	}


	for (auto& k : this->All_Cells)
	{
		/*for (int il = 0; il < 7; il++)
		{
			if (k->par[0].num_atoms > 0)
			{
				k->par[0].w_m[il] /= k->par[0].num_atoms;
			}
		}*/


		double no = (1.0 * AllNumber * k->Get_Volume_rotate(360.0));

		k->par[0].I_u = sum_s * k->par[0].I_u / no;
		k->par[0].I_v = sum_s * k->par[0].I_v / no;
		k->par[0].I_T = sum_s * k->par[0].I_T / no;

		k->par[0].II_u = sum_s * k->par[0].II_u / no;
		k->par[0].II_v = sum_s * k->par[0].II_v / no;
		k->par[0].II_T = sum_s * k->par[0].II_T / no;

		if (k->par[0].F_n > 0)
		{
			/*k->par[0].I_u = k->par[0].I_u / k->par[0].F_n;
			k->par[0].I_v = k->par[0].I_v / k->par[0].F_n;
			k->par[0].I_T = k->par[0].I_T / k->par[0].F_n;*/

			k->par[0].F_u = k->par[0].F_u / k->par[0].F_n;
			k->par[0].F_v = k->par[0].F_v / k->par[0].F_n;
			k->par[0].F_T = (2.0 / 3.0) * (k->par[0].F_T / k->par[0].F_n - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
		}
		else
		{
			/*k->par[0].I_u = 0.0;
			k->par[0].I_v = 0.0;
			k->par[0].I_T = 0.0;*/

			k->par[0].F_n = 0.0;
			k->par[0].F_u = 0.0;
			k->par[0].F_v = 0.0;
			k->par[0].F_T = 0.0;

		}

		if (k->par[0].F_n > 0)
		{
			for (int i = 0; i < 4; i++)
			{
				if (k->par[0].H_n[i] > 0)
				{
					k->par[0].H_u[i] = k->par[0].H_u[i] / k->par[0].H_n[i];
					k->par[0].H_v[i] = k->par[0].H_v[i] / k->par[0].H_n[i];
					//k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
					k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].H_u[i], k->par[0].H_v[i], 0.0));
				}
				else
				{
					k->par[0].H_u[i] = 0.0;
					k->par[0].H_v[i] = 0.0;
					k->par[0].H_T[i] = 0.0;
				}
			}
		}

		k->par[0].F_n = sum_s * k->par[0].F_n / no;
		for (int i = 0; i < 4; i++)
		{
			k->par[0].H_n[i] = sum_s * k->par[0].H_n[i] / no;
		}

		/*k->par[0].I_u = (n_p_LISM_ / Kn_) * (-k->par[0].ro * k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_ / Kn_) * (-k->par[0].ro * k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_ / Kn_) * (-(1.0 / 2.0) * k->par[0].ro * k->par[0].I_T);*/

		/*k->par[0].I_u = (n_p_LISM_ / Kn_) * (k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_ / Kn_) * (k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_ / Kn_) * (k->par[0].I_T);*/
		k->par[0].I_u = (n_p_LISM_) * (k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_) * (k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_) * (k->par[0].I_T);
		k->par[0].II_u = (n_p_LISM_) * (k->par[0].II_u);
		k->par[0].II_v = (n_p_LISM_) * (k->par[0].II_v);
		k->par[0].II_T = (n_p_LISM_) * (k->par[0].II_T);

		// Считаем мультифдюидные источники

		double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
		double U_H1, U_H2, U_H3, U_H4;
		double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
		double nu_H1, nu_H2, nu_H3, nu_H4;
		double q2_1, q2_2, q3;
		double u, v, ro, p, Q;

		double u_H1 = k->par[0].H_u[0], v_H1 = k->par[0].H_v[0], ro_H1 = k->par[0].H_n[0], p_H1 = 0.5 * k->par[0].H_T[0] * k->par[0].H_n[0];
		double u_H2 = k->par[0].H_u[1], v_H2 = k->par[0].H_v[1], ro_H2 = k->par[0].H_n[1], p_H2 = 0.5 * k->par[0].H_T[1] * k->par[0].H_n[1];
		double u_H3 = k->par[0].H_u[2], v_H3 = k->par[0].H_v[2], ro_H3 = k->par[0].H_n[2], p_H3 = 0.5 * k->par[0].H_T[2] * k->par[0].H_n[2];
		double u_H4 = k->par[0].H_u[3], v_H4 = k->par[0].H_v[3], ro_H4 = k->par[0].H_n[3], p_H4 = 0.5 * k->par[0].H_T[3] * k->par[0].H_n[3];
		u = k->par[0].u;
		v = k->par[0].v;
		ro = k->par[0].ro;
		p = k->par[0].p;

		if (ro <= 0.0)
		{
			ro = 0.0000001;
			p = 0.0;
		}
		if (ro_H1 <= 0.0)
		{
			ro_H1 = 0.0000001;
			p_H1 = 0.0;
		}
		if (ro_H2 <= 0.0)
		{
			ro_H2 = 0.0000001;
			p_H2 = 0.0;
		}
		if (ro_H3 <= 0.0)
		{
			ro_H3 = 0.0000001;
			p_H3 = 0.0;
		}
		if (ro_H4 <= 0.0)
		{
			ro_H4 = 0.0000001;
			p_H4 = 0.0;
		}

		U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H1 / ro_H1));
		U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H2 / ro_H2));
		U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H3 / ro_H3));
		U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H4 / ro_H4));

		U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H1 / ro_H1));
		U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H2 / ro_H2));
		U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H3 / ro_H3));
		U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H4 / ro_H4));

		sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
		sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
		sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
		sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

		nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
		nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
		nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
		nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

		k->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
			+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
		k->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
			+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));


		k->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
			(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
			nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
				(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
			nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
				(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
			nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
				(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));

		if (fabs(k->par[0].M_u) > 0.00001)
		{
			k->par[0].k_u = k->par[0].I_u / k->par[0].M_u;
			if (k->par[0].k_u > 3.0 || k->par[0].k_u < 0.33)
			{
				k->par[0].k_u = 1.0;
			}
		}
		else
		{
			k->par[0].k_u = 1.0;
		}

		if (fabs(k->par[0].M_v) > 0.00001)
		{
			k->par[0].k_v = k->par[0].I_v / k->par[0].M_v;
			if (k->par[0].k_v > 3.0 || k->par[0].k_v < 0.33)
			{
				k->par[0].k_v = 1.0;
			}
		}
		else
		{
			k->par[0].k_v = 1.0;
		}

		if (fabs(k->par[0].M_T) > 0.00001)
		{
			k->par[0].k_T = k->par[0].I_T / k->par[0].M_T;
			if (k->par[0].k_T > 3.0 || k->par[0].k_T < 0.33)
			{
				k->par[0].k_T = 1.0;
			}
		}
		else
		{
			k->par[0].k_T = 1.0;
		}
	}

	for (auto& k : this->Cell_other)
	{
		k->par[0].F_n = 1.0;
		k->par[0].F_u = Velosity_inf;
		k->par[0].F_v = 0.0;
		k->par[0].F_T = 1.0;
		k->par[0].H_n[0] = 0.0;
		k->par[0].H_u[0] = 0.0;
		k->par[0].H_v[0] = 0.0;
		k->par[0].H_T[0] = 0.0;
		k->par[0].H_n[1] = 0.0;
		k->par[0].H_u[1] = 0.0;
		k->par[0].H_v[1] = 0.0;
		k->par[0].H_T[1] = 0.0;
		k->par[0].H_n[2] = 0.0;
		k->par[0].H_u[2] = 0.0;
		k->par[0].H_v[2] = 0.0;
		k->par[0].H_T[2] = 0.0;
		k->par[0].H_n[3] = 1.0;
		k->par[0].H_u[3] = Velosity_inf;
		k->par[0].H_v[3] = 0.0;
		k->par[0].H_T[3] = 1.0;

		k->par[0].I_u = 0.0;
		k->par[0].I_v = 0.0;
		k->par[0].I_T = 0.0;
	}
}

void Setka::Fly_exchenge(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
	const double& mu_0, bool ExCh)
{
	//cout << "NEW" << endl;
	//cout  << x_0 << " " << sqrt(kvv(y_0, z_0, 0.0)) << " " << mu << endl;
	//double f1, f2;
	//now->Get_Center(f1, f2);
	//cout << "Cent = " << f1 << " " << f2 << endl;
	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = geo_accur * now->L / normV;
	//cout << "dt = " << dt << endl;


	// Нахождение времени до попадания в следующую ячейку
	//cout << "V = " << Vx << " " << Vy << " " << Vz << endl;
	while (now->belong(x_0 + 1000.0 * dt * Vx, sqrt(kv(y_0 + 1000.0 * dt * Vy) + kv(z_0 + 1000.0 * dt * Vz))) == false) // Слишком большой шаг
	{
		dt = dt / 2.0;
		//cout << "dt = " << dt << endl;
	}

	int a = 1000;
	int b = 10000;
	while (now->belong(x_0 + (1.0 * b) * dt * Vx, sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz))) == true)
	{
		b = b + 500;
		//cout << "b = " <<  x_0 + (1.0 * b) * dt * Vx << " " << sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz)) << " " << b <<  endl;
	}

	int k = (a + b) / 2;
	while (k != a)
	{
		if (now->belong(x_0 + (1.0 * k) * dt * Vx, sqrt(kv(y_0 + (1.0 * k) * dt * Vy) + kv(z_0 + (1.0 * k) * dt * Vz))) == true)
		{
			a = k;
		}
		else
		{
			b = k;
		}
		k = (a + b) / 2;
		//cout << a << " " << b << endl;
	}

	if (b != a + 1)  // Можно потом удалить проверку
	{
		cout << "ERRORIUEHUYCEUVDCC W" << endl;
		exit(-1);
	}

	double time = dt * (b); // Время нахождения в ячейке

	double uz, uz_M, uz_E, t1;// y, z, r;
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	bool change_was = false;
	double vx = now->par[0].u;
	double vy = now->par[0].v;
	double ro = now->par[0].ro;

	double t_ex = 0.0;
	double mu2 = mu;
	double nu_ex, u, alpha, mu_ex;
	double t2 = time;

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;
	double u1, u2, u3;

	if (ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		alpha = polar_angle(y_0, z_0);
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		nu_ex = (ro * uz * sigma2(uz, cp));
		double kappa = (nu_ex * time) / Kn_;
		t_ex = -(time / kappa) * log(1.0 - sens->MakeRandom() * (1.0 - exp(-kappa))); // Время до перезарядки
		t2 = t2 - t_ex;
		mu_ex = mu * (1.0 - exp(-kappa));
		mu2 = mu - mu_ex;
		alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);


		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			//cout << "TUT" << endl;
			x_ex = x_0 + (1.0 * a * dt) * Vx;  // Координаты перезарядки
			y_ex = y_0 + (1.0 * a * dt) * Vy;
			z_ex = z_0 + (1.0 * a * dt) * Vz;
			t_ex = (1.0 * a * dt);
			t2 = time - t_ex;
			alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
		}

		// Считаем источники для частицы до перезарядки
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
		uz_E = Velosity_3(u, cp);
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;


		now->mut.lock();
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].I_u += mu * t_ex * uz_M * uz * sigma(uz_M) * u1 / u;
		now->par[0].I_v += mu * t_ex * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		now->par[0].I_T += mu * t_ex * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
		now->mut.unlock();

		bool bb = false;
		if (mu_ex >= weight_ * mu_0)
		{
			//cout << mu_ex << " podhodit" << endl;
			bb = true;
		}
		else
		{
			if (mu_ex >= sens->MakeRandom() * weight_ * mu_0)
			{
				mu_ex = weight_ * mu_0;
				//cout << mu_ex << " ne podhodit = 0.25" << endl;
				bb = true;
			}
			else
			{
				bb = false;
				//cout << mu_ex << " zabili" << endl;
			}
		}

		if (bb)  // Иначе нужно забить на эти траектории
		{

			double alpha = polar_angle(y_ex, z_ex);

			double Vx2, Vy2, Vz2;

			Change_Velosity(sens, vx / cp, vy * cos(alpha) / cp, vy * sin(alpha) / cp,//
				Vx / cp, Vy / cp, Vz / cp, Vx2, Vy2, Vz2, cp);  // здесь подаются  r, theta, phi
			Vx2 = Vx2 * cp;
			Vy2 = Vy2 * cp;
			Vz2 = Vz2 * cp;

			//cout << "Change  " << endl;
			Fly_exchenge(sens, x_ex, y_ex, z_ex, Vx2, Vy2, Vz2, now, mu_ex, mu_0, true);
			//cout << "Vozvrat change" << endl;
		}
	}

	//cout << mu2 << " mu2" << endl;
	if (mu2 < weight_ * mu_0 && ExCh == false)
	{
		if (mu2 >= sens->MakeRandom() * weight_ * mu_0)
		{
			//cout << mu2 << " mu2 = 0.25" << endl;
			mu2 = weight_ * mu_0;
		}
		else
		{
			//cout << mu2 << " zabili" << endl;
			return;
		}
	}

	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
	uz = Velosity_1(u, cp);
	uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
	uz_E = Velosity_3(u, cp);
	u1 = vx - Vx;
	u2 = vy * cos(alpha) - Vy;
	u3 = vy * sin(alpha) - Vz;
	double skalar = Vx * u1 + Vy * u2 + Vz * u3;

	now->mut.lock();
	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].I_u += mu2 * t2 * uz_M * uz * sigma(uz_M) * u1 / u;
	now->par[0].I_v += mu2 * t2 * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
	now->par[0].I_T += mu2 * t2 * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
	now->mut.unlock();

	double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;


	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + (1.0 * b) * dt * Vx;
	yk = sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz));
	if (xk < Left_ || (sqrt(kv(xk) + kv(yk)) > R5_ && xk > 0) || (xk < 0 && yk > R5_))
	{
		next = nullptr;
	}
	else
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
		/*if (next == nullptr)
		{
			cout << xk << " " << yk << " " << "no sosed" << endl;
		}*/
	}

	if (next != nullptr)
	{
		//cout << "Next cell" << endl;
		Fly_exchenge(sens, X, Y, Z, Vx, Vy, Vz, next, mu2, mu_0, false);
		//cout << "Vozvrat next" << endl;
	}
	return;

	// Попытка сделать честно, удалять не стал пока
	//for (auto& i : now->Grans)
	//{
	//	a = i->a;
	//	b = i->b;
	//	A = -2.0 * a * b * Vx - 2.0 * kv(a) * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0;
	//	C = 2.0 * (kv(a) * kv(Vx) - kv(Vy) - kv(Vz));
	//	B = kv(2.0 * a * b * Vx + 2.0 * kv(a) * Vx * x_0 - 2.0 * Vy * y_0 - 2.0 * Vz * z_0) - 4.0 * (kv(a * Vx) - kv(Vy) - kv(Vz)) * //
	//		(kv(b) + 2.0 * a * b * x_0 + kv(a) * kv(x_0) - kv(y_0) - kv(z_0));
	//	if (fabs(C) > 0.000001)
	//	{
	//		if (B >= 0)
	//		{
	//			t1 = (A - sqrt(B)) / C;
	//			if (t1 > 0 && i->belong_gran(x_0 + t1 * Vx, sqrt(kv(y_0 + t1 * Vy) + kv(z_0 + t1 * Vz))))
	//			{
	//				if (t1 < t)
	//				{
	//					t = t1;

	//				}
	//			}
	//			t2 = (A + sqrt(B)) / C;
	//		}
	//	}
	//}
}

void Setka::Fly_exchenge_Imit_Korol(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
	int area, bool ExCh, const double& mu_start, int to_I = 0, int to_J = 0, bool georaschep = true)
// Новая версия полёта нейтрального атома. Функция без лишнего кода и насыщена комментариями
// Функция использует и геометрическое и физическое расщепление!!!
// Описание переменных:
// sens - датчик случайных чисел
// x_0 ... -  координаты начала полёта нейтрального атома в данной ячейке
// Vx ... -   скорость летящего нейтрального атома
// now -      текущая ячейка
// mu -       вес нейтрального атома
// area -     в какой области летит атом
// mu_start - начальный вес атома при инициализации (чтобы понимать, сколько от исходного осталось)
// to_I - в какую область предназначался атом
{
	Cell* next_mb = nullptr;
	// Алгоритм вырубания лишних траектори (траектории вырубаются при пересечении новой зоны, далее 
	// вес лмбо увеличивается, либо частица уничтожается)
	// Не уверен, что to_I бывает больше, чем now->zona, кажется не должно такого быть
	if (true && to_I < now->zona)
	{
		if (x_0 * Vx + y_0 * Vy + z_0 * Vz > 0.0 && to_I < now->zona) // Если летит вне солнца
		{
			double all = polar_angle(x_0, sqrt(kv(y_0) + kv(z_0)));
			to_I = now->zona;
			if (fabs(mu) < Mu[area][now->zona] * mu_start * 0.3 * sin(all))
			{
				if (fabs(mu) >= sens->MakeRandom() * Mu[area][now->zona] * mu_start * 0.3 * sin(all))
				{
					mu = Mu[area][now->zona] * mu_start * sign(mu) * 0.3 * sin(all);
				}
				else
				{
					return;
				}
			}
		}
	}

	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / normV;  // Время на один мини-шаг в ячейке (планируется сделать много шагов в ячейке)
	double ddt = 0.001 * dt;
	// Нахождение времени до попадания в следующую ячейку
	//cout << "V = " << Vx << " " << Vy << " " << Vz << endl;
	//cout << "A" << endl;
	double a = 0.0, b = dt;
	double t_per = 1000000000.0;
	// Блок проверки перигелия к оси
	/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
	if (t_per > 0.0000001)
	{
		if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
		{
			t_per = -1.0;
		}
	}
	else
	{
		t_per = -1.0;
	}*/

	//cout << "A" << endl;
	double tt1, tt2;
	double a1, a2, a3;

	// Пытаемся определить время выхода точно
	if (true)
	{
		double A, B, C;
		for (auto& i : now->Grans)
		{
			if (i->type == Axis)
			{
				continue;
			}
			A = i->aa;
			B = i->bb;
			C = i->cc;

			/*if (fabs(B) < 0.001)
			{
				tt1 = (-C / A - x_0) / Vx;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				continue;
			}*/
			/*cout << "1 - " << A << " " << B << " " << C << endl;
			cout << "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
			cout << "3 - " << Vx << " " << Vy << " " << Vz << endl;*/
			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
			//cout << "4 - " << a1 << " " << a2 << " " << a3 << endl;
			if (fabs(a1) < 0.0000001)
			{
				a1 = 0.0;
			}
			if (a1 >= 0.0)
			{
				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
				//cout << A << " " << B << " " << C << endl;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.01 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (t_per > tt2)
					{
						if (now->belong(x_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vx, sqrt(kv(y_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vy) + kv(z_0 + (tt2 + min(ddt, 0.01 * tt2)) * Vz))) == false)
						{
							t_per = tt2;
							next_mb = i->Sosed;
						}
					}
				}
			}
		}

		a = t_per - min(ddt, 0.005 * t_per);
		b = t_per + min(ddt, 0.005 * t_per);
	}

	bool hand = false;
	if (now->belong(x_0 + a * Vx, sqrt(kv(y_0 + a * Vy) + kv(z_0 + a * Vz))) == false)
	{
		hand = true;
	}
	else if (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true)
	{
		hand = true;
	}
	else if (t_per > 1000000.0)
	{
		//cout << "ERROR  8568  " << t_per << endl;
		hand = true;
	}

	double k;

	// Находим время  в ячейке
	int lk = 0;
	if (hand)
	{
		a = 0.0;
		b = dt;
		t_per = 1000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
		{
			lk++;
			if (lk > 1000)
			{
				cout << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 300.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; // Время нахождения в ячейке


	bool slay = false;
	double t1, t2;
	//Нужно проверить пересечение со сферой
	double D = kv(2.0 * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0) - 4.0 * kvv(Vx, Vy, Vz) * (kvv(x_0, y_0, z_0) - kv(Rmax_));
	if (D >= 0.0)
	{
		t1 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) - sqrt(D)) / (2.0 * kv(normV));
		t2 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) + sqrt(D)) / (2.0 * kv(normV));
		if (t1 > 0.0 && t1 < time)
		{
			if (x_0 + Vx * t1 >= 0.0)
			{
				time = t1;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
		if (t2 > 0.0 && t2 < time)
		{
			if (x_0 + Vx * t2 >= 0.0)  // Для правой полусферы
			{
				time = t2;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
	}
	// Теперь пересечение с правой границей
	if (fabs(Vx) > 0.000001)
	{
		t1 = -x_0 / Vx;
		if (t1 > 0.0 && t1 <= time && kvv(0.0, y_0 + t1 * Vy, z_0 + t1 * Vz) >= kv(Rmax_))
		{
			time = t1;
			a = time * 0.999;
			b = time * 1.001;
			slay = true;
		}
	}

	//Мини - проверка
	//for (double kk = 0.0; kk < a ; kk += time / 50.0)
	//{
	//	if (now->belong(x_0 + kk * Vx, sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz))) == false)
	//	{
	//		cout << "ERROR  8280 iruhfuerfefr  " << x_0 + kk * Vx << " " << sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz)) << " " << now->contour[0]->x << " " << now->contour[0]->y << " " << kk << " " << time << endl;
	//		cout << a << " " << b << endl;
	//		cout << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
	//		cout << x_0 + 0.01 * Vx << " " << sqrt(kv(y_0 + 0.01 * Vy) + kv(z_0 + 0.01 * Vz)) << endl;
	//		double xx, yy;
	//		now->Get_Center(xx, yy);
	//		cout << xx << " " << yy << endl;
	//		cout << " FFF " << endl;
	//		for (auto& i : now->Grans)
	//		{
	//			if (i->type == Axis)
	//			{
	//				continue;
	//			}
	//			double A = i->aa;
	//			double B = i->bb;
	//			double C = i->cc;
	//			cout << "1 - " << A << " " << B << " " << C << endl;
	//			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
	//				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
	//			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
	//			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
	//			cout << "4 - " << a1 << " " << a2 << " " << a3 << endl;
	//			if (fabs(a1) < 0.000000001)
	//			{
	//				a1 = 0.0;
	//			}
	//			if (a1 >= 0.0)
	//			{
	//				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
	//				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
	//				cout << "tt1 = " << tt1 << endl;
	//				cout << "tt2 = " << tt2 << endl;
	//			}
	//		}
	//		exit(-1);
	//	}
	//}
	
	// теперь время time, a, b уже определены
	// a - последний момент времени в ячейке, b - первый момент за пределами ячейки


	//cout << a << " " << b << endl;
	double uz, uz_M, uz_E;								// средние скорости в интегралах
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	double vx = now->par[0].u;							// Скорости плазмы в ячейке
	double vy = now->par[0].v;
	double ro = now->par[0].ro;

	double t_ex = 0.0;									// время до перезарядки
	// time - время нахождения атома в ячейке
	t2 = time;									// время после перезарядки (будет ниже)

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;		    // координаты перезарядки
	double u1, u2, u3;

	double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));                    // Расстояние до выхода из ячейки

	double alpha = 0.0;								    // угол 
	double u = 0.0;                                     // модуль относительной скорости атома и плазмы
	double nu_ex = 0.0;								    // частота перезарядки
	double mu_ex = 0.0;									// вес перезаряженного атома
	double mu2 = mu;									// вес не-перезаряженного атома

	if (true)//(ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		double kappa = 0.0;

		if (true)
		{
			for (double tt = 0.0; tt < time * 0.99; tt += time / 5.0)
			{
				alpha = polar_angle(y_0 + (tt + time / 10.0) * Vy, z_0 + (tt + time / 10.0) * Vz);
				u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
				uz = Velosity_1(u, cp);
				nu_ex = (ro * uz * sigma(uz)) / Kn_;
				kappa += (nu_ex * time / 5.0);
			}
		}
		else
		{
			alpha = polar_angle(y_0, z_0);
			u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
			uz = Velosity_1(u, cp);
			//double IS = MK.Int_cp_1(u / cp);
			nu_ex = (ro * uz * sigma(uz)) / Kn_;
			//nu_ex = (ro * IS) / Kn_;
			//kappa = (nu_ex * time) / Kn_;
			kappa = (nu_ex * time);
		}
		

		t_ex = -(time / kappa) * log(1.0 - sens->MakeRandom() * (1.0 - exp(-kappa))); // Время до перезарядки
		t2 = time - t_ex;   // Время сколько лететь после того, как атом перезарядился
		mu_ex = mu * (1.0 - exp(-kappa)); // вес перезаряженного атома
		mu2 = mu - mu_ex;                 // вес оставшегося неперезаряженного атома

		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;

		alpha = polar_angle(y_ex, z_ex);

		// проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			x_ex = x_0 + a * Vx;  // Координаты перезарядки
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			t2 = time - t_ex;
			alpha = polar_angle(y_ex, z_ex);
		}

		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		/*if (r > Rmax_ && x_ex > 0.0011)
		{
			cout << "8609 uhfuwvgfwdwdd" << endl;
			exit(-1);
			return;
		}*/
		double rr = sqrt(kvv(0.0, y_ex, z_ex));
		double all = polar_angle(x_ex, sqrt(kv(y_ex) + kv(z_ex)));

		// Считаем источники для частицы до перезарядки
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
		uz_E = Velosity_3(u, cp);
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		// определим область перезарядки (по геометрическим зонам)
		int i_z = 0, i_alp = 0;
		i_z = geo_zones(r);
		i_alp = alpha_zones(x_ex, rr);

		/*if (to_I == i_z && to_J == i_alp && georaschep == true)
		{
			this->m_m.lock();
			this->Mu_stat[area][i_z][i_alp] += (mu / max(sin(alpha), 0.3 * Sinus[i_alp])) / mu_start;
			georaschep = false;
			this->m_m.unlock();
		}*/

		now->mut.lock();

		//now->par[0].w_m[area] += mu / max(sin(alpha), 0.3 * Sinus[i_alp]);
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].H_n[area] += t_ex * mu;
		now->par[0].H_u[area] += t_ex * Vx * mu;
		now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

		now->par[0].I_u += -mu_ex * uz_M * u1 / u;
		now->par[0].I_v += -mu_ex * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
			uz_M * skalar / u);

		if (true)
		{
			for (double tt = 0.0; tt < t_ex * 0.99; tt += t_ex / 20.0)
			{
				alpha = polar_angle(y_0 + (tt + t_ex / 40.0) * Vy, z_0 + (tt + t_ex / 40.0) * Vz);

				now->par[0].F_v += (t_ex / 20.0) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
				now->par[0].H_v[area] += (t_ex / 20.0) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			}
		}
		else
		{
			now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		}

		//now->par[0].I_u += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu_ex * u1;
		now->par[0].II_v += -mu_ex * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu_ex * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/

		now->mut.unlock();

		int area2 = 0;
		// определим область в которой находится атом сейчас (это параметр самой ячейки)
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}
     
		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		double mu3;

		if (area2 == 0 || Ur/cp > 3.0)   // Без геометрического расщепления
		{
			bool kj = true;
			mu3 = mu_ex;
			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;  // Декартовы скорости после перезарядки
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);
			MK.Change_Velosity(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;
			dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
			double peregel, time_do_peregel;
			int ii_z;
			int ii_alp;

			if (Wr[0] < 0.0)
			{
				//wwt = sqrt(kv(Wphi[0]) + kv(Wthe[0]));
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex)/kvv(aa,bb,cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
			}
			else
			{
				ii_z = i_z;
				ii_alp = i_alp;
			}

			if (fabs(mu3) >= Mu[area2][ii_z] * mu_start * 0.3 *  sin(all))
			{
				kj = true;
			}
			else
			{
				if (fabs(mu3) >= sens->MakeRandom() * Mu[area2][ii_z] * mu_start * 0.3 * sin(all))
				{
					mu3 = Mu[area2][ii_z] * mu_start * sign(mu3) * 0.3 * sin(all);
					kj = true;
				}
				else
				{
					kj = false;
				}
			}

			if (kj == true)
			{
				now->mut.lock();
				now->par[0].II_u += mu3 * (Vx - aa);
				now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
				now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
				now->mut.unlock();
				Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true);
			}
			//this->Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, -log(1.0 - sens->MakeRandom()),//
			//	0.0, area2, mu_start, to_I, 0);
			
		}
		else  // делаем геометрическое расщепление
		{
			int I = geo_zones(r, 1.2);  // Число дополнительных траекторий
			if (I > to_I)  // Для того чтобы не расщеплять атом, который летел в to_I область на атомы из более крупных областей
			{
				I = to_I;
			}

			vector <double> Wr(I + 1);
			vector <double> Wthe(I + 1);
			vector <double> Wphi(I + 1);
			vector <double> mu_(I + 1);

			//cout << " C " << endl;
			bool bbb = MK.Change_Velosity(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
			//cout << " D " << endl;

			for (int i = 0; i <= I; i++)
			{
				Wr[i] = Wr[i] * cp;
				Wthe[i] = Wthe[i] * cp;
				Wphi[i] = Wphi[i] * cp;
			}

			double aa, bb, cc, time_do_peregel, peregel;
			int ii_z, ii_alp;
			bool kj = true;
			//now->par[0].num_atoms++;
			for (int i = 0; i < I; i++)
			{
				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
				mu3 = mu_[i] * mu_ex;
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));

				if (fabs(mu3) >= Mu[area2][ii_z] * mu_start * 0.3 * sin(all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= sens->MakeRandom() * Mu[area2][ii_z] * mu_start * 0.3 * sin(all))
					{
						mu3 = Mu[area2][ii_z] * mu_start * 0.3 * sin(all);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();
					Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true);
				}
			}

			if (bbb == true)  // Если вообще нужно запускать основной атом
			{
				mu3 = mu_[I] * mu_ex;

				// Нужно определить номер основного атома (от его перегелия)
				dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
				if (Wr[I] < 0.0)
				{
					time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
					peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
					ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
					ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				}
				else
				{
					ii_z = i_z;
					ii_alp = i_alp;
				}

				if (fabs(mu3) >= Mu[area2][ii_z] * mu_start * 0.3 * sin(all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= sens->MakeRandom() * Mu[area2][ii_z] * mu_start * 0.3 * sin(all))
					{
						mu3 = Mu[area2][ii_z] * mu_start * 0.3 * sin(all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}

				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();
					Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, ii_z, ii_alp, false);
				}
			}
			
		}


		// Здесь можно добавить рулетку для оставшегося неперезаряженного атома
		if (fabs(mu2) < Mu[area][to_I] * mu_start * 0.3 * sin(all))
		{
			if (fabs(mu2) >= sens->MakeRandom() * Mu[area][to_I] * mu_start * 0.3 * sin(all))
			{
				mu2 = Mu[area][to_I] * mu_start * 0.3 * sin(all) * sign(mu2);
			}
			else
			{
				return;
			}
		}

	}


	// этот блок кода будет работать и для частиц, которые не попали в предыдущий цикл
	// здесь t2, _ex будут иметь исходные значения, как при инициализации
	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	now->mut.lock();
	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[area] += t2 * mu2;
	now->par[0].H_u[area] += t2 * Vx * mu2;
	now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;
	if (true)
	{
		for (double tt = 0.0; tt < t2 * 0.99; tt += t2 / 20.0)
		{
			alpha = polar_angle(y_ex + (tt + t2 / 40.0) * Vy, z_ex + (tt + t2 / 40.0) * Vz);

			now->par[0].F_v += (t2 / 20.0) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
			now->par[0].H_v[area] += (t2 / 20.0) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		}
	}
	else
	{
		now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	}
	now->mut.unlock();

	double mu3 = mu2;
	/*double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;*/
	double X = x_0 + b * Vx;
	double Y = y_0 + b * Vy;
	double Z = z_0 + b * Vz;

	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	/*if (yk < 5.0)
	{
		cout << b << "  " << xk << endl;
	}*/

	if (slay == true)  //Если надо убить траекторию
	{
		return;
	}

	if (xk < Left_ )
	{
		return;
	}

	if (yk > R5_)
	{
		return;
	}

	if (r_k >= Rmax_ - 3.0 && xk >= -0.001)
	{
		return;
		next = nullptr;
	}

	if (next_mb != nullptr)
	{
		if (next_mb->belong(xk, yk) == true)
		{
			next = next_mb;
		}
	}

	if(next == nullptr)
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next == nullptr)
	{
		cout << "Long   find" << endl;
		//cout << xk << " " << yk << endl;
		//cout << x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << endl;
		//t_per = 100000000.0;
		//if (true)
		//{
		//	double A, B, C;
		//	for (auto& i : now->Grans)
		//	{
		//		double xx, yy;
		//		i->Get_Center(xx, yy);
		//		cout << "centr = " << xx << " " << yy << endl;
		//		cout << t_per << endl;
		//		if (i->type == Axis)
		//		{
		//			continue;
		//		}
		//		A = i->aa;
		//		B = i->bb;
		//		C = i->cc;
		//		if (fabs(B) < 0.00001)
		//		{
		//			tt1 = (-C / A - x_0) / Vx;
		//			cout << "TT1 = " << tt1 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.02 * tt1 * Vx, sqrt(kv(y_0 + 1.02 * tt1 * Vy) + kv(z_0 + 1.02 * tt1 * Vz))) == false)
		//					{
		//						t_per = tt1;
		//						next_mb = i->Sosed;
		//					}
		//				}
		//			}
		//			continue;
		//		}
		//		cout << "1 - " << A << " " << B << " " << C << endl;
		//		//cout << "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
		//		//cout << "3 - " << Vx << " " << Vy << " " << Vz << endl;
		//		a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
		//			B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
		//		a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
		//		a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
		//		cout << "4 - " << a1 << " " << a2 << " " << a3 << endl;
		//		if (fabs(a1) < 0.000001)
		//		{
		//			a1 = 0.0;
		//		}
		//		if (a1 >= 0.0)
		//		{
		//			tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
		//			tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
		//			//cout << A << " " << B << " " << C << endl;
		//			cout << tt1 << " " << tt2 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt1 * Vx, sqrt(kv(y_0 + 1.03 * tt1 * Vy) + kv(z_0 + 1.03 * tt1 * Vz))) == false)
		//					{
		//						cout << "A " << endl;
		//						t_per = tt1;
		//					}
		//				}
		//			}
		//			if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
		//			{
		//				if (t_per > tt2)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt2 * Vx, sqrt(kv(y_0 + 1.03 * tt2 * Vy) + kv(z_0 + 1.03 * tt2 * Vz))) == false)
		//					{
		//						cout << "B " << endl;
		//						t_per = tt2;
		//					}
		//				}
		//			}
		//		}
		//	}
		//	cout << t_per << endl;
		//	xk = x_0 + t_per * Vx;
		//	yk = sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz));
		//	cout << xk << " " << yk << endl;

		//	a = t_per * 0.99;
		//	b = t_per * 1.01;
		//}




		//exit(-1);
		for (auto& i : this->All_Cells)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next != nullptr)
	{
		Fly_exchenge_Imit_Korol(MK, sens, X, Y, Z, Vx, Vy, Vz, next, mu3, area, false, mu_start, to_I, to_J, georaschep);
	}
	else
	{

		if (xk > -1000.0 && xk < 1000.0 && yk < 1600.0)
		{
			cout << "ERROR  8670  poteryal atom" << endl;
			cout << xk << " " << yk << endl;
		}
	}

	return;
}

void Setka::Fly_exchenge_Imit(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, double KSI, //
	double I_do, int area, const double& mu_start, int to_I = 0, int iii = 0)
	// Имитационный метод
{
	//cout << "Start  " << iii <<  endl;
	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / (normV);  // Первое приближение времени шага

	// Алгоритм вырубания лишних траектори
	/*double Vr, Vthe, Vphi;
	spherical_skorost(y_0, z_0, x_0, Vy, Vz, Vx, Vr, Vphi, Vthe);
	if (Vr > 0.0 && to_I < now->zona)
	{
		to_I = now->zona;
		if (mu < Mu[area][now->zona])
		{
			if (mu >= sens->MakeRandom() * Mu[area][now->zona])
			{
				mu = Mu[area][now->zona];
			}
			else
			{
				return;
			}
		}
	}*/


	// Нахождение времени до попадания в следующую ячейку
	//cout << "V = " << Vx << " " << Vy << " " << Vz << endl;
	//cout << "A" << endl;
	double a = 0.0, b = dt;
	double t_per = 1000000000.0;
	// Блок проверки перигелия к оси
	/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
	if (t_per > 0.0000001)
	{
		if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
		{
			t_per = -1.0;
		}
	}
	else
	{
		t_per = -1.0;
	}*/

	double tt1, tt2;
	double a1, a2, a3;
	double k;
	// Находим время  в ячейке
	int lk = 0;
	if (true)
	{
		a = 0.0;
		b = dt;
		t_per = 1000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
		{
			lk++;
			if (lk > 1000)
			{
				cout << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 200.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; // Время нахождения в ячейке

	/*if (iii == 219)
	{
		cout << "a, b, time, dt,  a, now->L  =   " << a << "   " << b << " " << time << " " << dt << " " << now->L << " " << endl;
	}*/
	//cout << a << " " << b << " " << time << " " << dt << " " << now->L << " " << endl;

	double uz, uz_M, uz_E, t1;// y, z, r;
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	bool change_was = false;
	double vx = now->par[0].u;
	double vy = now->par[0].v;
	double ro = now->par[0].ro;

	double t_ex = 0.0;
	double mu2 = mu;
	double mu3 = 0.0;
	double nu_ex, u, alpha;
	double t2 = time;

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;
	double u1, u2, u3;

	int nk = 100;
	double sig = 0.0;
	double l = 0.0;
	double I = I_do;
	l = sqrt(kvv(time / nk * Vx, time / nk * Vy, time / nk * Vz));
	for (double tt = 0.0; tt < time * 0.9999; tt += time / nk)
	{
		alpha = polar_angle(y_0 + (tt + time / (2.0 * nk)) * Vy, z_0 + (tt + time / (2.0 * nk)) * Vz);
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		sig = sqrt(kvv(Vx, Vy, Vz)) / ((1.0 / Kn_) * ro * uz * sigma(uz));
		I += l / sig;
	}

	//double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));
	//alpha = polar_angle(y_0, z_0);
	//u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
	//uz = Velosity_1(u, cp);
	////double IS = MK.Int_cp_1(u / cp);
	////sig = sqrt(kvv(Vx, Vy, Vz)) / ((1.0 / Kn_) * ro * IS);
	//sig = sqrt(kvv(Vx, Vy, Vz)) / ((1.0 /Kn_) * ro * uz * sigma(uz));
	//double I = I_do + l / sig;



	if (I < KSI)  // Не произошла перезарядка
	{
		I_do = I;
	}
	else  // Перезарядка была
	{
		//cout << "C" << endl;

		double ksi = (KSI - I_do) * sig;
		t_ex = ksi / sqrt(kvv(Vx, Vy, Vz));
		I_do = 0.0;
		//KSI = -log(1.0 - sens->MakeRandom());
		alpha = polar_angle(y_0 + 0.5 * t_ex * Vy, z_0 + 0.5 * t_ex * Vz);  // Гарантирует расчёт угла в середине пути

		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			//cout << "TUT" << endl;
			x_ex = x_0 + a * Vx;  // Координаты перезарядки
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
		}


		// Считаем источники для частицы до перезарядки
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
		uz_E = Velosity_3(u, cp);
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		now->mut.lock();
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].H_n[area] += t_ex * mu;
		now->par[0].H_u[area] += t_ex * Vx * mu;
		now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

		now->par[0].I_u += -mu * uz_M * u1 / u;
		now->par[0].I_v += -mu * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		now->par[0].I_T += mu * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
			uz_M * skalar / u);

		//now->par[0].I_u += -mu * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu * u1;
		now->par[0].II_v += -mu * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/
		now->mut.unlock();

		//cout << "SDDD" << endl;
		int area2;
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}

		//cout << "SDDD 2" << endl;
		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		double rr = sqrt(kvv(0.0, y_ex, z_ex));
		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		//cout << "E" << endl;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		//double uu = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi)) / cp;
		//double UUr = (Ur - Vr) / cp;
		//double UUa = sqrt(kvv(0.0, Vthe - Uthe, Vphi - Uphi)) / cp;

		//f_way << x_ex << " " << sqrt(kvv(z_ex, y_ex, 0.0)) << " " << 15 << endl;
		//f_num++;


		if (r > Rmax_ && x_ex > 0.001)
		{
			return;
		}

		int i_z = geo_zones(r);

		// для uu > 5 значение интеграллов не вычисленно
		// для (UUr >= -5.0 && UUr <= 0.0)


		if (true)//(area2 == 0)   // Без расщепления
		{
			int num = 1;
			bool kj = true;
			mu3 = mu;

			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);

			MK.Change_Velosity(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;

			/*if (area2 == 1 && mu2 > 0.01)
			{
				cout << mu3 << " " << area << " " << Vr << " " << Wr[0] << " " << x_ex << " " << sqrt(kv(y_ex) + kv(z_ex)) << endl;
			}*/

			double wwt, peregel;
			int ii_z;
			bool rulet = true;

			if (Wr[0] < 0.0)
			{
				//cout << Wr[0] << " " << area << " " << mu << endl;
				rulet = false;
				wwt = sqrt(kv(Wphi[0]) + kv(Wthe[0]));
				peregel = r * wwt / (sqrt(kv(wwt) + kv(Wr[0])));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
			}
			else
			{
				rulet = true;
				ii_z = i_z;
				//if (ii_z == 2 || ii_z == 1)  // Слишком рано вырубать маленькие атомы не очень хорошо
				//{
				//	ii_z = 0;
				//}
			}

			if (true)
			{
				double ko = 1.0;

				/*if (area2 == 0 || area2 == 1)
				{
					ko = 0.001;
				}*/

				if (mu3 >= Mu[area2][ii_z] * ko)
				{
					kj = true;
				}
				else
				{
					if (mu3 >= sens->MakeRandom() * Mu[area2][ii_z] * ko)
					{
						mu3 = Mu[area2][ii_z] * ko;
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
			}

			if (kj == true)
			{
				dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
				KSI = -log(1.0 - sens->MakeRandom());
				Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, KSI, I_do, area2, mu_start, ii_z, iii);
			}
			
			return;
		}
		else   // С расщеплением
		{
			double alpha = polar_angle(y_ex, z_ex);

			int I = geo_zones(r, 1.5);  // Число дополнительных траекторий
			
			if (I > to_I + 1)   // Не расщепляем атомы, которые летели в свою область на более другие области, которые больше своей
			{
				I = to_I + 1;
			}

			if (I > 1)
			{
				I = 1;
			}

			int num = 1;
			/*if (rr < 50.0)
			{
				num = 3;
			}*/
			mu = mu / num;
			for (int kkk = 0; kkk < num; kkk++)
			{
				vector <double> Wr(I + 1);
				vector <double> Wthe(I + 1);
				vector <double> Wphi(I + 1);
				vector <double> mu_(I + 1);

				//bool bbb = true;
				bool bbb = MK.Change_Velosity(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I, x_ex, y_ex, z_ex);
				//Change_Velosity_Split(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I, x_ex, y_ex, z_ex);
				double aa, bb, cc;
				for (int i = 0; i <= I; i++)
				{
					Wr[i] = Wr[i] * cp;
					Wthe[i] = Wthe[i] * cp;
					Wphi[i] = Wphi[i] * cp;
				}


				// Удобный блок для рассмотрения траекторий расщеплённого атома в Техплоте
				/*if (area == 3 && I >= 5 && bbb == true)
				{
					cout << "to_I = " << to_I << endl;
					for (int i = 0; i < I; i++)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
						for (int ij = 0; ij < 27000; ij++)
						{
							cout << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << i << endl;
						}
					}
					if (bbb == true)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
						for (int ij = 0; ij < 30000; ij++)
						{
							cout << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << I << endl;
						}
					}
					exit(-1);
				}*/

				bool kj = true;
				double ko = 1.0;

				/*if (area2 == 0 || area2 == 1)
				{
					ko = 0.001;
				}*/

				/*if (area2 == 1)
				{
					for (int i = 0; i < I; i++)
					{
						cout << mu_[i] << endl;
					}
					exit(-1);
				}*/

				for (int i = 0; i < I; i++)
				{
					kj = true;
					mu3 = mu_[i] * mu;

					if (true)//(i != to_I)
					{
						if (fabs(mu3) >= Mu[area2][i] * ko / num)
						{
							kj = true;
						}
						else
						{
							if (fabs(mu3) >= sens->MakeRandom() * Mu[area2][i] * ko / num)
							{
								mu3 = Mu[area2][i] * ko * sign(mu3) / num;
								kj = true;
							}
							else
							{
								kj = false;
							}
						}
					}
					
					if (kj == true)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
						KSI = -log(1.0 - sens->MakeRandom());
						Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, KSI, 0.0, area2, mu_start / num, i, iii);
						//Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, i);
					}
				}

				double wwt, peregel;
				int ii_z;

				if (Wr[I] < 0.0)
				{
					wwt = sqrt(kv(Wphi[I]) + kv(Wthe[I]));
					peregel = r * wwt / (sqrt(kv(wwt) + kv(Wr[I])));
					ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				}
				else
				{
					ii_z = i_z;
				}

				kj = true;
				if (bbb == true)
				{

					mu3 = mu_[I] * mu;
					if (true)
					{
						if (fabs(mu3) >= Mu[area2][ii_z] * mu_start * ko / num)
						{
							kj = true;
						}
						else
						{
							if (fabs(mu3) >= sens->MakeRandom() * Mu[area2][ii_z] * mu_start * ko / num)
							{
								mu3 = Mu[area2][ii_z] * mu_start * ko * sign(mu3) / num;
								kj = true;
							}
							else
							{
								kj = false;
							}
						}
					}


					if (kj == true)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
						KSI = -log(1.0 - sens->MakeRandom());
						Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, KSI, 0.0, area2, mu_start / num, ii_z, iii);
						//Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, to_I);
					}
				}
			}

			//dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
			//KSI = -log(1.0 - sens->MakeRandom());
			//Fly_exchenge_Imit(sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu_[I] * mu, KSI, I_do, area, mu_start, to_I, iii);
			return;
		}
	}


	//cout << "B1" << endl;
	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	//cout << "B2" << endl;
	now->mut.lock();
	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[area] += t2 * mu2;
	now->par[0].H_u[area] += t2 * Vx * mu2;
	now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->mut.unlock();

	double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;


	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	//f_way << xk << " " << yk << " " << to_I << endl;
	//f_num++;

	if (r_0 <= Rmax_ && r_k >= Rmax_ && x_0 > 0.0)
	{
		//cout << x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << " " << xk << " " << yk << endl;
		next = nullptr;
		return;
	}
	else if (r_k >= Rmax_ - 0.001 && xk >= -0.001)
	{
		//cout << x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << " " << xk << " " << yk << endl;
		next = nullptr;
		return;
	}
	else
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next != nullptr)
	{
		Fly_exchenge_Imit(MK, sens, X, Y, Z, Vx, Vy, Vz, next, mu2, KSI, I_do, area, mu_start, to_I, iii);
	}
	else
	{
		if (xk > -1440.0 && xk < 1200.0 && yk < 200.0) // Проверка, если траектория пропадёт не на границе а в нуле скорее всего
		{
			cout << xk << " " << yk << endl;
		}
	}

	//cout << "Vihod" << endl;
	return;

}

//
//void Setka::Fly_exchenge_Split(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
//	const double& mu_0, bool ExCh, int zone)
//{
//	//cout << "NEW" << endl;
//	//cout  << x_0 << " " << sqrt(kvv(y_0, z_0, 0.0)) << " " << mu_0 << endl;
//	double f1, f2;
//	now->Get_Center(f1, f2);
//	//cout << "Cent = " << f1 << " " << f2 << endl;
//	double normV = sqrt(kvv(Vx, Vy, Vz));
//	double dt = geo_accur * now->L / normV;
//	//cout << "dt = " << dt << endl;
//
//
//	// Нахождение времени до попадания в следующую ячейку
//	//cout << "V = " << Vx << " " << Vy << " " << Vz << endl;
//	while (now->belong(x_0 + 100.0 * dt * Vx, sqrt(kv(y_0 + 100.0 * dt * Vy) + kv(z_0 + 100.0 * dt * Vz))) == false) // Слишком большой шаг
//	{
//		dt = dt / 2.0;
//		//cout << "dt = " << dt << endl;
//	}
//
//	int a = 100;
//	int b = 1000;
//	while (now->belong(x_0 + (1.0 * b) * dt * Vx, sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz))) == true)
//	{
//		b = b + 50;
//		//cout << "b = " <<  x_0 + (1.0 * b) * dt * Vx << " " << sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz)) << " " << b <<  endl;
//	}
//
//	int k = (a + b) / 2;
//	while (k != a)
//	{
//		if (now->belong(x_0 + (1.0 * k) * dt * Vx, sqrt(kv(y_0 + (1.0 * k) * dt * Vy) + kv(z_0 + (1.0 * k) * dt * Vz))) == true)
//		{
//			a = k;
//		}
//		else
//		{
//			b = k;
//		}
//		k = (a + b) / 2;
//		//cout << a << " " << b << endl;
//	}
//	
//	Smut.lock();
//	//cout << b << endl;
//	mmu1 += 1.0 * b;
//	mn1++;
//	Smut.unlock();
//
//	if (b != a + 1)  // Можно потом удалить проверку
//	{
//		cout << "ERRORIUEHUYCEUVDCC W" << endl;
//		exit(-1);
//	}
//
//	double time = dt * (b); // Время нахождения в ячейке
//
//	double uz, uz_M, uz_E, t1;// y, z, r;
//	double cp = sqrt(now->par[0].p / now->par[0].ro);
//	bool change_was = false;
//	double vx = now->par[0].u;
//	double vy = now->par[0].v;
//	double ro = now->par[0].ro;
//
//	double t_ex = 0.0;
//	double mu2 = mu;
//	double nu_ex, u, alpha, mu_ex;
//	double t2 = time;
//
//	double x_ex = x_0, y_ex = y_0, z_ex = z_0;
//	double u1, u2, u3;
//
//	if (ExCh == false) // Если перезарядки в этой ячейке ещё не было
//	{
//		alpha = polar_angle(y_0, z_0);
//		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
//		uz = Velosity_1(u, cp);
//		nu_ex = (ro * uz * sigma(uz));
//		double kappa = (nu_ex * time) / Kn_;
//		t_ex = -(time / kappa) * log(1.0 - sens->MakeRandom() * (1.0 - exp(-kappa))); // Время до перезарядки
//		t2 = t2 - t_ex;
//		mu_ex = mu * (1.0 - exp(-kappa));
//		mu2 = mu - mu_ex;
//		alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
//
//
//		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
//		y_ex = y_0 + t_ex * Vy;
//		z_ex = z_0 + t_ex * Vz;
//		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
//		{
//			x_ex = x_0 + (1.0 * a * dt) * Vx;  // Координаты перезарядки
//			y_ex = y_0 + (1.0 * a * dt) * Vy;
//			z_ex = z_0 + (1.0 * a * dt) * Vz;
//			t_ex = (1.0 * a * dt);
//			t2 = time - t_ex;
//			alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
//		}
//
//		// Считаем источники для частицы до перезарядки
//		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
//		uz = Velosity_1(u, cp);
//		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
//		uz_E = Velosity_3(u, cp);
//		u1 = vx - Vx;
//		u2 = vy * cos(alpha) - Vy;
//		u3 = vy * sin(alpha) - Vz;
//		double skalar = Vx * u1 + Vy * u2 + Vz * u3;
//
//
//		now->mut.lock();
//		now->par[0].F_n += t_ex * mu;
//		now->par[0].F_u += t_ex * Vx * mu;
//		now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
//		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
//		now->par[0].I_u += mu * t_ex * uz_M * uz * sigma(uz_M) * u1 / u;
//		now->par[0].I_v += mu * t_ex * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
//		now->par[0].I_T += mu * t_ex * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
//		now->mut.unlock();
//
//
//		double alpha = polar_angle(y_ex, z_ex);
//
//		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
//		double r = sqrt(kvv(x_ex, y_ex, z_ex));
//
//		spherical_skorost(y_ex, z_ex, x_ex, vy* cos(alpha), vy* sin(alpha), vx, Ur, Uphi, Uthe);
//		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
//
//		int I = 0;
//		//if (Ur <= 0)
//		//{
//		//	if (r > Ri[I_ - 2])  // Выбираем, сколько расщеплений атома будет.
//		//	{
//		//		I = 4;
//		//	}
//		//	else if (r > Ri[I_ - 3])
//		//	{
//		//		I = 3;
//		//	}
//		//	else if (r > Ri[I_ - 4])
//		//	{
//		//		I = 2;
//		//	}
//		//	else if (r > Ri[I_ - 5])
//		//	{
//		//		I = 1;
//		//	}
//		//	else
//		//	{
//		//		I = 0;
//		//	}
//		//}
//		//cout << "I = " << I << " " << Ur << endl;
//		vector <double> Wr(I + 1);
//		vector <double> Wthe(I + 1);
//		vector <double> Wphi(I + 1);
//		vector <double> mu_(I + 1);
//
//		//cout << 1 << endl;
//		Change_Velosity_Split(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
//		//cout << 2 << endl;
//
//		for (int i = 0; i <= I; i++)
//		{
//			Wr[i] = Wr[i] * cp;
//			Wthe[i] = Wthe[i] * cp;
//			Wphi[i] = Wphi[i] * cp;
//			//cout << mu_[i] << endl;
//		}
//		//exit(-1);
//
//		double aa, bb, cc;
//		bool bol = true;
//		for (int i = 0; i < I; i++)
//		{
//			if (mu_[i] >= 0.1)
//			{
//				bol = true;
//			}
//			else
//			{
//				if (mu_[i] >= sens->MakeRandom() * 0.1)
//				{
//					mu_[i] =0.1;
//					bol = true;
//				}
//				else
//				{
//					bol = false;
//				}
//			}
//			if (bol == true)
//			{
//				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
//				/*for (int k = 0; k < 30000; k++)
//				{
//					cout << x_ex + 0.3 * k * aa << " " << sqrt(kv(y_ex + 0.3 * k * bb) + kv(z_ex + 0.3 * k * cc)) << endl;
//				}*/
//				Fly_exchenge_Split(sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu_[i] * mu_ex, mu_ex, true, i);
//			}
//		}
//		//cout << 3 << endl;
//
//		// У основного атома ограничение на вес больше
//		if (mu_[I] >= Mu[I_ - 1])
//		{
//			bol = true;
//		}
//		else
//		{
//			if (mu_[I] >= sens->MakeRandom() * Mu[I_ - 1])
//			{
//				mu_[I] = Mu[I_ - 1];
//				bol = true;
//			}
//			else
//			{
//				bol = false;
//			}
//		}
//		//cout << 4 << endl;
//		if (bol == true)
//		{
//			dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
//			/*for (int k = 0; k < 30000; k++)
//			{
//				cout << x_ex + 0.3 * k * aa << " " << sqrt(kv(y_ex + 0.3 * k * bb) + kv(z_ex + 0.3 * k * cc)) << endl;
//			}*/
//			Fly_exchenge_Split(sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu_[I] * mu_ex, mu_ex, true, I_ - 1);
//		}
//		//exit(-1);
//
//		/*mmu5 += mu_[I];
//		mn5++;
//
//		if (I >= 1)
//		{
//			mmu1 += mu_[0];
//			mn1++;
//		}
//		if (I >= 2)
//		{
//			mmu2 += mu_[1];
//			mn2++;
//		}
//		if (I >= 3)
//		{
//			mmu3 += mu_[2];
//			mn3++;
//		}
//		if (I >= 4)
//		{
//			mmu4 += mu_[3];
//			mn4++;
//		}*/
//		
//	}
//
//	/*if (mu2 < weight_ * mu_0 && ExCh == false)
//	{
//		if (mu2 >= sens->MakeRandom() * weight_ * mu_0)
//		{
//			mu2 = weight_ * mu_0;
//		}
//		else
//		{
//			return;
//		}
//	}*/
//	bool bol = true;
//	//cout << 5 << endl;
//
//	if (mu2 < Mu[zone] && ExCh == false)
//	{
//		if (mu2 >= sens->MakeRandom() * Mu[zone])
//		{
//			mu2 = Mu[zone];
//		}
//		else
//		{
//			return;
//		}
//	}
//	//cout << 6 << endl;
//	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
//	u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
//	uz = Velosity_1(u, cp);
//	uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
//	uz_E = Velosity_3(u, cp);
//	u1 = vx - Vx;
//	u2 = vy * cos(alpha) - Vy;
//	u3 = vy * sin(alpha) - Vz;
//	double skalar = Vx * u1 + Vy * u2 + Vz * u3;
//
//	now->mut.lock();
//	now->par[0].F_n += t2 * mu2;
//	now->par[0].F_u += t2 * Vx * mu2;
//	now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
//	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
//	now->par[0].I_u += mu2 * t2 * uz_M * uz * sigma(uz_M) * u1 / u;
//	now->par[0].I_v += mu2 * t2 * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
//	now->par[0].I_T += mu2 * t2 * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
//	now->mut.unlock();
//
//	double X = x_ex + t2 * Vx;
//	double Y = y_ex + t2 * Vy;
//	double Z = z_ex + t2 * Vz;
//
//	//cout << 7 << endl;
//	Cell* next = nullptr;             // В какую ячейку попадаем следующую
//	double xk, yk;
//	xk = x_0 + (1.0 * b) * dt * Vx;
//	yk = sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz));
//	if (xk < Left_ || (sqrt(kv(xk) + kv(yk)) > R5_ && xk > 0) || (xk < 0 && yk > R5_))
//	{
//		next = nullptr;
//	}
//	else
//	{
//		for (auto& i : now->Grans)
//		{
//			if (i->Sosed != nullptr)
//			{
//				if (i->Sosed->belong(xk, yk) == true)
//				{
//					next = i->Sosed;
//				}
//			}
//		}
//		if (next == nullptr)
//		{
//			for (auto& k : now->Grans)
//			{
//				if (k->Sosed != nullptr)
//				{
//					for (auto& i : k->Sosed->Grans)
//					{
//						if (i->Sosed != nullptr)
//						{
//							if (i->Sosed->belong(xk, yk) == true)
//							{
//								next = i->Sosed;
//							}
//						}
//					}
//				}
//			}
//		}
//		/*if (next == nullptr)
//		{
//			cout << xk << " " << yk << " " << "no sosed" << endl;
//		}*/
//	}
//	//cout << 8 << endl;
//	if (next != nullptr)
//	{
//		//cout << "Next cell" << endl;
//		Fly_exchenge_Split(sens, X, Y, Z, Vx, Vy, Vz, next, mu2, mu_0, false, zone);
//		//cout << "Vozvrat next" << endl;
//	}
//	return;
//	
//}

void Setka::Change_Velosity(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, double& X, double& Y, double& Z, const double& cp)
{
	double x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi_ * x / (1.0 + 0.5 * sqrtpi_ * x);
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6;
	double om1, om2, om3, lo;
	double y1, y2, y3;
	double v1, v2, v3, u1, u2, u3;
	double uu, yy, vv, D, ko;
	do
	{
		ksi1 = s->MakeRandom();
		ksi2 = s->MakeRandom();
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();
		//cout << "sd " << endl;
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм
			/*do
			{
				om2 = 1.0 - 2.0 * s->MakeRandom();
				om3 = 1.0 - 2.0 * s->MakeRandom();
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
		uu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
	} while (((uu * sigma2(uu, cp)) / (sigma2(x, cp) * (x + yy))) <= ksi6);
	//cout << v2 << endl;
	X = v1;
	Y = v2;
	Z = v3;
}

void Setka::Change_Velosity_Split(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr, vector <double>& Wthe,//
	vector <double>& Wphi, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex, const double& y_ex, const double& z_ex)
	// I -  сколько дополнительных атомов вылетает (кроме основного)
{
	// I - в какой зоне мы сейчас находимся?
	// r - координата точки перезарядки
	double e1 = sqrtpi_ * kv(Ur);
	double e2 = 2.0 * fabs(Ur);
	double e3 = 0.5 * sqrtpi_;
	double p1 = e1 / (e1 + e2 + e3);
	double p2 = e2 / (e1 + e2 + e3);
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7, ksi8;
	double z = 0.0, h = 0.0;
	//cout << "A " << endl;
	for (int i = 0; i < I; i++)
	{
		do
		{
			ksi1 = s->MakeRandom();
			if (p1 > ksi1)
			{
				ksi2 = s->MakeRandom();
				ksi3 = s->MakeRandom();
				z = sqrt(-log(ksi2)) * cos(ksi3 * pi_);
			}
			else if (p1 + p2 > ksi1)
			{
				ksi2 = s->MakeRandom();
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
				ksi2 = s->MakeRandom();
				ksi3 = s->MakeRandom();
				ksi4 = s->MakeRandom();
				ksi5 = s->MakeRandom();
				z = sign(ksi5 - 0.5) * sqrt(-log(ksi2) - log(ksi3) * //
					kv(cos(pi_ * ksi4)));
			}

			Wr[i] = Ur + z;
			h = kv(Ur + z) / kv(fabs(Ur) + fabs(z));
			ksi6 = s->MakeRandom();
		} while (h <= ksi6 || z >= -Ur);
	}
	//cout << "B " << endl;
	//vector <double> w_(I + 1);
	vector <double> gam_(I + 1);
	vector <double> F_(I + 1);
	for (int i = 0; i < I; i++)
	{
		gam_[i] = 1.0 / (kv(r) / kv(Ri[i]) - 1.0);
		//cout << gam_[i] << endl;
		//w_[i] = sqrt(kv(Wr[i]) * gam_[i]);
		F_[i] = 0.5 * gam_[i] * ((0.5 + kv(Ur)) *//
			(1.0 - sign(Ur) * erf(fabs(Ur))) - //
			Ur * exp(-kv(Ur)) / sqrtpi_);
	}
	//cout << "C " << endl;
	double Val;
	double The;

	for (int i = 0; i < I; i++)
	{
		ksi7 = s->MakeRandom();
		ksi8 = s->MakeRandom();
		The = 2.0 * pi_ * ksi7;
		if (i > 0)
		{
			Val = sqrt(kv(w_c_v_s(r, Wr[i], i - 1)) + ksi8 * (kv(w_c_v_s(r, Wr[i], i)) - kv(w_c_v_s(r, Wr[i], i - 1))));
		}
		else
		{
			Val = sqrt(ksi8 * (kv(w_c_v_s(r, Wr[i], i))));
		}
		Wthe[i] = Val * cos(The);
		Wphi[i] = Val * sin(The);
	}
	//cout << "D " << endl;
	double x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double uu = exp(-kv(x)) / sqrtpi_ + (x + 1.0 / (2.0 * x)) * erf(x);

	for (int i = 0; i < I; i++)
	{
		double u = sqrt(kvv(Vr - Wr[i], Vthe - Wthe[i], Vphi - Wphi[i]));
		if (i > 0)
		{
			mu_[i] = (F_[i] - F_[i - 1]) * (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * exp(-kv(Wthe[i] - Uthe) - kv(Wphi[i]));
		}
		else
		{
			mu_[i] = F_[i] * (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * exp(-kv(Wthe[i] - Uthe) - kv(Wphi[i]));
		}
	}

	// Расчёт основного атома
	//x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi_ * x / (1.0 + 0.5 * sqrtpi_ * x);
	double om1, om2, om3, lo;
	double y1, y2, y3;
	double v1, v2, v3, u1, u2, u3;
	double yy, vv, D, ko, uuu;
	//cout << "E " << endl;
	do
	{
		ksi1 = s->MakeRandom();
		ksi2 = s->MakeRandom();
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();
		//cout << "sd " << endl;
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм
			/*do
			{
				om2 = 1.0 - 2.0 * s->MakeRandom();
				om3 = 1.0 - 2.0 * s->MakeRandom();
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
	} while (h <= ksi6 || (v1 <= 0 && kvv(v2, v3, 0.0) < kv(w_c_v_s(r, v1, I - 1))));
	//cout << v2 << endl;
	Wr[I] = v1;
	Wthe[I] = v2;
	Wphi[I] = v3;
	mu_[I] = 1.0;
	for (int i = 0; i < I; i++)
	{
		mu_[I] = mu_[I] - mu_[i];
	}

	/*if (mu_[I] < -5.0)
	{
		double aa, bb, cc;
		cout << "8742  MINUS  I = " << mu_[I] << "  "  << I << "  r = " << r << endl;
		cout << Ur << " " << Uthe << " " << Uphi << " " << Vr << " " << Vthe << " " << Vphi << endl;
		cout << "x = " << x << "   uu = " << uu << endl;
		cout << "erf(x) = " << erf(x) << endl;
		double u = sqrt(kvv(Vr - Wr[3], Vthe - Wthe[3], Vphi - Wphi[3]));
		cout << " u = " << u << endl;
		cout << (F_[3] - F_[2]) * (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * exp(-kv(Wthe[3] - Uthe) - kv(Wphi[3])) << endl;
		cout << (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) << endl;
		exit(-1);
		for (int i = 0; i < I; i++)
		{
			cout << "mu = " << mu_[i] << endl;
		}
		for (int i = 0; i < I; i++)
		{
			cout << "gam = " << gam_[i] << endl;
			cout << "F = " << F_[i] << endl;
		}
		for (int i = 0; i <= I; i++)
		{
			cout << "Wr = " << Wr[i] << "    Wthe = " << Wthe[i] << "    Wphi = " << Wphi[i] << endl;
		}
		for (int i = 0; i < I; i++)
		{
			dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
			for (int ij = 0; ij < 27000; ij++)
			{
				cout << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << i << endl;
			}
		}
		dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
		for (int ij = 0; ij < 30000; ij++)
		{
			cout << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << I << endl;
		}
		exit(-1);
	}*/
}

double Setka::w_c_v_s(const double& r, const double& v, int i)  // W для перезарядки
{
	double gamma = 1.0 / (kv(r) / kv(Ri[i]) - 1.0);
	return sqrt(gamma * kv(v));
}

double Setka::Velosity_1(const double& u, const double& cp)
{
	if (u < 0.00001)
	{
		return 2.0 * cp / sqrtpi_ + 2.0 * u * u / (3.0 * cp * sqrtpi_) - u * u * u * u / (15.0 * cp * cp * cp * sqrtpi_);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp / sqrtpi_ + (u + kv(cp) / (2.0 * u)) * erf(u / cp);
	}
}

double Setka::Velosity_2(const double& u, const double& cp)  // Считает на совсем скорость, а только её числитель (см. статью)
{
	if (u < 0.00001)
	{
		return (8.0 / 3.0) * kv(cp) * kv(cp) * pi * u + (8.0 / 15.0) * kv(cp) * pi * u * u * u - (4.0 / 105.0) * pi * kv(u) * kv(u) * u;
	}
	else
	{
		return  cp * cp * cp * pi * (exp(-u * u / kv(cp)) * cp * u * 2.0 * (kv(cp) + 2.0 * kv(u)) +//
			sqrtpi_ * (4.0 * kv(u) * kv(u) + 4.0 * cp * cp * kv(u) - kv(cp) * kv(cp)) * erf(u / cp)) / (4.0 * u * u);
	}
}

double Setka::Velosity_3(const double& u, const double& cp)
{
	if (u < 0.00001)
	{
		return 8.0 * cp / (3.0 * sqrtpi_) + 8.0 * u * u / (9.0 * cp * sqrtpi_) - 44.0 * u * u * u * u / (135.0 * cp * cp * cp * sqrtpi_);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp * (5.0 * kv(cp) + 2.0 * kv(u)) / (sqrtpi_ * (3.0 * kv(cp) + 2.0 * kv(u))) +//
			(4.0 * kv(u) * kv(u) + 12.0 * cp * cp * kv(u) + 3.0 * kv(cp) * kv(cp)) * erf(u / cp) / (2.0 * u * (3.0 * kv(cp) + 2.0 * kv(u)));
	}
}