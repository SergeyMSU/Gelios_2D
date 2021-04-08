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

}

Setka::Setka(int N1, int N2, int N3, int N4, int M1, int M2, int M3, int M4)
// Начальная инициализация сетки,  эту функцию лучше не трогать, она писалась поэтапно и результаты постоянно проверялись
// Изменение чего-либо может привести к ошибке
{
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
		double x = 1.0 * i / (N1 - 1);
		s = (0.35 * x * (1.0 - x) + kv(x)) * (pi_ / 2.0);  //i * (pi_ / 2.0) / (N1 - 1.0);
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
		for (int j = 0; j < M4; j++)
		{
			l = R4_ + (j + 1) * (R5_ - R4_) / M4;
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
			G->main_gran = true;
			K->Grans.push_back(G);
			G->Master = K;
			this->All_Gran.push_back(G);
		}

		for (int i = n * (N2 + N1 - 2) + (N3 - 2) * (N4 + M1 - 1); i < n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Axis);
			G->main_gran = true;
			K->Grans.push_back(G);
			G->Master = K;
			this->All_Gran.push_back(G);
		}


		for (int i = 0; i < N1 + N2 - 2; i++)
		{
			auto K = this->All_Cells[i * (M1 + M2 + M3 + M4 - 1)];
			auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
			auto G2 = new Gran(K->contour[0], K->contour[3], Usualy);
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
			G->main_gran = true;
			K->Grans.push_back(G);
			G->Master = K;
			this->All_Gran.push_back(G);
		}

		if (true)
		{
			auto K = this->All_Cells[this->All_Cells.size() - 1];
			auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
			G->main_gran = true;
			K->Grans.push_back(G);
			G->Master = K;
			this->All_Gran.push_back(G);
		}

	}

	int h = (N1 + N2 - 2) * (M1 + M2 + M3 + M4 - 1);
	for (int i = 0; i < N3 - 1; i++)
	{
		auto K = this->All_Cells[i * (M1 + N4 - 1) + h];
		auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
		auto G2 = new Gran(K->contour[0], K->contour[3], Usualy);
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
	auto G2 = new Gran(K->contour[0], K->contour[3], Usualy);
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
		G->main_gran = true;
		K->Grans.push_back(G);
		G->Master = K;
		this->All_Gran.push_back(G);
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
				return (this->polar_angle(x1, y1) < this->polar_angle(x2, y2));
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
		G2->main_gran = true;
		PO->Grans.push_back(G2);
		G2->Master = PO;
		this->All_Gran.push_back(G2);

		PO = centr_cell[centr_cell.size() - 1];
		G2 = new Gran(PO->contour[PO->contour.size() - 1], PO->contour[0], Axis);
		G2->main_gran = true;
		PO->Grans.push_back(G2);
		G2->Master = PO;
		this->All_Gran.push_back(G2);


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
				return (this->polar_angle(x1, y1) < this->polar_angle(x2, y2));
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
				return (this->polar_angle(x1, y1) < this->polar_angle(x2, y2));
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
				return (this->polar_angle(x1, y1) < this->polar_angle(x2, y2));
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
		"\"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->All_Cells)
	{
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
			" " << i->par[0].ro * ro_o << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o << " " << i->par[0].v * u_o << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " " //
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << //
			" " << i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << endl;

	}
}

void Setka::Proverka(void)
{
	cout << "Proverka  Start" << endl;
	for (auto& i : this->All_Cells)
	{
		if (i->Grans.size() != 4)
		{
			double x, y;
			i->Get_Center(x, y);
			cout << "Gran  82753   in  point:  " << x << " " <<  y << "   Have  sosed: " << i->Grans.size() << endl;
		}
	}

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

	cout << "Proverka  End" << endl;
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
		if (i->Sosed != nullptr)
		{
			fout << i->Master->number << " " << i->Sosed->number << " " << i->Gran_copy->number << endl;
		}
		else
		{
			fout << i->Master->number << " " << -1 << " " << -1 << endl;
		}
	}

	fout << this->All_Gran_copy.size() << endl;
	for (auto& i : this->All_Gran_copy)
	{
		fout << i->Master->number << " " << i->Sosed->number << " " << i->Gran_copy->number << endl;
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

void Setka::Move_Setka_Calculate(void)
{
	double x, y;
	double n1, n2, n;
	double Vs1 = 0.0, Vs2, Vs3 = 0.0, Vs;
	for (auto& i : this->A_Rails)
	{
		auto K1 = i->Key_point[0];
		auto K2 = i->Key_point[1];
		auto K3 = i->Key_point[2];
		auto K4 = i->Key_point[3];
		n = sqrt(kv(K1->x) + kv(K1->y));
		n1 = K1->x / n;
		n2 = K1->y / n;
		Vs1 = K1->Vx * n1 + K1->Vy * n2;
		Vs2 = K2->Vx * n1 + K2->Vy * n2;
		Vs3 = K3->Vx * n1 + K3->Vy * n2;

		for (int j = M1; j < M1 + M2 + 1; j++)
		{
			Vs = Vs1 + (j - M1) * (Vs2 - Vs1) / (M2);
			i->All_point[j]->Vx = n1 * Vs;
			i->All_point[j]->Vy = n2 * Vs;
		}

		for (int j = M1 + M2 + 1; j < M1 + M2 + 1 + M3; j++)
		{
			Vs = Vs2 + (j + 1 - (M1 + M2 + 1)) * (Vs3 - Vs2) / (M3);
			i->All_point[j]->Vx = n1 * Vs;
			i->All_point[j]->Vy = n2 * Vs;
		}

	}

	for (auto& i : this->B_Rails)
	{
		auto K1 = i->Key_point[0];
		auto K2 = i->Key_point[1];
		auto K3 = i->Key_point[2];
		auto K4 = i->Key_point[3];
		n1 = 0.0;
		n2 = 1.0;
		Vs1 = K1->Vx * n1 + K1->Vy * n2;
		Vs2 = K2->Vx * n1 + K2->Vy * n2;

		for (int j = M1; j < M1 + M2 + 1; j++)
		{
			Vs = Vs1 + (j - M1) * (Vs2 - Vs1) / (M2);
			i->All_point[j]->Vx =  n1 * Vs;
			i->All_point[j]->Vy =  n2 * Vs;
		}

		for (int j = M1 + M2 + 1; j < M1 + M2 + 1 + M3; j++)
		{
			Vs = Vs2 + (j + 1 - (M1 + M2 + 1)) * (Vs3 - Vs2) / (M3);
			i->All_point[j]->Vx =  n1 * Vs;
			i->All_point[j]->Vy =  n2 * Vs;
		}

	}

	for (auto& i : this->D_Rails)
	{
		auto K1 = i->Key_point[0];
		auto K2 = i->Key_point[1];
		auto K3 = i->Key_point[2];
		auto K4 = i->Key_point[3];
		n1 = 0.0;
		n2 = 1.0;
		Vs1 = K1->Vx * n1 + K1->Vy * n2;
		Vs2 = K2->Vx * n1 + K2->Vy * n2;

		for (int j = 0; j < M2 + 1; j++)
		{
			Vs = Vs1 + (j) * (Vs2 - Vs1) / (M2);
			i->All_point[j]->Vx =  n1 * Vs;
			i->All_point[j]->Vy =  n2 * Vs;
		}

		for (int j = M2 + 1; j < M2 + 1 + M3; j++)
		{
			Vs = Vs2 + (j + 1 - (M2 + 1)) * (Vs3 - Vs2) / (M3);
			i->All_point[j]->Vx =  n1 * Vs;
			i->All_point[j]->Vy =  n2 * Vs;
		}

	}
}

void Setka::Init_conditions(void)
{
	double x, y, dist;
	double ro = (389.988 * 389.988) / (chi_ * chi_);
	double P_E = ro * chi_ * chi_ / (ggg * 5.0 * 5.0);

	for (auto& i : this->All_Cells)
	{
		i->Get_Center(x, y);
		dist = sqrt(kv(x) + kv(y));
		if (dist < 150)
		{
			i->par[0] = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x / dist, chi_ * y / dist, ro / (dist * dist),//
		(2.0 * (sigma(chi_real) * 389.988) / (2.0 * Kn_ * chi_real))* (1.0 / dist - 0.1 / kv(dist)), 0.000001, chi_real * x / dist, chi_real * y / dist, 0.000001, 0.000001, 0.0, 0.0, 0.000001, 0.000001, 0.0, 0.0,//
			0.000001, 0.000001, 0.0, 0.0};
			i->par[1] = i->par[0];
		}
		else
		{
			i->par[0] = { 1.0, 1.0, Velosity_inf, 0.0, 100.0,//
			0.000001, 0.000001, 0.0, 0.0, 0.000001, 0.000001, 0.0, 0.0, 0.000001, 0.000001, 0.0, 0.0,//
			1.0, 0.5, Velosity_inf, 0.0};
			i->par[1] = i->par[0];
		}
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

			if (sqrt(kv(x) + kv(y)) < 100.0)
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
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc;
				int met = 0;
				/*if (step > 20000)
				{
					met = 1;
				}*/
				time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
							par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vc));
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
				printf("Problemsssss  par1.ro < 0!\n");
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

void Setka::Go_stationary_5_komponent(int step)
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
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc;
				int met = 0;

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro_H1, par1.Q, par1.p_H1, par1.u_H1, par1.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, 0.0, P, PQ, n1, n2, dist, met, Vc));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro_H2, par1.Q, par1.p_H2, par1.u_H2, par1.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, 0.0, P, PQ, n1, n2, dist, met, Vc));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro_H3, par1.Q, par1.p_H3, par1.u_H3, par1.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, 0.0, P, PQ, n1, n2, dist, met, Vc));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro_H4, par1.Q, par1.p_H4, par1.u_H4, par1.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, 0.0, P, PQ, n1, n2, dist, met, Vc));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}
				

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vc));
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
			if (false)//(par1.Q / par1.ro < 90)
			{
				u = par1.u * chi_real;
				v = par1.v * chi_real;
				ro = par1.ro / kv(chi_real);
				p = par1.p;
				Q = par1.Q / kv(chi_real);
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
			q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			/*u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;*/

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!\n");
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
			
			/*if (par1.Q / par1.ro < 90)
			{
				u = par1.u * chi_real;
				v = par1.v * chi_real;
				ro = par1.ro / kv(chi_real);
				p = par1.p;
				Q = par1.Q / kv(chi_real);
			}*/

			int k1 = 0, k2 = 0, k3 = 0, k4 = 0;

			if (kappa < 90)
			{
				k3 = 0;
				k4 = 0;
				if (((kv(u) + kv(v)) / (ggg * p / ro) > 1.0) || ((radius <= 50.0)))
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

void Setka::Go(int step)
{
	cout << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	vector <double> q1(5);
	vector <double> q2(5);
	vector <double> q(5);
	vector <double> n(3);
	double dsl, dsp;

	Solvers Sol;

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

		for (auto& i : this->Line_Contact)
		{
			i->A->Vx = 0.0;
			i->A->Vy = 0.0;
			i->B->Vx = 0.0;
			i->B->Vy = 0.0;
		}

		int bb = -1;
		for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
		{
			auto i = this->Line_Contact[j];
			Parametr par1 = i->Master->par[now1];
			Parametr par2 = i->Sosed->par[now1];
			double n1, n2;
			i->Get_normal(n1, n2);
			double VV;
			double P[4];
			double PQ;

			
			q1[0] = par1.ro;
			q1[1] = par1.p;
			q1[2] = par1.u;
			q1[3] = par1.v;
			q1[4] = 0.0;

			q2[0] = par2.ro;
			q2[1] = par2.p;
			q2[2] = par2.u;
			q2[3] = par2.v;
			q2[4] = 0.0;

			n[0] = n1;
			n[1] = n2;
			n[2] = 0.0;

			//Sol.Godunov_Solver_Alexashov(q1, q2, n, q, dsl, dsp, VV);
			this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
				par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, 1.0, 1, VV);

			double Max = sqrt((kv(par1.u) + kv(par1.v)) / (ggg * par1.p / par1.ro));
			VV = VV * 0.2;

			if (i->A->x < -200 && i->A->y < 300 && VV <= 0.0)
			{
				VV = 0.01;
			}
			//cout << i->A->x << " " << i->A->y << " " << VV << " " << n1 << " " << n2 <<  endl;
			if (Max > 1)
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
			}
		}

		this->Line_Contact[0]->A->Vx = this->Line_Contact[0]->B->Vx + (this->Line_Contact[0]->B->x - this->Line_Contact[0]->A->x);

		this->Line_Contact[this->Line_Contact.size() - 1]->B->Vy = this->Line_Contact[this->Line_Contact.size() - 1]->A->Vy + //
			(this->Line_Contact[this->Line_Contact.size() - 1]->A->y - this->Line_Contact[this->Line_Contact.size() - 1]->B->y);

		if(bb >= 0)
		{
			this->Line_Contact[bb]->A->Vx *= 2.0;
			this->Line_Contact[bb]->A->Vy *= 2.0;
		}

		this->Move_Setka_Calculate();

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
			//cout << "A3" << endl;

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc;
				double VV = (i->A->Vx + i->B->Vx) * 0.5 * n1 + (i->A->Vy + i->B->Vy) * 0.5 * n2;
				if (i->A->type == P_Contact && i->B->type == P_Contact)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, VV, P, PQ, n1, n2, dist, 1, Vc));
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, VV, P, PQ, n1, n2, dist, 0, Vc));
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
				mut.unlock();
			}

			double ro3, p3, u3, v3, Q33;
			double Volume = K->Get_Volume();
			double Volume2 = K->Get_Volume_posle(T[now1]);

			ro3 = par1.ro * Volume/Volume2 - T[now1] * (K->Potok[0] / Volume2 + par1.ro * par1.v / y);
			Q33 = par1.Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * par1.Q * par1.v / y;
			if (ro3 <= 0)
			{
				printf("Problemsssss  par1.ro < 0!\n");
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
			i->x += i->Vx * T[now1];
			i->y += i->Vy * T[now1];
		}

		//this->Line_Contact[0]->A->x = this->Line_Contact[0]->B->x;

		//this->Line_Contact[this->Line_Contact.size() - 1]->B->y = this->Line_Contact[this->Line_Contact.size() - 1]->A->y;
	}
}

double Setka::HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
	const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
	double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vc)
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


	vRR = UZ2 / UZ0;
	vLL = vRR;

	/*vLL = v1;
	vRR = v2;*/

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

double Setka::polar_angle(const double& x, const double& y)
{
	if (x < 0)
	{
		return atan(y / x) + 1.0 * pi_;
	}
	else if (x > 0 && y >= 0)
	{
		return atan(y / x);
	}
	else if (x > 0 && y < 0)
	{
		return atan(y / x) + 2.0 * pi_;
	}
	else if (y > 0 && x >= 0 && x <= 0)
	{
		return pi_ / 2.0;
	}
	else if (y < 0 && x >= 0 && x <= 0)
	{
		return  3.0 * pi_ / 2.0;
	}
	return 0.0;
}