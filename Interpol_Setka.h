#pragma once

#include "Help.h"
#include <vector>
#include <string>
#include <mutex>
#include "sensor2.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

class Setka;
class Point;
class Gran;
class Cell;

using namespace std;

class Interpol_Setka
{
public:
	vector <Point*> All_Points;   // ��� ����� (������ ����� �������� ����� + ��������������)
	vector <Gran*> All_Gran;           // ����� ��������
	vector <Gran*> All_Gran_copy;      // ����� ��������� (������� � ������ �������)
	vector <Cell*> All_Cells;          // ��� ������

	Interpol_Setka(Setka* SS);     // ���������� ���������������� ����� �� �������� SS

	void Print_Cell(string name);       // ������ �����
};

