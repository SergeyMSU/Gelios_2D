#pragma once

#include "Help.h"
#include <vector>

using namespace std;

enum Rail_type
{
    // ���� ��������� ������������� - ��� ��������� �������� ����� ���� ������
    // ������ ������������� ���������� ������� (�� ������ � �������)
    A,
    B,
    C,
    D,
    No,
};

class Point;

class Rail
{
public:
    Rail_type type;
    vector <Point*> All_point;
    vector <Point*> Key_point;
    int M1;
    int M2;
    int M3;
    int M4;

    double s; // ���� ���� ��� ���-�� ������ ��� ������� rail_type

	Rail(const double& s);

    void Init_start(Rail* T);
};
