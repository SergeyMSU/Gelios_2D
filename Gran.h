#pragma once
#include "Help.h"
#include <vector>

using namespace std;

enum Gran_type  // ��� ����� ����� ��� ��������� �������
{
    Usualy,                // �������
    Inner_sphere,          // ���������� �����
    Extern,                // �������� �������
    Input,                 // ������� �����
    Axis,                  // ��� ���������
    Upper_wall,            // ������� ������
};

class Point;
class Cell;
struct Parametr;

class Gran   // ������� � ����� ������� �������, ���� ���� �� ���� A � B. ������� ������ ��������� �� �������� ������.
{
public:
    int number;
    Gran_type type;
    Point* A;
    Point* B;
    bool main_gran;      // true - ��� �������� �����, false ��� ���������
    Cell* Master;        // ������, ������� ����������� �����.
    Cell* Sosed;         // ������ - ����� �� ���� �����
    Gran* Gran_copy;     // ����� - �����, �� � ������ ��������. ������� ��� ��������

    Gran(Point* A, Point* B, Gran_type type);
    void Get_Center(double& x, double& y);
    void Get_normal(double& n1, double& n2);
    double Get_square(void);

    void Get_par(Parametr& par, int i); // ����� �������� ��������� ������� �� �����


};

