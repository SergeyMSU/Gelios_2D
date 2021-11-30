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
    G_no,
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
    Cell* Master;        // ������, ������� ����������� �����.
    Cell* Sosed;         // ������ - ����� �� ���� �����
    Cell* Sosed_down;    // ��� ��� TVD
    Cell* Sosed_up;      // ����� ��� ��� �� ����������� ������� ������� � �����
    Gran* Gran_copy;     // ����� - �����, �� � ������ ��������. ������� ��� ��������
    double a;
    double b;
    double aa;           // aa * x + bb * y + cc = 0
    double bb;
    double cc;
    bool parallel;       // true ���� ����� ����������� ��� �
    bool main_gran;      // true - ��� �������� �����, false ��� ���������
    int koef;            // ����������� ��� ��������� ��������� ������ ��� ����������� � ����� ������� �����  (y - ax - b = 0)
                         // ������ ���� ������ ����

    Gran(Point* A, Point* B, Gran_type type = G_no);
    void Get_Center(double& x, double& y);
    void Get_Center_posle(double& x, double& y);
    void Get_normal(double& n1, double& n2);
    double Get_square(void);
    double Get_lenght(void);
    double Get_square_rotate(const double& angle);

    void Get_par(Parametr& par, int i); // ����� �������� ��������� ������� �� �����
    void Get_par_TVD(Parametr& par, int i);
    void Get_par_TVD_radial(Parametr& par, int i);

    bool belong(const double& x, const double& y);         // ����������� �� ����� ���������������� �����
    bool belong_gran(const double& x, const double& y);         // ����������� �� ����� �����

    void renew(void); // �������� �������� a � b
 

};

