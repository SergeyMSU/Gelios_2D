#pragma once
#include "Help.h"
#include <vector>

using namespace std;

enum Gran_type  // Тип грани нужен для граничных условий
{
    Usualy,                // Обычные
    Inner_sphere,          // Внутренняя сфера
    Extern,                // Выходные условия
    Input,                 // Входной поток
    Axis,                  // Ось симметрии
    Upper_wall,            // Верхняя стенка
    G_no,
};

class Point;
class Cell;
struct Parametr;

class Gran   // Нормаль к грани смотрит направо, если идти от узла A к B. Нормаль должна указывать на соседнюю ячейку.
{
public:
    int number;
    Gran_type type;
    Point* A;
    Point* B;
    bool main_gran;      // true - для основной грани, false для фантомной
    Cell* Master;        // Ячейка, которой принадлежит грань.
    Cell* Sosed;         // Ячейка - сосед по этой грани
    Cell* Sosed_down;
    Cell* Sosed_up;      // Сосед для ТВД по направлению внешней нормали к грани
    Gran* Gran_copy;     // Грань - копия, но с другой нормалью. Введена для удобства

    Gran(Point* A, Point* B, Gran_type type = G_no);
    void Get_Center(double& x, double& y);
    void Get_Center_posle(double& x, double& y);
    void Get_normal(double& n1, double& n2);
    double Get_square(void);

    void Get_par(Parametr& par, int i); // Здесь задаются граничные условия на грани
    void Get_par_TVD(Parametr& par, int i);
    void Get_par_TVD_radial(Parametr& par, int i);


};

