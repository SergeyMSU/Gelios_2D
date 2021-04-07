#include "Point.h"
#include <math.h>
using namespace std;

Point::Point(const double& x, const double& y)
{
	this->x = x;
	this->y = y;
	this->number = -1;
	this->type = Point_type::P_No;
	this->Vx = 0.0;
	this->Vy = 0.0;
}