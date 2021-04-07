#pragma once

#include "Setka.h"
#include "Rail.h"
#include "Point.h"
#include "Cell.h"
#include "Gran.h"
#include "Solvers.h"

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#define ga (5.0/3.0)          // Показатель адиабаты
#define ggg (5.0/3.0)
#define g1 (ga - 1.0)
#define kv(x) ( (x)*(x) )
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))

#define pi_ 3.14159265358979323846
#define sqrtpi_ 1.77245385

#define Velosity_inf -2.54186
#define chi_ 36.1059 // 36.1059
#define chi_real 36.1059
#define kurant  0.95
#define Kn_  242.944										// Число Кнудсена
#define a_2 0.102584
#define n_p_LISM_ (3.0) 
#define n_H_LISM_ (1.0)
#define sigma(x) (kv(1.0 - a_2 * log(x)))               // Диффернциальное сечение перезарядки
#define sigma2(x, y) (kv(1.0 - (a_2/(1.0 - a_2 * log(y))) * log(x)))  // Для другого обезразмеривания скорости на cp

#define R1_ 1.0
#define R2_ 90.0
#define R3_ 180.0
#define R4_ 400.0
#define R5_ 800.0
#define Left_ -800