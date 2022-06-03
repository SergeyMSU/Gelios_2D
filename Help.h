#pragma once

#include "Setka.h"
#include "Rail.h"
#include "Point.h"
#include "Cell.h"
#include "Gran.h"
#include "Solvers.h"
#include "sensor.h"
#include "sensor2.h"
#include "MKmethod.h"
#include <mutex>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <algorithm>


double sign(const double& x);
double minmod(const double& x, const double& y);
double polar_angle(const double& x, const double& y);
double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& x3, const double& t3, const double& y);
double linear(const double& x1, const double& t1, const double& x2, const double& t2, const double& y);
void polar_perenos(const double& x1, const double& y1, const double& x2, const double& y2, double& u, double& v);
void dekard_skorost(double x, double y, double z, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz);
void dekard_skorost2(double r2, double the_2, double phi_2, double Vr, double Vphi, double Vtheta, double& Vx, double& Vy, double& Vz);
void spherical_skorost(const double& x, const double& y, const double& z, const double& Vx,//
	const double& Vy, const double& Vz, double& Vr, double& Vphi, double& Vtheta);
double Godunov_squere_rad(const double& x1, const double& r1, const double& x2, const double& r2,//
	const double& x3, const double& r3, const double& x4, const double& r4);
void polar_provorot(const double& phi, double& u, double& v);
//double max(const double& x, const double& y);

#define ga (5.0/3.0)          // Показатель адиабаты
#define ggg (5.0/3.0)
#define g1 (ga - 1.0)
#define kv(x) ( (x)*(x) )
#define pow3(x) ( (x)*(x)*(x) )
#define pow4(x) ( (x)*(x)*(x)*(x) )
#define pow5(x) ( (x)*(x)*(x)*(x)*(x) )
#define pow6(x) ( (x)*(x)*(x)*(x)*(x)*(x) )
#define pow7(x) ( (x)*(x)*(x)*(x)*(x)*(x)*(x) )
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))

#define pi_ 3.14159265358979323846
#define sqrtpi_ 1.77245385

#define Velosity_inf -2.54338 //-4.54127 //-2.54127
#define M_inf 1.97009
#define chi_ 1.0 //36.1059 // 36.1059
#define chi_real 36.1275
#define kurant  0.1
//#define Kn_  0.622171
//#define Kn_  0.4326569808  // 0.622171 // 0.5  // 242.785									// Число Кнудсена
#define Kn_  6.0			
//#define Kn_  0.2	                                            // Число Кнудсена
//#define a_2 0.102578  // 0.10263
#define a_2 0.1307345665  // 0.102578  // 0.10263
#define n_p_LISM_ (3.0) 
#define n_H_LISM_ (1.0)
#define sigma(x) (kv(1.0 - a_2 * log(x)))               // Дифференциальное сечение перезарядки
//#define sigma2(x, y) (kv(1.0 - (a_2/(1.0 - a_2 * log(y))) * log(x)))  // Для другого обезразмеривания скорости на cp
#define sigma2(x, y) (kv(1.0 - a_2 * log((x) * (y))))  // Для другого обезразмеривания скорости на cp


#define n_inner 20  // Сколько ячеек во внутреннем слое
#define n_outer_shock 13  // До какой точки выделять ударную волну (зависит от разрешения сетки)
#define n_inner_shock 5  // До какой точки выделять ударную волну (считая от 90 градусов)
//#define m_inner 16  // Сколько ячеек во внутреннем слое  (до какой ячейки внутреений слой считается)
#define zazor 3.0   // Длина зазора между контактом и близжайшими точками 
#define zazor2 7.0   // Длина зазора между внешней волной и близжайшими точками 
#define s_k 0.025   // Сжатие области 4 к поверхностям
#define alpha_rot 0.1  // Угор вращения (полный сектор), то есть в одну сторону на половину угла и в другую
#define RR_ 390.222 // 389.8557              // Характерный размер задачи
#define R1_ (1.0/RR_)  // 1.0
#define R11_ (38.0/RR_)  // 38.0    // Действует! До какой внутренней границы не двигаем ячейки
#define R111_ (50.0/RR_)  // 50.0 
#define R2_ (100.0/RR_)  // 100.0
#define R3_ (150.0/RR_)  // 150.0
#define R4_ (600.0/RR_)  // 600.0
#define R5_ (2200.0/RR_)  // 2200.0
//#define Rmax_ (1300.0/RR_)  // 1300.0  // (1300.0)   // Максимальный радиус для Монте-Карло
#define Rmax_ (1800.0/RR_)  // 1390.0  // (1300.0)   // Максимальный радиус для Монте-Карло
#define Left_ (-1500.0/RR_)  // -1500
#define H_pow 0.5 //0.45  // Показатель убывания плотности для первой компоненты водорода
#define I_ 8 //5 // Количество зон по радиусу
#define J_ 11 //5 // Количество зон по углу
#define weight_ 0.25     // Веса для атомов
#define gam(i) ( 1.0/(kv(R5_ - 2.0/RR_)/kv(Ri[(i)]) - 1.0) )
#define R_m 4.0   // Для улучшения статистики на оси происходит запуск с диска этого радиуса
#define Ur_m 0.0   // Критерий расщепления
#define geo_step 100.0   // Сколько раз атом попадает в одну ячейку
#define geo_accur 100.0    // Точность геометрии вдоль траектории
#define Weyght 1000000.0  // Для весов чтобы не было потери точности
#define eta_ 0.5    // параметр точности
#define betta_ 1.0  // сжатие

#define k_sphere 1.0  // 0.9  с какого косинуса разыгрываем дополнительные траектории

#define polusum false  // Берём лт полусумму источников или только один
#define mu_statistic true  // считаем ли статистику весов по зонам
#define func_stat false     // считаем ли функцию распределения и моменты на сфере радиуса?
#define R_stat 80.0        // Сфера на которой считаем статистику
#define Al_stat 54         // Расчёт моментов по углу, на сколько дробим угол


// Предикат проверки прерывания веса атома

#define pred(i, j, al) (Mu[i][j] * mu_start * 0.3 * sin(al))

//#define pred(i, j, al) (Mu[i][j] * sin(al))