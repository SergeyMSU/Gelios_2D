#pragma once
#include "Help.h"

class MKmethod
{
public:
	// Переменные
	double R_[8];       // Радиусы начальных зон (и их количество)
	double alpha_[8];       // Коэффициенты зон при перезарядке(и их количество)
	double gam_[8];       // все гамма для запуска с границы
	double A0_;         // Главный интегралл сразу посчитанный 
	double A1_;         // Главный интегралл сразу посчитанный 
	double A2_;         // Главный интегралл сразу посчитанный 
	int num_area;      // Количество зон (именно количество дополнительных зон)
	double Int_[501];
	double Int_002[51][50];
	double Int_00625[51][50];
	double Int_02[51][50];
	double Int_055[51][50];

	MKmethod(void);


    // Здесь сами алгоритмы розыгрыша

	// розыгрыш начальных параметров запуска
	bool Init_Parametrs(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_);
	// Возвращает false, если не нужно запускать основной атом
	bool Init_Parametrs2(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_);
	int Init(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, double& X_);
	double FF(const double& gam, const double& Yr);

	// Розыгрыш скорости при перезарядке
	bool Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Второй алгоритм Маламы
	bool Change_Velosity2(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Самый первый алгоритм Маламы
	bool Change_Velosity3(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Мой алгоритм (по второму Маламы) с табличным вычислением весов
	bool Change_Velosity4(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Мой алгоритм (по первому Маламы) с табличным вычислением весов
	// Возвращает false, если не нужно запускать основной атом


	void TEST(void);   // Функция - тестирующая разные вспомогательные функции, нужна была для отладки
	double play_mho(Sensor* sens, const double& c);
	double play_mho2(Sensor* sens, const double& c);
	// Разыгрываем \mho 
	double norm_mho(const double& c);
	double h_mho(const double& x, const double& c);


	// Интерполяция интеграллов

	double Int_cp_1(const double& x);
	double Int_cp_2(const double& x);

	// Табличное вычисление интеграллов

	double Get_Int(const double& uu);
	double Get_Int002(const double& Ur, const double& uu);
	double Get_Int00625(const double& Ur, const double& uu);
	double Get_Int02(const double& Ur, const double& uu);
	double Get_Int055(const double& Ur, const double& uu);

	// Табличное вычисление интеграллов (новое)
	double int_1(const double& x, const double& cp);
	double int_2(const double& x, const double& cp);
	double int_3(const double& x, const double& cp);
	double int_1_f1(const double& x);
	double int_1_f2(const double& x);
	double int_1_f3(const double& x);
	double int_2_f1(const double& x);
	double int_2_f2(const double& x);
	double int_2_f3(const double& x);
	double int_3_f1(const double& x);
	double int_3_f2(const double& x);
	double int_3_f3(const double& x);

	double Lin_Interpolate(const double& x1, const double& y1, const double& x2, const double& y2, const double& x);

	// Для розыгрыша на сфере
	double A0(const double& Y);  // Норма функции распределения (название совпадает с файлом)    ПРОВЕРЕНО!
	double G(const double& gam, const double& Y);  // Норма плотности     ПРОВЕРЕНО!
	double F(const double& X, const double& gam, const double& Y);  // Функция распределения для X   // ПРОВЕРЕНО!
	double F0(const double& X, const double& Y);  // Функция распределения для X   // ПРОВЕРЕНО!
	double FI(const double& Z, const double& X, const double& gam, const double& Y);   // ПРОВЕРЕНО!
	double R(const double& X, const double& Y);              // ПРОВЕРЕНО!

	// Для розыгрыша перезарядки внутри области
	double f2(const double& V, const double& gam, const double& ur, const double& ut);
	double f2k(const double& V, const double& gam, const double& ur, const double& ut);  // учтено два члена ряда

	// Вспомогательные функции
	double Hx(const double& gam1, const double& gam2, const double& X, const double& Y, const double& ksi);
	double Hwr(const double& gam1, const double& gam2, const double& Z, const double& X, const double& Y, const double& ksi);

	double Hvr(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi);
	double Hvrk(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi);


private:
	

};

