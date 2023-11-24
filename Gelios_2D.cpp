﻿#include "Help.h"
#include <iomanip>
//#include <unistd.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iterator>
#include <omp.h>
#include "Help.h"


#if USEMPI
#include "mpi.h"
#endif

//extern inline void fff_velocity(const double& a, const double& b);

//template<typename Random_type, typename Distribution_type>
//auto Change(Random_type gen, Distribution_type dis)
//{
//    return (dis(gen));
//}

void fff_velocity_(const double& a, const double& b)
{
    double c = a + b + sin(a) + sin(b) + log(a * b);
}


int main(int argc, char** argv)
{
    std::cout << "Program prepared by Korolkov Sergey. All rights reserved!\n";


    Setka* SS;

    SS = new Setka(20, 15, 10, 30, 50, 20, 10, 10);    // n_inner 3
    //SS->TVD_prepare();
    SS->Proverka();
    SS->Print_cell2();


    // Заполняем ячейки, для которых будем считать магнитное поле
    for (auto& i : SS->All_Cells)
    {
        if (i->type == C_1 || i->type == C_2 || i->type == C_3)
        {
            i->par2[0] = new Parametr2;
            i->par2[1] = new Parametr2;

            i->par2[0]->psi = 0.0;
            i->par2[1]->psi = 0.0;

            double x, y, r, the;
            i->Get_Center(x, y);
            the = polar_angle(x, y);
            r = sqrt(x * x + y * y);
            i->par2[0]->psi = -BE * sin(the) / r;
            i->par2[1]->psi = -BE * sin(the) / r;

            SS->All_Cells_Inner.push_back(i);

            for (auto& j : i->Grans)
            {
                if (j->type == Inner_sphere)
                {
                    i->type = C_1_B;
                }
            }
        }
    }


    int now = 0;
    int now2 = 1;
    double newyzka = 0.0;
    double alpha = pi_ / 90.0;

    // Расчитываем уравнение Лапласа
    for (int iter = 0; iter <= 5000; iter++)
    {
        newyzka = 0.0;

        for (auto& i : SS->All_Cells_Inner)
        {
            double x, y, SS;
            double x2, y2;
            double x3, y3;
            double psi, psi2;
            double sivi = 0.0;
            double divS = 0.0;
            double d1, d2;
            i->Get_Center(x, y);
            psi = i->par2[now]->psi;


            if (i->type == C_1_B) continue;
            for (auto& j : i->Grans)
            {
                j->Get_Center(x2, y2);
                d1 = sqrt(kv(x2 - x) + kv(y2 - y));

                if (j->type == Extern)
                {
                    d2 = d1;
                    psi2 = psi;
                }
                else if (j->type == Axis)
                {
                    d2 = d1;
                    psi2 = 0.0;
                }
                else
                {
                    if (j->Sosed->type == C_4)
                    {
                        d2 = d1;
                        psi2 = psi;
                    }
                    else
                    {
                        j->Sosed->Get_Center(x3, y3);
                        d2 = sqrt(kv(x2 - x3) + kv(y2 - y3));
                        psi2 = (j->Sosed->par2[now]->psi * d1 + psi * d2)/(d1 + d2);
                    }
                }

                //SS = j->Get_square_rotate(alpha);
                SS = j->Get_square();
                sivi = sivi + SS / d1;
                divS = divS + psi2 * SS / d1;
            }
            //sivi = sivi + 2.0 * i->Get_Volume() / (y * alpha);
            //divS = divS + 2.0 * psi * i->Get_Volume() / (y * alpha);
            i->par2[now2]->psi = (divS / sivi);
            newyzka = max(newyzka, fabs(i->par2[now2]->psi - psi));
        }

        now = (now + 1) % 2; // Какие параметры сейчас берём
        now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем

        if (iter % 50 == 0)
        {
            cout << "Newyazka = " << newyzka << endl;
        }
    }

    for (auto& i : SS->All_Cells_Inner)
    {
        double x, y, SS;
        double x2, y2;
        double x3, y3;
        double psi, psi2;
        double sivix = 0.0;
        double siviy = 0.0;
        double n1, n2;
        double d1, d2;
        i->Get_Center(x, y);
        psi = i->par2[now]->psi;

        if (i->type == C_1_B) continue;
        for (auto& j : i->Grans)
        {
            j->Get_Center(x2, y2);
            d1 = sqrt(kv(x2 - x) + kv(y2 - y));

            if (j->type == Extern)
            {
                d2 = d1;
                psi2 = psi;
            }
            else if (j->type == Axis)
            {
                d2 = d1;
                psi2 = -psi;
            }
            else
            {
                if (j->Sosed->type == C_4)
                {
                    d2 = d1;
                    psi2 = psi;
                }
                else
                {
                    j->Sosed->Get_Center(x3, y3);
                    d2 = sqrt(kv(x2 - x3) + kv(y2 - y3));
                    psi2 = j->Sosed->par2[now]->psi;
                }
            }

            j->Get_normal(n1, n2);
            SS = j->Get_square();
            //SS = j->Get_square_rotate(alpha);
            sivix = sivix + n1 * (psi * d2 + psi2 * d1)/(d1 + d2) * SS;
            siviy = siviy + n2 * (psi * d2 + psi2 * d1)/(d1 + d2) * SS;
            //sivix = sivix + n1 * psi2 * SS;
            //siviy = siviy + n2 * psi2 * SS;
        }

        //siviy = siviy + 2.0 * sin(alpha) * psi * i->Get_Volume();

        i->par2[now2]->grad_psi_x = sivix / i->Get_Volume();
        i->par2[now2]->grad_psi_y = siviy / i->Get_Volume();

        //i->par2[now2]->grad_psi_x = sivix / i->Get_Volume_rotate(alpha);
        //i->par2[now2]->grad_psi_y = siviy / i->Get_Volume_rotate(alpha);
    }

    

    // Выводим всё в файл
    ofstream fout;
    fout.open("psi.txt");

    for (auto& i : SS->All_Cells_Inner)
    {
        double x, y;
        i->Get_Center(x, y);
        fout << x << " " << y << " " << i->par2[now2]->psi <<
            " " << i->par2[now2]->grad_psi_x << " " << i->par2[now2]->grad_psi_y << 
            " " << (kv(i->par2[now2]->grad_psi_x) + kv(i->par2[now2]->grad_psi_y))/(8.0 * pi_) << endl;
    }
    fout.close();

    // БЛОК для расчёта газовой динамики (на CPU)
    //if (parameter_4)
    //{
    //    for (auto& i : SS->All_Cells)
    //    {
    //        i->par[0].I_u = 0.0;
    //        i->par[0].I_v = 0.0;
    //        i->par[0].I_T = 0.0;
    //        i->par[0].II_u = 0.0;
    //        i->par[0].II_v = 0.0;
    //        i->par[0].II_T = 0.0;
    //        for (int j = 0; j < pop_atom; j++)
    //        {
    //            i->par[0].H_n[j] = 0.0;
    //            i->par[0].H_u[j] = 0.0;
    //            i->par[0].H_v[j] = 0.0;
    //            i->par[0].H_T[j] = 0.0;
    //            i->par[0].H_uu[j] = 0.0;
    //            i->par[0].H_uv[j] = 0.0;
    //            i->par[0].H_vv[j] = 0.0;
    //            i->par[0].H_uuu[j] = 0.0;
    //        }

    //        i->par[0].k_u = 0.0;
    //        i->par[0].k_v = 0.0;
    //        i->par[0].k_T = 0.0;
    //    }
    //    //SS->Download_Source_MK("source_vers7_7.txt");

    //    SS->Download_Source_MK(parameter_22); // ТУТ ЗАНУЛЯЮТСЯ ИСТОЧНИКИ

    //    //SS->Download_Source_MK("source_vers7_9.txt");
    //    double norm_istok = 1.0;
    //    for (auto& i : SS->All_Cells)
    //    {
    //        i->par[0].I_u /= norm_istok;
    //        i->par[0].I_v /= norm_istok;
    //        i->par[0].I_T /= norm_istok;
    //        i->par[0].II_u /= norm_istok;
    //        i->par[0].II_v /= norm_istok;
    //        i->par[0].II_T /= norm_istok;
    //        for (int j = 0; j < pop_atom; j++)
    //        {
    //            i->par[0].H_n[j] /= norm_istok;
    //            i->par[0].H_u[j] /= norm_istok;
    //            i->par[0].H_v[j] /= norm_istok;
    //            i->par[0].H_T[j] /= norm_istok;
    //            i->par[0].H_uu[j] /= norm_istok;
    //            i->par[0].H_uv[j] /= norm_istok;
    //            i->par[0].H_vv[j] /= norm_istok;
    //            i->par[0].H_uuu[j] /= norm_istok;
    //        }

    //        i->par[0].k_u /= norm_istok;
    //        i->par[0].k_v /= norm_istok;
    //        i->par[0].k_T /= norm_istok;
    //    }

    //    //SS->Print_Gran("gran_vers7_7.txt");
    //    double start;
    //    double end = 0.0;
    //    start = omp_get_wtime();

    //    //SS->normir(0);
    //    //SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);
    //    cout << "SSSSS   " << endl;
    //    cout << SS->B_Rails[SS->B_Rails.size() - 1]->Key_point[0]->x << " " << SS->B_Rails[SS->B_Rails.size() - 1]->Key_point[0]->y << endl;
    //    cout << SS->C_Rails[0]->Key_point[0]->x << " " << SS->C_Rails[0]->Key_point[0]->y << endl;
    //    //SS->Print_cell();
    //    SS->Print_Gran("DO_gran_" + name_gd);
    //    SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);

    //    for (int k = 0; k < 0; k++)  // 10
    //    {
    //        SS->Go_stationary_5_komponent_inner_MK(15000);
    //        SS->Go_5_komponent_MK(30000, false);
    //        SS->Print_Tecplot_MK();
    //    }

    //    int max_k = 100;  // 150  считаются 30 минут, просить 40 минут   900
    //    //SS->Init_conditions();
    //    cout << "FFFFF = " << 0.0 << endl;
    //    for (int k = 1; k <= max_k; k++)  // 10
    //    {
    //        if (true)//(k % 5 == 0)
    //        {
    //            cout << "Global step = " << k << "  is " << max_k << endl;
    //        }
    //        //SS->Go_stationary_5_komponent_inner_2(50000);
    //        //SS->Go_5_komponent_2(50000);
    //        SS->Go_stationary_5_komponent_inner_MK(300);
    //        cout << "Outer" << endl;
    //        //SS->Go_5_komponent__MK2(5000);
    //        SS->Go_5_komponent_MK(1500);
    //        //SS->Init_conditions();
    //        if (k % 1 == 0 && k > 1)
    //        {
    //            SS->Print_Gran("gran_" + name_gd);
    //        }

    //        if ((k % 5 == 0 && k > 1)||(k == 5))
    //        {
    //            SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);
    //        }

    //        if (k % 100 == 0 && k > 1)
    //        {
    //            SS->Save_Setka_ALL_ALPHA(name_gd);
    //        }
    //    }

    //    //SS->normir(1);

    //    double seconds = difftime(end, start);
    //    end = omp_get_wtime();
    //    printf("Work took %f seconds\n", end - start);
    //    SS->Save_Setka_ALL_ALPHA(name_gd);
    //    //SS->Print_Gran("gran_" + name_gd);
    //    //SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);
    //    SS->Print_Tecplot_MK_2d("tecplot_MK_" + name_gd);
    //    return 0;
    //}

    return 0;

}
