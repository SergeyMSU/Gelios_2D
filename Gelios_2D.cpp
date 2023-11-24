#include "Help.h"
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

    SS = new Setka(20, 10, 10, 10, 20, 20, 10, 10);    // n_inner 3
    //SS->TVD_prepare();
    SS->Proverka();
    SS->Print_cell2();

    // Нужно отдельно тут заполнять массивы ячеек (внутренние и т.д.)

    // Подготовим функции распределения


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
