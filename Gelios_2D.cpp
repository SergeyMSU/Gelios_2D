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


    if (false)
    {
        Sensor* s = new Sensor(1, 4, 8);
        Sensor& adv = *s;

        //std::random_device rd;
        //std::mt19937 gen(rd());  // mt19937    mt19937
        //std::uniform_real_distribution<double> dis(0.0000001, 1.0);

        /*std::mt19937 gen2(rd());
        std::uniform_real_distribution<double> dis2(0.0, 1.0);*/
        //float aa;
        //float bb = 1.23;
        //float cc = 2324.2323;

        double start;
        double end;
        start = omp_get_wtime();
        

//#pragma omp parallel for
        for (unsigned long int n = 0; n < 2000; ++n) {
            //std::mt19937 gen2(rd());
            //std::uniform_real_distribution<double> dis2(0.0, 1.0);

            for (unsigned long int n2 = 0; n2 < 1000; ++n2) {
                for (unsigned long int n3 = 0; n3 < 1000; ++n3) {
                    //aa = exp(sin(sqrt(bb * bb + cc * cc + n))) + log(bb * cc);
                    //aa = Change(gen2, dis2);
                    //dis(gen);
                    //double a = 1.23423432;
                    //double b = 2.23424123;
                    //Change(a, b);
                    adv.MakeRandom();
                    //aa = s->MakeRandom();
                }
            }
        }

        //cout << s->a1_ << " " << s->a2_ << " " << s->a3_ << endl;
        //cout << adv.a1_ << " " << adv.a2_ << " " << adv.a3_ << endl;

        //cout << SIZE_MAX << endl;
        end = omp_get_wtime();
        printf("Work took %f seconds\n", end - start);

        exit(-1);
    }

   

    //CC->TVD_prepare();
    //CC->Proverka();

    //CC->Init_conditions();

    //CC->M_K_prepare();     // Нужно комментить, если не считается монте-карло, там удаляются источники
    //CC->MK_start_new();
  
    ////CC->Print_Gran("surface6_test.txt");
    //CC->Print_Tecplot_MK();
    ////CC->Print_cell2();
    //CC->Save_Setka_ALL_ALPHA("vers_test1.txt");

    //exit(-1);

    // КОНЕЦ ТЕСТОВОЙ ЗАДАЧИ ---------------------------------------------------------------


    Setka* SS, * K, *SS2, * SS3;
    if (false)
    {
        SS = new Setka(80, 20, 15, 20, 50, 30, 40, 25);    // n_inner 30
        K = new Setka();
        K->Download_Setka_ALL_ALPHA_2_0("vers7_19_gaz_dinamic.txt");
        SS->Copy(K);
    }

    SS = new Setka();

    //SS->Download_Setka_ALL_ALPHA_2_0("vers6_106.txt");  // 17    IPROBE
    //cout << "SLEEP  sec" << endl;
    //sleep(7800);
    //cout << "NOT SLEEP" << endl;
    //SS->Download_Setka_ALL_ALPHA_2_0("vers6_100.txt");  // 17       IEX

    SS->Download_Setka_ALL_ALPHA_2_0(parameter_1);  // 16_7  5  10      IEX    vers18_21.txt  vers7_9
    string name_gd = parameter_21;
    SS->TVD_prepare();
    SS->Proverka();

    // Подготовим функции распределения


    int ijk = 0;
    string nmk;
    // Считаем функцию распределения для Насти
    if (false)
    {
        for (auto& KL : SS->All_Cells)
        {
            double xx, yy;
            KL->Get_Center(xx, yy);
            //cout << "35 point   " << xx << " " << yy << endl;

            double rrr = sqrt(xx * xx + yy * yy);
            if (rrr < 78.0 / RR_ || rrr > 82.0 / RR_)
            {
                continue;
            }

            ijk++;
            nmk = "S4_" + to_string(ijk);
            auto asd = new Dist_func(30, 30, 30, -3.0, 4.0, 0.0, 3.5, 0.0, 5.0);
            asd->xxx = xx;
            asd->yyy = yy;
            asd->call_name(nmk);
            SS->Dist_func_all.push_back(asd);
            KL->df_s4 = asd;
            KL->df_s4_bool = true;

            asd = new Dist_func(30, 30, 30, -2.0, 3.5, 0.0, 4.1, 0.0, 4.1);
            asd->xxx = xx;
            asd->yyy = yy;
            nmk = "S3_" + to_string(ijk);
            asd->call_name(nmk);
            SS->Dist_func_all.push_back(asd);
            KL->df_s3 = asd;
            KL->df_s3_bool = true;
        }
    }
    
    
    // PUI
    //SS->Print_for_Igor();
    //exit(-1);

    //SS->Init_conditions();
    //SS->Print_Gran("gran_vers18_4.txt");
    /*SS->Print_Tecplot_MK();
    exit(-1);*/

    //SS->Print_cell2();
    //SS->Print_Gran();
    //SS->Print_Gran("sur7_101.txt");

    /*double xx, yy;
    SS->All_Cells[34]->Get_Center(xx, yy);
    cout << xx << " " << yy << endl;
    cout << SS->All_Cells[34]->number << endl;
    exit(-1);*/


    // Двигаем узлы вручную
    if (true)
    {

        for (auto& i : SS->All_Points)
        {
            i->Vx = 0.0;
            i->Vy = 0.0;
        }

        //SS->Line_Inner[47]->B->Vx = -0.022;
        //SS->Line_Inner[48]->B->Vx = -0.005;

        SS->Move_Setka_Calculate_2(0.0);

        for (auto& i : SS->All_Points)
        {
            i->x = i->x2;
            i->y = i->y2;
            i->Vx = 0.0;
            i->Vy = 0.0;
            i->count = 0;
        }
    }


    // Подготовка массивов для внутренней области счёта и перенормировка параметров, если надо. НЕ УДАЛЯТЬ
    for (auto i : SS->All_Cells)
    {
        if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
        {
            //Для счёта монте-карло
           //i->par[0].u = i->par[0].u * (chi_real / chi_);       // Перенормировка
           //i->par[0].v = i->par[0].v * (chi_real / chi_);
           //i->par[0].ro = i->par[0].ro / kv(chi_real / chi_);
           //i->par[0].Q = i->par[0].Q / kv(chi_real / chi_);

           //i->par[1].u = i->par[1].u * (chi_real / chi_);       // Перенормировка
           //i->par[1].v = i->par[1].v * (chi_real / chi_);
           //i->par[1].ro = i->par[1].ro / kv(chi_real / chi_);
           //i->par[1].Q = i->par[1].Q / kv(chi_real / chi_);


           // для счёта плазмы
           //i->par[1].u = i->par[0].u = i->par[0].u / (chi_real / chi_);       // Перенормировка
           //i->par[1].v = i->par[0].v = i->par[0].v / (chi_real / chi_);
           //i->par[1].ro = i->par[0].ro = i->par[0].ro * kv(chi_real / chi_);
           //i->par[1].Q = i->par[0].Q = i->par[0].Q * kv(chi_real / chi_);
        }


        if (false)//(i->type == C_centr)
        {
            i->par[0].ro_H1 = i->par[1].ro_H1 = 0.00000001;
            i->par[0].p_H1 = i->par[1].p_H1 = 0.00000001;
            i->par[0].u_H1 = i->par[1].u_H1 = 0.0;
            i->par[0].v_H1 = i->par[1].v_H1 = 0.0;

            i->par[0].ro_H2 = i->par[1].ro_H2 = 0.00000001;
            i->par[0].p_H2 = i->par[1].p_H2 = 0.00000001;
            i->par[0].u_H2 = i->par[1].u_H2 = 0.0;
            i->par[0].v_H2 = i->par[1].v_H2 = 0.0;

            i->par[0].ro_H3 = i->par[1].ro_H3 = 0.00000001;
            i->par[0].p_H3 = i->par[1].p_H3 = 0.00000001;
            i->par[0].u_H3 = i->par[1].u_H3 = 0.0;
            i->par[0].v_H3 = i->par[1].v_H3 = 0.0;

            i->par[0].ro_H4 = i->par[1].ro_H4 = 0.00000001;
            i->par[0].p_H4 = i->par[1].p_H4 = 0.00000001;
            i->par[0].u_H4 = i->par[1].u_H4 = 0.0;
            i->par[0].v_H4 = i->par[1].v_H4 = 0.0;
        }


        double x, y;
        i->Get_Center(x, y);
        if (sqrt(x * x + y * y) <= R111_)
        {
            SS->All_Cells_Inner.push_back(i);
        }


        if (i->type == C_centr)
        {
            SS->All_Cells_zero.push_back(i);
        }
    }


    // Блок расчёта с Kn -> infty
    if (false)
    {
        for (auto i : SS->All_Cells)
        {
            if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
            {
                //Для счёта монте-карло
               i->par[0].u = i->par[0].u / (chi_real / chi_);       // Перенормировка
               i->par[0].v = i->par[0].v / (chi_real / chi_);
               i->par[0].ro = i->par[0].ro * kv(chi_real / chi_);
               i->par[0].Q = i->par[0].Q * kv(chi_real / chi_);

               i->par[1].u = i->par[1].u / (chi_real / chi_);       // Перенормировка
               i->par[1].v = i->par[1].v / (chi_real / chi_);
               i->par[1].ro = i->par[1].ro * kv(chi_real / chi_);
               i->par[1].Q = i->par[1].Q * kv(chi_real / chi_);
            }
        }
        //SS->Print_Gran("gran_vers7_7.txt");
        double start;
        double end = 0.0;
        start = omp_get_wtime();

        int max_k = 1000;  // 150  считаются 30 минут, просить 40 минут

        for (int k = 0; k < max_k; k++)  // 10
        {
            cout << "Global step = " << k + 1 << "  is " << max_k << endl;
            SS->Go_stationary_inner_infty(1500);
            SS->Go_5_komponent_infty(3000, true);
            if (k % 200 == 0 && k > 1)
            {
                SS->Print_Gran("gran_vers18_2.txt");
                SS->Save_Setka_ALL_ALPHA("vers18_2.txt");
                SS->Print_Tecplot_MK();
            }
        }

        for (auto i : SS->All_Cells)
        {
            if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
            {
                //Для счёта монте-карло
                i->par[0].u = i->par[0].u * (chi_real / chi_);       // Перенормировка
                i->par[0].v = i->par[0].v * (chi_real / chi_);
                i->par[0].ro = i->par[0].ro / kv(chi_real / chi_);
                i->par[0].Q = i->par[0].Q / kv(chi_real / chi_);

                i->par[1].u = i->par[1].u * (chi_real / chi_);       // Перенормировка
                i->par[1].v = i->par[1].v * (chi_real / chi_);
                i->par[1].ro = i->par[1].ro / kv(chi_real / chi_);
                i->par[1].Q = i->par[1].Q / kv(chi_real / chi_);
            }
        }

        double seconds = difftime(end, start);
        end = omp_get_wtime();
        printf("Work took %f seconds\n", end - start);
        SS->Print_Gran("gran_vers18_1.txt");
        SS->Save_Setka_ALL_ALPHA("vers18_1.txt");
        SS->Print_Tecplot_MK();
        return 0;
    }

    // БЛОК для расчёта газовой динамики (на CPU)
    if (parameter_4)
    {
        for (auto& i : SS->All_Cells)
        {
            i->par[0].I_u = 0.0;
            i->par[0].I_v = 0.0;
            i->par[0].I_T = 0.0;
            i->par[0].II_u = 0.0;
            i->par[0].II_v = 0.0;
            i->par[0].II_T = 0.0;
            for (int j = 0; j < pop_atom; j++)
            {
                i->par[0].H_n[j] = 0.0;
                i->par[0].H_u[j] = 0.0;
                i->par[0].H_v[j] = 0.0;
                i->par[0].H_T[j] = 0.0;
                i->par[0].H_uu[j] = 0.0;
                i->par[0].H_uv[j] = 0.0;
                i->par[0].H_vv[j] = 0.0;
                i->par[0].H_uuu[j] = 0.0;
            }

            i->par[0].k_u = 0.0;
            i->par[0].k_v = 0.0;
            i->par[0].k_T = 0.0;
        }
        //SS->Download_Source_MK("source_vers7_7.txt");
        
        SS->Download_Source_MK(parameter_22); // ТУТ ЗАНУЛЯЮТСЯ ИСТОЧНИКИ
        
        //SS->Download_Source_MK("source_vers7_9.txt");
        double norm_istok = 1.0;
        for (auto& i : SS->All_Cells)
        {
            i->par[0].I_u /= norm_istok;
            i->par[0].I_v /= norm_istok;
            i->par[0].I_T /= norm_istok;
            i->par[0].II_u /= norm_istok;
            i->par[0].II_v /= norm_istok;
            i->par[0].II_T /= norm_istok;
            for (int j = 0; j < pop_atom; j++)
            {
                i->par[0].H_n[j] /= norm_istok;
                i->par[0].H_u[j] /= norm_istok;
                i->par[0].H_v[j] /= norm_istok;
                i->par[0].H_T[j] /= norm_istok;
                i->par[0].H_uu[j] /= norm_istok;
                i->par[0].H_uv[j] /= norm_istok;
                i->par[0].H_vv[j] /= norm_istok;
                i->par[0].H_uuu[j] /= norm_istok;
            }

            i->par[0].k_u /= norm_istok;
            i->par[0].k_v /= norm_istok;
            i->par[0].k_T /= norm_istok;
        }

        //SS->Print_Gran("gran_vers7_7.txt");
        double start;
        double end = 0.0;
        start = omp_get_wtime();

        //SS->normir(0);
        //SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);
        cout << "SSSSS   " << endl;
        cout << SS->B_Rails[SS->B_Rails.size() - 1]->Key_point[0]->x << " " << SS->B_Rails[SS->B_Rails.size() - 1]->Key_point[0]->y << endl;
        cout << SS->C_Rails[0]->Key_point[0]->x << " " << SS->C_Rails[0]->Key_point[0]->y << endl;
        //SS->Print_cell();
        SS->Print_Gran("DO_gran_" + name_gd);
        SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);
        return 0;

        for (int k = 0; k < 0; k++)  // 10
        {
            SS->Go_stationary_5_komponent_inner_MK(15000);
            SS->Go_5_komponent_MK(30000, false);
            SS->Print_Tecplot_MK();
        }

        int max_k = 900;  // 150  считаются 30 минут, просить 40 минут   900
        //SS->Init_conditions();
        for (int k = 0; k < max_k; k++)  // 10
        {
            if (k % 10 == 0 && k > 1)
            {
                cout << "Global step = " << k + 1 << "  is " << max_k << endl;
            }
            //SS->Go_stationary_5_komponent_inner_2(50000);
            //SS->Go_5_komponent_2(50000);
            SS->Go_stationary_5_komponent_inner_MK(1500);
            //SS->Go_5_komponent__MK2(5000);
            SS->Go_5_komponent_MK(3000);
            //SS->Init_conditions();
            if (k % 10 == 0 && k > 1)
            {
                SS->Print_Gran("gran_" + name_gd);
            }

            if (k % 170 == 0 && k > 1)
            {
                SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);
            }

            if (k % 230 == 0 && k > 1)
            {
                SS->Save_Setka_ALL_ALPHA(name_gd);
            }
        }

        //SS->normir(1);

        double seconds = difftime(end, start);
        end = omp_get_wtime();
        printf("Work took %f seconds\n", end - start);
        SS->Save_Setka_ALL_ALPHA(name_gd);
        SS->Print_Gran("gran_" + name_gd);
        SS->Print_Tecplot_MK("tecplot_MK_" + name_gd);
        return 0;
    }

    //SS->culc_PUI();

    //exit(-1);
    SS->M_K_prepare();     // Нужно комментить, если не считается монте-карло, там удаляются источники
    double start;
    double end = 0.0;
    start = omp_get_wtime();

    SS->MK_start_new();  // Была эта

    
    //SS->Download_Source_MK(parameter_22);
    SS->func_pogloshenie();
    //SS->Print_Tecplot_MK();

    //SS->MK_start_2_0();

    //SS->MPI_MK_start(argc, argv);


    int rank = 0, size = 0;

#if USEMPI 
    MPI_Comm_size(MPI_COMM_WORLD, &size);               // Получить общее число процессов - компов
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);               // Получить номер текущего процесса - компа

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "Barrer " << rank << endl;
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (rank == 0)
    {
        double seconds = difftime(end, start);
        end = omp_get_wtime();
        printf("Work took %f seconds\n", end - start);
        //SS->Print_for_Igor();

       //SS->GD_prepare();

        SS->Print_cell2();
        //SS2 = new Setka();
        //SS3 = new Setka();
        //SS2->Download_Setka_ALL_ALPHA_2_0("vers6_118.txt");  // 84 до добавления источников 
        //SS3->Download_Setka_ALL_ALPHA_2_0("vers6_121.txt");  // 84 до добавления источников 
        if (false)
        {
            for (int i = 0; i < SS->All_Cells.size(); i++)
            {
                auto A1 = SS->All_Cells[i];
                auto A2 = SS2->All_Cells[i];
                auto A3 = SS3->All_Cells[i];

                A1->par[0].H_n2[0] = A2->par[0].H_n[0];
                A1->par[0].H_n2[1] = A2->par[0].H_n[1];
                A1->par[0].H_n2[2] = A2->par[0].H_n[2];
                A1->par[0].H_n2[3] = A2->par[0].H_n[3];

                A1->par[0].H_u2[0] = A2->par[0].H_u[0];
                A1->par[0].H_u2[1] = A2->par[0].H_u[1];
                A1->par[0].H_u2[2] = A2->par[0].H_u[2];
                A1->par[0].H_u2[3] = A2->par[0].H_u[3];

                A1->par[0].H_v2[0] = A2->par[0].H_v[0];
                A1->par[0].H_v2[1] = A2->par[0].H_v[1];
                A1->par[0].H_v2[2] = A2->par[0].H_v[2];
                A1->par[0].H_v2[3] = A2->par[0].H_v[3];

                A1->par[0].H_T2[0] = A2->par[0].H_T[0];
                A1->par[0].H_T2[1] = A2->par[0].H_T[1];
                A1->par[0].H_T2[2] = A2->par[0].H_T[2];
                A1->par[0].H_T2[3] = A2->par[0].H_T[3];

                A1->par[0].k_u2 = A2->par[0].k_u;
                A1->par[0].k_v2 = A2->par[0].k_v;
                A1->par[0].k_T2 = A2->par[0].k_T;

                A1->par[0].H_n3[0] = A3->par[0].H_n[0];
                A1->par[0].H_n3[1] = A3->par[0].H_n[1];
                A1->par[0].H_n3[2] = A3->par[0].H_n[2];
                A1->par[0].H_n3[3] = A3->par[0].H_n[3];

                A1->par[0].H_u3[0] = A3->par[0].H_u[0];
                A1->par[0].H_u3[1] = A3->par[0].H_u[1];
                A1->par[0].H_u3[2] = A3->par[0].H_u[2];
                A1->par[0].H_u3[3] = A3->par[0].H_u[3];

                A1->par[0].H_v3[0] = A3->par[0].H_v[0];
                A1->par[0].H_v3[1] = A3->par[0].H_v[1];
                A1->par[0].H_v3[2] = A3->par[0].H_v[2];
                A1->par[0].H_v3[3] = A3->par[0].H_v[3];

                A1->par[0].H_T3[0] = A3->par[0].H_T[0];
                A1->par[0].H_T3[1] = A3->par[0].H_T[1];
                A1->par[0].H_T3[2] = A3->par[0].H_T[2];
                A1->par[0].H_T3[3] = A3->par[0].H_T[3];

                A1->par[0].k_u3 = A3->par[0].k_u;
                A1->par[0].k_v3 = A3->par[0].k_v;
                A1->par[0].k_T3 = A3->par[0].k_T;
            }
        }
        //delete SS2;
        //delete SS3;

        // SS->Download_Source_MK("source_vers19_1.txt");
        //SS->Print_for_Igor();
        //SS->culc_K_Istok();
        //SS->Save_Source_MK("source_vers18_16.txt");
        //SS->Print_Tecplot_MK();
        //SS->Print_Gran("gran_vers19_1_0.txt");
        //SS->Print_Tecplot_MK("21_2");

        for (int k = 0; k < 0; k++)  // 10
        {
            start = omp_get_wtime();
            cout << "Global step = " << k + 1 << endl;
            //SS->Go_stationary_5_komponent_inner_2(50000);
            //SS->Go_5_komponent_2(50000);
            SS->Go_stationary_5_komponent_inner_MK(15000);
            //SS->Go_5_komponent__MK2(5000);
            SS->Go_5_komponent_MK(30000, false);
            end = omp_get_wtime();
            printf("Work took %f seconds\n", end - start);
        }

        for (int k = 0; k < 0; k++)  // 10
        {
            start = omp_get_wtime();
            cout << "Global step = " << k + 1 << endl;
            //SS->Go_stationary_5_komponent_inner_2(50000);
            //SS->Go_5_komponent_2(50000);
            SS->Go_stationary_5_komponent_inner_MK(10000);
            //SS->Go_5_komponent__MK2(5000);
            SS->Go_5_komponent_MK(30000);
            end = omp_get_wtime();
            printf("Work took %f seconds\n", end - start);
            if (k % 30 == 0 && k > 1)
            {
                string nm;
                nm = "vers19_1_" + to_string(k) + ".txt";
                SS->Print_Gran("gran_" + nm);
                SS->Save_Setka_ALL_ALPHA(nm);
            }
        }

        //SS->Save_Setka_ALL_ALPHA("vers7_6.txt");
        SS->Save_Source_MK(parameter_22);
        SS->Print_Tecplot_MK();
        //SS->Print_cell2();
       // SS->Print_Gran("gran_vers19_1.txt");
        //SS->Print_Tecplot_MK();
        //SS->Print_Sourse();
        //SS->Save_Setka_ALL_ALPHA("vers6_107.txt");
        //SS->Save_Setka_ALL_ALPHA("vers19_1.txt");

        for (auto& i : SS->Dist_func_all)
        {
            //i->print_1d(1);
            //i->print_1d(2);
            //i->print_1d(3);
            //i->print_3d();
            break;
        }

        ofstream fout_cr;
        fout_cr.open("Pogloshenie.txt");

        for (int k = 3; k < 4; k++)
        {
            for (int i = 0; i < 1; i++)
            {
               // double alf = (i + 1) * (pi_ / 2.0) / (pogl_alf_ + 1.0);
                for (int j = 0; j < pogl_rad_; j++)
                {
                  //  fout_cr << (SS->pogVmin + (j + 0.5) * (SS->pogVmax - SS->pogVmin) / pogl_rad_) << " "
                  //      << SS->pogloshenie[k][i][j] * 3.0/(2.0 * pi_ * sin(alf) * ((SS->pogVmax - SS->pogVmin) / pogl_rad_) ) << endl;
                }
            }
        }

        fout_cr.close();

    }

#if USEMPI 
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "EXIT MPI  " << rank << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return 0;
    exit(-1);
    // ----------------------------------------------------------------------------------------------------
    delete SS;
    SS = new Setka();
    SS->Download_Setka_ALL_ALPHA_2_0("vers11_2.txt");
    SS->TVD_prepare();
    SS->Proverka();
    for (auto i : SS->All_Cells)
    {
        if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
        {
            //Для счёта монте-карло
           //i->par[0].u = i->par[0].u * (chi_real / chi_);       // Перенормировка
           //i->par[0].v = i->par[0].v * (chi_real / chi_);
           //i->par[0].ro = i->par[0].ro / kv(chi_real / chi_);
           //i->par[0].Q = i->par[0].Q / kv(chi_real / chi_);

           //i->par[1].u = i->par[1].u * (chi_real / chi_);       // Перенормировка
           //i->par[1].v = i->par[1].v * (chi_real / chi_);
           //i->par[1].ro = i->par[1].ro / kv(chi_real / chi_);
           //i->par[1].Q = i->par[1].Q / kv(chi_real / chi_);


           // для счёта плазмы
           //i->par[1].u = i->par[0].u = i->par[0].u / (chi_real / chi_);       // Перенормировка
           //i->par[1].v = i->par[0].v / (chi_real / chi_);
           //i->par[1].ro = i->par[0].ro = i->par[0].ro * kv(chi_real / chi_);
           //i->par[1].Q = i->par[0].Q = i->par[0].Q * kv(chi_real / chi_);
        }


        if (i->type == C_centr)
        {
            i->par[0].ro_H1 = i->par[1].ro_H1 = 0.00000001;
            i->par[0].p_H1 = i->par[1].p_H1 = 0.00000001;
            i->par[0].u_H1 = i->par[1].u_H1 = 0.0;
            i->par[0].v_H1 = i->par[1].v_H1 = 0.0;

            i->par[0].ro_H2 = i->par[1].ro_H2 = 0.00000001;
            i->par[0].p_H2 = i->par[1].p_H2 = 0.00000001;
            i->par[0].u_H2 = i->par[1].u_H2 = 0.0;
            i->par[0].v_H2 = i->par[1].v_H2 = 0.0;

            i->par[0].ro_H3 = i->par[1].ro_H3 = 0.00000001;
            i->par[0].p_H3 = i->par[1].p_H3 = 0.00000001;
            i->par[0].u_H3 = i->par[1].u_H3 = 0.0;
            i->par[0].v_H3 = i->par[1].v_H3 = 0.0;

            i->par[0].ro_H4 = i->par[1].ro_H4 = 0.00000001;
            i->par[0].p_H4 = i->par[1].p_H4 = 0.00000001;
            i->par[0].u_H4 = i->par[1].u_H4 = 0.0;
            i->par[0].v_H4 = i->par[1].v_H4 = 0.0;
        }


        double x, y;
        i->Get_Center(x, y);
        if (sqrt(x * x + y * y) <= R111_)
        {
            SS->All_Cells_Inner.push_back(i);
        }

        //if (x < -700 && y < 200)
        //{
        //    i->par[1].u = i->par[0].u = Velosity_inf;       // Перенормировка
        //    i->par[1].v = i->par[0].v = 0.0;
        //}

        //if (sqrt(x * x + y * y) <= 80.0)
        //{
        //    i->par[1].u = i->par[0].u = 36.12 / (chi_real / chi_) * x / sqrt(x * x + y * y);       // Перенормировка
        //    i->par[1].v = i->par[0].v = 36.12 / (chi_real / chi_) * y / sqrt(x * x + y * y);
        //    i->par[1].ro = i->par[0].ro = 116.667 * kv(chi_real / chi_) / (x * x + y * y);
        //    i->par[1].p = i->par[0].p = kv(36.12 / (chi_real / chi_)) * (116.667 * kv(chi_real / chi_)) / (ggg * kv(10.0)) * pow(1.0 / sqrt(x * x + y * y), 2.0 * ggg);
        //}

        if (i->type == C_centr)
        {
            SS->All_Cells_zero.push_back(i);
        }
    }
    SS->M_K_prepare();     // Нужно комментить, если не считается монте-карло, там удаляются источники
    SS->MK_start_new();
    for (int k = 0; k < 40; k++)  // 10
    {
        cout << "Global step 2 = " << k + 1 << endl;
        SS->Go_stationary_5_komponent_inner_MK2(40000);
        SS->Go_5_komponent__MK2(50000);
    }

    SS->Print_Gran("surface11_3.txt");
    SS->Print_Tecplot_MK();
    SS->Save_Setka_ALL_ALPHA("vers11_3.txt");

    // ----------------------------------------------------------------------------------------------------
    delete SS;
    SS = new Setka();
    SS->Download_Setka_ALL_ALPHA_2_0("vers11_3.txt");
    SS->TVD_prepare();
    SS->Proverka();
    for (auto i : SS->All_Cells)
    {

        if (i->type == C_centr)
        {
            i->par[0].ro_H1 = i->par[1].ro_H1 = 0.00000001;
            i->par[0].p_H1 = i->par[1].p_H1 = 0.00000001;
            i->par[0].u_H1 = i->par[1].u_H1 = 0.0;
            i->par[0].v_H1 = i->par[1].v_H1 = 0.0;

            i->par[0].ro_H2 = i->par[1].ro_H2 = 0.00000001;
            i->par[0].p_H2 = i->par[1].p_H2 = 0.00000001;
            i->par[0].u_H2 = i->par[1].u_H2 = 0.0;
            i->par[0].v_H2 = i->par[1].v_H2 = 0.0;

            i->par[0].ro_H3 = i->par[1].ro_H3 = 0.00000001;
            i->par[0].p_H3 = i->par[1].p_H3 = 0.00000001;
            i->par[0].u_H3 = i->par[1].u_H3 = 0.0;
            i->par[0].v_H3 = i->par[1].v_H3 = 0.0;

            i->par[0].ro_H4 = i->par[1].ro_H4 = 0.00000001;
            i->par[0].p_H4 = i->par[1].p_H4 = 0.00000001;
            i->par[0].u_H4 = i->par[1].u_H4 = 0.0;
            i->par[0].v_H4 = i->par[1].v_H4 = 0.0;
        }


        double x, y;
        i->Get_Center(x, y);
        if (sqrt(x * x + y * y) <= R111_)
        {
            SS->All_Cells_Inner.push_back(i);
        }


        if (i->type == C_centr)
        {
            SS->All_Cells_zero.push_back(i);
        }
    }
    SS->M_K_prepare();     // Нужно комментить, если не считается монте-карло, там удаляются источники
    SS->MK_start_new();
    for (int k = 0; k < 40; k++)  // 10
    {
        cout << "Global step 2 = " << k + 1 << endl;
        SS->Go_stationary_5_komponent_inner_MK2(40000);
        SS->Go_5_komponent__MK2(50000);
    }

    SS->Print_Gran("surface11_4.txt");
    SS->Print_Tecplot_MK();
    SS->Save_Setka_ALL_ALPHA("vers11_4.txt");

    exit(-1);

    
    Setka S;


    //S.Download_Setka_ALL_ALPHA("all_save_3_122.txt");  
    /// 122 - начальные параметры газовой динамики до счёта монте-карло    ALPHA


    S.Download_Setka_ALL_ALPHA_2_0("all_save_4_177.txt"); // 125 до движения поверхностей
    S.Print_Gran();
    //S.Print_Tecplot_MK();
    //exit(-1);

    //S.Print_Sourse();
    //return 0;


    // Для выгрузки данных в равномерную 2д сетку
    if (false)
    {
        ofstream file_save;
        file_save.open("area_point_.txt");

        int N = 2048;  //1792 //1792                 // Количество ячеек по x
        int M = 1280; //1280 //1280                 // Количество ячеек по y
        int K = (N * M);                // Количество ячеек в сетке
        double x_min = -1500.0; // -2500.0 // -1300  //-2000                // -1500.0
        double x_max = 1500.0; // 450.0
        double y_max = 2000.0; // 1600.0 //1840.0
        double y_min = (y_max / (2.0 * M));
        double dx = ((x_max - x_min) / (N - 1));   // Величина грани по dx
        double dy = ((y_max) / (M));     // Величина грани по dy
        bool bk = false;
        for (int k = 0; k < K; k++)  // Заполняем начальные условия  до K
        {
            bk = false;
            if (k % 1000 == 0)
            {
                cout << k << " " << K << endl;
            }
            int n = k % N;                                   // номер ячейки по x (от 0)
            int m = (k - n) / N;                             // номер ячейки по y (от 0)
            double y = y_min + m * dy;
            double x = x_min + n * dx;
            double dist = sqrt(kvv(0.0, x, y));
            if (dist <= 44.0)
            {
                for (auto& i : S.All_Cells)
                {
                    if (i->belong(x, y) == true)
                    {
                        double xx, yy;
                        double ro, p, u, v, ro1, p1, u1, v1, Q;
                        i->Get_Center(xx, yy);
                        double dist2 = sqrt(kvv(0.0, xx, yy));

                        if (dist > 1.0)
                        {
                            ro = i->par[0].ro * kv(dist2) / kv(dist);
                            p = i->par[0].p * pow(dist2 / dist, H_pow);
                            ro1 = i->par[0].ro_H1 * pow(dist2 / dist, H_pow);
                            ////par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
                            ////par2.p_H1 = par2.p_H1 * pow(dd / dis, 2 * ggg);
                            Q = i->par[0].Q * kv(dist2) / kv(dist);
                            //u = i->par[0].u;
                            //v = i->par[0].v;
                            //u1 = i->par[0].u_H1;
                            //v1 = i->par[0].v_H1;
                            //polar_perenos(x, y, xx, yy, u, v);
                            //polar_perenos(x, y, xx, yy, u1, v1);
                        }
                        else
                        {
                            ro = i->par[0].ro;
                            p = i->par[0].p;
                            ro1 = i->par[0].ro_H1;
                            Q = i->par[0].Q;
                        }

                        file_save << ro << " " << p << " " << i->par[0].u << " " << i->par[0].v << " " //
                            << ro1 << " " << i->par[0].p_H1 << " " << i->par[0].u_H1 << " " << i->par[0].v_H1 << " " //
                            << i->par[0].ro_H2 << " " << i->par[0].p_H2 << " " << i->par[0].u_H2 << " " << i->par[0].v_H2 << " " //
                            << i->par[0].ro_H3 << " " << i->par[0].p_H3 << " " << i->par[0].u_H3 << " " << i->par[0].v_H3 << " " //
                            << i->par[0].ro_H4 << " " << i->par[0].p_H4 << " " << i->par[0].u_H4 << " " << i->par[0].v_H4 << " " //
                            << Q << endl;
                        bk = true;
                        break;
                    }
                }
            }

            if (bk == false)
            {
                file_save << 1.0 << " " << 1.0 << " " << Velosity_inf << " " << 0.0 << " " //
                    << 0.0000001 << " " << 0.0000001 << " " << 0.0 << " " << 0.0 << " " //
                    << 0.0000001 << " " << 0.0000001 << " " << 0.0 << " " << 0.0 << " " //
                    << 0.0000001 << " " << 0.0000001 << " " << 0.0 << " " << 0.0 << " " //
                    << 1.0 << " " << 0.5 << " " << Velosity_inf << " " << 0.0 << " " << 100.0 << endl;
            }
            
        }
        file_save.close();
        return 0;
    }

    /// 122 - начальные параметры газовой динамики до счёта монте-карло    ALPHA
    /// 126 - счёт монте-карло с большим числом частиц    ALPHA_2_0

    // Если нужно перенормировать параметры в сетке
    if (false)
    {
        for (auto& i : S.All_Cells)
        {
            /*i->par[0].ro_H2 = i->par[0].ro_H2 * 0.5;
            i->par[0].p_H2 = i->par[0].p_H2 * 0.5;
            i->par[0].u_H2 = i->par[0].u_H2 * 0.7;
            i->par[0].v_H2 = i->par[0].v_H2 * 0.7;*/

            if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
            {
                //Для счёта монте-карло
               //i->par[0].u = i->par[0].u * (chi_real / chi_);       // Перенормировка
               //i->par[0].v = i->par[0].v * (chi_real / chi_);
               //i->par[0].ro = i->par[0].ro / kv(chi_real / chi_);
               //i->par[0].Q = i->par[0].Q / kv(chi_real / chi_);

               //i->par[1].u = i->par[1].u * (chi_real / chi_);       // Перенормировка
               //i->par[1].v = i->par[1].v * (chi_real / chi_);
               //i->par[1].ro = i->par[1].ro / kv(chi_real / chi_);
               //i->par[1].Q = i->par[1].Q / kv(chi_real / chi_);


               // для счёта плазмы
               //i->par[0].u = i->par[0].u / (chi_real / chi_);       // Перенормировка
               //i->par[0].v = i->par[0].v / (chi_real / chi_);
               //i->par[0].ro = i->par[0].ro * kv(chi_real / chi_);
               //i->par[0].Q = i->par[0].Q * kv(chi_real / chi_);
            }

            //double x, y;
            //i->Get_Center(x, y);
            //double dist = sqrt(x * x + y * y);

            //if (dist < 30.0)//(i->type == C_centr)
            //{
            //    i->par[0].ro_H1 = 0.000001066666;
            //    i->par[0].p_H1 = (0.000001066666 * chi_real * chi_real / (ggg * 14.1344 * 14.1344));
            //    i->par[0].u_H1 = chi_real * x / dist;
            //    i->par[0].v_H1 = chi_real * y / dist;
            //    i->par[1].ro_H1 = 0.000001066666;
            //    i->par[1].p_H1 = (0.000001066666 * chi_real * chi_real / (ggg * 14.1344 * 14.1344));
            //    i->par[1].u_H1 = chi_real * x / dist;
            //    i->par[1].v_H1 = chi_real * y / dist;
            //}


        }
    }

    //cout << min << " " << max << endl;


    //Setka S = Setka(14, 5, 5, 5, 7, 10, 7, 8);
    //Setka S = Setka(30, 12, 13, 30, 40, 20, 20, 15);  // Нужно чтобы количиство ячеек по углу делилось на 8 или 10 (не было простым)


    //S.Print_point();
    //S.Print_cell2();
    /*S.Print_Gran();
    S.Print_cell_type();
    S.Print_connect();
    S.Proverka();
    S.Print_Gran_type();*/
    //S.Save_Setka_ALL_ALPHA("aaa.txt");
    //S.Print_point_connect();
    //S.Print_cell2();
    //S.Init_conditions();
    //S.Go_stationary(350000);
    //S.Save_G_D();
    //S.Download_G_D_5_komponent();
    //S.Go_stationary_5_komponent(200000);
    //S.Go_5_komponent(1);
    //S.Go_stationary(100000);

    /*for (auto& i : S.All_Cells)
    {
        double x, y;
        i->Get_Center(x, y);
        if (x < -600 && y < 150)
        {
            i->par[1].u = i->par[0].u = -3.0;
        }
    }*/
    /*S.Move_Setka_Calculate(0);
    for (auto& i : S.All_Points)
    {
        i->x = i->x2;
        i->y = i->y2;
        i->Vx = 0.0;
        i->Vy = 0.0;
    }*/

    S.TVD_prepare();
    S.M_K_prepare();   // Нужно комментить, если не считается монте-карло, т.к. там удаляются источники
    //S.Print_TVD();


    S.Proverka();
    //S.Print_cell_type();

    // Блок просмотра
    if (false)
    {
        ofstream fout;
        string name_f = "111.txt";
        fout.open(name_f);

        //int jk = 0;
        //for (auto& i : S.All_Gran)
        //{
        //    if (i->Sosed_down != nullptr)
        //    {
        //        jk++;
        //    }
        //}

        //fout << "TITLE = \"HP\" ";
        //fout << " VARIABLES = \"X\", \"Y\", \"Z\"  ZONE T= \"HP\", N=  " << jk * 2;
        //fout << " , E= " << jk;
        //fout << " , F=FEPOINT, ET=LINESEG  " << endl;
        //int jkk = 0;
        //for (auto& i : S.All_Gran)
        //{
        //    if (i->Sosed_down != nullptr)
        //    {
        //        double x, y, n1, n2;
        //        i->Get_Center(x, y);
        //        i->Sosed_down->Get_Center(n1, n2);
        //        fout << x << " " << y << " " << jkk << endl;
        //        fout << n1 << " " << n2 << " " << jkk << endl;
        //        jkk++;
        //        /*if (x < 37.2 && x > 36.8 && y < 1.4 && y > 1.0)
        //        {
        //            cout << x << " " << y << endl;
        //            cout << n1 << " " << n2 << endl;
        //        }*/
        //    }
        //}
        //for (int i = 0; i < jk; i++)
        //{
        //    fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
        //}

            for (auto& i : S.All_Gran)
        {
            if (i->Sosed_down == nullptr)
            {
                double x, y, n1, n2;
                i->Get_Center(x, y);
                i->Get_normal(n1, n2);
                fout << x << " " << y << " " << 2 << endl;
                fout << x + 0.1 * n1 << " " << y + 0.1 * n2 << " " << 2 << endl;
            }
        }
        exit(-1);
    }

    


    // Если мы хотим подвинуть сетку до начала счёта то помогает следующий блок кода
   /* S.Move_Setka_Calculate(0);
    for (auto& i : S.All_Points)
    {
        i->x = i->x2;
        i->y = i->y2;
        i->Vx = 0.0;
        i->Vy = 0.0;
    }*/


    // Подготовка массивов для внутренней области счёта. НЕ УДАЛЯТЬ
    for (auto i : S.All_Cells)
    {
        double x, y;
        i->Get_Center(x, y);
        if (sqrt(x * x + y * y) <= R111_)
        {
            S.All_Cells_Inner.push_back(i);
        }

        if (i->type == C_centr)
        {
            S.All_Cells_zero.push_back(i);
        }
    }

    cout << "V nule yacheek = " << S.All_Cells_zero.size() << endl;

    //S.Go_stationary_5_komponent_inner(10000);

    //S.Go_stationary_5_komponent_inner(30000);
    //S.Go_stationary_5_komponent_inner(250000);
    //S.Print_Gran();
    //S.Go_stationary_5_komponent_inner(250000);
    //S.Move_Setka_Calculate(0);

    //S.Init_conditions();

    // Монте-карло блок
    S.MK_start_new();
    /*cout << 1.0 * S.mmu1 / (1.0 * S.mn1) << " " << S.mn1 << endl;
    cout << S.mmu2 / (1.0 * S.mn2) << " " << S.mn2 << endl;
    cout << S.mmu3 / (1.0 * S.mn3) << " " << S.mn3 << endl;
    cout << S.mmu4 / (1.0 * S.mn4) << " " << S.mn4 << endl;
    cout << S.mmu5 / (1.0 * S.mn5) << " " << S.mn5 << endl;
    for (int i = 0; i < 270; i++)
    {
        cout << S.Sensors2[i]->num << endl;
    }*/
    //S.Go_stationary_5_komponent_inner(200000);

    ofstream ffout;
    ffout.open("outer_.txt");
    for (int i = 1; i <= 0; i++)
    {
        cout << "ITER = " << i << "\\" << 60 << endl;
        //S.Go_stationary_5_komponent_inner_MK(50000);
        S.Go_stationary_5_komponent_inner(70000);
        //S.Go_stationary_5_komponent_inner(200000);
        //S.Go_5_komponent_MK(50000);
        S.Go_5_komponent(50000);
        S.Print_Gran();
        //S.Print_Tecplot_MK();
        ffout << "Inner = " << S.Line_Inner[0]->A->x << endl;
        ffout << "Contact = " << S.Line_Contact[0]->A->x << endl;
        ffout << "Outer = " << S.Line_Outer[0]->A->x << endl;
    }
    ffout.close();
    //S.Go_5_komponent(150000);
    //S.Go_5_komponent(100000);
    //S.Go_5_komponent(300000);
    //S.Go_5_komponent(5000000);
    //S.Go_5_komponent(500000);

    S.Print_Tecplot_MK();
    S.Print_cell2();
    S.Print_Gran();
    S.Print_connect();
    S.Save_Setka_ALL_ALPHA("all_save_4_178.txt");
    S.Print_Sourse();
}
