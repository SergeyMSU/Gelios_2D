#include "Help.h"
#include <iomanip>


int main()
{
    std::cout << "Program prepared by Korolkov Sergey. All rights reserved!\n";

    //Проверка распадника
    //double P[4];
    //double PQ;
    //auto RR = new Setka();
    //double a, b, c;
    //double ro1 = 1.0;
    //double p1 = 1.0;
    //double u1 = 1.0;
    //double ro2 = 2.0;
    //double p2 = 1.0;
    //double u2 = 1.0;
    //double n1 = 1.0;
    //double w = 3.1;


    //RR->HLLC_2d_Korolkov_b_s(ro1, 1.0, p1, u1, 0.0,// 
    //    ro2, 1.0, p2, u2, 0.0,//
    //    w, P, PQ, n1, 0.0, 10.0, 1, a, b, c);
    //for (int i = 0; i < 4; i++)
    //{
    //    cout << P[i] << endl;
    //}

    //exit(-1);

    //ofstream fout;
    //fout.open("Source_3_cp12.txt");

    //Setka* CC;
    //CC = new Setka(4, 3, 3, 10, 13, 6, 7, 7);

   /* MKmethod* MK;
    MK = new MKmethod();
    double c = 1.2;
    double x = 0.7;
    cout << MK->for_Wr_2(0.0, 0.4, 1.1, 4.0) - MK->for_Wr_2(0.0, 0.1, 1.1, 4.0) << endl;*/

    //double u = 0.0, cp = 12.0;
    //for (double u = 0.05; u/cp <= 7.0; u = u + 0.05)
    //{
    //    double skalar = 20.0;
    //    double u1 = 0.8;
    //    double uz = CC->Velosity_1(u, cp);
    //    double uz_M = CC->Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
    //    double uz_E = CC->Velosity_3(u, cp);
    //    double AA = (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz));
    //    double BB = (-0.5 * MK->int_3(u, cp) / MK->int_1(u, cp));
    //    cout << u << " " << u/MK->int_1(u, cp) << endl;
    //    //cout << u << " " << BB << endl;
    //}


    //exit(-1);

    // SS = new Setka(15, 5, 7, 15, 30, 10, 20, 20);   //  vers2_15   n_inner 15
    // SS = new Setka(30, 7, 11, 20, 60, 20, 40, 40);    // vers3_   n_inner 30
    // SS = new Setka(34, 7, 11, 20, 60, 20, 50, 60);    // vers4_7    n_inner 30   результаты одинаковыс с 3
    // SS = new Setka(40, 10, 15, 25, 60, 30, 60, 60);    // vers5_7   n_inner 30   не досчитал, вроде какие-то проблемы с устойчивостью
    // Привет
   /* double u = 0.0;
    double v = 1.0;
    polar_provorot(pi_ * alpha_rot / 360.0, u, v);
    cout << u << " " << v << endl;
    exit(-1);*/

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
    
    //SS->Download_Setka_ALL_ALPHA_2_0("vers6_100.txt");  // 17       IEX
    SS->Download_Setka_ALL_ALPHA_2_0("vers6_130.txt");  // 17       IEX
    SS->Init_conditions();

    SS->TVD_prepare();
    SS->Proverka();
    //SS->Print_cell2();
    //SS->Print_Gran();
    //SS->Print_Gran("sur7_101.txt");
    //exit(-2);

    /*double xx, yy;
    SS->All_Cells[34]->Get_Center(xx, yy);
    cout << xx << " " << yy << endl;
    cout << SS->All_Cells[34]->number << endl;
    exit(-1);*/

    /*SS->Line_Outer[2]->B->Vx = -0.0028;

    SS->Move_Setka_Calculate_2(1.0);
    for (auto& i : SS->All_Points)
    {
        i->x = i->x2;
        i->y = i->y2;
        i->Vx = 0.0;
        i->Vy = 0.0;
        i->count = 0;
    }*/

    // Подготовка массивов для внутренней области счёта. НЕ УДАЛЯТЬ
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


    for (int k = 0; k < 0; k++)  // 10
    {
        cout << "Global step = " << k + 1 << endl;
        //SS->Go_stationary_5_komponent_inner_2(50000);
        //SS->Go_5_komponent_2(50000);
        SS->Go_stationary_5_komponent_inner_MK2(40000);
        SS->Go_5_komponent__MK2(50000);
    }

    //SS->Print_cell2();
    SS->Print_Gran("surface6_131.txt");
    SS->Print_Tecplot_MK();
    //SS->Print_Sourse();
    //SS->Save_Setka_ALL_ALPHA("vers6_107.txt");
    SS->Save_Setka_ALL_ALPHA("vers6_131.txt");

    exit(-1);
    delete SS;
    SS = new Setka();
    SS->Download_Setka_ALL_ALPHA_2_0("vers6_124.txt");
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
    for (int k = 0; k < 50; k++)  // 10
    {
        cout << "Global step 2 = " << k + 1 << endl;
        SS->Go_stationary_5_komponent_inner_MK2(30000);
        SS->Go_5_komponent__MK2(50000);
    }

    SS->Print_Gran("surface6_125.txt");
    SS->Print_Tecplot_MK();
    SS->Save_Setka_ALL_ALPHA("vers6_125.txt");

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
