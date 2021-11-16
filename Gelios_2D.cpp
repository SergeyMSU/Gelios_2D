#include "Help.h"
#include <iomanip>


int main()
{
    std::cout << "Hello World!\n";
    /*double Vx, Vy, Vz;
    dekard_skorost(y, z, x, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);
    dekard_skorost(100.0, 100.0, 100.0, -1.0, 0.0, 0.1, Vy, Vz, Vx);*/
    //cout << Vx << " " << Vy << " " << Vz << endl;

    /*Setka A;
    Sensor* sens = new Sensor(1, 2, 3);
    MKmethod MK = MKmethod();
    cout << A.Velosity_1(0.0000001, 3.0) << endl;
    cout << sigma(A.Velosity_1(0.0000001, 3.0)) << endl;
    exit(-1);*/

    /*double Sum5 = 0.0;
    for (int i = 0; i < 10000000; i++)
    {
        Sum5 += MK.play_mho(sens, 3.346573);
    } 
    cout << Sum5 / (1.0 * 10000000) << endl;*/ 
    //cout << MK.play_mho(sens, 0.0) << endl;
    /*double a = 0.00000000000001;
    double b = 0.00000000001;
    double c = 1.0;

    cout << scientific << setprecision(16) << 2.0 * a + b + c << endl;

    exit(-10);*/

    // Делаем датчики

    //ofstream fout;
    //string name_f = "datchik3.txt";
    //fout.open(name_f);

    //Sensor* sens1 = new Sensor(15236, 12353, 96755);
    ///*Sensor* sens2 = new Sensor(0, 1, 0);
    //Sensor* sens3 = new Sensor(15236, 12353, 96755);
    //fout << sens1->a1_ << " " << sens1->a2_ << " " << sens1->a3_ << endl;
    //fout << sens2->a1_ << " " << sens2->a2_ << " " << sens2->a3_ << endl;
    //fout << sens3->a1_ << " " << sens3->a2_ << " " << sens3->a3_ << endl;*/
    //fout << sens1->a1_ << " " << sens1->a2_ << " " << sens1->a3_ << endl;
    //int a, b, c;
    //for (int i = 0; i < 272; i++)
    //{
    //    cout << i << endl;
    //    for (int j = 0; j < 1000000000; j++)
    //    {
    //        sens1->MakeRandom();
    //        if (sens1->a1_ == 15236 && sens1->a2_ == 12353 && sens1->a3_ == 96755)
    //        {
    //            cout << "STOP" << endl;
    //            exit(-1);
    //        }
    //        if (sens1->a1_ == 1 && sens1->a2_ == 0 && sens1->a3_ == 123)
    //        {
    //            cout << "STOP" << endl;
    //            exit(-1);
    //        }
    //        if (sens1->a1_ == 0 && sens1->a2_ == 1 && sens1->a3_ == 0)
    //        {
    //            cout << "STOP" << endl;
    //            exit(-1);
    //        }

    //    }
    //    fout << sens1->a1_ << " " << sens1->a2_ << " " << sens1->a3_ << endl;
    //    //fout << sens2->a1_ << " " << sens2->a2_ << " " << sens2->a3_ << endl;
    //    //fout << sens3->a1_ << " " << sens3->a2_ << " " << sens3->a3_ << endl;
    //}
    //fout.close();
    //exit(-1);

    //vector<double> v1(6);
    //vector<double> v2(6);
    //vector<double> v3(6);
    //vector<double> v4(6);
    //vector<double> v5(6);
    //double Sum1 = 0.0;
    //double Sum2 = 0.0;
    //double Sum3 = 0.0;
    //double Sum4 = 0.0;
    //double Sum5 = 0.0;
    //double Sum6 = 0.0;
    //double m1 = 0.0;
    //double m2 = 0.0;
    //double m3 = 0.0;
    //double m4 = 0.0;
    //double m5 = 0.0;
    //double m6 = 0.0;

    //cout << MK.f2(-2.0, 0.0, -1.0, 0.1) << endl;

    /*for (int i = 0; i < 100000; i++)
    {
        MK.Init_Parametrs(sens, v1, v2, v3, v4, v5);
        Sum6 += v5[5];
        m6 += 1.0;
        Sum5 += v5[4];
        m5 += 1.0;
        Sum1 += v5[3];
        m1 += 1.0;
        Sum2 += v5[2];
        m2 += 1.0;
        Sum3 += v5[1];
        m3 += 1.0;
        Sum4 += v5[0];
        m4 += 1.0;
    }

    cout << Sum6 / m6 << endl;
    cout << Sum5 / m5 << endl;
    cout << Sum1 / m1 << endl;
    cout << Sum2 / m2 << endl;
    cout << Sum3 / m3 << endl;
    cout << Sum4 / m4 << endl;*/

    //exit(-1);

   // cout << MK.Get_Int002(-4.47265e-05, 1.23034) << endl;
    /*double a, b, c;
    double a2, b2, c2;
    spherical_skorost(1, 1, 1, 1, 2, 3, a, b, c);
    spherical_skorost(1, 1, 1, 3, 4, 1, a2, b2, c2);
    cout << sqrt(kvv(1, 2, 3)) << endl;
    cout << sqrt(kvv(a, b, c)) << endl;
    cout << sqrt(kvv(1 - 3, 2 - 4, 3 - 1)) << endl;
    cout << sqrt(kvv(a - a2, b - b2, c - c2)) << endl;
    exit(-1);*/
    /*MK.Change_Velosity2(sens, -2.3, 0.3, 0.0, -1.1, 0.2, 0.03, v1, v2, v3, v4, 1.0, 800.0, 5);
    for (int i = 0; i < 6; i++)
    {
        cout << v4[i] << endl;
    }*/
    

    Setka S;

    // Блок для проверки вычисления интеграллов (лучше не удалять, полезный)
    //double Vx = -3.0;
    //double Vy = 1.1;
    //double Vz = -0.2;
    //double cp = 2.8;
    //double vx = -1.9;
    //double vy = 0.1;
    //double vz = -0.2;
    //double alpha = 0.0;
    //double u = sqrt(kvv(Vx - vx, Vy - vy, Vz - vz));
    //double uz = S.Velosity_1(u, cp);
    //double uz_M = S.Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
    //double uz_E = S.Velosity_3(u, cp);
    //double u1 = vx - Vx;
    //double u2 = vy - Vy;
    //double u3 = vz - Vz;
    //double skalar = Vx * u1 + Vy * u2 + Vz * u3;

    //cout << -uz_M * (sigma(uz_M) / sigma(uz)) * u1 / u << " " << -uz_M * u1 / u <<  endl;
    //cout << (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma(uz_E) / sigma(uz)) - //
    //    uz_M * (sigma(uz_M) / sigma(uz)) * skalar / u) << " " << (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz)  - //
    //        uz_M * skalar / u) << endl;
    //cout << u << " " << uz << " " << uz_M << " " << uz_E << endl;

    //exit(-1);

    //S.Download_Setka_ALL_ALPHA("all_save_3_122.txt");  
    /// 122 - начальные параметры газовой динамики до счёта монте-карло    ALPHA


    S.Download_Setka_ALL_ALPHA_2_0("all_save_4_138.txt"); // 125 до движения поверхностей
    //S.Print_Tecplot_MK();
    //exit(-1);

    /// 122 - начальные параметры газовой динамики до счёта монте-карло    ALPHA
    /// 126 - счёт монте-карло с большим числом частиц    ALPHA_2_0

    //S.Print_Tecplot_MK();
    // Нужно перенормировать параметры плазмы
    //double min = 100.0, max = 0.0, cp;
    //for (int i = 0; i < 95 ; i++)
    //{
    //    //cout << i << " " << S.All_Cells[i]->contour[0]->x << endl;
    //    S.All_Cells[i]->par[0].H_n[0] = S.All_Cells[i + 95]->par[0].H_n[0];
    //    S.All_Cells[i]->par[0].H_n[1] = S.All_Cells[i + 95]->par[0].H_n[1];
    //    S.All_Cells[i]->par[0].H_n[2] = S.All_Cells[i + 95]->par[0].H_n[2];
    //    S.All_Cells[i]->par[0].H_n[3] = S.All_Cells[i + 95]->par[0].H_n[3];

    //    S.All_Cells[i]->par[0].H_u[0] = S.All_Cells[i + 95]->par[0].H_u[0];
    //    S.All_Cells[i]->par[0].H_u[1] = S.All_Cells[i + 95]->par[0].H_u[1];
    //    S.All_Cells[i]->par[0].H_u[2] = S.All_Cells[i + 95]->par[0].H_u[2];
    //    S.All_Cells[i]->par[0].H_u[3] = S.All_Cells[i + 95]->par[0].H_u[3];

    //    S.All_Cells[i]->par[0].H_v[0] = S.All_Cells[i + 95]->par[0].H_v[0];
    //    S.All_Cells[i]->par[0].H_v[1] = S.All_Cells[i + 95]->par[0].H_v[1];
    //    S.All_Cells[i]->par[0].H_v[2] = S.All_Cells[i + 95]->par[0].H_v[2];
    //    S.All_Cells[i]->par[0].H_v[3] = S.All_Cells[i + 95]->par[0].H_v[3];

    //    S.All_Cells[i]->par[0].H_T[0] = S.All_Cells[i + 95]->par[0].H_T[0];
    //    S.All_Cells[i]->par[0].H_T[1] = S.All_Cells[i + 95]->par[0].H_T[1];
    //    S.All_Cells[i]->par[0].H_T[2] = S.All_Cells[i + 95]->par[0].H_T[2];
    //    S.All_Cells[i]->par[0].H_T[3] = S.All_Cells[i + 95]->par[0].H_T[3];

    //}   

    for (auto& i : S.All_Cells)
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
            //i->par[0].u = i->par[0].u / (chi_real / chi_);       // Перенормировка
            //i->par[0].v = i->par[0].v / (chi_real / chi_);
            //i->par[0].ro = i->par[0].ro * kv(chi_real / chi_);
            //i->par[0].Q = i->par[0].Q * kv(chi_real / chi_);
        }
        //Поправим коэффициенты в источниках
        /*if (fabs(i->par[0].k_u) > 2.0 || (i->par[0].k_u) < 0.5)
        {
            i->par[0].k_u = 1.0;
        }
        if (fabs(i->par[0].k_v) > 2.0 || (i->par[0].k_v) < 0.5)
        {
            i->par[0].k_v = 1.0;
        }
        if (fabs(i->par[0].k_T) > 2.0 || (i->par[0].k_T) < 0.5)
        {
            i->par[0].k_T = 1.0;
        }*/
        /*cp = sqrt(i->par[0].p / i->par[0].ro);
        cout << i->contour[0]->x << " " << i->contour[0]->y << " " << cp << endl;
        if (cp < min)
        {
            min = cp;
        }
        if (cp > max)
        {
            max = cp;
        }*/
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
    //S.M_K_prepare();   // Нужно комментить, если не считается монте-карло, т.к. там удаляются источники
    //S.Print_TVD();


    S.Proverka();
    //S.Print_cell_type();
    ofstream fout;
    string name_f = "111.txt";
    fout.open(name_f);
    for (auto& i : S.All_Gran)  
    {
        if (i->Gran_copy == nullptr)
        {
            double x, y, n1, n2;
            i->Get_Center(x, y);
            i->Get_normal(n1, n2);
            fout << x << " " << y << " " << 1 << endl;
            fout << x + 0.1 * n1 << " " << y + 0.1 * n2 << " " << 1 << endl;
        }
    }
    for (auto& i : S.All_Gran_copy)
    {
        if (i->Gran_copy == nullptr)
        {
            double x, y, n1, n2;
            i->Get_Center(x, y);
            i->Get_normal(n1, n2);
            fout << x << " " << y << " " << 2 << endl;
            fout << x + 0.1 * n1 << " " << y + 0.1 * n2 << " " << 2 << endl;
        }
    }
    exit(-1);


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
    //S.MK_start_new();
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

    for (int i = 0; i < 1; i++)
    {
        //S.Go_stationary_5_komponent_inner_MK(20000);
        //S.Go_stationary_5_komponent_inner(100000);
        //S.Go_stationary_5_komponent_inner(200000);
        //S.Go_5_komponent_MK(200000);
        S.Go_5_komponent(50000);
    }
    //S.Go_5_komponent(150000);
    //S.Go_5_komponent(100000);
    //S.Go_5_komponent(300000);
    //S.Go_5_komponent(5000000);
    //S.Go_5_komponent(500000);

    S.Print_Tecplot_MK();
    S.Print_cell2();
    S.Print_Gran();
    S.Print_connect();
    S.Save_Setka_ALL_ALPHA("all_save_4_139.txt");
    

}
