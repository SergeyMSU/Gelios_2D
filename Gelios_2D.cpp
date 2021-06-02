#include "Help.h"



int main()
{
    std::cout << "Hello World!\n";
    Setka S = Setka();
    S.Download_Setka_ALL_ALPHA("all_save_3_123.txt");   // 96 

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
    //S.Init_conditions();
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
    //S.Print_TVD();

    /*auto A = S.All_Cells[2]->Grans[1];
    double x1, x2, x3, x4, y1, y2, y3, y4;
    double n1, n2;
    A->Get_Center(x1, y1);
    A->Get_normal(n1, n2);
    cout << x1 << " " << y1 << endl;
    cout << n1 << " " << n2 << endl;
    A->Sosed_down->Get_Center(x2, y2);
    cout << x2 << " " << y2 << endl;
    A->Gran_copy->Sosed_down->Get_Center(x3, y3);
    cout << x3 << " " << y3 << endl;*/

    S.Proverka();
    //S.Print_cell_type();



    // Если мы хотим подвинуть сетку до начала счёта то помогает следующий блок кода
   /* S.Move_Setka_Calculate(0);
    for (auto& i : S.All_Points)
    {
        i->x = i->x2;
        i->y = i->y2;
        i->Vx = 0.0;
        i->Vy = 0.0;
    }*/


    for (auto i : S.All_Cells)
    {
        double x, y;
        i->Get_Center(x, y);
        if (sqrt(x * x + y * y) <= R111_)
        {
            S.All_Cells_Inner.push_back(i);
        }

        /*if (i->type == C_1 || i->type == C_centr)
        {
            i->par[0].u_H2 = 0.0;
            i->par[0].v_H2 = 0.0;
            i->par[0].ro_H2 = 0.0003;
            i->par[1] = i->par[0];
        }*/
    }

    //S.Go_stationary_5_komponent_inner(10000);

    //S.Go_stationary_5_komponent_inner(30000);
    //S.Go_stationary_5_komponent_inner(250000);
    //S.Print_Gran();
    //S.Go_stationary_5_komponent_inner(250000);
    //S.Move_Setka_Calculate(0);
    //S.Init_conditions();
    //S.Go_stationary_5_komponent_inner(200000);
    cout << "sin = " << sin(0.5 * alpha_rot * pi_ / 180.0) << endl;
    for (int i = 0; i < 1; i++)
    {
        //S.Go_stationary_5_komponent_inner(100000);
        S.Go_5_komponent(400000);
        //S.Go_stationary_5_komponent_inner(50000);
    }
    //S.Go_5_komponent(100000);
    //S.Go_5_komponent(300000);
    //S.Go_5_komponent(5000000);
    //S.Go_5_komponent(500000);
    S.Print_Tecplot();
    S.Print_cell2();
    S.Print_Gran();
    S.Print_connect();
    S.Save_Setka_ALL_ALPHA("all_save_3_124.txt");
    

}
