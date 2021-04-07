#include "Help.h"



int main()
{
    std::cout << "Hello World!\n";
    Setka S = Setka(14, 5, 5, 5, 7, 10, 7, 8);
    //Setka S = Setka(35, 12, 12, 45, 30, 20, 20, 15);

    S.Print_point();
    S.Print_cell2();
    S.Print_Gran();
    S.Print_cell_type();
    S.Print_connect();
    S.Proverka();
    //S.Init_conditions();
    //S.Go_stationary(250000);
    //S.Save_G_D();
    //S.Download_G_D_5_komponent();
    /*for (auto& i : S.All_Cells)
    {
        if (i->par[0].Q / i->par[0].ro < 90)
        {
            i->par[0].u = i->par[0].u * chi_real;
            i->par[0].v = i->par[0].v * chi_real;
            i->par[0].ro = i->par[0].ro / kv(chi_real);
            i->par[0].p = i->par[0].p;
            i->par[0].Q = i->par[0].Q / kv(chi_real);
        }
    }*/
    /*S.Go_stationary_5_komponent(100000);
    S.Print_Tecplot();
    S.Save_G_D_5_komponent();*/
    

}
