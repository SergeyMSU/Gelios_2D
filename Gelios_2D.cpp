#include "Help.h"



int main()
{
    std::cout << "Hello World!\n";
    Setka S = Setka();
    S.Download_Setka_ALL_ALPHA("aaa.txt");

    //Setka S = Setka(14, 5, 5, 5, 7, 10, 7, 8);
    //Setka S = Setka(35, 12, 13, 45, 30, 20, 20, 15);  // Нужно чтобы колучиство ячеек по углу делилось на 8 или 10 (не было простым)

    S.Print_point();
    S.Print_cell2();
    S.Print_Gran();
    S.Print_cell_type();
    S.Print_connect();
    S.Proverka();
    S.Print_Gran_type();
    //S.Save_Setka_ALL_ALPHA("aaa.txt");
    //S.Print_point_connect();

    //S.Init_conditions();
    //S.Go_stationary(350000);
    //S.Save_G_D();
    ////S.Download_G_D_5_komponent();
    //S.Go_stationary_5_komponent(600000);
    //S.Print_Tecplot();
    //S.Save_G_D_5_komponent();
    

}
