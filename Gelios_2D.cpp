#include "Help.h"



int main()
{
    std::cout << "Hello World!\n";
    Setka S = Setka();
    S.Download_Setka_ALL_ALPHA("all_save_9.txt");

    //Setka S = Setka(14, 5, 5, 5, 7, 10, 7, 8);
    //Setka S = Setka(35, 12, 13, 45, 30, 20, 20, 15);  // Нужно чтобы колучиство ячеек по углу делилось на 8 или 10 (не было простым)


    S.Contact.push_back(new Point(200, 0));
    S.Contact.push_back(new Point(197, 58));
    S.Contact.push_back(new Point(174, 138));
    S.Contact.push_back(new Point(129, 217));
    S.Contact.push_back(new Point(28, 315));
    S.Contact.push_back(new Point(15, 326));
    S.Contact.push_back(new Point(-86, 380));
    S.Contact.push_back(new Point(-150, 407));
    S.Contact.push_back(new Point(-286, 450));
    S.Contact.push_back(new Point(-400, 450));
    S.Contact.push_back(new Point(-600, 450));


    S.Inner.push_back(new Point(109, 0));
    S.Inner.push_back(new Point(102, 45));
    S.Inner.push_back(new Point(82, 87));
    S.Inner.push_back(new Point(44, 125));
    S.Inner.push_back(new Point(-21, 158));
    S.Inner.push_back(new Point(-110, 164));
    S.Inner.push_back(new Point(-177, 130));
    S.Inner.push_back(new Point(-224, 83));
    S.Inner.push_back(new Point(-254, 0));
    
    S.Outer.push_back(new Point(370, 0));
    S.Outer.push_back(new Point(364, 84));
    S.Outer.push_back(new Point(313, 277));
    S.Outer.push_back(new Point(290, 334));
    S.Outer.push_back(new Point(214, 468));
    S.Outer.push_back(new Point(73, 566));
    S.Outer.push_back(new Point(0, 595));

    //S.Print_point();
    //S.Print_cell2();
    /*S.Print_Gran();
    S.Print_cell_type();
    S.Print_connect();
    S.Proverka();
    S.Print_Gran_type();*/
    //S.Save_Setka_ALL_ALPHA("aaa.txt");
    //S.Print_point_connect();

    //S.Init_conditions();
    //S.Go_stationary(350000);
    //S.Save_G_D();
    //S.Download_G_D_5_komponent();
    S.Go_stationary_5_komponent(300000);
    S.Go_5_komponent(25000);
    S.Print_Tecplot();
    S.Print_cell2();
    S.Save_Setka_ALL_ALPHA("all_save_10.txt");
    //S.Save_G_D_5_komponent();
    

}
