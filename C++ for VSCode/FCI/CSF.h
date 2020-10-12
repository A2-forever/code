#ifndef CSF_H_
#define CSF_H_

#include "FCI.h"

class CSF{
    //构建所有符合条件的CI组态，输入指定电子数，轨道数与自旋z分量，CI组态存于AI_Array中
    friend int CSF_new(int nelec_ex, int nOrb_ex, const double &MS, std::vector<int> &Orbital_ex, std::vector<CSF> &CSF_Array);

    private:
        int nOrb = 0;
        int nelec = 0;
        double MS = 0;
        double S = 0;
        std::vector<int> CSF;//一个四进制数，表示
        std::vector<Slater_det> Slater_CI;
        std::vector<double> coefficient;

        bool Slater_create(const double &MS);

    public:
        CSF();                                   //默认构造函数
        CSF(std::vector<int> &Orbital_ex, const double &MS);                                   //默认构造函数
        ~CSF();                                  //默认析构函数
        bool CSF2Slater(std::vector<int> &Orbital_ex);
        int vector2Slater(std::vector<int> &Orbital_ex, std::vector<int> &tmp_Slater, std::vector<int> &single);
};


class CI{
    private:
        int nOrb = 0;          //分子轨道数
        double h_nuc = 0;      //核积分项
        std::vector<double> h; //单电子积分
        std::vector<double> g; //双电子积分

        int dim = 0;
        int dim2 = 0;
        int dim3 = 0;

    public:
        CI(); //默认构造函数
        CI(const double &h_nuc_ex, const std::vector<double> &h_ex, const std::vector<double> &g_ex, const int nOrb_ex);
        ~CI(); //默认析构函数
        double get(const int &i, const int &j, const int &k, const int &l);

        //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
        double H_ij(Slater_det &k1, Slater_det &k2);
        double H(Slater_det &k1, Slater_det &k2);

        //计算分子轨道i与分子轨道j关于hamilton算符的耦合项单电子部分
        //sigma表示电子自旋，0为alpha，1为beta
        double F_ij(Slater_det &k1, Slater_det &k2, const int sigma, const std::vector<int> &Num, const std::vector<int> &index);

        //计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
        //sigma表示电子自旋，0为alpha，1为beta
        double G_ij(Slater_det &k1, Slater_det &k2, const std::vector<int> &Num, const std::vector<int> &index);
        double G1_ij(Slater_det &k1, Slater_det &k2, const int sigma, const std::vector<int> &Num, const std::vector<int> &index); //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
        double G2_ij(Slater_det &k1, Slater_det &k2, const std::vector<int> &Num, const std::vector<int> &index);                  //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
};



int sgn(const int &n)//奇数返回-1，偶数返回1
{
    return n%2==0 ? 1:-1; 
};


#endif