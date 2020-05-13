#ifndef FCI_H_
#define FCI_H_

#include <cmath>
#include <vector>
#include <string>
#include "Array.cpp"

class Slater_det{
    //判断分子轨道i与分子轨道j之间alpha或beta电子的差别,index存放位置，N表示alpha和beta电子的数目
    friend bool find(Slater_det &k1, Slater_det &k2, int *Num, int ***index);
    friend std::ostream &operator<<(std::ostream &os, const Slater_det &k);

private:
    int nOrb=0;                             //分子轨道数
    int nelec=0;                            //占据电子数
    int **Orbital=NULL;                     //1表示alpha轨道，2表示beta轨道,此数组表示分子轨道的电子占据，0代表空，1代表占据电子
    int **nelec_occ=NULL;                   //1表示alpha轨道，2表示beta轨道，此数组表示分子轨道i前的轨道的电子占据数
                                            //eg:nelec[i][1]代表分子轨道i前的占据的beta电子数

    public:
    Slater_det();                                                //默认构造函数
    Slater_det(const int nelec_ex, const int nOrb_ex, int **Orbital_ex);//将轨道信息写入k轨道中,nelec_ex表示占据电子数
    ~Slater_det();                                               //默认析构函数
    int Orb(const int I, const int sigma);                      //用于调用分子轨道的占据情况
    bool occ_cal();                                              //计算轨道i前的轨道的电子占据数
    int gamma(const int I,const int sigma);                      //组态，I_n位置的Gamma

};


class CI{
    private:
    int nOrb=0;                            //分子轨道数
    double h_nuc=0;                     //核积分项
    double **h=NULL;                    //指向单电子积分的指针
    double ****g=NULL;                  //指向双电子积分的指针

    public:
    CI();                                                //默认构造函数
    CI(double h_nuc_ex, double **h_ex, double ****g_ex, int nOrb_ex);
    ~CI();                                               //默认析构函数
    bool get_Int(double **h,double ****g);
    double get(int i, int j, int k, int l);


    //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
    double H_ij(Slater_det &k1, Slater_det &k2);

    //计算分子轨道i与分子轨道j关于hamilton算符的耦合项单电子部分
    //sigma表示电子自旋，0为alpha，1为beta
    double F_ij(Slater_det &k1, Slater_det &k2, const int sigma, int *Num, int ***index);

    //计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
    //sigma表示电子自旋，0为alpha，1为beta
    double G_ij(Slater_det &k1, Slater_det &k2, int *Num, int ***index);
    double G1_ij(Slater_det &k1, Slater_det &k2, const int sigma, int *Num, int ***index);//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
    double G2_ij(Slater_det &k1, Slater_det &k2, int *Num, int ***index);//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
};

int sgn(const int n);
bool CI_new(const int nelec, const int nOrb, const int nelec_ex, const int nOrb_ex, int *Orbital_ex, std::vector<Slater_det> CI_Array);


#endif