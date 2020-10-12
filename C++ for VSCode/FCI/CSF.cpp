#ifndef CSF_CPP
#define CSF_CPP

#include "CSF.h"


#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using std::vector;
using std::string;
using std::endl;


CSF::CSF()                                   //默认构造函数
{
}

CSF(std::vector<int> &Orbital_ex, const double &MS_ex)                                   //默认构造函数
{
    nOrb = Orbital_ex.size();
    MS = MS_ex;

    for (int i = 0; i < nOrb; i++)
        CSF[i] = Orbital_ex[i];

    for (int i = 0; i < nOrb; i++)
    {
        if(CSF[i] == 1)
        {
            S += 0.5;
            nelec++;
        }
        else if(CSF[i] == 2)
        {
            S -= 0.5;
            nelec++;
        }
        else if(CSF[i] == 3)
            nelec += 2;
    }

    this->CSF2Slater(&Orbital_ex);
}

~CSF()                                  //默认析构函数
{
}

bool CSF2Slater(std::vector<int> &Orbital_ex)
{
    int num = 1;
    vector<int> single;
    vector<int> tmp_Slater(2 * nOrb);
    for (int i = 0; i < nOrb; i++)
    {
        if (CSF[i] == 1 || CSF[i] == 2)
        {
            num *= 2;
            single.push_back(i);
        }
        else if(CSF[i] == 3)
        {
            tmp_Slater[0 * nOrb + i] = 1;
            tmp_Slater[1 * nOrb + i] = 1;
        }
    }

    this->vector2Slater(Orbital_ex, tmp_Slater, single);
}

CI::CI()                                                 //默认构造函数
{
}

CI::CI(const double &h_nuc_ex, const vector<double> &h_ex, const vector<double> &g_ex, const int nOrb_ex)
{
    h_nuc = h_nuc_ex;
    h = h_ex;
    g = g_ex;
    nOrb = nOrb_ex;

    dim =nOrb_ex;
    dim2 = pow(nOrb_ex, 2);
    dim3 = pow(nOrb_ex, 3);

}

CI::~CI()                                                   //默认析构函数
{
}

double CI::get(const int &i = -1, const int &j = -1, const int &k = -1, const int &l = -1)
{
    if(k!=-1)
        return g[i * dim3 + j * dim2 + k * dim + l];
    else if (i != -1){
        return h[i * dim + j];
    }
    else{
        return h_nuc;
    }
}

//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
double CI::H_ij(Slater_det &k1, Slater_det &k2)
{
    vector<int> Num(2);//位置不同的电子的数量
    vector<int> index(2 * 2 * nOrb);//2个组态，2种自旋，最多nOrb个轨道
    bool flag = 1;
    flag = find(k1, k2, Num, index);
    if(!flag)
        return 0;

    if (Num[0] + Num[1] > 2)
        return 0;
    /*
	for (int i = 0; i < 2;i++)
		cout << Num[i] << "\t";
    cout << endl;
    for (int i = 0; i < 2;i++)
    {
        for (int j = 0; j < 2;j++)
        {
            for (int k = 0; k < Num[j];k++)
                cout << index[i * 2 * nOrb + k * 2 + j] << "\t";
            cout << endl;
        }
        cout << endl;
    }
    */

    double f_ij_alpha = this->F_ij(k1, k2, 0, Num, index);      //单电子耦合项的alpha部分
    //cout << f_ij_alpha << endl;
    double f_ij_beta = this->F_ij(k1, k2, 1, Num, index);       //单电子耦合项的beta部分
    //cout << f_ij_beta << endl;

    double g_ij = this->G_ij(k1, k2, Num, index);               //双电子耦合项
    //cout << g_ij << endl;

    //cout << f_ij_alpha + f_ij_beta + g_ij << endl;
    return f_ij_alpha+f_ij_beta+g_ij;

}


//计算分子轨道i与分子轨道j关于hamilton算符的耦合项单电子部分
//sigma表示电子自旋，0为alpha，1为beta
double CI::F_ij(Slater_det &k1, Slater_det &k2, const int sigma, const vector<int> &Num, const vector<int> &index)
{
    if (Num[1 - sigma] != 0) //在进行alpha部分的耦合时，beta部分必须完全一样，反之亦然
        return 0;

    double sum_int=0;
    if (Num[sigma] >= 2)                                        //k1与k2的n自旋轨道相差两个及以上的电子
    {
        sum_int = 0;
    }
    else if(Num[sigma] == 1)                                    //k1与k2的n自旋轨道相差一个电子
    {
        int I = index[1 * 2 * dim + sigma * dim + 0];           //k2的I_sigma分子轨道的n自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + sigma * dim + 0];           //k1的J_sigma分子轨道的n自旋轨道，位置占据电子，k2未占据
        sum_int = k2.gamma(I, sigma) * k1.gamma(J, sigma) * h[I * dim + J];
    }
    else if(Num[sigma] == 0)                                    //k1与k2的n自旋轨道相同 
    {
        for (int P = 0; P < nOrb; P++)
            sum_int += h[P * dim + P] * k1.Orb(P, sigma);
    }

    return sum_int;
}

//计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
//sigma表示电子自旋，0为alpha，1为beta
double CI::G_ij(Slater_det &k1, Slater_det &k2, const vector<int> &Num, const vector<int> &index)
{
    double g1_ij_alpha = this->G1_ij(k1, k2, 0, Num, index);
    double g1_ij_beta = this->G1_ij(k1, k2, 1, Num, index);

    double g2_ij = this->G2_ij(k1, k2, Num, index);

    return g1_ij_alpha + g1_ij_beta + g2_ij;

}
//计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
double CI::G1_ij(Slater_det &k1, Slater_det &k2, const int sigma, const vector<int> &Num, const vector<int> &index)
{
    if (Num[1 - sigma] != 0) //在进行alpha部分的耦合时，beta部分必须完全一样，反之亦然
        return 0;

    double sum_int = 0;
    if(Num[sigma] >= 3)                                         //k1与k2的n自旋轨道相差三个及以上的电子
    {
        sum_int=0;
    }
    else if(Num[sigma] == 2)                                    //k1与k2的n自旋轨道相差两个电子
    {
        int I = index[1 * 2 * dim + sigma * dim + 0];           //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[1 * 2 * dim + sigma * dim + 1];           //k2的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据
        int K = index[0 * 2 * dim + sigma * dim + 0];           //k1的K分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int L = index[0 * 2 * dim + sigma * dim + 1];           //k1的L分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        sum_int = k2.gamma(I, sigma) * k2.gamma(J, sigma) * k1.gamma(K, sigma) * k1.gamma(L, sigma) * (g[I * dim3 + K * dim2 + J * dim + L] - g[I * dim3 + L * dim2 + J * dim + K]);
    }
    else if(Num[sigma] == 1)                                    //k1与k2的n自旋轨道相差一个电子
    {
        int I = index[1 * 2 * dim + sigma * dim + 0];           //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + sigma * dim + 0];           //k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for (int R = 0; R < nOrb; R++)
            sum_int += k1.Orb(R, sigma) * (g[I * dim3 + J * dim2 + R * dim + R] - g[I * dim3 + R * dim2 + R * dim + J]);
        
        sum_int = k2.gamma(I, sigma) * k1.gamma(J, sigma) * sum_int;
    }
    else if(Num[sigma] == 0)                               //k1与k2的sigma自旋轨道相同
    {
        for (int P = 0; P < nOrb; P++)
            for (int R = 0; R < nOrb; R++)
                sum_int += k1.Orb(P, sigma) * k1.Orb(R, sigma) * (g[P * dim3 + P * dim2 + R * dim + R] - g[P * dim3 + R * dim2 + R * dim + P]);

        sum_int *= 0.5;
    }
    

    return sum_int;
    
}

double CI::G2_ij(Slater_det &k1, Slater_det &k2, const vector<int> &Num, const vector<int> &index)//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
{
    double sum_int=0;

    if(Num[0] == 0 && Num[1] == 0)                            //k1与k2的自旋轨道相同
    {
        for (int P = 0; P < nOrb; P++)
            for (int R = 0; R < nOrb; R++)
                sum_int += k1.Orb(P, 0) * k1.Orb(R, 1) * g[P * dim3 + P * dim2 + R * dim + R];
    }
    else if(Num[0] == 1 && Num[1] == 0)                         //k1与k2的alpha自旋轨道相差一个电子,beta自旋轨道相同
    {
        int I = index[1 * 2 * dim + 0 * dim + 0];               //1*nOrb*nOrb+0*nOrb+0,k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + 0 * dim + 0];               //0*nOrb*nOrb+0*nOrb+0,k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for (int R = 0; R < nOrb; R++)
            sum_int += k1.Orb(R, 1) * g[I * dim3 + J * dim2 + R * dim + R];

        sum_int = k2.gamma(I, 0) * k1.gamma(J, 0) * sum_int;
    }
    else if(Num[0] == 0 && Num[1] == 1)                         //k1与k2的beta自旋轨道相差一个电子,alpha自旋轨道相同
    {
        int I = index[1 * 2 * dim + 1 * dim + 0];               //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + 1 * dim + 0];               //k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for (int P = 0; P < nOrb; P++)
            sum_int += k1.Orb(P, 0) * g[P * dim3 + P * dim2 + I * dim + J];

        sum_int = k2.gamma(I, 1) * k1.gamma(J, 1) * sum_int;
    }
    else if(Num[0] == 1 && Num[1] == 1)                      //k1与k2的alpha与beta自旋轨道均相差一个电子
    {
        int I_alpha = index[1 * 2 * dim + 0 * dim + 0];      //k2的I分子轨道的alpha自旋轨道，位置占据电子，k1未占据
        int J_alpha = index[0 * 2 * dim + 0 * dim + 0];      //k1的J分子轨道的alpha自旋轨道，位置占据电子，k2未占据
        int I_beta = index[1 * 2 * dim + 1 * dim + 0];       //k2的I分子轨道的beta自旋轨道，位置占据电子，k1未占据
        int J_beta = index[0 * 2 * dim + 1 * dim + 0];       //k1的J分子轨道的beta自旋轨道，位置占据电子，k2未占据

        sum_int = k2.gamma(I_alpha, 0) * k1.gamma(J_alpha, 0) * k2.gamma(I_beta, 1) * k1.gamma(J_beta, 1) * g[I_alpha * dim3 + J_alpha * dim2 + I_beta * dim + J_beta];
    }
    else
    {
        sum_int=0;
    }

    return sum_int;
}


#endif