#include "FCI.h"

Slater_det::Slater_det()                                //默认构造函数
{
}

Slater_det::Slater_det(const int nelec_ex, const int nOrb_ex, int **Orbital_ex)                                //将轨道信息写入k轨道中
{
    nelec=nelec_ex;
    nOrb=nOrb_ex;
    Orbital=int2_new(nOrb,2);
    nelec_occ=int2_new(nOrb,2);
    for(int i=0;i<nOrb;i++)
        for(int sigma=0;sigma<2;sigma++)
            Orbital[i][sigma]=Orbital_ex[i][sigma];

    this->occ_cal();
}

Slater_det::~Slater_det()                                //默认析构函数
{
}

int Slater_det::Orb(const int I, const int sigma)
{
    return Orbital[I][sigma];
}

bool Slater_det::occ_cal()                               //计算轨道i前的轨道的电子占据数
{
    bool flag=0;
    for(int sigma=0;sigma<2;sigma++){
        nelec_occ[0][sigma] = 0;
        for(int I=1;I<nOrb;I++)
            nelec_occ[I][sigma]=nelec_occ[I-1][sigma]+Orbital[I-1][sigma];//利用动态规划完成计算，sigma表示电子自旋
    }
    flag=1;
    return flag;
}

int Slater_det::gamma(const int I,const int n)
{
    int i=0;
    i=sgn(nelec_occ[I][n]);//轨道的I_n的gamma值

    return i;
}

//判断分子轨道i与分子轨道j之间占据情况不同的轨道与数量
//sigma表示电子自旋，0为alpha，1为beta
//index用于存储电子不同的位置，Num用于存储不同的占据情况不同的轨道数量
//eg:index[0][I][1]表示k1的I_beta位置占据电子，而k2的该位置未占据
//eg:Num[0]表示占据情况不同的alpha轨道的数量
//返回0表示出现问题
bool find(Slater_det &k1, Slater_det &k2, int *Num, int ***index)
{
    int count1=0;
    int count2=0;
    bool flag=1;

    for(int sigma=0;sigma<2;sigma++){
        count1=0;
        count2=0;
        for(int i=0;i<k1.nOrb;i++){
            if(k1.Orbital[i][sigma]!=k2.Orbital[i][sigma]){         //出现占据情况不同
                if(k1.Orbital[i][sigma]==1){                //k1该位置占据电子，k2未占据
                    index[0][count1][sigma]=i;
                    count1++;
                }
                else{                               //k2该位置占据电子，k1未占据
                    index[1][count2][sigma]=i;
                    count2++;
                }
                
            }
        }
        if(count1!=count2)
            flag=0;
        Num[sigma]=count1;
    }

    return flag;
}

//输出Slater_det类的轨道占据情况
std::ostream &operator<<(std::ostream &os, const Slater_det &k)
{
    string E_sigma[2];
    E_sigma[0]="alpha: ";
    E_sigma[1]="beta:  ";
    os << endl;
    os << "nelec: " << k.nelec << endl;
    os << "nOrb:  " << k.nOrb << endl;
    for (int sigma = 0; sigma < 2; sigma++){
        os << E_sigma[sigma];
        for (int i = 0; i < k.nOrb; i++){
            os << k.Orbital[i][sigma] << " ";
        }
        os << endl;
    }

    return os;
}



CI::CI()                                                 //默认构造函数
{
}

CI::CI(double h_nuc_ex, double **h_ex, double ****g_ex, int nOrb_ex)
{
    h=h_ex;
    g=g_ex;
    h_nuc=h_nuc_ex;
    nOrb=nOrb_ex;
}

CI::~CI()                                                   //默认析构函数
{
}

double CI::get(int i = -1, int j = -1, int k = -1, int l = -1)
{
    if(k!=-1)
        return g[i][j][k][l];
    else if (i != -1){
        return h[i][j];
    }
    else{
        return h_nuc;
    }
}

//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
double CI::H_ij(Slater_det &k1, Slater_det &k2)
{
    int *Num = new int[2];
    int ***index = int3_new(2, nOrb, 2);
    bool flag = 1;
    flag = find(k1, k2, Num, index);
    if(!flag){
        delete[] Num;
        int3_delete(index, 2, nOrb);
        return 0;
    }

    double f_ij_alpha=this->F_ij(k1,k2,0,Num,index);          //单电子耦合项的alpha部分
    double f_ij_beta=this->F_ij(k1,k2,1,Num,index);           //单电子耦合项的beta部分

    double g_ij=this->G_ij(k1,k2,Num,index);                  //双电子耦合项

    delete[] Num;
    int3_delete(index, 2, nOrb);  
    return f_ij_alpha+f_ij_beta+g_ij;

}


//计算分子轨道i与分子轨道j关于hamilton算符的耦合项单电子部分
//sigma表示电子自旋，0为alpha，1为beta
double CI::F_ij(Slater_det &k1, Slater_det &k2, const int sigma, int *Num, int ***index)
{
    double sum_int=0;

    if(Num[sigma]==0){                                  //k1与k2的n自旋轨道相同
        for(int P=0;P<nOrb;P++)
            sum_int+=h[P][P]*k1.Orb(P,sigma);
    }
    else if(Num[sigma]==1){                                  //k1与k2的n自旋轨道相差一个电子
        int I=index[1][0][sigma];                //k2的I_sigma分子轨道的n自旋轨道，位置占据电子，k1未占据
        int J=index[0][0][sigma];                //k1的J_sigma分子轨道的n自旋轨道，位置占据电子，k2未占据
        sum_int=k2.gamma(I,sigma)*k1.gamma(I,sigma)*h[I][J];
    }
    else{                                               //k1与k2的n自旋轨道相差两个及以上的电子
        sum_int=0;
    }

    return sum_int;
    
}

//计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
//sigma表示电子自旋，0为alpha，1为beta
double CI::G_ij(Slater_det &k1, Slater_det &k2, int *Num, int ***index)
{
    double g1_ij_alpha=this->G1_ij(k1,k2,0,Num,index);
    double g1_ij_beta=this->G1_ij(k1,k2,1,Num,index);

    double g2_ij=this->G2_ij(k1,k2,Num,index);

    return g1_ij_alpha+g1_ij_beta+g2_ij;

}
//计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
double CI::G1_ij(Slater_det &k1, Slater_det &k2, const int sigma, int *Num, int ***index)
{
    double sum_int=0;

    if(Num[sigma]==0){                                  //k1与k2的sigma自旋轨道相同
        for(int P=0;P<nOrb;P++)
            for(int R=0;R<nOrb;R++)
                sum_int+=k1.Orb(P,sigma)*k1.Orb(R,sigma)*(g[P][P][R][R]-g[P][R][R][P]);
        
        sum_int=0.5*sum_int;
    }
    else if(Num[sigma]==1){                             //k1与k2的n自旋轨道相差一个电子
        int I=index[1][0][sigma];                    //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J=index[0][0][sigma];                    //k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for(int R=0;R<nOrb;R++)
            sum_int+=k1.Orb(R,sigma)*(g[I][J][R][R]-g[I][R][R][J]);

        sum_int=k2.gamma(I,sigma)*k1.gamma(J,sigma)*sum_int;
    }
    else if(Num[sigma]==2){                             //k1与k2的n自旋轨道相差两个及以上的电子
        int I=index[1][0][sigma];                    //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J=index[1][1][sigma];                    //k2的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据
        int K=index[0][0][sigma];                    //k1的K分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int L=index[0][1][sigma];                    //k1的L分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        sum_int=k2.gamma(I,sigma)*k2.gamma(J,sigma)*k1.gamma(K,sigma)*k1.gamma(L,sigma)*(g[I][K][J][L]-g[I][L][J][K]);
    }
    else{
        sum_int=0;
    }

    return sum_int;
    
}

double CI::G2_ij(Slater_det &k1, Slater_det &k2, int *Num, int ***index)//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
{
    double sum_int=0;

    if(Num[0]==0&&Num[1]==0){                           //k1与k2的自旋轨道相同
        for(int P=0;P<nOrb;P++)
            for(int R=0;R<nOrb;R++)
                sum_int+=k1.Orb(P,0)*k1.Orb(R,1)*g[P][P][R][R];
        
    }
    else if(Num[0]==1&&Num[1]==0){                   //k1与k2的alpha自旋轨道相差一个电子,beta自旋轨道相同
        int I=index[1][0][0];                        //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J=index[0][0][0];                        //k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for(int R=0;R<nOrb;R++)
            sum_int+=k1.Orb(R,1)*g[I][J][R][R];

        sum_int=k2.gamma(I,0)*k1.gamma(J,0)*sum_int;
    }
    else if(Num[0]==0&&Num[1]==1){                   //k1与k2的beta自旋轨道相差一个电子,alpha自旋轨道相同
        int I=index[1][0][1];                        //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J=index[0][0][1];                        //k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for(int P=0;P<nOrb;P++)
            sum_int+=k1.Orb(P,0)*g[P][P][I][J];

        sum_int=k2.gamma(I,1)*k1.gamma(J,1)*sum_int;
    }
    else if(Num[0]==1&&Num[1]==1){                   //k1与k2的alpha与beta自旋轨道均相差一个电子
        int I_alpha=index[1][0][0];                  //k2的I分子轨道的alpha自旋轨道，位置占据电子，k1未占据
        int J_alpha=index[1][1][0];                  //k2的J分子轨道的alpha自旋轨道，位置占据电子，k2未占据
        int I_beta=index[0][0][1];                   //k1的I分子轨道的beta自旋轨道，位置占据电子，k1未占据
        int J_beta=index[0][1][1];                   //k1的J分子轨道的beta自旋轨道，位置占据电子，k2未占据
        
        sum_int=k2.gamma(I_alpha,0)*k2.gamma(J_alpha,0)*k1.gamma(I_beta,1)*k1.gamma(J_beta,1)*g[I_alpha][J_alpha][I_beta][J_beta];
    }
    else{
        sum_int=0;
    }

    return sum_int;
}


int sgn(const int n)//奇数返回-1，偶数返回1
{
    if(n%2==1)
        return -1;
    
    return 1;
}

//ex表示现存的轨道或电子数，未加ex表示实际的总电子数
//Orbital_ex用于表示暂时表示组态的数组
//CI_Array用于存储组态
bool CI_new(const int nelec, const int nOrb, const int nelec_ex, const int nOrb_ex, int **Orbital_ex, std::vector<Slater_det> &CI_Array)
{ 
    //cout << nelec_ex <<" "<< nOrb_ex<< endl;
    if(nelec_ex==0){
        Slater_det new_det(nelec,nOrb,Orbital_ex);
        CI_Array.push_back(new_det);
        return 1;
    }

    //本次使用的i用于计数，即标记占据最高能级的电子e_m的几种排列方式
    //现有nelec_ex个电子，n_Orb_ex个分子轨道，那么电子e_m可以从第nelec_ex号轨道排列到第n_Orb_ex号轨道
    //Oribital_ex[I][sigma]表示第i=I+1+n_Orb*sigma号轨道
    //sigma=(i-1) / n_Orb,I=i-sigma*n_Orb-1
    bool flag=0;
    int I=-1;
    int sigma=-1;
    for(int i=nelec_ex;i<=nOrb_ex;i++){
	    sigma = (i-1) / nOrb;
        I=i-sigma*nOrb-1;
        Orbital_ex[I][sigma]=1;
        flag = CI_new(nelec, nOrb, nelec_ex - 1, i - 1, Orbital_ex, CI_Array);
        Orbital_ex[I][sigma]=0;
    }

    return flag;
}
