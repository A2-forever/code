#include "Eigen.h"

bool eigenh(double **Matrix, double **dbVectors, double dbEigenvalues[], const int ndim, const int nJt)
{
    double **M=resemble(Matrix,ndim,ndim);
    
    int count=0;
    int p=-1;
    int q=-1;
    bool flag=1;

    while(1)
    {
        max(M,ndim,p,q);//寻找到非对角线的绝对值最大的数的位置

        if(fabs(M[p][q])<eps){
            flag=1;
            break;
        }//非对角线元素尽量逼近0
            
        if(count>nJt){
            flag=0;
            break;
        }//超过迭代次数

        count++;//进行一次迭代计数
        jacobi(M,dbVectors,ndim,p,q);//进行一次jacobi迭代
        
        /*
        std::cout<<"Matrix is: "<<std::endl;
        output(M,ndim);
        std::cout<<"dbVectors is: "<<std::endl;
        output(dbVectors,ndim);
        */

    }

    for(int i=0;i<ndim;i++)
        dbEigenvalues[i]=M[i][i];

    double2_delete(M,ndim);

    return flag;
}

void jacobi(double **M, double **dbVectors, const int n, const int p, const int q)//进行一次jacobi分解的循环过程
{
    double tan2phi=-2*M[p][q]/(M[q][q]-M[p][p]);
    double phi=0.5*atan(tan2phi);

    //对下标为pp,qq,pq,qp的数组成员进行旋转操作
    double b_pp=M[p][p]*cos(phi)*cos(phi)+M[q][q]*sin(phi)*sin(phi)+2*M[p][q]*cos(phi)*sin(phi);
    double b_qq=M[p][p]*sin(phi)*sin(phi)+M[q][q]*cos(phi)*cos(phi)-2*M[p][q]*cos(phi)*sin(phi);
    double b_pq=0.5*(M[q][q]-M[p][p])*sin(2*phi)+M[p][q]*cos(2*phi);
    double b_p;
    double b_q;

    M[p][p]=b_pp;
    M[q][q]=b_qq;
    M[p][q]=M[q][p]=b_pq;

    //对下标为pi,qi,ip,iq的数组成员进行旋转操作=
    for(int i=0;i<n;i++){
        if(i==p||i==q)
            continue;
        else{
            b_p=M[p][i]*cos(phi)+M[q][i]*sin(phi);
            b_q=-M[p][i]*sin(phi)+M[q][i]*cos(phi);
            M[p][i]=b_p;
            M[q][i]=b_q;

            b_p=M[i][p]*cos(phi)+M[i][q]*sin(phi);
            b_q=-M[i][p]*sin(phi)+M[i][q]*cos(phi);
            M[i][p]=b_p;
            M[i][q]=b_q;
        }
    }

    rotate(dbVectors,n,p,q,phi);//计算特征向量，对特征向量进行旋转操作
}

bool max(double **M, const int n, int &p,int &q)//返回非对角线绝对值最大值的坐标p,q，要求p<q
{
    double temp=0;
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            if(temp<fabs(M[i][j])){
                temp=fabs(M[i][j]);
                p=i;
                q=j;
            }
        }
    }

    return 1;
}

double **multiply(double **M, double **N, const int n)//对于两个矩阵相乘
{
    double **ans=double2_new(n,n);                   //新建一个n x n的二维数组

    double temp = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            temp = 0;
            for (int k = 0; k < n; k++)
                temp += M[i][k] * N[k][j];
            ans[i][j]=temp;
        }
    
    return ans;
}

bool rotate(double **V, const int n, const int p,const int q,const double phi)//右旋转矩阵V（指VU），大小为n x n, 旋转角为phi,位置为pq 
{
    double b_p;
    double b_q;

    for(int i=0;i<n;i++){
        b_p=V[i][p]*cos(phi)+V[i][q]*sin(phi);
        b_q=-V[i][p]*sin(phi)+V[i][q]*cos(phi);
        V[i][p]=b_p;
        V[i][q]=b_q;
    }
    
    return 1;
}

double **ones(const int n)                        //构建单位矩阵M,大小为n x n
{
    double **M=double2_new(n,n);                   //新建一个n x n的二维数组

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++){
            if(i==j)
                M[i][j]=1.0f;
            else
                M[i][j]=0.0f;
        }

    return M;
}


double **resemble(double **M, const int row, const int column)//复制二维数组M,大小为row x column
{
    double **N=double2_new(row,column);                   //新建一个row x column的二维数组

    for(int i=0;i<row;i++)
        for(int j=0;j<column;j++)
            N[i][j]=M[i][j];

    return N;
}
