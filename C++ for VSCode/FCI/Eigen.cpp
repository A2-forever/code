#ifndef EIGEN_CPP_
#define EIGEN_CPP_

#include "Eigen.h"

bool eigenh(vector<double> Matrix, vector<double> &dbVectors, vector<double> &dbEigenvalues, const int nJt)
{
    int count=0;
    bool flag = 1;
    int ndim = sqrt(Matrix.size());
    vector<int> index(2);
    int p=-1;
    int q=-1;
    
    for(int i=0;i<ndim;i++)//将dbVectors变为单位矩阵，dbVectors的输入必须为一个全为0的矩阵或对角矩阵
        dbVectors[i * ndim + i] = 1;
/*
    std::ofstream logfile;
	logfile.open("test\\test.log");
	logfile.setf(std::ios::showpoint); //设置为始终输出小数点后的数字，就是说 a = 3，它也输出 3.00000 这样
	logfile.precision(6);
	logfile.setf(std::ios::fixed); //设置为小数位始终有 6 位，没有这个的话就会像上面那个代码那样固定的不是小数点后面的数字了。
*/
    while(1)
    {
        //logfile<<"Matrix is: "<<std::endl;
        //output(logfile,Matrix,ndim);
        //logfile<<"dbVectors is: "<<std::endl;
        //output(logfile,dbVectors,ndim);

        index=max(Matrix);//寻找到非对角线的绝对值最大的数的位置
        p = index[0];//p<q
        q = index[1];

        if(fabs(Matrix[p*ndim+q])<eps){
            flag=1;
            break;
        }//非对角线元素尽量逼近0
            
        if(count>nJt){
            flag=0;
            break;
        }//超过迭代次数

        count++;//进行一次迭代计数
        jacobi(Matrix, dbVectors, index); //进行一次jacobi迭代
    }

    for(int i=0;i<ndim;i++)
        dbEigenvalues[i*ndim+i]=Matrix[i*ndim+i];


    return flag;
}

void jacobi(vector<double> &Matrix, vector<double> &dbVectors, const vector<int> &index)//进行一次jacobi分解的循环过程
{
    int p = index[0];
    int q = index[1];
    int ndim = sqrt(Matrix.size());
    double tan2phi=-2*Matrix[p*ndim+q]/(Matrix[q*ndim+q]-Matrix[p*ndim+p]);
    double phi=0.5*atan(tan2phi);

    //对下标为pp,qq,pq,qp的数组成员进行旋转操作
    double b_pp=Matrix[p*ndim+p]*cos(phi)*cos(phi)+Matrix[q*ndim+q]*sin(phi)*sin(phi)+2*Matrix[p*ndim+q]*cos(phi)*sin(phi);
    double b_qq=Matrix[p*ndim+p]*sin(phi)*sin(phi)+Matrix[q*ndim+q]*cos(phi)*cos(phi)-2*Matrix[p*ndim+q]*cos(phi)*sin(phi);
    double b_pq=0.5*(Matrix[q*ndim+q]-Matrix[p*ndim+p])*sin(2*phi)+Matrix[p*ndim+q]*cos(2*phi);
    double b_p;
    double b_q;

    Matrix[p*ndim+p]=b_pp;
    Matrix[q*ndim+q]=b_qq;
    Matrix[p*ndim+q]=Matrix[q*ndim+p]=b_pq;

    //对下标为pi,qi,ip,iq的数组成员进行旋转操作=
    for(int i=0;i<ndim;i++){
        if(i==p||i==q)
            continue;
        else{
            b_p=Matrix[p*ndim+i]*cos(phi)+Matrix[q*ndim+i]*sin(phi);
            b_q=-Matrix[p*ndim+i]*sin(phi)+Matrix[q*ndim+i]*cos(phi);
            Matrix[p*ndim+i]=b_p;
            Matrix[q*ndim+i]=b_q;

            b_p=Matrix[i*ndim+p]*cos(phi)+Matrix[i*ndim+q]*sin(phi);
            b_q=-Matrix[i*ndim+p]*sin(phi)+Matrix[i*ndim+q]*cos(phi);
            Matrix[i*ndim+p]=b_p;
            Matrix[i*ndim+q]=b_q;
        }
    }

    rotate(dbVectors,index,phi);//计算特征向量，对特征向量进行旋转操作
}


bool rotate(vector<double> &dbVectors, const vector<int> &index, const double phi)//右旋转矩阵V（指VU），大小为n x n, 旋转角为phi,位置为pq 
{
    int p = index[0];
    int q = index[1];
    int ndim = sqrt(dbVectors.size());
    double b_p;
    double b_q;

    for(int i=0;i<ndim;i++){
        b_p=dbVectors[i*ndim+p]*cos(phi)+dbVectors[i*ndim+q]*sin(phi);
        b_q=-dbVectors[i*ndim+p]*sin(phi)+dbVectors[i*ndim+q]*cos(phi);
        dbVectors[i*ndim+p]=b_p;
        dbVectors[i*ndim+q]=b_q;
    }
    
    return 1;
}

vector<int> max(const vector<double> &Matrix)//返回非对角线绝对值最大值的坐标p,q，要求p<q
{
    vector<int> index(2);
    int ndim = sqrt(Matrix.size());
    double temp = 0;

    for (int i = 0; i < ndim; i++){
        for (int j = i + 1; j < ndim; j++){
            if(temp < fabs(Matrix[i * ndim + j])){
                temp = fabs(Matrix[i * ndim + j]);
                index[0] = i;
                index[1] = j;
            }
        }
    }

    return index;
}


std::ostream &operator<<(std::ostream &os, const vector<double> &Matrix)
{
    int ndim = sqrt(Matrix.size());

    os << std::endl;
    for (int i = 0; i < ndim;i++){
        for (int j = 0; j < ndim;j++)
        {
            os << Matrix[i * ndim + j];
            if(Matrix[i * ndim + j]<0)
                os << "       ";
            else
                os << "        ";
        }
        os << std::endl;
    }
    return os;
}


vector<double> multiply(const vector<double> &Matrix1, const vector<double> &Matrix2)//对于两个方阵相乘
{
    int ndim = sqrt(Matrix1.size());
    vector<double> ans(Matrix1.size());
    double temp = 0;
    for (int i = 0; i < ndim; i++)
        for (int j = 0; j < ndim; j++)
        {
            temp = 0;
            for (int k = 0; k < ndim; k++)
                temp += Matrix1[i*ndim+k] * Matrix2[k*ndim+j];
            ans[i*ndim+j]=temp;
        }
    
    return ans;
}

#endif