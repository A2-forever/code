#ifndef EIGEN_H_
#define EIGEN_H_

#include <cmath>
#include "Array.cpp"

#define eps 1e-6


bool eigenh(double **Matrix, double **dbVectors, double dbEigenvalues[], const int ndim, const int nJt);//对角化实对称矩阵,返回
void jacobi(double **M, double **dbVectors, const int n, const int p, const int q);//进行一次jacobi分解的循环过程

double **multiply(double **M, double **N, const int n);//对于两个矩阵相乘
bool rotate(double **V,const int n, const int p,const int q,const double phi);//旋转矩阵V，大小为n x n, 旋转角为phi,位置为pq

bool max(double **M, const int n, int &p,int &q); //返回非对角线绝对值最大值的坐标p,q，要求p<q
double **ones(const int n);//构建单位矩阵M,大小为n x n
double **resemble(double **M, const int row, const int column);//复制二维数组M,大小为row x column

#endif