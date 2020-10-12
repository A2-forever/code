#ifndef EIGEN_H_
#define EIGEN_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

using std::vector;


bool eigenh(vector<double> Matrix, vector<double> &dbVectors, vector<double> &dbEigenvalues, const int nJt);//对角化实对称矩阵,返回
void jacobi(vector<double> &Matrix, vector<double> &dbVectors, const vector<int> &index);//进行一次jacobi分解的循环过程

bool rotate(vector<double> &dbVectors, const vector<int> &index, const double phi);//旋转矩阵V，大小为n x n, 旋转角为phi,位置为pq
vector<int> max(const vector<double> &Matrix); //返回非对角线绝对值最大值的坐标p,q，要求p<q

std::ostream &operator<<(std::ostream &os, const vector<double> &Matrix);
vector<double> multiply(const vector<double> &Matrix1, const vector<double> &Matrix2);          //对于两个方阵相乘

void Maqr(vector<double> &Matrix, vector<double> &Q, int m, int n);                          //进行一般实矩阵QR分解的函数
#endif