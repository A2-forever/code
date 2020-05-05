#ifndef ARRAY_H_
#define ARRAY_H_

#include <iostream>
#include <cstdlib>

double **double2_new(const int row, const int column);//创建二维double数组
double **double2_new(const long long row, const long long column); //创建二维double数组,long long
double ***double3_new(const int r1, const int r2, const int r3);//创建三维double数组
double ****double4_new(const int r1, const int r2, const int r3, const int r4);//创建四维double数组

void double2_delete(double **M, const int row);//释放二维double数组
void double2_delete(double **M, const long long row);//释放二维double数组,long long
void double3_delete(double ***M, const int r1, const int r2);//释放三维double数组
void double4_delete(double ****M, const int r1, const int r2, const int r3);//释放四维double数组


int **int2_new(const int row, const int column);//创建二维int数组
int ***int3_new(const int r1, const int r2, const int r3);//创建三维int数组
int ****int4_new(const int r1, const int r2, const int r3, const int r4);//创建四维int数组

void int2_delete(int **M, const int row);//释放二维int数组
void int3_delete(int ***M, const int r1, const int r2);//释放三维int数组
void int4_delete(int ****M, const int r1, const int r2, const int r3);//释放四维int数组

void output(double **M,const int n);

#endif