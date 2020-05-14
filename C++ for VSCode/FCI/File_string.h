#ifndef FILE_STRING_H_
#define FILE_STRING_H_

#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <cmath>
#include "Array.cpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;

bool read_int(const string &file_name, int &nelec, int &n_Orb, double &h_nuc, double **h, double ****g); //读取文件中的积分，电子数与分子轨道数
bool read_line_int(const string &line, double &h_nuc, double **h, double ****g);//读取轨道积分
double string2int(const string &s);//将字符串读成数字，形式x.xxxxxExxx
bool judge_start(const string &line, const string &delim);//判断是否可以开始读取积分
vector<string> split(const string &s, const string &delim);

#endif