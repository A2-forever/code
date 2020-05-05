#ifndef FILE_CHAR_H_
#define FILE_CHAR_H_

#include <fstream>
#include <iostream>
#include <string.h>
#include <cmath>


bool read_int(const char *file_name, double &nelec, int &n_Orb, double **h, double ****g);//读取文件中的积分，电子数与分子轨道数
double string2int(char *s);//将字符串读成数字，形式x.xxxxxExxx
bool read_line_int(char *line, double **h, double ****g);//读取轨道积分

#endif