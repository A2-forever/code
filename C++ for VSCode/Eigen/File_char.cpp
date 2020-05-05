#include "File_char.h"


bool read_int(const char *file_name, double &nelec, int &n_Orb, double **h, double ****g)//读取文件中的积分，电子数与分子轨道数
{
    bool flag=1;
    char tmp[256];
    
    return flag;
}

double string2int(char *s)//将字符串读成数字，形式x.xxxxxExxx
{
    char delim[]="E";
    char *p=0;
    double res;
    int t;

    p=strtok(s,delim);
    res=atof(p);

    p=strtok(NULL,delim);
    if(!p)
        return res;

    t=atof(p);

    res=res*pow(10,t);

    return res;

}

bool read_line_int(char *line, double **h, double ****g)//读取轨道积分
{
    bool flag=1;
    char delim[]=" ";

    char *t=strtok(line,delim);
	char *p = strtok(NULL, delim);

    int index[4]={0};
    int count=0;//计数四次

    while(p){
        index[count]=atof(p);
        p=strtok(NULL,delim);
        count++;
    }

    double INT=string2int(t);

    if(index[2]==0)
        h[index[0]-1][index[1]-1]=INT;
    else{
        g[index[0]-1][index[1]-1][index[2]-1][index[3]-1]=INT;
    }

    return flag;
}