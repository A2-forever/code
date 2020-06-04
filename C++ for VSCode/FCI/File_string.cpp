#ifndef FILE_STRING_CPP_
#define FILE_STRING_CPP_

#include "File_string.h"


bool read_int(const string &file_name, int &nelec, int &nOrb, double &MS, double &h_nuc, vector<double> &h, vector<double> &g)//读取文件中的积分，电子数与分子轨道数
{
    bool flag=0;
    std::ifstream file;
    file.open(file_name.c_str());

    if(!file.is_open()){
        cout << "未成功打开积分文件" << endl;
        return 0;
    }

    string line;
    string delim;
    
    //第一行格式如右，用于读取轨道个数与电子个数:"&FCI NORB=114,NELEC=42,MS2=0,"
    getline(file, line);
    delim = "=";
    vector<string> v = split(line, delim);

    delim = ",";
    vector<string> v1 = split(v[1], delim);
    nOrb = atof(v1[0].c_str());
    v1 = split(v[2], delim);
    nelec = atof(v1[0].c_str());
    v1 = split(v[3], delim);
    MS = atof(v1[0].c_str()) ;
    /*
    cout << nelec << endl
         << nOrb << endl;
    */

    //右为读取开始的标志:" &END"
    delim = "&END";

    h.resize(nOrb * nOrb);
    g.resize(nOrb * nOrb * nOrb * nOrb);

    while(getline(file,line)){
        if(!flag)//用于判断是否可以开始读取积分,flag=1表示可以开始读取积分
            flag = judge_start(line, delim);
        else{
            flag = read_line_int(line, h_nuc, h, g);
            if(flag==0)//读到空行，读取结束
                break;
            }
    }

    file.close();

    return 1;
}

bool read_line_int(const string &line, double &h_nuc, vector<double> &h, vector<double> &g)//读取轨道积分
{
    string delim=" ";
    vector<string> v = split(line, delim);


    double INT=string2int(v[0]);

    int i = atof(v[1].c_str()) - 1;
    int j = atof(v[2].c_str()) - 1;
    int k = atof(v[3].c_str()) - 1;
    int l = atof(v[4].c_str()) - 1;

    int dim2 = h.size();
    int dim = sqrt(dim2);       //数组的维度
    int dim3 = dim * dim2;

    //cout << i << " " << j << " " << k << " " << l << "\t" << INT << std::endl;
    if(k!=-1){//双电子积分
        g[i * dim3 + j * dim2 + k * dim + l] = INT;              //ijkl
        g[i * dim3 + j * dim2 + l * dim + k] = INT;              //jikl
        g[j * dim3 + i * dim2 + k * dim + l] = INT;              //ijlk
        g[j * dim3 + i * dim2 + l * dim + k] = INT;              //jilk
    
        g[k * dim3 + l * dim2 + i * dim + j] = INT;              //klij
        g[k * dim3 + l * dim2 + j * dim + i] = INT;              //klji
        g[l * dim3 + k * dim2 + i * dim + j] = INT;              //lkij
        g[l * dim3 + k * dim2 + j * dim + i] = INT;              //lkji
        //cout << i << " " << j << " " << k << " " << l << "\t" << g[i-1][j-1][k-1][l-1] << std::endl;
    }
    else if (i != -1){//单电子积分
        h[i * dim + j] = INT;
        h[j * dim + i] = INT;
        //cout << i << " " << j << "\t" << h[i - 1][j - 1] << std::endl;
    }
    else{//核积分项
        h_nuc = INT;
        //cout << h_nuc << endl;
        return 0;//表示读取结束，读到文件末尾
    }

    return 1;//表示读取正常，读取仍未结束
}

double string2int(const string &s)//将字符串读成数字，形式x.xxxxxExxx
{
    string delim="E";
    vector<string> v=split(s,delim);
    int n=v.size();
    double tmp=0;

    if(n==1){
        tmp=atof(v[0].c_str());
    }
    else{
        double t1=atof(v[0].c_str());
        int t2=atof(v[1].c_str());
        tmp=t1*pow(10,t2);
    }
    return tmp;

}

bool judge_start(const string &line, const string &delim)//判断是否可以开始读取积分
{
    bool flag = 0;
    vector<string> v1 = split(line, delim);
    if(v1.size()==0){
        flag = 1;
    }
    else{
        string delim2 = " ";
        vector<string> v2 = split(v1[0], delim2);
        if(v2.size()==0)
            flag = 1;
    }

    return flag;
}

vector<string> split(const string& s, const string& delim)
{
    vector<string> v;
    string::size_type pos1, pos2;
    pos2 = s.find(delim);
    pos1 = 0;
    while (string::npos != pos2){
        if(pos1!=pos2)//防止出现连续的分隔符delim
            v.push_back(s.substr(pos1, pos2 - pos1)); //开始位置与获取字符串的长度，取出来放到vector中

        pos1 = pos2 + delim.size();
        pos2 = s.find(delim, pos1);

    }
        
    if(pos1 != s.length())
        v.push_back(s.substr(pos1));
    
    return v;
}

#endif