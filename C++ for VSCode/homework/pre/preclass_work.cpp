#include<iostream>
using namespace std;

struct Love{
    long i;
    Love* rt;
};

int main()
{

    unsigned long n;//行数
    long i1,i2;
    cin>>n;
    Love* num=0;

    for(unsigned long i=0;i<n;i++)
    {
        cin>>i1>>i2;//输入
        num->i=i1+i2;
    }

    for(unsigned long i=0;i<n;i++)
        cout<<Nin[i]<<endl;//输出

    return 0;
}