#include <iostream>
#include<ctime>
#include<stdlib.h>
#include "skiplist2.h"
#include "skiplist2.cpp"
using namespace std;



int main()
{  
    Index *Beginning[21]={NULL};
    int n_level=-1;
    int p=0;
    int k=2000000;
    long i=0;
    clock_t startTime,endTime;
    
    startTime = clock();//计时开始
    srand(time(NULL));

    

    for(i=1;i<=k;i++)
        Request(Beginning,n_level,i*50,i);

    for(i=1;i<=k;i++)
    {
        p=Query(Beginning,n_level,i*50);
        if(i%10000==0)
            cout<<p<<endl;
    }
    cout<<endl;
    
    for(i=1;i<=k;i++)
    {
        p=Depart(Beginning,n_level);
        if(i%100000==0)
            cout<<p<<endl;
        if(i>k-10)
            cout<<p<<endl;
    }

    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    endTime = clock();//计时结束
    cout << "The run time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    system("pause");
    return 0;
}  
