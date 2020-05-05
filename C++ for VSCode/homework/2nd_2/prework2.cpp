#include <iostream>
#include<ctime>
#include<stdlib.h>

using namespace std;

typedef struct Data
{
    long t;
    int c;
}Flight;


bool Random()
{
    bool r;
    r=rand()%2;
    return r;
};

int main()
{
    int n;                                                                        //the number of lines we will have

    int op, c;                                                                    //operation,passengers
    long t;                                                                       //departure time
    int n_level=0;
    bool flag=0;                                                                  //title wether our need is satisfied

    Index *Beginning[16];                                                        //define levels of index,15 levels
    Index *result=NULL;                                                          //the result that we will output

    cin >> n;                                                                     //the number of lines we will have

    srand(time(NULL));
    for (int i = 0; i < n; i++)
    {
        cin >> op;
        if (op == 0)                                                              //insert a new order leaving at time t with c passsengers
            cin >> t >> c;
        else if(op==2)                                                            //query the number of passengers leaving at time t
            cin>>t;

        switch (op)
        {
        case 0:
            flag=Request(Beginning, n_level, t,c);
            if(!flag)
                cout<<"0"<<endl;
            break;
        case 1:
            result=Depart(Beginning,n_level);
            Output(result);
            break;
        case 2:
            result=Query(Beginning[n_level],t);
            Output(result);
            break;
        }
    }

    return 0;
}