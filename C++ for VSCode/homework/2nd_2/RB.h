#ifndef RB_H_
#define RB_H_
#include <iostream>
#include<ctime>
#include<stdlib.h>
using namespace std;

typedef struct Data
{
    long t;
    int c;
    Index *left=NULL;
    Index *down=NULL;
    Index *parent=NULL;
}Flight;

bool Random()
{
    bool r;
    r=rand()%2;
    return r;
};

Flight *Create_RBNODE(long T,int C)                                       //create an index pointing to flight
{
    struct Flight *p = new struct Flight;
    t=T;
    c=C;
    return p;
};

void Output(Flight *f)                                                             //output the information of the node or the flight
{
    if(f==NULL)
        cout<<"0"<<endl;
    else
        cout<<f->c<<endl;
};

void
#endif