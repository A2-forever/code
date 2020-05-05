#ifndef SKIPLIST_H_
#define SKIPLIST_H_

#include <iostream>
#include<ctime>
#include<stdlib.h>
using namespace std;

typedef struct Data
{
    long t;
    int c;
}Flight;

typedef struct Index
{
    struct Data *point=NULL;
    Index *right=NULL;
    Index *down=NULL;
}Index;


bool Random()
{
    bool r;
    r=rand()%2;
    return r;
};

struct Data *Create_Flight(long T,int C)                                           //create a new node and define it with the input data
{
    struct Data *p = new struct Data;
    p->t=T;
    p->c=C;
    return p;
};

Index *Create_Index(Flight *f)                                       //create an index pointing to flight
{
    struct Index *p = new struct Index;
    p->point=f;
    return p;
};

void Output(Index *index)                                                             //output the information of the node or the flight
{
    if(index==NULL)
        cout<<"0"<<endl;
    else
        cout<<(index->point)->c<<endl;
};

void Request(Index **Beginning,int &n_level,long T,int C);                             //insert the node in the skiplist
Index *Depart(Index **Beginning,int &n_level);                                         //delete the earlist node and report the number of passengers 
Index *Query(Index **Beginning,int &n_level,long T);                                                     //search the needed time in the skiplist

bool Check(Index **Beginning,int n_level,long T);                                      //check whether the conflict flight exists
void Insert(Index *Index_Inserted,Index *Index_Near);                                  //insert the node p near temp
Index *Find(Index *index,long T);                                                      //search the needed time in one level

#endif