#ifndef SKIPLIST2_H_
#define SKIPLIST2_H_

#include <iostream>
#include<ctime>
#include<stdlib.h>
using namespace std;


typedef struct Node
{
    long t;
    int c;
    struct Node *right=NULL;
    struct Node *down=NULL;
}Index;


Index *Create_Index(long T=-100,int C=-100)                                       //create an index pointing to flight
{
    struct Node *p = new struct Node;
    p->t=T;
    p->c=C;
    p->right=NULL;
    p->down=NULL;
    return p;
};

void Output(Index *index)                                                             //output the information of the node or the flight
{
    if(index==NULL)
        cout<<"0"<<endl;
    else
        cout<<index->c<<endl;
};

bool Request(Index *Beginning[],int &n_level,long T,int C);                             //insert the node in the skiplist
int Depart(Index *Beginning[],int &n_level);                                            //delete the earlist node and report the number of passengers 
int Query(Index *Beginning[], int &n_level,long T);                                    //search the needed time in the skiplist

bool Check(Index **Beginning,int n_level,long T,Index **point);                                      //check whether the conflict flight exists
void Insert(Index *Index_Inserted,Index *Index_Near);                                  //insert the node p near temp
Index *Find(Index *index,long T);                                                      //search the needed time in one level




#endif