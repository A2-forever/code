#include <iostream>
#include<ctime>
#include<stdlib.h>
#include "skiplist.h"

using namespace std;

void Request(Index **Beginning, int &n_level,long T,int C)
{
    bool flag=0;
    flag=Check(Beginning,n_level,T);                                      //check whether the conflict flight exists
    if(flag)
        return;                                                      //means that this request can't be accepted
    

    Flight *f=Create_Flight(T,C);
    Index *down_index=NULL;                                                //reserve this index for the connection between near level
    Index *temp_index=NULL;                                                //create a new index pointing to p
    Index *near_index=NULL;                                                //the node will be near the inserted node
    int i=-1;
    
    flag=1;
    while(flag)
    {
        i++;                                                                      //the level positioned presently

        temp_index=Create_Index(f);                                               //create a new index pointing to p
        temp_index->down=down_index;
        
        if(Beginning[i]==NULL)
            Beginning[i]=temp_index;
        else
        {
            near_index=Find(Beginning[i],T);                                                //the node will be near the inserted node
            if(near_index==NULL)
            {
                temp_index->right=Beginning[i]->right;
                Beginning[i]=temp_index;
            }
            else
                Insert(temp_index,near_index);                                                  //insert new index into this skiplist in the ith level
        }

        down_index=temp_index;                                                    //reserve this index for the connection between near level
        flag=Random();                                                            //determine whether to come to next floor
        if(i==15)                                                                 //in the 15th level, pause this loop
            break;                                              
    }
    if(i>n_level)
        n_level=i;

    return;                                                                  //means that this request has been accepted
}

Index *Depart(Index **Beginning,int &n_level)
{
    Index *temp=NULL;
    int n=n_level;
    for(int i=0;i<=n;i++)
    {
        temp=Beginning[i];
        Beginning[i]=Beginning[i]->right;

        if(Beginning[i]==NULL)
            n_level=n_level-1;

        if(Beginning[i+1]==NULL||i==n)
            break;
    }
    return temp;
}

Index *Query(Index **Beginning, int &n_level,long T)
{
    Index *p=Beginning[n_level];
    int i=n_level;
    while(p!=NULL&&i>=0)
    {
        if(Beginning[i]->point->t>T)
        {
            i--;
            p=Beginning[i];
            continue;
        }

        p=Find(p,T);

        if((p->point)->t==T)                                                            //we have found the flight
            return p;
        p=p->down;
    }

    p=NULL;
    return p;                                                                           //we doesn't find the flight
}

bool Check(Index **Beginning,int n_level,long T)
{
    Index *temp=Beginning[n_level];                                                                  //find the time which is near T
    int i=n_level;

    while(temp!=NULL)
    {
        if(Beginning[i]->point->t>T+30)
        {
            if(i==0)
                return false;
            i--;
            temp=Beginning[i];
            continue;
        }
        else if(Beginning[i]->point->t>=T-30)
            return true;


        temp=Find(temp,T);

        if((temp->point->t)>=T-30)
            return true;
        else
        {
            if(temp->right==NULL||(temp->right->point->t)>T+30)
                temp=temp->down;
            else
                return true;
        }
        
    }
    
    return false;
}

void Insert(Index *Index_Inserted,Index *Index_Near)
{
    Index *temp_index;

    temp_index=Index_Near->right;                                                      //Near,inserted,temp
    Index_Near->right=Index_Inserted;                                               //the nodes near inserted node point to it
    Index_Inserted->right=temp_index;
    
}

Index *Find(Index *index,long T)
{
    Index *p=index;

    if(p->point->t>T)
        return NULL;
    
    while(p->right!=NULL)
    {
        if(p->right->point->t<=T)
            p=p->right;
        else
            break;
    }

    return p;
}