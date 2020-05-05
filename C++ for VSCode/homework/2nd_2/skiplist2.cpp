#include <iostream>
#include<ctime>
#include<stdlib.h>
#include "skiplist2.h"

using namespace std;


bool Request(Index *Beginning[], int &n_level,long T,int C)
{
    Index *down_index=NULL;                                                //reserve this index for the connection between near level
    Index *temp_index=NULL;                                                //create a new index pointing to p
    Index *point[16]={NULL};
    int i=-1;
    bool flag=1;
    
    if(Beginning[0]!=NULL)
    {
        flag=Check(Beginning,n_level,T,point);                                      //check whether the conflict flight exists
        if(!flag)
            return false;                                                      //means that this request can't be accepted
    }

    while(flag)
    {
        i++;                                                                      //the level positioned presently

        temp_index=Create_Index(T,C);                                               //create a new index pointing to p
        temp_index->down=down_index;                                               //connect two nodes vertically
        down_index=temp_index;                                                    //reserve this index for the connection between near level

        if(Beginning[i]==NULL)
        {
            Beginning[i]=Create_Index();                                      //the head node is empty to ensure each node is in its right side
            if(i>0)
                Beginning[i]->down=Beginning[i-1];
            point[i]=Beginning[i];
        }
        
        Insert(temp_index,point[i]);                                                  //insert new index into this skiplist in the ith level

        flag=rand()%2;;                                                            //determine whether to come to next floor
        
        /*if(temp_index->right!=NULL)
            cout<<point[i]->t<<"|"<<temp_index->t<<"|"<<temp_index->right->t<<"       ";
        else
            cout<<point[i]->t<<"|"<<temp_index->t<<"|-100       ";*/
        if(i==15)                                                                 //in the 15th level, pause this loop
            break;
    }

    if(i>n_level)
        n_level=i;
    
    //cout<<endl;

    return true;                                                                  //means that this request has been accepted
}

int Depart(Index **Beginning,int &n_level)
{
    if(n_level==-1)
        return 0;
    
    int n=n_level;
    long T=Beginning[0]->right->t;
    int C=Beginning[0]->right->c;

    for(int i=0;i<=n;i++)
    {
        if(Beginning[i]->right->t!=T)
            break;
        
        Beginning[i]->right=Beginning[i]->right->right;
 
        if(Beginning[i]->right==NULL)
        {
            n_level=n_level-1;
            Beginning[i]=NULL;
        }
        
    }
    
    return C;
}

int Query(Index *Beginning[], int &n_level,long T)
{
    if(n_level==-1)
        return 0;

    Index *p=Beginning[n_level];
    while(p!=NULL)
    {
        p=Find(p,T);

        if(p->t==T)                                                            //we have found the flight
            return p->c;
        p=p->down;
    }

    return 0;                                                                           //we doesn't find the flight
}

bool Check(Index **Beginning,int n_level,long T,Index **point)
{
    Index *temp=Beginning[n_level];                                                                  //find the time which is near T
    int i=n_level;

    while(temp!=NULL)
    {
        temp=Find(temp,T);

        if(temp->t>=T-30)
            return false;
        else
        {
            if(temp->right==NULL||(temp->right->t)>T+30)
            {
                point[i]=temp;
                i--;
                temp=temp->down;
            }
            else
                return false;
        }
        
    }
    
    return true;
}

void Insert(Index *Index_Inserted,Index *Index_Near)
{

    Index_Inserted->right=Index_Near->right;                                               //the nodes near inserted node point to it
    Index_Near->right=Index_Inserted;
}

Index *Find(Index *index,long T)
{
    Index *p=index;
    if(p==NULL)
        return p;
    
    while(p->right!=NULL)
    {
        if(p->right->t<=T)
            p=p->right;
        else
            break;
    }

    return p;
}