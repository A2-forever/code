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

bool Random()
{
    bool r;
    r=rand()%2;
    return r;
};

Index *Create_Index(long T,int C)                                       //create an index pointing to flight
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

void Request(Index *Beginning[],int &n_level,long T,int C);                             //insert the node in the skiplist
Index *Depart(Index **Beginning,int &n_level);                                         //delete the earlist node and report the number of passengers 
Index *Query(Index *index, long T);                                                     //search the needed time in the skiplist

void Check(Index **Beginning,int n_level,long T,Index **point);                                      //check whether the conflict flight exists
void Insert(Index *Index_Inserted,Index *Index_Near);                                  //insert the node p near temp
Index *Find(Index *index,long T);                                                      //search the needed time in one level

int main()
{  
    Index *Beginning[20]={NULL};
    Index *p=NULL;
    int n_level=0;
    int n=0;
    long i=0;
    srand(time(NULL));

    int k=1;
    for(i=1;i<=500000;i++)
    {
        Request(Beginning,n_level,60*i,i);
    }
    for(i=1;i<=500000;i++)
    {
        Request(Beginning,n_level,60*(1000000-i+1),i+500000);
    }

    n=n_level;
    cout<<n_level<<endl;
    
    for(i=1;i<=100;i++)
    {
        k=i*60*10000;
        p=Query(Beginning[n_level],k);
        Output(p);
    }

    cout<<endl;

    k=1;
    for(i=1;i<=1000000;i++)
    {
        p=Depart(Beginning,n_level);
        if(i==k*10000)
        {
            k++;
            Output(p);
        }
    }

    return 0;
}  

void Request(Index *Beginning[], int &n_level,long T,int C)
{
    Index *down_index=NULL;                                                //reserve this index for the connection between near level
    Index *temp_index=NULL;                                                //create a new index pointing to p
    Index *point[21]={NULL};
    int i=-1;
    bool flag=1;
    
    if(Beginning[0]==NULL)
    {
        Beginning[0]=Create_Index(-100,-100);
        temp_index=Create_Index(T,C);
        Insert(temp_index,Beginning[0]);
    }

    Check(Beginning,n_level,T,point);                                      //check whether the conflict flight exists
    if(point[0]==NULL)
        return;                                                      //means that this request can't be accepted
    
    while(flag)
    {
        i++;                                                                      //the level positioned presently

        temp_index=Create_Index(T,C);                                               //create a new index pointing to p
        temp_index->down=down_index;                                               //connect two nodes vertically
        down_index=temp_index;                                                    //reserve this index for the connection between near level

        if(Beginning[i]==NULL)
        {
            Beginning[i]=Create_Index(-100,-100);                                      //the head node is empty to ensure each node is in its right side
            if(i>0)
                Beginning[i]->down=Beginning[i-1];
            point[i]=Beginning[i];
        }
        
        Insert(temp_index,point[i]);                                                  //insert new index into this skiplist in the ith level
        
        flag=Random();                                                            //determine whether to come to next floor
        
        /*if(temp_index->right!=NULL)
            cout<<point[i]->c<<"|"<<temp_index->c<<"|"<<temp_index->right->c<<" ";
        else
            cout<<point[i]->c<<"|"<<temp_index->c<<"|-1 ";
        */

        if(i==20)                                                                 //in the 15th level, pause this loop
            break;
    }

    if(i>n_level)
        n_level=i;
    
    //cout<<endl;

    return;                                                                  //means that this request has been accepted
}

Index *Depart(Index **Beginning,int &n_level)
{
    Index *temp=NULL;
    int n=n_level;
    long T=Beginning[0]->right->t;
    for(int i=0;i<=n;i++)
    {
        if(Beginning[i]->right->t!=T)
            break;
        
        temp=Beginning[i]->right;
        Beginning[i]->right=temp->right;

        if(Beginning[i]->right==NULL)
        {
            n_level=n_level-1;
            Beginning[i]=NULL;
        }
        
    }
    
    if(n_level==-1)
        n_level=0;
    
    return temp;
}

Index *Query(Index *index, long T)
{
    Index *p=index;
    while(p!=NULL)
    {
        p=Find(p,T);

        if(p->t==T)                                                            //we have found the flight
            return p;
        p=p->down;
    }

    return p;                                                                           //we doesn't find the flight
}

void Check(Index **Beginning,int n_level,long T,Index **point)
{
    Index *temp=Beginning[n_level];                                                                  //find the time which is near T
    int i=n_level;

    while(temp!=NULL)
    {
        temp=Find(temp,T);

        if(temp->t>=T-30)
            return;
        else
        {
            if(temp->right==NULL||(temp->right->t)>T+30)
            {
                point[i]=temp;
                i--;
                temp=temp->down;
            }
            else
                return;
        }
        
    }
    
    return;
}

void Insert(Index *Index_Inserted,Index *Index_Near)
{
    Index *temp_index=NULL;

    temp_index=Index_Near->right;                                                      //Near,inserted,temp
    Index_Near->right=Index_Inserted;                                               //the nodes near inserted node point to it
    Index_Inserted->right=temp_index;
    
}

Index *Find(Index *index,long T)
{
    Index *p=index;
    
    while(p->right!=NULL)
    {
        if(p->right->t<=T)
            p=p->right;
        else
            break;
    }

    return p;
}