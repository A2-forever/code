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


Index *Create_Index(long T=-100,int C=-100)
{
    struct Node *p = new struct Node;
    p->t=T;
    p->c=C;
    p->right=NULL;
    p->down=NULL;
    return p;
};

void Request(Index *Beginning[],int &n_level,long T,int C);
void Depart(Index *Beginning[],int &n_level);
void Query(Index *Beginning[], int &n_level,long T);

bool Check(Index **Beginning,int n_level,long T,Index **point);
void Insert(Index *Index_Inserted,Index *Index_Near);
Index *Find(Index *index,long T);

int main()
{
    long n;

    int op, c;
    long t;
    int n_level=-1;

    Index *Beginning[16]={NULL};

    scanf("%ld",&n);

    srand(time(NULL));
    for (long i = 0; i < n; i++)
    {
        scanf("%d",&op);
        //cin >> op;
        if (op == 0)
            scanf("%ld%d",&t,&c);
            //cin >> t >> c;
        else if(op==2)
            scanf("%ld",&t);
            //cin>>t;

        switch (op)
        {
        case 0:
            Request(Beginning,n_level,t,c);
            break;
        case 1:
            Depart(Beginning,n_level);
            break;
        case 2:
            Query(Beginning,n_level,t);
            break;
        }
    }

    return 0;
}



void Request(Index *Beginning[], int &n_level,long T,int C)
{
    Index *down_index=NULL;
    Index *temp_index=NULL;
    Index *point[16]={NULL};
    int i=-1;
    bool flag=1;
    
    if(Beginning[0]!=NULL)
    {
        flag=Check(Beginning,n_level,T,point);
        if(!flag)
            return;
    }

    while(flag)
    {
        i++;

        temp_index=Create_Index(T,C);
        temp_index->down=down_index;
        down_index=temp_index;

        if(Beginning[i]==NULL)
        {
            Beginning[i]=Create_Index();
            if(i>0)
                Beginning[i]->down=Beginning[i-1];
            point[i]=Beginning[i];
        }
        
        Insert(temp_index,point[i]);
        
        flag=rand()%2;
        
        if(i==15)
            break;
    }

    if(i>n_level)
        n_level=i;
    

    return;
}

void Depart(Index **Beginning,int &n_level)
{
    if(n_level==-1)
    {
        cout<<"0"<<endl;
        return;
    }

    
    printf("%d\n",Beginning[0]->right->c);

    int n=n_level;
    long T=Beginning[0]->right->t;

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
    
    return;
}

void Query(Index *Beginning[], int &n_level,long T)
{
    if(n_level==-1)
    {
        printf("0\n");
        return;
    }

    Index *p=Beginning[n_level];
    while(p!=NULL)
    {
        p=Find(p,T);

        if(p->t==T)
        {
            printf("%d\n",p->c);
            return;
        }
        p=p->down;
    }

    printf("0\n");
    return;
}

bool Check(Index **Beginning,int n_level,long T,Index **point)
{
    Index *temp=Beginning[n_level];
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

    Index_Inserted->right=Index_Near->right;
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