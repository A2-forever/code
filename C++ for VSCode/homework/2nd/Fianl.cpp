#include <iostream>
#include<stdlib.h>

using namespace std;

struct NODE{
    int key;
    int order;
    NODE* left=NULL;
    NODE* right=NULL;
};

int FIND(NODE *root,int key);
NODE* BUILD(int *A,int n);
void INSERT(NODE *root,NODE *P);

int main()
{
    int n;
    long long count=0;
    NODE* head;

	int i;
	int j;
    int flag;
    bool flag1;

    
    cin>>n;

	int* A1 = new int[n]; 
	int* A2 = new int[n];
    

    for(i=0;i<n;i++)
    {
        cin>>A1[i];
    }

    for(i=0;i<n;i++)
    {
        cin>>A2[i];
    }
    
    head=BUILD(A2,n);

    for(i=0;i<n;i++)
    {
        flag1=1;
        for(j=0;j<n;j++)
        {
            flag=FIND(head,A1[j]);
            if(i==flag)
            {
                flag1=0;
            }
            if(flag>=i&&flag1==1)
            {
                count++;
            }
            
        }


    }

    cout<<count;

}

int FIND(NODE *root,int key)
{
    NODE *p=root;

    while(p!=NULL)
    {
        if(p->key==key)
        {
            return p->order;
        }

        if(key<p->key)
        {
            p=p->left;
        }
        else
        {
            p=p->right;
        }
        
    }
    return -1;
}


NODE* BUILD(int *A,int n)
{
    NODE *root=new NODE;
    root->key=A[0];
    root->order=0;

    for(int i=1;i<n;i++)
    {
        NODE *p=new NODE;
        p->key=A[i];
        p->order=i;
        INSERT(root,p);
    }
    return root;
}


void INSERT(NODE *root,NODE *p)
{
    NODE* x=root;
    NODE* y=NULL;

    while(x!=NULL)
    {
        y=x;
        if(p->key<x->key)
        {
            x=x->left;
        }
        else
        {
            x=x->right;
        }
    }

    if(y==NULL)
    {
        root=p;
    }
    else if(p->key<y->key)
    {
        y->left=p;
    }
    else
    {
        y->right=p;
    }
    
}