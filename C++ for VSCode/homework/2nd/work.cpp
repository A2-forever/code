#include <iostream>
#include<stdlib.h>

using namespace std;

long find(long *A,long key,long n);
void sort(long *A,long p,long r);
long RANDOM_PARTITION(long *A,long p,long r);
long PARTITION(long *A,long p,long r);
void swap(long *A,long i,long j);

int main()
{
    long n;
    long i=0;
    long j=0;
    long flag;
    bool flag1;
    bool flag2;
    long long count=0;

    
    cin>>n;

    long A1[n];
    long A2[n];
    long A3[n];
    bool A11[n];
    bool A21[n];
    

    for(i=0;i<n;i++)
    {
        cin>>A1[i];
    }

    for(i=0;i<n;i++)
    {
        cin>>A2[i];
    }

    
    for(i=0;i<n;i++)
    {
        A3[i]=A1[i];
    }

    sort(A3,0,n-1);

    for(i=0;i<n;i++)
    {
        flag1=1;
        flag2=1;
        for(j=0;j<n;j++)
        {
            flag=find(A3,A1[j],n);
            if(A1[j]==A3[i])
            {
                flag1=0;
            }
            
            A11[flag]=flag1;

            flag=find(A3,A2[j],n);
            if(A2[j]==A3[i])
            {
                flag2=0;
            }
            A21[flag]=flag2;
        }
        
        for(j=i;j<n;j++)
        {
            if(A11[j]!=A21[j])
            {
                count++;
            }
        }

    }

    cout<<count;
}

long find(long *A,long key,long n)
{
    long p=0;
    long r=n-1;
    long flag=(p+r)/2;

    while(p!=r)
    {
        if(A[flag]==key)
        {
            return flag;
        }
        else if(A[flag]>key)
        {
            r=flag-1;
            flag=(r+p)/2;
        }
        else
        {
            p=flag+1;
            flag=(r+p)/2;
        }
        
    }

    return p;
}

void sort(long *A,long p,long r)
{
    if(p<r)
    {
        long q;
        q=RANDOM_PARTITION(A,p,r);
        sort(A,p,q-1);
        sort(A,q+1,r);
    }
}

long RANDOM_PARTITION(long *A,long p,long r)
{
        long flag;
        flag=(rand() % (r-p+1))+ p;
        swap(A,flag,r);

        return PARTITION(A,p,r);
}

long PARTITION(long *A,long p,long r)
{
    long x=A[r];
    long key=p-1;

    for(long i=p;i<r;i++)
    {
        if(A[i]<=x)
        {
            key++;
            swap(A,i,key);
        }
    }

    swap(A,key+1,r);
    return key+1;
}


void swap(long *A,long i,long j)
{
    long temp=A[i];
    A[i]=A[j];
    A[j]=temp;
}
