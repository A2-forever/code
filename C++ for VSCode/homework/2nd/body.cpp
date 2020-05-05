#include <iostream>

using namespace std;
long long CountSort(int *A,int p,int r);
long long CountMerge(int *A,int p,int med,int r);

int main()
{
    int n;
	int i;
    long long count;

    cin>>n;
    
    int *A1=new int[n];
    int *A2=new int[n];

    for(i=0;i<n;i++)
    {
        A1[i]=i+1;
        A2[i]=n-i;
    }
   
    int *order=new int[n];
    for(i=0;i<n;i++)
    {
        order[A2[i]-1]=i;
    }
    delete [] A2;

    int *order1=new int[n];
    for(i=0;i<n;i++)
    {
        order1[i]=order[A1[i]-1];
    }
    delete [] order;
    delete [] A1;
    for(i=0;i<n;i++)
    {
        cout<<order1[i]<<' ';
    }
    cout<<endl;
    count=CountSort(order1,0,n-1);
    
    cout<<count;
}


long long CountSort(int *A,int p,int r)
{
    int l=r-p+1;
    if(l==1)
    {
        return 0;
    }
    long long n1,n2,n3,n;

    n1=CountSort(A,p,(p+r)/2);
    n2=CountSort(A,(p+r)/2+1,r);
    n3=CountMerge(A,p,(r-p)/2+1,r);

    n=n1+n2+n3;
    return n;
}


long long CountMerge(int *A,int p,int mid,int r)
{
    int i=p;
    int j=mid;
    int k=0;
    long long n=0;
    int *temp=new int[r-p+1];

    while(i<mid&&j<=r)
    {
        if(A[i] < A[j])
        {
            temp[k] = A[i];
            i++;
        }
        else
        {
            temp[k] = A[j];
            n = n + (mid-i);
            j++;
        }
        k++;
    }
    
    while(i<mid)
    {
        temp[k] = A[i];
        k++;
        i++;
    }

    while(j<=r)
    {
        temp[k] = A[j];
        k++;
        j++;
    }

    for(i=p;i<=r;i++)
    {
        A[i]=temp[i-p];
    }
    delete []temp;

    return n;
}