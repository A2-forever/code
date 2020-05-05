#include <iostream>
struct naive
{
    char al;
    struct naive *p[2];
};

struct naive *Create() //新建结点并初始化
{
    struct naive *point = new struct naive;
    point->p[0] = NULL;
    point->p[1] = NULL;
    return point;
}

using namespace std;
void Insert(const char ch, naive *point[2], int s);
char Delete(naive *point[2], int s);
char Query(naive *point);

int main()
{
    int n;         //the number of lines we will have
    int l, s, p;   //lines,start or end,operation
    char ch = '0'; //insert lowcase letter 
    char result = '0';

    cin >> n;

    naive *Naive[100][2] = {NULL}; //it stores the points of the first and the last elements of each line

    for (int i = 0; i < n; i++)
    {
        cin >> l >> s >> p;

        if (p == 0)
            cin >> ch;
        
        result = '0';

        if (l >= 100)
            continue;

        switch (p)
        {
        case 0:
            Insert(ch, Naive[l], s);
            break;
        case 1:
            result = Delete(Naive[l], s);
            break;
        case 2:
            result = Query(Naive[l][s]);
            break;
        }

        if (result != '0')
            cout << result << endl;
    }

    return 0;
}

void Insert(const char ch, naive *point[2], int s)
{
    if (int(ch) < 97 || int(ch) > 122)
    {
        cout << "Wrong" << endl;
        return;
    }

    naive *New = Create();
    New->al = ch;

    if (point[s] == NULL)
    {
        point[0] = New;
        point[1] = New;
        return;
    }
    else
    {
        New->p[s] = point[s];
        point[s]->p[1 - s] = New;
        point[s] = New;
        return;
    }
}

char Delete(naive *point[2], int s)
{
    if (point[s] == NULL)
        return '#';
    char result = point[s]->al;

    point[s] = point[s]->p[s];
    if (point[s] != NULL)
        point[s]->p[1 - s] = NULL;
    else
        point[1 - s] = NULL;

    return result;
}

char Query(naive *point)
{
    if (point == NULL)
        return '#';

    return point->al;
}