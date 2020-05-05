#include <iostream>
struct naive
{
    char al;
    struct naive *p[2];
};

struct naive *Create() //新建结点并初始化
{
    struct naive *point = new struct naive;
    point->al = NULL;
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
    /*if(n>100)
    {
        cout<<"Wrong";
        return 0;
    }*/

    naive *Naive[100][2] = {NULL}; //it stores the points of the first and the last elements of each line

    for (int i = 0; i < n; i++)
    {        cin >> l >> s >> p;

        result = '0';
        cout << l << s << p << endl;

        if (p == 0)
        {
            cin >> ch;
            cout << ch << endl;
        }

        if (l >= 100)
            continue;

        switch (p)
        {
        case 0:
            cout << &Naive[l][s] << endl;
            cout << &Naive[l][1 - s] << endl;
            if (Naive[l][s] != NULL)
            {
                cout << Naive[l][s]->al << endl;
                cout << Naive[l][1 - s]->al << endl;
            }
            Insert(ch, Naive[l], s);
            break;
        case 1:
            result = Delete(Naive[l], s);
            break;
        case 2:
            result = Query(Naive[l][s]);
            break;
        }

        //if (result != '0')
        cout << "result=" << result << endl;

        cout << &Naive[l][s] << endl;
        cout << &Naive[l][1 - s] << endl;
        cout << Naive[l][s] << endl;
        cout << Naive[l][1 - s] << endl;

        if (Naive[l][s] != NULL)
        {
            cout << Naive[l][s]->al << endl;
            cout << Naive[l][1 - s]->al << endl;
        }

        cout << endl;
    }

    return 0;
}

void Insert(const char ch, naive *point[2], int s)
{

    cout << &point[s] << endl;
    cout << &point[1 - s] << endl;

    if (int(ch) < 97 || int(ch) > 122)
    {
        cout << "Wrong" << endl;
        return;
    }

    cout << "s=" << s << endl;

    if (point[s] != NULL)
        cout << point[s]->al;
    naive *New = Create();
    New->al = ch;

    cout << New << endl;

    if (point[s] == NULL)
    {
        point[0] = New;
        point[1] = New;
        cout << "we have create a new line" << endl;
        cout << point[s]->al << endl;
        cout << point[1 - s]->al << endl;
        cout << point[s] << endl;
        cout << point[1 - s] << endl;
        return;
    }
    else
    {
        cout << point[s] << " " << point[1 - s] << endl;
        cout << point[s]->al << endl;
        New->p[s] = point[s];
        cout << New->p[s]->al << endl;
        point[s]->p[1 - s] = New;
        cout << (point[s]->p[1 - s])->al << endl;
        point[s] = New;
        cout << point[s]->al << endl;
        cout << point[1 - s]->al << endl;
        return;
    }
}

char Delete(naive *point[2], int s)
{
    if (point == NULL)
        return '#';
    char result = point[s]->al;

    cout << point[s]->al << endl;
    cout << point[s]->p[s]->al << endl;
    point[s] = point[s]->p[s];
    //cout << point[s]->al << endl;
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
    char result = point->al;
    return result;
}