from numpy import sqrt, ones, zeros
from numpy.linalg import matrix_rank
from sympy import symbols, simplify, expand
eps = 1e-6
x = symbols("x")


def Normalize(M):
    (m, n) = M.shape
    for i in range(0, n):
        temp_Normalization = 0
        for j in range(0, n):
            temp_Normalization += M[i, j] * M[i, j]
        temp_Normalization = sqrt(temp_Normalization)
        if M[i, 0] < -eps:
            temp_Normalization = -temp_Normalization
        for k in range(0, n):
            M[i, k] = M[i, k] / temp_Normalization


def Gauss(M):
    (m, n) = M.shape
    rank = matrix_rank(M)
    cout = rank
    n0 = 0
    record = -1 * ones(n, int)
    for i1 in range(0, n):
        if abs(M[i1, i1]) <= eps and i1 >= cout:
            record[n0] = i1
            n0 = n0 + 1
            continue
        elif abs(M[i1, i1]) <= eps and i1 < cout:
            flag = 0
            for position in range(i1 + 1, n):
                if abs(M[position, i1]) > eps:
                    swap(M, i1, position)
                    flag = 1
                    break
            if flag == 0:
                cout = cout + 1
                record[n0] = i1
                n0 = n0 + 1
                continue
        for i2 in range(0, n):
            if abs(M[i2, i1]) <= eps or i2 == i1:
                continue
            co = M[i2, i1] / M[i1, i1]
            for i3 in range(0, n):
                M[i2, i3] = M[i2, i3] - co * M[i1, i3]

    return record


def swap(M, i, j, str='row'):
    (n, n) = M.shape
    if str == 'row':
        for num in range(0, n):
            temp = M[i, num]
            M[i, num] = M[j, num]
            M[j, num] = temp
    elif str == 'column':
        for num in range(0, n):
            temp = M[num, i]
            M[num, i] = M[num, j]
            M[num, j] = temp
    else:
        print("Please enter 'row' or 'column'")


def det(M):
    M_size = M.shape
    n = M_size[0]
    if n == 1:
        return M[0, 0]

    sum = x
    for i in range(0, n):
        if M[0, i] == 0.0:
            continue
        N = create(M, i)
        sum += ((-1)**(i)) * M[0, i] * det(N)

    return expand(sum - x)


def create(M, k):
    size = M.shape
    n = size[0]
    N = zeros([n - 1, n - 1]) * x

    j1 = 0
    for j in range(0, n):
        if j == k:
            continue
        for i in range(1, n):
            N[i - 1, j1] = M[i, j]
        j1 += 1

    return N