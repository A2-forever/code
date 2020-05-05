from sys import exit
from time import process_time

from numpy import empty, eye, ones, sqrt, zeros
from numpy.linalg import matrix_rank
from sympy import solveset, symbols, expand, simplify_logic

import function

t1 = process_time()
eps = 1e-6
fa = "φ"
x = symbols("x")

with open("Input_Conjugate.txt", 'r') as f:
    cout = 0
    for line in f:
        if cout == 0:
            line.count('\n')
            (n, m) = line.split()
            try:
                n = int(n)
                m = int(m)
            except:
                with open('Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            if n <= 0 or m < 0 or m < n:
                with open('Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            con_M = zeros([n, n])
            cout = 1
        elif cout == 1:
            if line.count('\n') == len(line):
                cout = 2
                continue
            ne = line.split()
            try:
                ne = int(ne[0])
            except:
                with open('Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            if ne > 2 * n or ne <= 0:
                with open('Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            cout = 2
        else:
            if line.count('\n') == len(line):
                continue
            (i, j) = line.split()
            try:
                i = int(i)
                j = int(j)
            except:
                with open('Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            if i > n or j > n or i < 1 or j < 1:
                with open('Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            con_M[i - 1, j - 1] = 1
            con_M[j - 1, i - 1] = 1

E = eye(n)
for i in range(0, m - 1):
    con_M[i, i + 1] = 1
    con_M[i + 1, i] = 1

p = con_M + E * x
d = function.det(p)
print(d)