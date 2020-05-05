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

with open("input\\Input_Conjugate.txt", 'r') as f:
    cout = 0
    for line in f:
        if cout == 0:
            line.count('\n')
            (n, m) = line.split()
            try:
                n = int(n)
                m = int(m)
            except:
                with open('input\\Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            if n <= 0 or m < 0 or m < n:
                with open('input\\Output_Conjugate.txt', 'w') as f:
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
                with open('input\\Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            if ne > 2 * n or ne <= 0:
                with open('input\\Output_Conjugate.txt', 'w') as f:
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
                with open('input\\Output_Conjugate.txt', 'w') as f:
                    f.write(
                        "Please enter the right serial numbers of the carbons!"
                    )
                    exit()
            if i > n or j > n or i < 1 or j < 1:
                with open('input\\Output_Conjugate.txt', 'w') as f:
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
solution = solveset(d, x)
print(d)
l = len(solution)
solution = list(solution)

k = 0
store_c = zeros([n, n])
energy = empty(n)
for i in range(0, l):
    if i > 0:
        if abs(solution[i] - solution[i - 1]) < eps:
            continue

    M = con_M + E * float(solution[i])

    rank = matrix_rank(M)
    n_zero = n - rank
    for j in range(k, k + n_zero):
        energy[j] = solution[i]

    record_p = function.Gauss(M)

    flag = 0
    for I in range(k, k + n_zero):
        for J in range(0, n):
            if J == record_p[flag]:
                store_c[I, J] = 1
            elif J not in record_p:
                store_c[I, J] = -M[J, record_p[flag]] / M[J, J]
        flag += 1
    k += n_zero

for i in range(0, n - 1):
    if abs(energy[i] - energy[i + 1]) < eps:
        for j in range(0, n):
            store_c[i, j] = store_c[i, j] - store_c[i + 1, j]
            store_c[i + 1, j] = store_c[i, j] + 2*store_c[i + 1, j]

function.Normalize(store_c)

t1 = int(ne / 2) - 1
if abs(energy[t1] - energy[t1 + 1]) > eps:
    t2 = t1
else:
    t2 = t1 + 1

e_density = zeros(n)
e_bond = zeros([n, n])

for i in range(0, n):
    for j in range(0, t1):
        e_density[i] += 2 * store_c[j, i] * store_c[j, i]
    e_density[
        i] += store_c[t1, i] * store_c[t1, i] + store_c[t2, i] * store_c[t2, i]

for i in range(0, n):
    for j in range(i, n):
        if con_M[i, j] > eps:
            for k in range(0, t1):
                e_bond[i, j] += 2 * store_c[k, i] * store_c[k, j]
            e_bond[i, j] += store_c[t1, i] * store_c[t1, j] + store_c[
                t2, i] * store_c[t2, j]

e_free = ones(n) * sqrt(3)
for i in range(0, n):
    for j in range(0, n):
        if con_M[i, j] > eps:
            e_free[i] = e_free[i] - e_bond[i, j] - e_bond[j, i]

with open('input\\Output_Conjugate.txt', 'w') as f:
    f.write("Molecular Energy Level and Orbital: ")
    f.write("\n")
    for i in range(0, n):
        if energy[i] > eps:
            f.write("Energy: α" + str('%.4f' % (-energy[i])) + "β ：")
        elif energy[i] < -eps:
            f.write("Energy: α+" + str('%.4f' % (-energy[i])) + "β ：")
        else:
            f.write("Energy: α " + '\t' + "：")

        for j in range(0, n):
            if store_c[i, j] > eps and j != 0:
                f.write("+")
            f.write(str('%.4f' % store_c[i, j]) + fa + str(j + 1) + '    ')
        f.write("\n")
    f.write("\n")

    f.write("Charge Density: ")
    f.write("\n")
    for i in range(0, n):
        f.write("Carbon" + str(i + 1) + ": ")
        f.write(str('%.4f' % e_density[i]))
        f.write("\n")
    f.write("\n")

    f.write("Bond Level: ")
    f.write("\n")
    for i in range(0, n):
        for j in range(i, n):
            if con_M[i, j] > eps:
                f.write("Carbon" + str(i + 1) + " and Carbon " + str(j + 1) +
                        " is: ")
                f.write(str('%.4f' % e_bond[i, j]))
                f.write("    ")
        f.write("\n")
    f.write("\n")

    f.write("Free Valence: ")
    f.write("\n")
    for i in range(0, n):
        f.write("Carbon" + str(i + 1) + ": ")
        f.write(str('%.4f' % e_free[i]))
        f.write("\n")
    f.write("\n")

    t2 = process_time()
    t = t1 - t2
    f.write("cost " + str(t) + "s")
