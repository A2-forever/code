from numpy import eye,zeros,empty
from numpy.linalg import matrix_rank
from sympy import Matrix,symbols,det,solveset
from re import findall
import function

eps=1e-6
fa="φ"
x = symbols("x")

n = int(input("Please enter the number of carbons of your compound:"))
m = int(
    input(
        "Please enter the number of carbons which are in the main chain:"))

E = eye(n)
con_M = zeros([n,n])
for i in range(0, m - 1):
    con_M[i, i + 1] = 1
    con_M[i + 1, i] = 1

print(
    'Please enter the serial numbers of the two connected carbons which are not adjacent in serial number(enter any key to exit):'
)
while 1:
    i, j = map(str, input('i,j:').split())
    if not findall('^[0-9]+$', i):
        break
    if not findall('^[0-9]+$', j):
        break
    i = int(i)
    j = int(j)
    con_M[i - 1, j - 1] = 1
    con_M[j - 1, i - 1] = 1

p = Matrix(con_M) + Matrix(E) * x
d = det(p)
solution = solveset(d, x)
l = len(solution)
solution=list(solution)

k=0
store_c = zeros([n,n])
energy=empty(n)
for i in range(0,l):
    M = con_M + E * float(solution[i])

    rank=matrix_rank(M)
    n_zero=n-rank
    for j in range(k,k+n_zero):
        energy[j]=solution[i]
    
    record_p=function.Gauss(M)

    flag=0
    for I in range(k,k+n_zero):
        for J in range(0,n):
            if J==record_p[flag]:
                store_c[I,J]=1
            elif J not in record_p:
                store_c[I,J]=-M[J,record_p[flag]]/M[J,J]
        flag+=1

    
    k+=n_zero


function.Normalize(store_c)
print(store_c) 
print(energy)


with open('Output_Conjugate_Keyboard.txt','w') as f:
    for i in range(0,n):
        if energy[i]>eps:
            f.write("Energy: α+"+str('%.4f'%energy[i])+"β ： ")
        elif energy[i]<-eps:
            f.write("Energy: α"+str('%.4f'%energy[i])+"β ： ")
        else:
            f.write("Energy: α ： ")

        for j in range(0,n):
            if store_c[i,j]<-eps or j==0:
                f.write(str('%.4f'%store_c[i,j])+fa+str(j+1))
            else:
                f.write("+"+str('%.4f'%store_c[i,j])+fa+str(j+1))
        f.write("\n")
