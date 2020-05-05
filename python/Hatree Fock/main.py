from time import process_time

import numpy as np
from numpy import linalg as LA
from FUNCTION_HF1 import *


eps = 1e-6
t1 = process_time()

with open("Input\input.json","r") as f:
    line=f.readlines()
file=line[0].split(".")
file_name=file[0]+"."+file[1]

n,n_alpha,n_beta=Get_n(file_name)                       #基函数数，alpha与beta电子数
S,H,TE_INT=Get_INT(file_name,n)                         #获取重叠积分，核哈密顿积分以及双电子积分


s,U=LA.eigh(S)                                          #U是正交矩阵或者说酉矩阵。s是其特征值,S是overlap矩阵，是一个厄米矩阵（或者实对称矩阵）
                                                        #eigh针对对称矩阵，对称矩阵的分解有单独的算法
s1=1./np.sqrt(s)                                        #将所有特征值开根号，取倒数
S_sqrt=np.diag(s1)                                      #张成对应矩阵
X=U*S_sqrt                                              #X=U*S_sqrt(s)


#HF_SCF方法的具体迭代
Energy=0                                                #初始能量设为0
C =np.mat(np.zeros([n,n]))                              #初始密度矩阵设为0
P=C[:,0:n_alpha]*C[:,0:n_alpha].T                       #设置初始密度矩阵，n个alpha电子的闭壳层

Var=1
while(Var>=eps):                                        #误差在接受范围以内
    Energy_old=Energy                                   #上次迭代的能量
    G=Fock_G(P,TE_INT)                                  #构建G矩阵
    F=H+G                                               #Fock矩阵的构建

    F_U=X.T*F*X                                         #F`=X+FX，F'矩阵的构建
                                                        #X为正交化C的矩阵，目的在于消除交叉

    Eigenvalue,C_U=LA.eig(F_U)
    sorted_indices = np.argsort(Eigenvalue)             #从小到大排序，存储下标
    Eigenvalue=Eigenvalue[sorted_indices]               #能量排序
    C_U=C_U[:,sorted_indices]                           #系数排序

    C=X*C_U                                             #新的基函数系数
    
    
    P=C[:,0:n_alpha]*C[:,0:n_alpha].T                   #新的密度矩阵
    Energy=np.trace(P*(H+F))                            #P*(H+F)的对角线元素之和(矩阵的迹)

    Var=np.abs(Energy_old-Energy)                       #本次计算的能量与上次的差
    print(Var)


t2 = process_time()
with open("input\\output "+file[0]+".json ", 'w') as f:
    f.write(file[0]+"\n")                               #体系名称
    f.write(str(n)+" basis functions \t")               #基函数数目
    f.write(str(n_alpha)+" alpha electrons \t")         #alpha电子数目
    f.write(str(n_beta)+" beta electrons ")             #beta电子数目
    f.write("\n\n")

    f.write("Energy:\n"+str(Energy))                    #体系基态能量
    f.write("\n\n")

    f.write("Eigenvalue:\n")                            #本征值
    for i in range(0,len(Eigenvalue)):
        f.write(str(Eigenvalue[i])+"\t")
    f.write("\n\n")

    f.write("C:\n")
    for i in range(0,C.shape[0]):
        for j in range(0,C.shape[1]):
            f.write(str(C[i,j])+"\t")
        f.write("\n")
    
    f.write("\n")

    f.write("X:\n")
    for i in range(0,X.shape[0]):
        for j in range(0,X.shape[1]):
            f.write(str(X[i,j])+"\t")
        f.write("\n")

    f.write("\n")

    f.write("Process time: "+str(t2-t1)+"s")            #运算时间
