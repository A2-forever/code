import numpy as np
from FUNCTION_HF2 import *


eps = 1e-6

def Get_n(file_name):                                           #获取基函数数，alpha与beta电子数
    flag=0
    with open("Input\\"+file_name, 'r') as f:
        for line in f:
            temp=line.split()
            l=len(temp)
            if(l<6):
                continue

            if(temp[1]+" "+temp[2]=="basis functions,"):
                n=int(temp[0])                                  #基函数数
                flag=flag+1
            elif(temp[1]+" "+temp[2]=="alpha electrons"):
                n_alpha=int(temp[0])                            #alpha电子数
                n_beta=int(temp[3])                             #beta电子数
                flag=flag+1
            
            if(flag==2):
                break
        

    return n,n_alpha,n_beta                                     #基函数数，alpha与beta电子数

    
def Get_INT(file_name,n):                                       #获取重叠积分，核哈密顿积分以及双电子积分
    
    S=np.mat(np.zeros([n,n]))
    H=np.mat(np.zeros([n,n]))
    TE_INT=np.zeros([n,n,n,n])

    index=[0]*5                                                 #Gaussian一般一行五个
    count_completed = 0                                         #用于完成表示矩阵读取的进度,本次只需完成三个矩阵读取
    count_processing=0                                          #用于表示正在矩阵读取的进度,本次只需完成三个矩阵读取
    flag=0                                                      #0表示本行不读取，下一个；1表示读取结束；2表示正常读取



    with open("Input\\"+file_name, 'r') as f:
        for  line in  f:
            if(line==" *** Overlap *** \n" ):
                count_processing=1
                continue
            elif(line==" ****** Core Hamiltonian ****** \n"):
                count_processing=2
                continue
            elif(line==" *** Dumping Two-Electron integrals ***\n"):
                count_processing=3
                continue


            if(count_processing==count_completed):              #表示已经读取完一个矩阵，但还未开始下一个矩阵的读取
                continue                                        #一般完成进度与正在读取进度不一样，表示正在读取，但是在每次读取的开头会有一些特殊情况
            elif(count_processing==3):                          #正在读取双电子积分
                flag=Read_TE_INT(line,TE_INT)
                if(flag==1):
                    count_completed=3
                    break
            elif(count_processing==1):                          #正在读取交叉矩阵
                flag=Read_Matrix_INT(line,S,index)
                if(flag==1):
                    count_completed=1
            elif(count_processing==2):                          #正在读取核哈密顿矩阵
                flag=Read_Matrix_INT(line,H,index)
                if(flag==1):
                    count_completed=2

    return S,H,TE_INT



def Fock_G(P,TE_INT):                                           #基函数参数矩阵，与双电子积分数组（4维）
    n=P.shape[0]
    G = np.zeros_like(P)
    #a=1-3j
    #G=G+a

    for i in range(0,n):
        for j in range(0,n):                                    #Fock矩阵G矩阵的单位元构建
            for k in range(0,n):
                for l in range(0,n):
                    G[i,j]=G[i,j]+P[l,k]*(2*TE_INT[i][j][k][l]-TE_INT[i][l][k][j])
            #print(G[i][j])
    
    return G


'''
def Write_Matrix(M,file_name):                                  #在文件中写入二维矩阵
 
    for i in range(0,M.shape[0]):
        for j in range(0,M.shape[1]):
            f.write(str(S[i,j].real)+"\t")
        f.write("\n")
'''