import numpy as np


eps = 1e-6


def Read_INT(INT):                              #在Gaussian中读取积分值
    (left_string,index_string)=INT.split("D")
    Left_float=float(left_string)
    index=int(index_string)
        
    return Left_float*(10**index)


def Read_TE_INT(line,TE_INT):                   #m = [[[[0] * 2] * 3] * 4]*5  定义一个5*4*3*2的四维数组的方法
    A=line.split()
    if(len(A)==0):                              #双电子积分读取还未开始，但进入了读取区域
        #print(0)
        return 0
    elif(A[0]!="I=" and A[0]!="Leave"):         #双电子积分读取还未开始，但进入了读取区域
        #print(0)
        return 0
    elif(A[0]=="Leave"):                        #双电子积分读取结束
        #print(1)
        return 1

    i=int(A[1])
    j=int(A[3])
    k=int(A[5])
    l=int(A[7])
    INT=Read_INT(A[9])

    TE_INT[i-1][j-1][k-1][l-1]=INT              #ijkl
    TE_INT[j-1][i-1][k-1][l-1]=INT              #jikl
    TE_INT[i-1][j-1][l-1][k-1]=INT              #ijlk
    TE_INT[j-1][i-1][l-1][k-1]=INT              #jilk
    
    TE_INT[k-1][l-1][i-1][j-1]=INT              #klij
    TE_INT[k-1][l-1][j-1][i-1]=INT              #klji
    TE_INT[l-1][k-1][i-1][j-1]=INT              #lkij
    TE_INT[l-1][k-1][j-1][i-1]=INT              #lkji
    
    return 2                                    #正常读取


def Read_Matrix_line(line):                     #在Gaussian中读取矩阵的单行的积分值
    temp=line.split()
    l=len(temp)

    line_Matrix=[0] * l
    if(l==0):
        return [],True                          #返回空数组，空行下一个
    elif(temp[0]=="SSO" or temp[0]=="***"):     #前者是核哈密顿矩阵读取结束，后者是overlap矩阵的读取结束
        return [],False                         #返回空数组，本次矩阵结束
    elif(l==1):                                 #开启一个新的序列，存储新的第二坐标,不过只有一个坐标
        line_Matrix[0]=int(temp[0])
        flag=True                               #开启一个新的序列，存储新的第二坐标,只是新的坐标只有一个
        return line_Matrix,flag

    if(len(temp[1])>10):
        for i in range(1,l):
            line_Matrix[i]=Read_INT(temp[i])
        line_Matrix[0]=int(temp[0])
        flag=False                              #没有开启一个新的序列
    else:
        for i in range(0,l):
            line_Matrix[i]=int(temp[i])
        flag=True                               #开启一个新的序列，存储新的第二坐标

    return line_Matrix,flag


def Read_Matrix_INT(line,M,index):              #在Gaussian中读取矩阵，并存储
    (temp,flag)=Read_Matrix_line(line)
    l=len(temp)
    if(l==0 and flag==True):                    #空行，下一个
        return 0
    elif(l==0 and flag==False):                 #本次矩阵读取结束
        return 1

    #数组存在的情况
    if(flag==True):
        for i in range(0,l):
            index[i]=temp[i]
    else:
        t=temp[0]
        for i in range(1,l):
            M[t-1,index[i-1]-1]=temp[i]
            M[index[i-1]-1,t-1]=temp[i]

    return 2                                    #正常读取


