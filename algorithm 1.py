import numpy as np
import math
def PSNE(mat):
    rowmin1=1
    if mat[0][0]>mat[0][1]:
        rowmin1=2
    rowmin2=1
    if mat[1][0]>mat[1][1]:
        rowmin2=2
    rowmax=[1,rowmin1]
    if mat[0][rowmin1-1]<mat[1][rowmin2-1]:
        rowmax=[2,rowmin2] #rowmaxmin
    colmax1=1
    if mat[0][0]<mat[1][0]:
        colmax1=2
    colmax2=1
    if mat[0][1]<mat[1][1]:
        colmax2=2
    colmin=[colmax1,1]
    if mat[colmax1-1][0]>mat[colmax2-1][1]:
        colmin=[colmax2,2] #colminmax
    if colmin==rowmax:
        return colmin #the position of saddle point
    else:
        return [-1,-1]
def MSNE(mat):
    a=mat[0][0]
    b=mat[0][1]
    c=mat[1][0]
    d=mat[1][1]
    q=(d-b)/(a-b-c+d) #by solving equations formed by equating gains
    p=(d-c)/(a-b-c+d)
    if q<0 or q>1 or p<0 or p>1:
        return PSNE(mat)
    return [[p,1-p],[q,1-q]]
def matsum(mat1,mat2,a,b):    #weighted sum of matrices
    zero = [[0 for j in range(2)] for i in range(2)]
    for i in range(2):
        for j in range(2):
            mat1[i][j] = a*mat1[i][j]
    for i in range(2):
        for j in range(2):
            mat2[i][j] = b*mat2[i][j]
    for i in range(2):
        for j in range(2):
            zero[i][j] = mat1[i][j] + mat2[i][j]
    return zero
def sampling(mat):    #inserting noise inside the game matrix
    k=np.random.normal(0,1)
    for i in range(2):
        for j in range(2):
            mat[i][j]+=k
    return mat
def algo1(e,d,mat):      #finds epsilon good solution
    T=int(8*math.log(16/d)/(e*e))
    mean=[[0,0],[0,0]]
    n=0              #no. of samples
    for t in range(1,T+1):
        mean=matsum(mean,sampling(mat),n/(n+1),1/(n+1))
        n+=1
        delta=(2*math.log(16*T/d)/t)**(0.5)
        delmin=min(abs(mean[0][0]-mean[0][1]),abs(mean[1][0]-mean[1][1]),abs(mean[0][0]-mean[1][0]),abs(mean[0][1]-mean[1][1]))
        delbar=abs(mean[0][0]-mean[0][1]-mean[1][0]+mean[1][1])
        check=(delmin+2*delta)/(delmin-2*delta)
        if check>=1 and check<=1.5 and PSNE(mean)[0]>-1:
            return PSNE(mean)
        elif check>=1 and check<=1.5 and delbar<10*e:
            for i in range(T-t):
                mean=matsum(mean,sampling(mat),n/(n+1),1/(n+1))
                n+=1
            return MSNE(mean)
        elif check>=1 and check<=1.5 and delbar>=10*e:
            N=(80*math.log(16*T/d)/(e*delbar))
            if N>T-t:
                for i in range(T-t):
                    mean=matsum(mean,sampling(mat),n/(n+1),1/(n+1))
                    n+=1
                return MSNE(mean)
            for i in range(N):
                mean=matsum(mean,sampling(mat),n/(n+1),1/(n+1))
                n+=1
            return MSNE(mean)
    return MSNE(mean)