import numpy as np
import math

def red(mat):   #reducing matrix by removing rows with less priority
    remove = []
    for i in range(len(mat)):
        for j in range(len(mat)):
            if mat[j][0]>mat[i][1] and mat[j][0]>mat[i][1]:
                #rmat.pop(i)
                remove.append(mat[i])

    rmat = [ele for ele in mat if ele not in remove]
    mat = rmat
    return mat


def Supp(ls):   #to return Support of a probability vector
    suppls = []
    for i in ls:
        if i>0:
            suppls.append(i)
    return suppls

def PSNE(mat):  #finding Pure Strategy Nash Equilibrium for n*2 matrix
    # invoke minimax and maximin
    cminmax=min(max(mat[i][j] for i in range(len(mat))) for j in range(len(mat[0])))
    rmaxmin=max(min(mat[i][j] for j in range(len(mat[0]))) for i in range(len(mat)))
    position = [-1,-1]  #if PSNE does not exist
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j]==cminmax and mat[i][j]==rmaxmin:
                position=[i+1,j+1]
    return position     #returns position where PSNE occurs

def matsum(mat1,mat2,a,b): #used to update empirical means for n*2 matrix with noise
    zero = [[0 for j in range(len(mat1[0]))] for i in range(len(mat1))]
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            mat1[i][j] = a*mat1[i][j]
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            mat2[i][j] = b*mat2[i][j]
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            zero[i][j] = mat1[i][j] + mat2[i][j]
    return zero

def sampling(mat):  #Introducing 1-sub-gaussian noise and querying values
    k=np.random.normal(0,1)
    for i in range(n):
        for j in range(2):
            mat[i][j]+=k

def dotp(ls1,ls2): #dot product of two vectors
    if len(ls1) == len(ls2):
        s = 0
        for i in len(ls1):
            s += ls1[i]*ls2[i]
        return s

def rtilda(i,mat,i1,i2):    #used to define r_i^tilda's as in algorithm 3 
    numerator = abs(mat[i1][0]-mat[i1][1])+abs(mat[i2][0]-mat[i2][1])
    denominator = abs(mat[i1][0]-mat[i1][1])+abs(mat[i2][0]-mat[i2][1])+abs(mat[i][0]-mat[i][1])
    frac = numerator/denominator
    return frac
def algo3(e,d,mat): #defining algorithm 3
    T = int(8*math.log(8*len(mat)/d)/(e*e))
    mean = [[0,0] for i in range(len(mat))] 
    n=0
    for t in range(1,T+1):
        mean=matsum(mean,sampling(mat),n/(n+1),1/(n+1))
        n+=1
        delta=(2*math.log(8*len(mat)*T/d)/t)**(0.5)
        dp = min(abs(mean[i][0] - mean[i][1]) for i in range(len(mat)))
        dq = min((abs(mean[j][0] - mean[k][0]) for j in range(k)) for k in range(1,len(mat)+1)) 
        dr = min((abs(mean[j][0] - mean[k][0]) for j in range(k)) for k in range(1,len(mat)+1))
        delmin = min(dp, dq, dr)
        check=(delmin+2*delta)/(delmin-2*delta)
        if check>=1 and check<=1.5 and PSNE(mean)[0]>-1:
            #print(n,mean)
            return PSNE(mean)
        elif check>=1 and check<=1.5 and PSNE(mean)[0]==-1:
            red(mean)
            for i in range(t,T+1):
                mean=matsum(mean,sampling(mat),n/(n+1),1/(n+1))
                n+=1
                deld = (2*math.log(8*n*T/d)/i)**(0.5)
                # Return nash equilibrium of A_bar, that is mean, say (x,y)
                # Return if i = T, nash equilibrium of A_bar
                # if len(Supp(x)) == 2:
                #       i1, i2 = Supp(x)[0], Supp(x)[1]
                #       for j in len(mean):
                #           if j not in supp(x):
                #               delg = min(rtilda(j,mean,i1,i2)*(value_of_mean - dotp(y,mean[j]))
                #       if delg >= 4*deld:
                #           supportnash = [[i1+1,i2+1],[1,2]]
                #           return supportnash
    # return nash equilibrium of A_bar, i.e. mean




