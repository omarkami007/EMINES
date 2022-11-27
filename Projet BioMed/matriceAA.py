import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math
import scipy
from scipy import sparse as sp
from scipy.sparse import dia_matrix

np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)


def matriceAA(M,N): #N>M
    Lx=2
    Ly=N*Lx/M
    a=(M/Lx)**2
    c=-4*a
    eta=10*(-8)
    b=M/(eta*Lx)
    #A5
    A5=np.identity((N+1)*(M+1))
    L=[i for i in range(M+2,N*(M+1)-1)]
    for i in range(2,N):
        L.remove(i*(M+1))
        L.remove(i*(M+1)-1)
    for x in L:
        A5[x,x]=c
        A5[x,x-1],A5[x,x+1]=a,a
        A5[x,x-(M+1)]=a
        A5[x,x+(M+1)]=a
    A9=A5.copy()
    A6=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    A8=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    A1=np.identity((N+1)*(M+1))
    for x in L:
        A1[x,x]=0
    for i in range(1,N):
        A1[i*(M+1)][i*(M+1)+1]=-1
        A1[(i+1)*(M+1)-1][(i+1)*(M+1)-2]=-1
    A3=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    for x in L:
        A3[x,x]=1
        A3[x,x-1]=-1
    '''for i in range(1,N):
        A3[i*(M+1)]=1
        A3[(i+1)*(M+1)-1]=1'''
    A7=-b*A3.copy()
    A2=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    for x in L:
        A2[x,x]=1
        A2[x,x-M-1]=-1
    A4=-b*A2.copy()
    A= np.bmat([[A1,A2,A3],[A4,A5,A6],[A7,A8,A9]])
    return(A)

def matriceBB(M,N):
    P2=10
    P1=5
    B= np.zeros(3*(N+1)*(M+1))
    B[:M+1]=P2
    B[N*(M+1):(N+1)*(M+1)]=P1
    return(B)

def inconnuX(M,N):
    A = matriceAA(M,N)
    B = matriceBB(M,N)
    X=np.linalg.solve(A,B)
    P = X[:(N+1)*(M+1)]
    P = P.reshape((M+1,N+1))
    P = np.transpose(P)

    U = X[(N+1)*(M+1):2*(N+1)*(M+1)]
    U = U.reshape((M+1,N+1))
    U = np.transpose(U)

    V = X[2*(N+1)*(M+1):3*(N+1)*(M+1)]
    V = V.reshape((M+1,N+1))
    V = np.transpose(V)
    return(P,U,V)