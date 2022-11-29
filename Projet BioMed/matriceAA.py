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
    Ly=(6.7)*1e-3
    Lx=(47.6)*1e-3
    a=(N/Lx)**2
    b=(M/Ly)**2
    c=-2*(a+b)
    eta=(18.5)*10**(-6)
    d=M/(2*eta*Lx)
    #b=2
    L=[i for i in range(M+2,N*(M+1)-1)]
    for i in range(2,N):
        L.remove(i*(M+1))
        L.remove(i*(M+1)-1)

    A1=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    A1[N*(M+1):,N*(M+1):]=np.eye(M+1)
    A1[:M+1,:M+1]=np.eye(M+1)
    '''for i in range(1,N):
        A1[i*(M+1)][i*(M+1)]=1
        A1[i*(M+1)][i*(M+1)+2]=1
        A1[i*(M+1)][i*(M+1)+1]=-4
        A1[i*(M+1)][i*(M+1)+1-M-1]=1
        A1[i*(M+1)][i*(M+1)+1+M+1]=1


        A1[(i+1)*(M+1)-1][(i+1)*(M+1)-1]=1
        A1[(i+1)*(M+1)-1][(i+1)*(M+1)-3]=1
        A1[(i+1)*(M+1)-1][(i+1)*(M+1)-2]=-4
        A1[(i+1)*(M+1)-1][(i+1)*(M+1)-2+M+1]=1
        A1[(i+1)*(M+1)-1][(i+1)*(M+1)-2-M-1]=1'''



    A2=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    for x in L:
        A2[x,x]=1
        A2[x,x-M-1]=-1
    A3=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    for x in L:
        A3[x,x+1]=1
        A3[x,x]=-1
    for i in range(1,N):
        A3[i*(M+1)][i*(M+1)+1]=1
        A3[i*(M+1)][i*(M+1)]=-1
        A3[(i+1)*(M+1)-1][(i+1)*(M+1)-2]=1
        A3[(i+1)*(M+1)-1][(i+1)*(M+1)-1]=-1


    A4=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    for x in L:
        A4[x,x]=d
        A4[x,x+M+1]=-d


    A5=np.identity((N+1)*(M+1))
    for x in L:
        A5[x,x]=c
        A5[x,x-1],A5[x,x+1]=b,b
        A5[x,x-(M+1)]=a
        A5[x,x+(M+1)]=a
    for j in range(M+1):
        A5[:M+1][j,j+N*(M+1)]=-1   #periodicité sur u
        A5[N*(M+1):][j,j]=-1
        A5[N*(M+1):][j,j+M+1]=1
        A5[N*(M+1):][j,j+(N-1)*(M+1)]=1
        A5[N*(M+1):][j,j+N*(M+1)]=-1
    A5[0][N*(M+1)]=0
    A5[M][-1]=0
    A5[N*(M+1)][0]=0
    A5[N*(M+1)][M+1]=0
    A5[N*(M+1)][(N-1)*(M+1)]=0

    A5[(N+1)*(M+1)-1][M]=0
    A5[(N+1)*(M+1)-1][2*(M+1)-1]=0
    A5[(N+1)*(M+1)-1][N*(M+1)-1]=0


    A6=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))

    A7=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    for x in L:
        A7[x,x]=d
        A7[x,x+1]=-d

    A8=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))

    A9=A5.copy()

    dA2=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    dA3=np.zeros(((N+1)*(M+1),(N+1)*(M+1)))
    dA1=dPmatriceA(M,N,a,b)


    A= np.block([[dA1,dA2,dA3],[A4,A5,A6],[A7,A8,A9]])


    #return(A1,A2,A3,A4,A5,A6,A7,A8,A9)
    return(A)


def matriceBB(M,N):
    P2=(1.037)*np.ones(M+1)
    P1=1*np.ones(M+1)
    B= np.zeros(3*(N+1)*(M+1))
    B[:M+1]=P2
    B[N*(M+1):(N+1)*(M+1)]=P1

    return(B)

def inconnuX(M,N):
    A = matriceAA(M,N)
    B = matriceBB(M,N)
    X=np.linalg.solve(A,B)
    P = X[:(N+1)*(M+1)]
    P = P.reshape((N+1,M+1))
    P = np.transpose(P)

    U = X[(N+1)*(M+1):2*(N+1)*(M+1)]
    U = U.reshape((N+1,M+1))
    U = np.transpose(U)

    V = X[2*(N+1)*(M+1):3*(N+1)*(M+1)]
    V = V.reshape((N+1,M+1))
    V = np.transpose(V)

    return(P,U,V)

def tracage(M,N):
    Ly=(6.7)*1e-3
    Lx=(47.6)*1e-3
    x=np.linspace(0,Lx,N+1)
    y=np.linspace(0,Ly,M+1)
    X,Y = np.meshgrid(x,y)
    fig = plt.figure(figsize = (11,7), dpi=100)
    plt.contourf(X, Y, inconnuX(M,N)[1],levels=50)
    plt.colorbar()
    plt.quiver(X, Y, inconnuX(M,N)[1], inconnuX(M,N)[2])
    #ax = sns.heatmap(champsU(M,N))
    plt.show()

