import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
import seaborn as sns
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)



def matriceA(M,a,b):
    c=-2*(a+b)
    A=np.identity((M+1)*(M+1))
    L=[i for i in range(M+2,M*(M+1)-1)]
    for i in range(2,M):
        L.remove(i*(M+1))
        L.remove(i*(M+1)-1)
    for x in L:
        A[x,x]=c
        A[x,x-1],A[x,x+1]=b,b
        A[x,x-(M+1)]=a
        A[x,x+(M+1)]=a
    return A

def matriceB(M,Lx,Ly,CD,CG,CB,CH,g):
    ''' CD=[f(O,y)]
    CG=[f(Lx,y)] y variant dans (0,Ly)
    CB=[f(x,0)]
    CH=[f(x,Ly)]'''
    x = np.linspace(0, Lx , M+1)
    y = np.linspace(0, Ly, M+1)
    X,Y = np.meshgrid(x, y)
    G=np.zeros((M+1,M+1))
    for i in range(M+1):
        for j in range(M+1):
            G[i,j]=g(x[i],y[j])
    B=G.reshape(((M+1)**2,1))
    B[:M+1,0]=CD
    B[M*(M+1):,0]=CG
    for i in range(1,M):
        B[i*(M+1)]=CB[i]
        B[(i+1)*(M+1)-1]=CH[i]
    return(B)


def reshape(f,M):
    F=f.reshape((M+1,M+1))
    return(F)


def D2(M,Lx,Ly,CD,CG,CB,CH,g):
    a = (M/Lx)**2
    b = (M/Ly)**2
    A = matriceA(M,a,b)
    B = matriceB(M,Lx,Ly,CD,CG,CB,CH,g)
    f = np.linalg.solve(A,B)
    F = reshape(f,M)
    return(A,B,f,np.transpose(F))

def tracage(F,M):
    x=np.linspace(0,1,M+1)
    y=np.linspace(0,1,M+1)
    X,Y = np.meshgrid(x,y)

    fig = plt.figure(figsize = (11,7), dpi=100)
    plt.quiver(X, Y, F, np.zeros((M+1,M+1)))
    #ax = sns.heatmap(F)
    plt.show()
