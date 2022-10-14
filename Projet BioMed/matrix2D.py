import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import math
import seaborn as sns; sns.set()


np.set_printoptions(linewidth=999999)

def matriceM(mx,my,a,b):
    c=-2*(a+b)
    M=np.identity((mx+1)*(my+1))
    L=[i for i in range(my+2,mx*(my+1)-1)]
    for i in range(2,mx):
        L.remove(i*(my+1))
        L.remove(i*(my+1)-1)
    for x in L:
        M[x,x]=c
        M[x,x-1],M[x,x+1]=b,b
        M[x,x-(my+1)]=a
        M[x,x+(my+1)]=a
    return M

def B(u0,M,Lx,Ly):
    X=np.linspace(0,Lx,M+1)
    Y=np.linspace(0,Ly,M+1)
    U0=-(Y*(Y-Ly)) #u0 profil de vitesse initial
    B=[-10**6]*(M+1)*(M+1)
    mx=M
    my=M
    for i in range(1,mx+1):
        B[i*(my+1)]=0
        B[i*(my+1)-1]=0
    for k in range(M+1):
        B[k]=U0[k]
    for j in range(mx*(my+1),(mx+1)*(my+1)):
        B[j]=B[j-(mx*(my+1))]
    print("U0 is:",U0)
    B=np.array(B)
    B = B.reshape(len(B),1)
    return(B)


def reshap(U,M):
    u= [[0]*(M+1)]*(M)
    u=U.reshape(M+1,M+1)
    return(u)

# def u(x,y):
#     for i in range(len(X)):
#         for j in range(len(Y)):
#             if list(X)[i] == x and list(Y)[j] == y:
#                 return u[i][j]

def D2(M,L):
    a = (M/L)**2
    b = a
    A = matriceM(M,M,a,b)
    C = B(lambda y: -y*(y-L),M,L,L)
    U=np.linalg.solve(A,C)
    u = np.transpose(reshap(U,M))
    ax = sns.heatmap(u)
    plt.show()

    return(U,u)