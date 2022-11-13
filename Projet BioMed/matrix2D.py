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

def matriceB(M,Lx,Ly,CD,CG,CB,CH,G):
    ''' CG=[f(O,y)]
    CD=[f(Lx,y)] y variant dans (0,Ly)
    CB=[f(x,0)]
    CH=[f(x,Ly)]'''
    x = np.linspace(0, Lx , M+1)
    y = np.linspace(0, Ly, M+1)
    X,Y = np.meshgrid(x, y)
    B=np.zeros((M+1,M+1))
    for i in range(M+1):
        for j in range(M+1):
            B[i,j]=G[i,j]
    B = B.reshape(((M+1)**2,1))
    B[0:M+1,0]=CG
    B[M*(M+1):,0]=CD
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
    F2 = np.transpose(F.reshape((M+1,M+1)))
    fig = plt.figure(figsize = (11,7), dpi=100)
    #plt.quiver(X, Y, F, np.zeros((M+1,M+1)))
    ax = sns.heatmap(F2)
    plt.show()




#Variables:
M = 15
Lx = 2
Ly = 1
a = (M/Lx)**2
b = (M/Ly)**2
x=np.linspace(0,Lx,M+1)
y=np.linspace(0,Ly,M+1)
X,Y = np.meshgrid(x,y)

    #Calcul du champs P:
        #Conditions auxbords:
P1=10
P2=5
#dp/dy = 0 a gauche et droite

def PmatriceA(M,a,b):
    A = matriceA(M,a,b)
    for i in range(1,M):
        A[i*(M+1)][i*(M+1)+1]=-1
        A[(i+1)*(M+1)-1][(i+1)*(M+1)-2]=-1
    return(A)

def PmatriceB(M,a,b):
    B = matriceB(M,Lx,Ly,P2,P1,np.zeros(M+1),np.zeros(M+1),np.zeros((M+1,M+1)))
    return(B)

def champsP(M):
    A = PmatriceA(M,a,b)
    B = PmatriceB(M,a,b)
    P = np.linalg.solve(A,B)
    P = P.reshape((M+1,M+1))
    return(np.transpose(P))



##Calcul du champs de vitesse:
U1 = 0
U2 = 0

eta = 1



#Calcul de U= ux:

def UmatriceA(M,a,b):
    A = matriceA(M,a,b)
    A[M*(M+1):]=np.zeros((M+1,(M+1)**2))
    for j in range(M+1):
        A[:M+1][j][j+M*(M+1)]=-1

        A[M*(M+1):,j][j]=-1
        A[M*(M+1):][j][j+M+1]=1
        A[M*(M+1):][j][j+(M+1)*(M-1)]=1
        A[M*(M+1):][j][j+(M+1)*M]=-1

    return(A)

def dPdx(M):
    P =  champsP(M)
    dP = np.zeros((M+1,M+1))
    d = eta*((a)**-0.5)
    for i in range(M+1):
        for j in range(M):
            dP[i,j]=(P[i,j+1]-P[i,j])/d
        dP[i,M]=(P[i,M]-P[i,M-1])/d
    return(dP)


def UmatriceB(M,a,b):
    B = matriceB(M,Lx,Ly,U2,U1,np.zeros(M+1),np.zeros(M+1),dPdx(M))
    return(B)


def champsU(M):
    A = UmatriceA(M,a,b)
    B = UmatriceB(M,a,b)
    U = np.linalg.solve(A,B)
    return(U)


#Calcul de V = uy:

V2 = 0
V1 = 0

def VmatriceA(M,a,b):
    A = matriceA(M,a,b)
    A[M*(M+1):]=np.zeros((M+1,(M+1)**2))
    for j in range(M+1):
        A[:M+1][j][j+M*(M+1)]=-1
        if j !=0 and j != M:
            A[j*(M+1)+1]=np.zeros((M+1)**2)
            A[j*(M+1)+1][j*(M+1)+1]=1
            A[j*(M+1)+1][j*(M+1)]=-1
            A[(j+1)*(M+1)-2]=np.zeros((M+1)**2)
            A[(j+1)*(M+1)-2][(j+1)*(M+1)-2]=1
            A[(j+1)*(M+1)-2][(j+1)*(M+1)-1]=-1

        A[M*(M+1):,j][j]=-1
        A[M*(M+1):][j][j+M+1]=1
        A[M*(M+1):][j][j+(M+1)*(M-1)]=1
        A[M*(M+1):][j][j+(M+1)*M]=-1

    return(A)

def dPdy(M):
    P =  champsP(M)
    dP = np.zeros((M+1,M+1))
    d = eta*((b)**-0.5)
    for j in range(M+1):
        for i in range(M):
            dP[i,j]=(P[i+1,j]-P[i,j])/d
        dP[M,j]=(P[M,j]-P[M-1,j])/d
    return(dP)


def VmatriceB(M,a,b):
    B = matriceB(M,Lx,Ly,V2,V1,np.zeros(M+1),np.zeros(M+1),dPdy(M))
    for i in range(1,M):
        B[i*(M+1)+1]=0
        B[(i+1)*(M+1)-2]=0
    return(B)


def champsV(M):
    A = VmatriceA(M,a,b)
    B = VmatriceB(M,a,b)
    V = np.linalg.solve(A,B)
    return(V)


fig = plt.figure(figsize = (11,7), dpi=100)
plt.quiver(Y, X, champsU(M) , champsV(M))
plt.show()

