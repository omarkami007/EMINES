import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)




M = 15
N = 15
Ly=(6.7)*1e-3
Lx=(47.6)*1e-3
a = (N/Lx)**2
b = (M/Ly)**2
c=-2*(a+b)

x=np.linspace(0,Lx,N+1)
y=np.linspace(0,Ly,M+1)
X,Y = np.meshgrid(x,y)


def matriceA(M,N,a,b):
    c=-2*(a+b)
    A=np.identity((N+1)*(M+1))
    L=[i for i in range(M+2,N*(M+1)-1)]
    for i in range(2,N):
        L.remove(i*(M+1))
        L.remove(i*(M+1)-1)
    for x in L:
        A[x,x]=c
        A[x,x-1],A[x,x+1]=b,b
        A[x,x-(M+1)]=a
        A[x,x+(M+1)]=a
    return(A)

def matriceB(M,N,Lx,Ly,CD,CG,CB,CH,G):
    ''' CG=[f(O,y)]
    CD=[f(Lx,y)] y variant dans (0,Ly)
    CB=[f(x,0)]
    CH=[f(x,Ly)]'''
    x = np.linspace(0, Lx , N+1)
    y = np.linspace(0, Ly, M+1)
    X,Y = np.meshgrid(x, y)
    B=np.zeros((M+1,N+1))
    for i in range(M+1):
        for j in range(N+1):
            B[i,j]=G[i,j]
    B = B.reshape(((N+1)*(M+1),1))
    B[0:M+1,0]=CG
    B[N*(M+1):,0]=CD
    for i in range(1,N):
        B[i*(M+1)]=CB[i]
        B[(i+1)*(M+1)-1]=CH[i]
    return(B)



##Calcul du champs P:
        #Conditions auxbords:
P1=(1.037)*np.ones(M+1)
P2=1*np.ones(M+1)

#dp/dy = 0 en haut et bas


def dPmatriceA(M,N,a,b):
    A = matriceA(M,N,a,b)
    for i in range(1,N):
        A[i*(M+1)][i*(M+1)+1]=-1
        A[(i+1)*(M+1)-1][(i+1)*(M+1)-2]=-1


    return(A)





def dPmatriceA2(M,N,a,b):
    A = matriceA(M,N,a,b)
    for i in range(1,N):
        A[i*(M+1)][i*(M+1)+1]=-1             #dp/dy=0 en bas
        A[(i+1)*(M+1)-1][(i+1)*(M+1)-2]=-1   # **** en haut

    for i in range(3,5):
        A[i*(M+1)+1]=np.zeros((M+1)*(N+1))
        A[i*(M+1)+1][i*(M+1)+1]=1         #dp/dy=0 au bord de la defor
        A[i*(M+1)+1][i*(M+1)+2]=-1
    A[2*(M+1)+1]=np.zeros((M+1)*(N+1))  #dp/dx=0 au bords verticaux de la deform
    A[2*(M+1)+1][2*(M+1)+1]=1
    A[2*(M+1)+1][3*(M+1)+1]=-1
    A[5*(M+1)+1]=np.zeros((M+1)*(N+1))
    A[5*(M+1)+1][5*(M+1)+1]=1
    A[5*(M+1)+1][4*(M+1)+1]=-1


    return(A)




def dPmatriceB(M,N,a,b):
    #dP/dy=0 aux bords
    B = matriceB(M,N,Lx,Ly,P2,P1,np.zeros(N+1),np.zeros(N+1),np.zeros((M+1,N+1)))

    return(B)

def dPmatriceB2(M,N,a,b):
    #dP/dy=0 aux bords
    B = matriceB(M,N,Lx,Ly,P2,P1,np.zeros(N+1),np.zeros(N+1),np.zeros((M+1,N+1)))
    for i in range(3,5):
        B[i*(M+1)+1]=0
    B[2*(M+1)+1]=0
    B[5*(M+1)+1]=0

    return(B)



def dchampsP(M,N):
    A = dPmatriceA(M,N,a,b)
    B = dPmatriceB(M,N,a,b)
    P = np.linalg.solve(A,B)
    P = P.reshape((N+1,M+1))
    return(np.transpose(P))

def dchampsP2(M,N):
    A = dPmatriceA2(M,N,a,b)
    B = dPmatriceB2(M,N,a,b)
    P = np.linalg.solve(A,B)
    P = P.reshape((N+1,M+1))
    return(np.transpose(P))


## Calcul de Ux:
U1 = 0
U2 = 0

eta=(18.5)*10**(-6)


def UmatriceA(M,N,a,b):
    A = matriceA(M,N,a,b)
    A[N*(M+1):]=0
    for j in range(M+1):
        A[:M+1][j,j+N*(M+1)]=-1   #periodicité sur u
        A[N*(M+1):][j,j]=-1
        A[N*(M+1):][j,j+M+1]=1
        A[N*(M+1):][j,j+(N-1)*(M+1)]=1
        A[N*(M+1):][j,j+N*(M+1)]=-1
    for i in range(3,5):
        A[i*(M+1)+1]=np.zeros((M+1)*(N+1))
        A[i*(M+1)+1][i*(M+1)+1]=1
    A[2*(M+1)+1]=np.zeros((M+1)*(N+1))
    A[5*(M+1)+1]=np.zeros((M+1)*(N+1))
    A[2*(M+1)+1][2*(M+1)+1]=1
    A[2*(M+1)+1][3*(M+1)+1]=-1
    A[5*(M+1)+1][5*(M+1)+1]=1
    A[5*(M+1)+1][4*(M+1)+1]=-1
    return(A)

def dPdx(M,N):
    P =  dchampsP2(M,N)
    dP = np.zeros((M+1,N+1))
    d = eta*((a)**-0.5)
    for j in range(M+1):
        for i in range(N):
            dP[j,i]=(P[j,i+1]-P[j,i])/d
        dP[j,N]=(P[j,N]-P[j,N-1])/d
    return(dP)


def UmatriceB(M,N,a,b):
    B = matriceB(M,N,Lx,Ly,U2,U1,np.zeros(N+1),np.zeros(N+1),dPdx(M,N))
    for i in range(3,5):
        B[i*(M+1)+1]=0
    B[2*(M+1)+1]=0
    B[5*(M+1)+1]=0
    return(B)


def champsU(M,N):
    A = UmatriceA(M,N,a,b)
    B = UmatriceB(M,N,a,b)
    U = np.linalg.solve(A,B)
    U = np.transpose(U.reshape((N+1,M+1)))
    return(U)

##Calcul de V = uy:

V2 = 0
V1 = 0

def VmatriceA(M,N,a,b):
    A = matriceA(M,N,a,b)
    A[N*(M+1):]=np.zeros((M+1,(N+1)*(M+1)))
    for j in range(M+1):
        A[:M+1][j][j+N*(M+1)]=-1
    for i in range(1,N):
        A[i*(M+1)+1]=np.zeros((N+1)*(M+1))
        A[i*(M+1)+1][i*(M+1)+1]=1
        A[i*(M+1)+1][i*(M+1)]=-1
        A[(i+1)*(M+1)-2]=np.zeros((N+1)*(M+1))
        A[(i+1)*(M+1)-2][(i+1)*(M+1)-2]=1
        A[(i+1)*(M+1)-2][(i+1)*(M+1)-1]=-1
    for j in range(M+1):
        A[N*(M+1):,j][j]=-1
        A[N*(M+1):][j][j+M+1]=1
        A[N*(M+1):][j][j+(M+1)*(N-1)]=1
        A[N*(M+1):][j][j+(M+1)*N]=-1
    for i in range(3,5):
        A[i*(M+1)+1]=np.zeros((M+1)*(N+1))
        A[i*(M+1)+1][i*(M+1)+1]=1
        A[i*(M+1)+2]=np.zeros((M+1)*(N+1))
        A[i*(M+1)+2][i*(M+1)+2]=1
        A[i*(M+1)+2][i*(M+1)+1]=-1


    return(A)

def dPdy(M,N):
    P =  dchampsP2(M,N)
    dP = np.zeros((M+1,N+1))
    d = eta*((b)**-0.5)
    for i in range(N+1):
        for j in range(M):
            dP[j,i]=(P[j+1,i]-P[j,i])/d
        dP[M,i]=(P[M,i]-P[M-1,i])/d
    return(dP)


def VmatriceB(M,N,a,b):
    B = matriceB(M,N,Lx,Ly,V2,V1,np.zeros(N+1),np.zeros(N+1),dPdy(M,N))
    for i in range(1,N):
        B[i*(M+1)+1]=0
        B[(i+1)*(M+1)-2]=0
    for i in range(3,5):
        B[i*(M+1)+1]=0
        B[i*(M+1)+2]=0

    return(B)


def champsV(M,N):
    A = VmatriceA(M,N,a,b)
    B = VmatriceB(M,N,a,b)
    V = np.linalg.solve(A,B)
    V = np.transpose(V.reshape((N+1,M+1)))
    return(V)


##KOLCHI:











fig = plt.figure(figsize = (11,7), dpi=100)
plt.contourf(X, Y, np.sqrt((champsU(M,N))**2+(champsV(M,N))**2),levels=50)
clb=plt.colorbar()
clb.ax.set_title('Vitesse U (en m/s')
#ax = sns.heatmap(champsU(M,N))
plt.xlabel('Longueur en(mm)')
plt.ylabel('Largeur en (mm)')
plt.title("Ecoulement dans une bronche deformée de génération 4")
plt.quiver(X, Y, champsU(M,N), champsV(M,N))
#ax = sns.heatmap(dchampsP2(M,N))
plt.show()