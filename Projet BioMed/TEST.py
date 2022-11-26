import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
import seaborn as sns
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)



#Variables:
M = 20
Lx = 2
Ly = 1
a = (M/Lx)**2
b = (M/Ly)**2
x=np.linspace(0,Lx,M+1)
y=np.linspace(0,Ly,M+1)
X,Y = np.meshgrid(x,y)





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

    #Calcul du champs P:
        #Conditions auxbords:
P1=np.zeros(M+1)+10
P2=np.zeros(M+1)+5
P1[0],P2[0]=0,0
P1[-1],P2[-1]=0,0
#dp/dy = 0 a gauche et droite

def dPmatriceA(M,a,b):
    A = matriceA(M,a,b)
    for i in range(1,M):
        '''A[i*(M+1)][i*(M+1)+1]=-1
        #A[(i+1)*(M+1)-1][(i+1)*(M+1)-2]=-1'''
    return(A)

def dPmatriceB(M,a,b):
    #dP/dx=0 aux bords
    B = matriceB(M,Lx,Ly,P2,P1,np.zeros(M+1),np.zeros(M+1),np.zeros((M+1,M+1)))
    return(B)

def dchampsP(M):
    A = PmatriceA(M,a,b)
    B = PmatriceB(M,a,b)
    P = np.linalg.solve(A,B)
    P = P.reshape((M+1,M+1))
    return(np.transpose(P))