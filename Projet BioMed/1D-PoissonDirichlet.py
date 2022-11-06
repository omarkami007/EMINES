import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
import seaborn as sns
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)


def matriceA_1D(N):
    a=[1,-2,1]+[0]*(N-2)
    A=[]
    for i in range(1,N):
        A.append(a)
        a=[0]+a
        a.pop(-1)
    A.append([0]*N+[1])
    return np.array(A)

def D1PoissonDirichlet(N,L,C0,CL,g):
    Dx=L/N #pas
    X = np.linspace(0,L,N+1) #[0,L] discretisé
    G = [g(x) for x in X] #g([0,L]) discretisé
    A = np.array([[1]+[0]*N] + list(matriceA_1D(N))) #la derivee deuxieme discretisé
    b=[C0]
    for i in range(1,len(G)-1):
        b+=[G[i]*Dx**2]
    b+=[CL]
    B=np.array(b)
    F = np.linalg.solve(A, B)
    return(F)

###Solution 1ere etape: Uy=0 partout
def Etape1(g):
    M=10
    u = D1PoissonDirichlet(M,10,0,0,g)
    X,Y= np.meshgrid(np.linspace(0,10,M+1),np.linspace(0,10,M+1))

    U = np.zeros((len(u),len(u)))
    for j in range(len(U)):
        U[j,:]=u[j]
    print(U)
    #ax = sns.heatmap(U)
    fig = plt.figure(figsize = (11,7), dpi=100)
    plt.quiver(X, Y, U, np.zeros((M+1,M+1)))
    plt.show()
    return()

