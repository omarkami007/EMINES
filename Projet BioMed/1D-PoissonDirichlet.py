import numpy as np
import matplotlib.pyplot as plt
import math
def matriceA_1D(N):
    a=[1,-2,1]+[0]*(N-2)
    A=[]
    for i in range(1,N):
        A.append(a)
        a=[0]+a
        a.pop(-1)
    A.append([0]*N+[1])
    return np.array(A)

def D1PoissonDirichlet(N,C0,CL,f,g):
    L=10
    Dx=L/N #pas
    X = np.linspace(0,L,N+1) #[0,L] discretisé
    T = np.linspace(0, L, 10*N) # [0,L] discretisé 10 fois plus
    G = [g(x) for x in X] #g([0,L]) discretisé
    F = [f(x) for x in X] #l'expression de f
    A = np.array([[1]+[0]*N] + list(matriceA_1D(N))) #la derivee deuxieme discretisé
    b=[C0]
    for i in range(1,len(G)-1):
        b+=[G[i]*Dx**2]
    b+=[CL]
    B=np.array(b)
    Y = np.linalg.solve(A, B)
    plt.plot(X,Y,label="approx")
    plt.plot(X,F,label = 'solution')
    plt.legend()
    plt.show()
