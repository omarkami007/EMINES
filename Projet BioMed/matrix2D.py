import numpy as np
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

def B(u0,M):
    Lx=10
    Ly=10
    X=np.linspace(0,Lx,M+1)
    Y=np.linspace(0,Ly,M+1)
    U0=[u0(y) for y in Y] #u0 profil de vitesse initial
    B=[0]*(M+1)*(M+1)
    for k in range(M+1):
        B[k]=U0[k]
    for k in range(M*(M-1)+1,M*M+1):
        B[k]=U0[k-(M*(M-1)+1)]
    return(np.array(B))

def couplage (i,j) :
    return (((i+j)*(i+j+1))/2+j)


def inverse (n) : #la bijection réciproque de la fct couplage
    p=math.floor((-1+math.sqrt(1+8*n))/2)
    return (p*(p+3)/2-n,n-p*(p+1)/2)

def reshap(U,M):
    u= [[0]*(M+1)]*(M)
    for i in range(M+1):
        for j in range(M+1):
            u[i][j]=U[i*(M+1)+j+1]
    return(u)

A=matriceM(3,3,0.09,0.09)
B=B(lambda y: y*(y-10),3)
U=np.linalg.solve(A,B)