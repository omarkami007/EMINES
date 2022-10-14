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


