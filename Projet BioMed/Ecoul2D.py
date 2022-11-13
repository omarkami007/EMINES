import numpy
import numpy as np
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
import seaborn as sns
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)


##Conditions aux bords:

##variable declarations

M,Lx,Ly = 10,2,1
nx = M+1
ny = M+1
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)




##plotting aids
x = numpy.linspace(0, Lx, nx)
y = numpy.linspace(0, Ly, ny)
X,Y = np.meshgrid(x,y)






##calcul de P



def laplace2d(dx, dy, l1norm_target):
    ##initial conditions
    p = numpy.zeros((ny, nx))  # create a XxY vector of 0's
    ##boundary conditions
    p[:, 0] = 0  # p = P1 @ x = 0
    p[:, -1] = 1  # p = P2 @ x = Lx
    p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
    p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = Ly


    l1norm = 1
    pn = numpy.empty_like(p)

    while l1norm > l1norm_target:
        pn = p.copy()
        p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) +
                         dx**2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) /
                        (2 * (dx**2 + dy**2)))

        p[:, 0] = 0  # p = P1 @ x = 0
        p[:, -1] = 1  # p = P2 @ x = Lx
        p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
        p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = Ly
        l1norm = (numpy.sum(numpy.abs(p[:]) - numpy.abs(pn[:])) /
                numpy.sum(numpy.abs(pn[:])))

    return p







##calcul de Ux:

##initial conditions
u = numpy.zeros((ny, nx))  # create a XxY vector of 0's

u[:, 0] = 1  # u = u1 @ x = 0
u[:, -1] = 1
u[0, :] = 0  # u = 0 @ y = 0
u[-1, :] = 0  # u = 0 @ y = Ly
u[:, -2] = u[:, 0] - u[:, 1] + u[:, -1]

def Upoisson2d(dx, dy, l1norm_target):
    ##initial conditions
    u = numpy.zeros((ny, nx))  # create a XxY vector of 0's

    u[:, 0] = 1  # u = u1 @ x = 0
    u[:, -1] = 1
    u[0, :] = 0  # u = 0 @ y = 0
    u[-1, :] = 0  # u = 0 @ y = Ly
    u[:, -2] = u[:, 0] - u[:, 1] + u[:, -1]


    l1norm = 1
    eta  = 1 #viscosité
    un = numpy.empty_like(u)
    b = np.zeros((ny,nx))
    P = laplace2d(dx, dy, 1e-4)
    b[1:-1, 1:-1]= (P[1:-1, 2:]-P[1:-1, 1:-1])/(eta*dx)
    while l1norm > l1norm_target:
        un = u.copy()

        u[1:-1,1:-1] = (((un[1:-1, 2:] + un[1:-1, :-2]) * dy**2 +
                    (un[2:, 1:-1] + un[:-2, 1:-1]) * dx**2 -
                    b[1:-1, 1:-1] * dx**2 * dy**2) /
                    (2 * (dx**2 + dy**2)))

        u[:, 0] = 1  # u = u1 @ x = 0
        u[:, -1] = 1
        u[0, :] = 0  # u = 0 @ y = 0
        u[-1, :] = 0  # u = 0 @ y = Ly
        u[:, -2] = u[:, 0] - u[:, 1] + u[:, -1]

        l1norm = (numpy.sum(numpy.abs(u[:]) - numpy.abs(un[:])) /
                numpy.sum(numpy.abs(un[:])))

    return u



##Calcul de Uy:



def Vpoisson2d( dx, dy, l1norm_target):
    ##initial conditions
    v = numpy.zeros((ny, nx))  # create a XxY vector of 0's




    ##boundary conditions
    v[:, 0] = 0  # v = v1 @ x = 0
    v[:, -1] = 0  # v = v2 @ x = Lx
    v[0, :] = v[1, :]  # dv/dy = 0 @ y = 0
    v[-1, :] = v[-2, :]  # dv/dy = 0 @ y = Ly
    l1norm = 1
    vn = numpy.empty_like(v)

    while l1norm > l1norm_target:
        vn = v.copy()
        v[1:-1, 1:-1] = ((dy**2 * (vn[1:-1, 2:] + vn[1:-1, 0:-2]) +
                         dx**2 * (vn[2:, 1:-1] + vn[0:-2, 1:-1])) /
                        (2 * (dx**2 + dy**2)))

        v[:, 0] = 0  # v = v1 @ x = 0
        v[:, -1] = 0  # v = v2 @ x = Lx
        v[0, :] = v[1, :]  # dv/dy = 0 @ y = 0
        v[-1, :] = v[-2, :]  # dv/dy = 0 @ y = Ly

        if numpy.sum(numpy.abs(vn[:])) == 0 : break

        l1norm = (numpy.sum(numpy.abs(v[:]) - numpy.abs(vn[:])) /
                numpy.sum(numpy.abs(vn[:])))

    return v


'''#A supprimer:
def plot2D(x, y, p):
    fig = pyplot.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = numpy.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    pyplot.show()'''




U = Upoisson2d(dx, dy, 1e-4)
V = Vpoisson2d(dx, dy, 1e-4)
fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X[::3, ::3], Y[::3, ::3], U[::3, ::3], V[::3, ::3])
pyplot.show()
