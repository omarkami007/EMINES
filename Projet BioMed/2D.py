import numpy                       #here we load numpy
from matplotlib import pyplot
from matplotlib import cm #here we load matplotlib
import time, sys
from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots
import sys
import seaborn as sns
numpy.set_printoptions(threshold=sys.maxsize)



x=numpy.linspace(0,10,11)
y=numpy.linspace(0,10,11)
nx=ny=11
dx= dy = 1
X,Y = numpy.meshgrid(x,y)
u = numpy.zeros((ny,nx))
for i in range(nx):
    u[:,i]=y
u[0, :] = 0
u[-1, :] = 0

for j in range(1,ny-1):
    u[j,:]=0.5+(u[j+1,:]+u[j-1,:])/2

print(u)
def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm = 1
    pn = numpy.empty_like(p)

    while l1norm > l1norm_target:
        pn = p.copy()
        for j in range(1,ny-1):
            p[j,:]=0.5+(p[j+1,:]+p[j-1,:])/2

        p[0, :] = 0
        p[-1, :] = 0
        print(p)
        l1norm = (numpy.sum(numpy.abs(p[:]) - numpy.abs(pn[:])) /
                numpy.sum(numpy.abs(pn[:])))

    return p

u = laplace2d(u, y, dx, dy, 1e-4)
p = numpy.ones((ny,nx))

print(u)
fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X, Y, u, numpy.zeros((ny,nx)))
ax = sns.heatmap(u)

pyplot.show()
