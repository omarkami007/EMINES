import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math
import scipy
from scipy import sparse as sp
from scipy.sparse import dia_matrix

np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)


data = np.array([[1, 2, 3, 4],[0, 1, 2, 0]])
offsets = np.array([1, 0])
print(dia_matrix((data, offsets), shape=(4, 4)).toarray())


def matriceA(M,N,a):
    Lx=2
    Ly=N*Lx/M
    a=(M/Lx)**2
    c=-4*a
    eta=10*(-8)
    b=M/(eta*Lx)
    for i in range()