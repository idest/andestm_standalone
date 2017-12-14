import numpy as np
import scipy as sp
from scipy.interpolate import griddata, RectBivariateSpline
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt


def func(x, y):
    return np.cos(np.pi*x) + np.sin(np.pi*y)

X = np.linspace(0,1,num=100)
Y = np.linspace(0,1,num=200)
print(X)
print(Y)
XX, YY = np.meshgrid(X,Y)
Z = func(XX,YY).reshape(len(X), len(Y))

x_obs = np.random.rand(100)
y_obs = np.random.rand(100)
z_obs = np.random.rand(100)

plt.pcolor(X, Y, Z.T)

interpolator = RectBivariateSpline(X, Y, Z)

z_interpolated = interpolator.ev(x_obs, y_obs)

plt.scatter(x_obs, z_obs, s=4, c=z_interpolated)
plt.show()

"""
#Suppose we want to interpolate the 2-D function
def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

# on a grid in [0, 1]x[0, 1]
grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

#but we only know its values at 1000 data points:
points = np.random.rand(1000, 2)

values = func(points[:,0], points[:,1])

grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
"""
"""
plt.subplot(221)
plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()
"""
