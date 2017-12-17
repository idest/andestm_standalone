import numpy as np
np.set_printoptions(threshold=np.nan)
import scipy as sp
from scipy.interpolate import griddata, RectBivariateSpline
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

def func(x, y):
    return np.cos(x) + np.sin(y)

X = np.linspace(10,20,num=200)
Y = np.linspace(-10,-20,num=100)
XX, YY = np.meshgrid(X,Y)

Z = func(XX,YY)
Z[:,0:100] = np.nan
Z[np.isnan(Z)] = -9999

plt.pcolor(X,Y,Z,vmin=-1.5,vmax=1.5)
plt.colorbar()

x_obs = np.random.rand(100)*(9)+10.5
y_obs = np.random.rand(100)*(9)-19.5
z_obs = np.random.rand(100)*(3)-1.5

interpolator = RectBivariateSpline(X, Y[::-1], Z[::-1].T)

z_interpolated = interpolator.ev(x_obs, y_obs)

plt.scatter(x_obs, y_obs, c=z_interpolated, vmin=-1.5, vmax=1.5)

print(max(z_interpolated))

plt.show()
