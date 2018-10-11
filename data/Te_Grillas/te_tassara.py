import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
from utils import module_from_file
colormaps = module_from_file('src', '../../src/colormaps.py')
eet_tassara_07 = colormaps.eet_tassara_07

# Grilla EET Tassara et al., 2007
te_tassara = np.loadtxt('Tassara/te.SA.xyz')
x_values = te_tassara[:,0]
y_values = te_tassara[:,1]
eet = te_tassara[:,2]

# Lons
lon_max= np.nanmax(x_values)
lon_min = np.nanmin(x_values)
lon_step = 1/12
# Lats
lat_max = np.nanmax(y_values)
lat_min = np.nanmin(y_values)
lat_step = 1/12

# Axes
x_axis = np.linspace(lon_min, lon_max,
    num=int(round(abs(lon_max-lon_min)/lon_step+1)),
    endpoint=True)
y_axis = np.linspace(lat_max, lat_min,
    num=int(round(abs(lat_min-lat_max)/lat_step+1)),
    endpoint=True)

# Grids
x_grid = x_values.reshape(len(y_axis), len(x_axis))
y_grid = y_values.reshape(len(y_axis), len(x_axis))
eet_grid = eet.reshape(len(y_axis), len(x_axis))
xx, yy = np.meshgrid(x_axis, y_axis)


# Interpolation
y_axis_model = np.linspace(-10, -45,
    num=int(round(abs(-45+10)/0.2+1)),
    endpoint=True)
x_axis_model = np.linspace(-80, -60,
    num=int(round(abs(-60+80)/0.2+1)),
    endpoint=True)
xx_model, yy_model = np.meshgrid(x_axis_model,y_axis_model)
eet_interpolator = RectBivariateSpline(x_axis, y_axis[::-1], eet_grid[::-1].T)
eet_interpolated = eet_interpolator(x_axis_model, y_axis_model[::-1])
eet_interpolated = eet_interpolated[:,::-1].T

np.savetxt('Te_Tassara_07.txt', eet_interpolated.T)
eet_interpolated = np.loadtxt('Te_Tassara_07.txt').T

# Plot
#fig, ax = plt.subplots()
fig = plt.figure()
gs = gridspec.GridSpec(1,2)
# Original Map
ax1 = fig.add_subplot(gs[0,0])
map1 = Basemap(
    llcrnrlon=-85, llcrnrlat=-58,
    urcrnrlon=-30, urcrnrlat=20,
    epsg=4326, resolution = 'l', suppress_ticks=False)
map1.drawcoastlines(linewidth=0.5)
map1.drawlsmask(land_color='0.8', ocean_color='0.8', resolution='l')
heatmap1 = map1.pcolormesh(xx, yy, eet_grid, shading='gouraud',
    vmin=0, vmax=100, cmap=eet_tassara_07)
rect = Rectangle((-80,-45), 20, 30, facecolor='none', edgecolor='red')
ax1.add_patch(rect)
# Colorbar
divider = make_axes_locatable(ax1)
cbar_ax1 = divider.append_axes('right', '5%', pad='12%')
plt.sca(ax1)
cbar = plt.colorbar(heatmap1, cax=cbar_ax1)
# Interpolated Map
ax2 = fig.add_subplot(gs[0,1])
map2 = Basemap(
    llcrnrlon=-80, llcrnrlat=-45,
    urcrnrlon=-60, urcrnrlat=-10,
    epsg=4326, resolution = 'l', suppress_ticks=False)
map2.drawcoastlines(linewidth=0.5)
map2.drawlsmask(land_color='0.8', ocean_color='0.8', resolution='l')
heatmap2 = map2.pcolormesh(xx_model, yy_model, eet_interpolated, shading='gouraud',
    vmin=0, vmax=100, cmap=eet_tassara_07)
# Colorbar
divider = make_axes_locatable(ax2)
cbar_ax2 = divider.append_axes('right', '5%', pad='12%')
plt.sca(ax2)
cbar = plt.colorbar(heatmap2, cax=cbar_ax2)

plt.tight_layout()
plt.show()
