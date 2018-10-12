import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
from utils import module_from_file
from scipy.io import netcdf

colormaps = module_from_file('src', 'src/colormaps.py')
eet_pg_07 = colormaps.eet_pg_07

# Grilla EET Perez-Gussinye et al., 2007
te_pg = netcdf.netcdf_file('data/Te_invertido/Perez-Gussinye/Te_400.grd',
    mmap=False)

# Axes
x_axis = te_pg.variables['x'][:]
y_axis = te_pg.variables['y'][::-1]

# Grids
eet_grid = te_pg.variables['z'][::-1, :]
xx, yy = np.meshgrid(x_axis, y_axis)

# Interpolation
y_axis_model = np.linspace(-10, -45,
    num=int(round(abs(-45+10)/0.2+1)),
    endpoint=True)
x_axis_model = np.linspace(-80, -60,
    num=int(round(abs(-60+80)/0.2+1)),
    endpoint=True)
xx_model, yy_model = np.meshgrid(x_axis_model,y_axis_model)
eet_grid_masked = eet_grid.copy()
eet_grid[np.isnan(eet_grid)] = 1.e-10000000000000
eet_interpolator = RectBivariateSpline(x_axis, y_axis[::-1], eet_grid[::-1].T)
eet_interpolated = eet_interpolator(x_axis_model, y_axis_model[::-1])
eet_interpolated = eet_interpolated[:,::-1].T

np.savetxt(
    'data/Te_invertido/Interpolados/Te_PG_07_400.txt', eet_interpolated.T)
eet_interpolated = np.loadtxt(
    'data/Te_invertido/Interpolados/Te_PG_07_400.txt').T

# Plot
#fig, ax = plt.subplots()
fig = plt.figure()
gs = gridspec.GridSpec(1,2)
# Original Map
ax1 = fig.add_subplot(gs[0,0])
map1 = Basemap(
    llcrnrlon=-90, llcrnrlat=-57,
    urcrnrlon=-30, urcrnrlat=15,
    epsg=4326, resolution = 'l', suppress_ticks=False)
map1.drawcoastlines(linewidth=0.5)
map1.drawlsmask(land_color='0.8', ocean_color='0.8', resolution='l')
heatmap1 = map1.pcolormesh(xx, yy, eet_grid, shading='gouraud',
    vmin=0, vmax=100, cmap=eet_pg_07)
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
    vmin=0, vmax=100, cmap=eet_pg_07)
# Colorbar
divider = make_axes_locatable(ax2)
cbar_ax2 = divider.append_axes('right', '5%', pad='12%')
plt.sca(ax2)
cbar = plt.colorbar(heatmap2, cax=cbar_ax2)

plt.tight_layout()
plt.show()
#plt.close()

