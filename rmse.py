import numpy as np
np.set_printoptions(threshold=np.nan)
import scipy as sp
from scipy.interpolate import RectBivariateSpline, interp2d, RegularGridInterpolator
from termomecanico import CS, TM
from mpl_toolkits.basemap import Basemap
import numpy.ma as ma
import matplotlib.pyplot as plt
import os

x_axis = CS.get_x_axis()
y_axis = CS.get_y_axis()[::-1]

datos_q = np.loadtxt('datos_Q/QsObs.txt')
datos_q = datos_q[datos_q[:,0] > -80.]
datos_q = datos_q[datos_q[:,0] < -60.]
datos_q = datos_q[datos_q[:,1] > -45.]
datos_q = datos_q[datos_q[:,1] < -10.]

datos_q_x = datos_q[:,0].T
datos_q_y = datos_q[:,1].T
datos_q_shf = -datos_q[:,2].T * 1.e-3# TODO: multiplicar por escala

weight = datos_q[:,-1]
weight[weight == 4] = 0.2 # Marine Geophysics
weight[weight == 3] = 0.4 # Geochemical
weight[weight == 2] = 0.8 # Land Borehole
weight[weight == 1] = 1.0 # ODP Borehole




###
#Interpolar modelo en puntos de datos Q con RegularGridInterpolator
###
surface_heat_flow = TM.get_surface_heat_flow()[::-1]
surface_heat_flow[np.isnan(surface_heat_flow)] = 0

interpolator = RegularGridInterpolator((x_axis, y_axis),surface_heat_flow)
q_x = datos_q_x[:, np.newaxis]
q_y = datos_q_y[:, np.newaxis]
obs_pts = np.append(q_x, q_y, axis=1)
print(obs_pts.shape)

interpolated_data = interpolator(obs_pts)

rmse_rgi=np.sqrt(sum(((interpolated_data - datos_q_shf)**2)*(weight/sum(weight))))
print("RMSE (RGI Interpolation):", rmse_rgi)



# Cargar modelo surface heat flow
surface_heat_flow = TM.get_surface_heat_flow()[::-1]

"""
valid_shf_indexes = np.where(np.isfinite(surface_heat_flow))
valid_x = x_axis[valid_shf_indexes[0]]
valid_y = y_axis[::-1][valid_shf_indexes[1]]
valid_shf = surface_heat_flow[valid_shf_indexes]
shf_interpolator_i2d = interp2d(valid_x, valid_y, valid_shf)
interpolated_shf_i2d = shf_interpolator_i2d(datos_q_x, datos_q_y)

# Calcular rmse
rmse_i2d=np.sqrt(sum(((interpolated_shf_i2d - datos_q_shf)**2)*(weight/sum(weight))))
print("RMSE (Interp2d Interpolation):", rmse_i2d)
"""




###
# Interpolar modelo en puntos de datos Q con Bivariate Spline Interpolation #
###
surface_heat_flow[np.isnan(surface_heat_flow)] = -9999 # nans = -9999
shf_interpolator_bsi = RectBivariateSpline(x_axis, y_axis, surface_heat_flow)
interpolated_shf_bsi = shf_interpolator_bsi.ev(datos_q_x, datos_q_y)
#print(np.nanmax(interpolated_shf_bsi))
#print(np.nanmin(interpolated_shf_bsi))

# Reemplazar valores interpolados negativos por nan y calcular sus indices
interpolated_shf_bsi[interpolated_shf_bsi < -1000] = np.nan
valid_interp_shf_idxs = np.where(np.isfinite(interpolated_shf_bsi))

# Seleccionar solo datos que no corresponden a un valor nan de shf interpolado
valid_interp_shf = interpolated_shf_bsi[valid_interp_shf_idxs]
valid_datos_q_shf = datos_q_shf[valid_interp_shf_idxs]
valid_datos_q_x = datos_q_x[valid_interp_shf_idxs]
valid_datos_q_y = datos_q_y[valid_interp_shf_idxs]
valid_weight = weight[valid_interp_shf_idxs]
#print(valid_interp_shf)
#print(valid_datos_q_shf)
#print(valid_weight)

print(valid_interp_shf)

# Calcular rmse
rmse_bsi=np.sqrt(sum(((valid_interp_shf - valid_datos_q_shf)**2)*(valid_weight/sum(valid_weight))))
print("RMSE (Bivariate Spline Interpolation):", rmse_bsi)

<<<<<<< Updated upstream
#Mapa que compara modelo con valores interpolados
map = Basemap(llcrnrlon= -79.8, llcrnrlat= -44.8, urcrnrlon= -58.0, urcrnrlat= -10.0, epsg= 4326, resolution = 'f')
#map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
map.drawparallels(np.arange(-90,90,3), labels=[1,0,0,0])
map.drawmeridians(np.arange(-180,180,4), labels=[0,0,1,0])
map.etopo()
map.drawcoastlines(linewidth=0.5)
#hacer grid y cargar los datos para la paleta de colores del mapa
mlon, mlat = map(valid_datos_q_x,valid_datos_q_y)
x = np.linspace(map.llcrnrx, map.urcrnrx, CS.get_x_axis().shape[0])
y = np.linspace(map.llcrnry, map.urcrnry, CS.get_y_axis().shape[0])
xx, yy = np.meshgrid(x, y)
interpolate_q = ma.masked_invalid(valid_interp_shf)
datam = ma.masked_invalid(TM.get_surface_heat_flow())
M = map.pcolormesh(xx, yy[::-1], datam.T, cmap='afmhot_r', shading= 'gouraud')
#Graficar datos de Q y barra de color
plt.clim(0,-0.2)
plot_q = map.scatter(mlon, mlat, c = interpolate_q, cmap = 'afmhot_r')
#plt.clim(0,-0.2)
cbar = plt.colorbar(M)
cbar.set_label('Flujo de Calor (W/m2)', rotation=90, labelpad=-60)
#nombre = "Mapa_Q_interpolate"
if not os.path.exists('Mapas'):
    os.makedirs('Mapas')
os.chdir('Mapas')
plt.savefig('Mapa_Q_interpolate')
os.chdir('../')
plt.close()



=======
"""
>>>>>>> Stashed changes
###
# Interpolar modelo en puntos de datos Q con Interp2d
###
"""
valid_shf_indexes = np.where(np.isfinite(surface_heat_flow))
valid_x = x_axis[valid_shf_indexes[0]]
valid_y = y_axis[::-1][valid_shf_indexes[1]]
valid_shf = surface_heat_flow[valid_shf_indexes]
shf_interpolator_i2d = interp2d(valid_x, valid_y, valid_shf)
interpolated_shf_i2d = shf_interpolator_i2d(datos_q_x, datos_q_y)

# Calcular rmse
rmse_i2d=np.sqrt(sum(((interpolated_shf_i2d - datos_q_shf)**2)*(weight/sum(weight))))
print("RMSE (Interp2d Interpolation):", rmse_i2d)
"""
"""
directory = os.fsencode('surface_heat_flow_dir')

for file in os.listdir(directory):
    filename = os.fsdecode(file)

    # Cargar modelo surface heat flow e interpolar en puntos de datos_Q
    surface_heat_flow = np.loadtxt(filename)[::-1]
    surface_heat_flow[np.isnan(surface_heat_flow)] = -9999 # nans = -9999
    shf_interpolator = RectBivariateSpline(x_axis, y_axis, surface_heat_flow)
    interpolated_shf = shf_interpolator.ev(datos_q_x, datos_q_y)

    # Reemplazar valores interpolados negativos por nan y calcular sus indices
    interpolated_shf[interpolated_shf < 0] = np.nan
    valid_interp_shf_idxs = np.where(np.isfinite(interpolated_shf))

    # Seleccionar solo datos que no corresponden a un valor nan de shf interpolado
    valid_interp_shf = interpolated_shf[valid_interp_shf_idxs]
    valid_datos_q_shf = datos_q_shf[valid_interp_shf_idxs]
    valid_weight = weight[valid_interp_shf_idxs]

    rmse=np.sqrt(sum(((valid_datos_q_shf - valid_interp_shf)**2)*valid_weight)
                 /sum(valid_weight))
"""

"""
surface_heat_flow = TM.get_surface_heat_flow()
valid_shf_indexes = np.where(np.isfinite(surface_heat_flow))
valid_x = x_axis[valid_shf_indexes[0]]
valid_y = y_axis[::-1][valid_shf_indexes[1]]
valid_shf = surface_heat_flow[valid_shf_indexes]
shf_interpolator = interp2d(valid_x, valid_y, valid_shf)
a = shf_interpolator(datos_q_x, datos_q_y)
"""
