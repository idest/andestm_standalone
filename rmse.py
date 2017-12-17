import numpy as np
np.set_printoptions(threshold=np.nan)
import scipy as sp
from scipy.interpolate import RectBivariateSpline, interp2d, RegularGridInterpolator
from termomecanico import CS, TM
from mpl_toolkits.basemap import Basemap
import numpy.ma as ma
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os

x_axis = CS.get_x_axis()
y_axis = CS.get_y_axis()

datos_q = np.loadtxt('datos_Q/QsObs.txt')
datos_q = datos_q[datos_q[:,0] > -80.]
datos_q = datos_q[datos_q[:,0] < -60.]
datos_q = datos_q[datos_q[:,1] > -45.]
datos_q = datos_q[datos_q[:,1] < -10.]

datos_q_x = datos_q[:,0]
datos_q_y = datos_q[:,1]
datos_q_shf = -datos_q[:,2] * 1.e-3# TODO: multiplicar por escala

weight = datos_q[:,-1]
weight[weight == 4] = 0.2 # Marine Geophysics
weight[weight == 3] = 0.4 # Geochemical
weight[weight == 2] = 0.8 # Land Borehole
weight[weight == 1] = 1.0 # ODP Borehole

variable_models_directory = 'Output/var_models/'

for cdir in next(os.walk(variable_models_directory))[1]:
    var = cdir
    full_cdir = variable_models_directory + cdir
    full_cdir_encoded = os.fsencode(full_cdir)
    var_values = []
    rmses = []
    for file in os.listdir(full_cdir_encoded):
        filename = os.fsdecode(file)
        if not filename.endswith('.png'):
            print(filename)
            ext = ".txt"
            var_value = filename[:filename.find(ext)].split('_')[-1]
            var_values.append(var_value)
            # Cargar datos modelo
            surface_heat_flow = np.loadtxt(full_cdir+'/'+filename)[::-1]
            # Interpolar modelo en puntos de datos Q con Bivariate Spline
            # Interpolation
            shf_min = np.nanmin(surface_heat_flow)
            shf_max = np.nanmax(surface_heat_flow)
            surface_heat_flow[np.isnan(surface_heat_flow)] = 0 # nans = -9999
            shf_interpolator_bsi = RectBivariateSpline(x_axis, y_axis[::-1],
                                                    surface_heat_flow[:,::-1])
            interpolated_shf_bsi = shf_interpolator_bsi.ev(datos_q_x,
                                                           datos_q_y)
            #print(np.nanmax(interpolated_shf_bsi))
            #print(np.nanmin(interpolated_shf_bsi))
            """
            # Descartar ciertos valores de interpolated_shf_bsi y de datos_q
            interpolated_shf_bsi[abs(interpolated_shf_bsi) > 1] = np.nan
            valid_interp_shf_idxs = np.where(np.isfinite(interpolated_shf_bsi))
            interpolated_shf_bsi = interpolated_shf_bsi[valid_interp_shf_idxs]
            datos_q_shf = datos_q_shf[valid_interp_shf_idxs]
            datos_q_x = datos_q_x[valid_interp_shf_idxs]
            datos_q_y = datos_q_y[valid_interp_shf_idxs]
            weight = weight[valid_interp_shf_idxs]
            """
            # Calcular rmse
            rmse=np.sqrt(sum(((interpolated_shf_bsi - datos_q_shf)**2)*weight)
                                   /sum(weight))
            #print("RMSE (Bivariate Spline Interpolation):", rmse_bsi)
            rmses.append(rmse)
    # Ordenar var_values de forma ascendente
    var_values, rmses = (list(n) for n in zip(*sorted(zip(var_values, rmses))))
    # Graficar resultados RMSE
    index = np.arange(len(var_values))
    plt.figure(figsize=(10,5))
    plt.bar(index, rmses)
    diff = max(rmses) - min(rmses)
    plt.ylim(min(rmses)-0.2*diff, max(rmses)+0.2*diff)
    plt.xticks(index, var_values)
    plt.title('Modeled Surface Heat Flow RMSE')
    plt.ylabel('RMSE')
    plt.xlabel('{} value'.format(var))
    plt.tight_layout()
    plt.savefig(full_cdir + '/RMSE.png')

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
# Interpolar datos modelo y calcular rmse
surface_heat_flow = np.loadtxt(TM.get_surface_heat_flow())[::-1]
# Interpolar modelo en puntos de datos Q con Bivariate Spline
# Interpolation
surface_heat_flow = TM.get_surface_heat_flow()
shf_min = np.nanmin(surface_heat_flow)
shf_max = np.nanmax(surface_heat_flow)
surface_heat_flow[np.isnan(surface_heat_flow)] = 0 # nans = -9999
shf_interpolator_bsi = RectBivariateSpline(x_axis, y_axis[::-1],
                                           surface_heat_flow[:,::-1])
interpolated_shf_bsi = shf_interpolator_bsi.ev(datos_q_x, datos_q_y)
#print(np.nanmax(interpolated_shf_bsi))
#print(np.nanmin(interpolated_shf_bsi))
# Descartar ciertos valores de interpolated_shf_bsi y de datos_q
interpolated_shf_bsi[abs(interpolated_shf_bsi) > 1] = np.nan
valid_interp_shf_idxs = np.where(np.isfinite(interpolated_shf_bsi))
interpolated_shf_bsi = interpolated_shf_bsi[valid_interp_shf_idxs]
datos_q_shf = datos_q_shf[valid_interp_shf_idxs]
datos_q_x = datos_q_x[valid_interp_shf_idxs]
datos_q_y = datos_q_y[valid_interp_shf_idxs]
weight = weight[valid_interp_shf_idxs]
# Calcular rmse
rmse=np.sqrt(sum(((interpolated_shf_bsi - datos_q_shf)**2)*weight)
                       /sum(weight))
#print("RMSE (Bivariate Spline Interpolation):", rmse_bsi)
rmses.append(rmse)
# Graficar valores modelados + valores interpolados en un plot de basemap
map = Basemap(projection='merc', resolution='l', epsg=4326,
              llcrnrlon=-80.0, llcrnrlat=-45.0,
              urcrnrlon=-60.0, urcrnrlat=-10.0)
#map.etopo()
#map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000,
                 verbose= True)
map.drawparallels(np.arange(-90,90,3), labels=[1,0,0,0])
map.drawmeridians(np.arange(-180,180,4), labels=[0,0,1,0])
map.drawcoastlines(linewidth=0.5)
mlon, mlat = map(datos_q_x, datos_q_y)
x1, y1 = map(x_axis[0], y_axis[-1])
x2, y2 = map(x_axis[-1], y_axis[0])
x = np.linspace(x1, x2, len(x_axis))
y = np.linspace(y1, y2, len(y_axis))
xx, yy = np.meshgrid(x,y)
datam = ma.masked_invalid(surface_heat_flow)
map.pcolormesh(xx, yy[::-1], datam.T, cmap='afmhot_r', shading='gouraud',
                   vmax=0, vmin=shf_min)
map.scatter(mlon, mlat, c = interpolated_shf_bsi, cmap='afmhot_r',
                     vmax=0, vmin=shf_min)
map.plot(mlon, mlat, 'bx', ms=1)
cbar = plt.colorbar()
cbar.set_label('Flujo de Calor (W/m2)', rotation=90, labelpad=-60)
plt.show()
#plt.savefig('interpolated_shf_map')
#if not os.path.exists('Mapas'):
#    os.makedirs('Mapas')
#os.chdir('Mapas')
#plt.savefig('Mapa_Q_interpolate')
#os.chdir('../')
#plt.close()
"""
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
###
#Interpolar modelo en puntos de datos Q con RegularGridInterpolator
###
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
surface_heat_flow = TM.get_surface_heat_flow()
valid_shf_indexes = np.where(np.isfinite(surface_heat_flow))
valid_x = x_axis[valid_shf_indexes[0]]
valid_y = y_axis[::-1][valid_shf_indexes[1]]
valid_shf = surface_heat_flow[valid_shf_indexes]
shf_interpolator = interp2d(valid_x, valid_y, valid_shf)
a = shf_interpolator(datos_q_x, datos_q_y)
"""
