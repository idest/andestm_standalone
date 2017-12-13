import numpy as np
np.set_printoptions(threshold=np.nan)
import scipy as sp
from scipy.interpolate import RectBivariateSpline, interp2d
from termomecanico import CS, TM

x_axis = CS.get_x_axis()
y_axis = CS.get_y_axis()[::-1]

datos_q = np.loadtxt('datos_Q/QsObs.txt')
datos_q = datos_q[datos_q[:,0] > -80.]
datos_q = datos_q[datos_q[:,0] < -60.]
datos_q = datos_q[datos_q[:,1] > -45.]
datos_q = datos_q[datos_q[:,1] < -10.]

datos_q_x = datos_q[:,0].T
datos_q_y = datos_q[:,1].T
datos_q_shf = datos_q[:,2].T * 0.0001# TODO: multiplicar por escala

weight = datos_q[:,-1]
weight[weight == 4] = 0.2 # Marine Geophysics
weight[weight == 3] = 0.4 # Geochemical
weight[weight == 2] = 0.8 # Land Borehole
weight[weight == 1] = 1.0 # ODP Borehole

# Cargar modelo surface heat flow e interpolar en puntos de datos_Q
surface_heat_flow = TM.get_surface_heat_flow()[::-1]
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

# Calcular rmse
valid_weight = np.ones(valid_weight.shape)
rmse=np.sqrt(sum(((valid_datos_q_shf - valid_interp_shf)**2)*valid_weight)
             /sum(valid_weight))

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
valid_x = valid_shf_indexes[0]
valid_y = valid_shf_indexes[1]
valid_shf = surface_heat_flow[valid_shf_indexes]
shf_interpolator = interp2d(valid_x, valid_y, valid_shf)
a = shf_interpolator(datos_q_x, datos_q_y)
"""
