import numpy as np
from src.setup import data_setup, input_setup
from termomecanico import termomecanico
from src.datos_q_pandas import shf_data
from scipy.interpolate import RectBivariateSpline
from src.compute import SpatialArray2D
import pandas as pd
#from src.plot import heatmap_map

def interpolate_2D_array(array_2D, x, y):
    array_2D_masked = array_2D.copy()
    array_2D_masked[np.isnan(array_2D)] = 1.e-1000000000
    x_axis = array_2D.cs.get_x_axis()
    y_axis = array_2D.cs.get_y_axis()
    array_2D_interpolator = RectBivariateSpline(
            x_axis,
            y_axis[::-1],
            array_2D_masked[:,::-1]
            )
    array_2D_interpolated = array_2D_interpolator.ev(x, y)
    return array_2D_interpolated

def __calc_Ho_M1(shf, k, delta, Zm, Zb, Tb):
    Zm = -Zm*1.e3
    Zb = -Zb*1.e3
    shf = shf/1.e3
    Ho = (shf*Zb - Tb*k)*np.exp(Zm/delta)/(delta*(Zb*np.exp(Zm/delta) - Zb + Zm - delta*np.exp(Zm/delta) + delta))
    return Ho

def __calc_k_M1(shf, Ho, delta, Zm, Zb, Tb):
    Zm = -Zm*1.e3
    Zb = -Zb*1.e3
    shf = shf/1.e3
    k = (Ho*delta*(Zb - Zm - delta) + (-Ho*Zb*delta + Ho*delta**2 + shf*Zb)*np.exp(Zm/delta))*np.exp(-Zm/delta)/Tb
    return k

def __calc_Huc_M2(shf, Hlc, k, Zi, Zm, Zb, Tb):
    Zi = -Zi*1.e3
    Zm = -Zm*1.e3
    Zb = -Zb*1.e3
    shf = shf/1.e3
    Huc = (2*Hlc*Zb*Zi - 2*Hlc*Zb*Zm - Hlc*Zi**2 + Hlc*Zm**2 + 2*shf*Zb - 2*Tb*k)/(Zi*(2*Zb - Zi))
    return Huc

def __calc_k_M2(shf, Hlc, Huc, Zi, Zm, Zb, Tb):
    Zi = -Zi*1.e3
    Zm = -Zm*1.e3
    Zb = -Zb*1.e3
    shf = shf/1.e3
    k = (Hlc*Zb*Zi - Hlc*Zb*Zm - Hlc*Zi**2/2 + Hlc*Zm**2/2 - Huc*Zb*Zi + Huc*Zi**2/2 + shf*Zb)/Tb
    return k

def __calc_Huc_M2_c(shf, Hlc, k_uc, k_lcm, Zi, Zm, Zb, Tb):
    Zi = -Zi*1.e3
    Zm = -Zm*1.e3
    Zb = -Zb*1.e3
    shf = shf/1.e3
    Huc = (2*Hlc*Zb*Zi*k_uc - 2*Hlc*Zb*Zm*k_uc - Hlc*Zi**2*k_uc + Hlc*Zm**2*k_uc + 2*shf*Zb*k_uc + 2*shf*Zi*k_lcm - 2*shf*Zi*k_uc - 2*Tb*k_lcm*k_uc)/(Zi*(2*Zb*k_uc + Zi*k_lcm - 2*Zi*k_uc))
    return Huc

def __calc_k_uc_M2_c(shf, Huc, Hlc, k_lcm, Zi, Zm, Zb, Tb):
    Zi = -Zi*1.e3
    Zm = -Zm*1.e3
    Zb = -Zb*1.e3
    shf = shf/1.e3
    k_uc = Zi*k_lcm*(-Huc*Zi + 2*shf)/(-2*Hlc*Zb*Zi + 2*Hlc*Zb*Zm + Hlc*Zi**2 - Hlc*Zm**2 + 2*Huc*Zb*Zi - 2*Huc*Zi**2 - 2*shf*Zb + 2*shf*Zi + 2*Tb*k_lcm)
    return k_uc

Ho = 3e-6
Huc = 1.65e-6
Hlc = 4e-7
k = 2.0
k_uc = 3.0
k_lcm = 1.0
Tb = 1300
delta = 10.e3

_, m_input = input_setup()
t_input = {
        'delta_icd': False,
        't_lat': True,
        'k': 3.0,
        'k_cs': 3.2,
        'k_ci': 3.2,
        'k_ml': 3.2,
        'H_0': 3e-06,
        'H_cs': 1.65e-06,
        'H_ci': 4e-07,
        'H_ml': 0.0,
        'kappa': 1e-06,
        'Tp': 1375.0,
        'G': 0.0004,
        'V': 66000.0,
        'b': 1.0,
        'dip': 20.0,
        'D': 0.0015,
        'delta': 12,
        't': 39.13
        }
model = termomecanico(t_input, m_input)

geometries_interpolated = []
for geometry in model.gm.get_boundaries():
    geometry_interpolated = interpolate_2D_array(
            geometry,
            shf_data['lons'],
            shf_data['lats']
            )
    geometries_interpolated.append(geometry_interpolated)

Zi = geometries_interpolated[1]
Zm = geometries_interpolated[2]
Zb = geometries_interpolated[3]

q_geom = pd.read_csv('data/Q_GEOM.txt', sep='\s+', header=None)

shf_tassara = q_geom.iloc[:, 2]
Zi_tassara = q_geom.iloc[:, 6]
Zm_tassara = q_geom.iloc[:, 5]
Zb_tassara = q_geom.iloc[:, 4]

Tb = interpolate_2D_array(SpatialArray2D(model.tm.slab_lab_temp,model.cs), shf_data['lons'], shf_data['lats'])

#Ho_M1 = __calc_Ho_M1(shf_data['data_values'], k, delta, Zm, Zb, Tb)
k_M1 = __calc_k_M1(shf_data['data_values'], Ho, delta, Zm, Zb, Tb)
k_M1_tassara = __calc_k_M1(shf_tassara, Ho, delta, Zm_tassara, Zb_tassara, Tb)
#Huc_M2 = __calc_Huc_M2(shf_data['data_values'], Hlc, k, Zi, Zm, Zb, Tb)
#k_M2 = __calc_k_M2(shf_data['data_values'], Huc, Hlc, Zi, Zm, Zb, Tb)

#print(k_M1.to_string())

output = pd.DataFrame({
    'lons': q_geom.iloc[:, 0],
    'lats': q_geom.iloc[:, 1],
    'k_M1_tassara': k_M1_tassara,
    'k_M1_mio': k_M1,
    'Tb': Tb
    })

print(output.to_string())
