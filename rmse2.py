import numpy as np
from datos_q import shf_data_formatted as shf_data, shf_data_x, shf_data_y


def rmse(surface_heat_flow):
    shf_interpolated = interpolate_surface_heat_flow(surface_heat_flow,
                                                     shf_data_x,
                                                     shf_data_y)
    return

def calc_rmse(model, data):
    rmse = np.sqrt((model - data))**2).mean()

def calc_rmse_error(model, data, data_min, data_max, error):
    diff = model - data
    data_salida = np.zeros(len(diff))
    for i in range(len(diff)):
        if abs(diff[i]) < abs(error[i]*1e-3):
            data_salida[i] = model[i]
        elif diff[i] > 0:
            data_salida[i] = data_min[i]
        else:
            data_salida[i] = data_max[i]
    rmse = calc_rmse(model, data_salida)
    return rmse, data_salida
     
def interpolate_surface_heat_flow(surface_heat_flow, x, y):
    surface_heat_flow_masked = surface_heat_flow.copy()
    surface_heat_flow_masked[np.isnan(surface_heat_flow)] = 1.e-1000000000
    x_axis = surface_heat_flow.cs.get_x_axis()
    y_axis = surface_heat_flow.cs.get_y_axis()
    shf_interpolator = RectBivariateSpline(x_axis, y_axis[::-1],
                                           surface_heat_flow_masked[:,::-1])
    shf_interpolated = shf_interpolator.ev(shf_data_x, shf_data_y)
    return shf_interpolated
