import numpy as np
from scipy.interpolate import RectBivariateSpline
from dotmap import DotMap
import datos_q2 as dq

def rmse(surface_heat_flow, weigh_error=False, return_ishf=False):
    shf_interpolated = interpolate_surface_heat_flow(surface_heat_flow,
                                                     dq.shf_data_x,
                                                     dq.shf_data_y)
    if weigh_error is True:
        rmse, diff, data = calc_rmse_error(shf_interpolated, dq.shf_data,
                                           dq.shf_data_min, dq.shf_data_max,
                                           dq.error)
        dic = {'rmse': rmse, 'diff': diff, 'shf_data_weighted': data}
    else:
        rmse, diff = calc_rmse(shf_interpolated, dq.shf_data)
        dic = {'rmse': rmse, 'diff': diff}
    return_tuple = []
    return_tuple.append(DotMap(dic))
    if return_ishf:
        return_tuple.append(shf_interpolated)
    return return_tuple

def calc_rmse(model, data):
    diff = model-data
    rmse = np.sqrt((diff**2).mean())
    return rmse, diff

def calc_rmse_error(model, data, data_min, data_max, error):
    diff = model - data
    data_salida = np.zeros(len(diff))
    for i in range(len(diff)):
        if abs(diff[i]) < abs(error[i]):
            data_salida[i] = model[i]
        elif diff[i] > 0:
            data_salida[i] = data_min[i]
        else:
            data_salida[i] = data_max[i]
    rmse, diff = calc_rmse(model, data_salida)
    np.savetxt('shf_data.txt', data)
    np.savetxt('ishf.txt', model)
    np.savetxt('shf_data_error.txt', data_salida)
    np.savetxt('diff.txt', diff)
    return rmse, diff, data_salida

def interpolate_surface_heat_flow(surface_heat_flow, x, y):
    surface_heat_flow_masked = surface_heat_flow.copy()
    surface_heat_flow_masked[np.isnan(surface_heat_flow)] = 1.e-1000000000
    x_axis = surface_heat_flow.cs.get_x_axis()
    y_axis = surface_heat_flow.cs.get_y_axis()
    shf_interpolator = RectBivariateSpline(x_axis, y_axis[::-1],
                                           surface_heat_flow_masked[:,::-1])
    shf_interpolated = shf_interpolator.ev(x, y)
    return shf_interpolated
