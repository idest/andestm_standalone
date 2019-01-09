import numpy as np
from scipy.interpolate import RectBivariateSpline
#from sklearn.metrics import mean_absolute_error
from dotmap import DotMap
#import src.datos_q as dq
#from src.datos_q import shf_df
import pandas as pd

def rmse(
    surface_heat_flow, shf_df, weigh_error=False, return_ishf=False,
    save_dir=None):
    if save_dir is None:
        save_dir = 'Output/'
    # Surface Heat Flow Model Interpolation
    shf_interpolated = interpolate_surface_heat_flow(
        surface_heat_flow, shf_df['lon'], shf_df['lat'])
    shf_df = shf_df.assign(model=shf_interpolated)
    diff = shf_df['model'] - shf_df['data']
    #RMSE
    if weigh_error is True:
        shf_df = shf_df.dropna(subset=['data_error'])
        rmse = calc_rmse_weighted(
            shf_df['model'], shf_df['data'], shf_df['data_error'])
        e_prom, sigmas = sigma_weighted(
            shf_df['model'], shf_df['data'], shf_df['data_error'])
    else:
        #shf_df = shf_df.drop(columns=['data_error'])
        rmse = calc_rmse(shf_df['model'], shf_df['data'])
        e_prom, sigmas = sigma(shf_df['model'], shf_df['data'])
    estimators = {'rmse': rmse, 'e_prom': e_prom, 'sigmas': sigmas}
    print_estimators_table(estimators, 'estimadores')
    dic = {
        'diff': diff,
        'rmse': rmse, 'e_prom': e_prom, 'sigmas': sigmas}
    #Standard deviation
    return_tuple = []
    return_tuple.append(DotMap(dic))
    if return_ishf:
        return_tuple.append(shf_df['model'])
    return return_tuple

def calc_rmse(model, data):
    diff = model - data
    rmse = np.sqrt((diff**2).mean())
    return rmse

def calc_rmse_weighted(model, data, data_error):
    data_weight = 1 / data_error
    diff = model - data
    rmse = np.sqrt(sum((data_weight/sum(data_weight))*(diff**2)))
    return rmse

def calc_rmse_weighted_aggressive(model, data, data_error):
    diff = model - data
    data_min = data + data_error
    data_max = data - data_error
    data_salida = np.zeros(len(diff))
    for i in range(len(diff)):
        if abs(diff[i]) < abs(data_error[i]):
            data_salida[i] = model[i]
        elif diff[i] > 0:
            data_salida[i] = data_max[i]
        else:
            data_salida[i] = data_min[i]
    rmse, diff2 = calc_rmse(model, data_salida)
    #np.savetxt('shf_data.txt', data)
    #np.savetxt('ishf.txt', model)
    #np.savetxt('shf_data_error.txt', data_salida)
    #np.savetxt('diff.txt', diff)
    return rmse, data_salida

def sigma(shf_interpolated, data):
    diff = shf_interpolated - data
    sigma = np.std(diff) #np.sqrt(((diff-diff.mean())**2).mean())
    n_1_sigma = diff.mean() - sigma
    p_1_sigma = diff.mean() + sigma
    n_2_sigma = diff.mean() - 2*sigma
    p_2_sigma = diff.mean() + 2*sigma
    #e_prom = mean_absolute_error(data, shf_interpolated)
    e_prom = diff.mean()
    sigmas = {
        'p_1_sigma': p_1_sigma, 'n_1_sigma': n_1_sigma,
        'p_2_sigma': p_2_sigma, 'n_2_sigma': n_2_sigma}
    sigmas = DotMap(sigmas)
    #moda = np.nanmax(abs(diff))
    return e_prom, sigmas

def sigma_weighted(shf_interpolated, data, data_error):
    diff = shf_interpolated - data
    data_weight = 1 / data_error
    sigma = np.sqrt(sum((data_weight/sum(data_weight))*(diff-diff.mean())**2))
    n_1_sigma = diff.mean() - sigma
    p_1_sigma = diff.mean() + sigma
    n_2_sigma = diff.mean() - 2*sigma
    p_2_sigma = diff.mean() + 2*sigma
    e_prom = sum((data_weight/sum(data_weight))*diff)
    sigmas = {
        'p_1_sigma': p_1_sigma, 'n_1_sigma': n_1_sigma,
        'p_2_sigma': p_2_sigma, 'n_2_sigma': n_2_sigma}
    sigmas = DotMap(sigmas)
    #moda = np.nanmax(abs(diff))
    return e_prom, sigmas

def interpolate_surface_heat_flow(surface_heat_flow, x, y):
    surface_heat_flow_masked = surface_heat_flow.copy()
    surface_heat_flow_masked[np.isnan(surface_heat_flow)] = 1.e-1000000000
    x_axis = surface_heat_flow.cs.get_x_axis()
    y_axis = surface_heat_flow.cs.get_y_axis()
    shf_interpolator = RectBivariateSpline(x_axis, y_axis[::-1],
                                           surface_heat_flow_masked[:,::-1])
    shf_interpolated = shf_interpolator.ev(x, y)
    return shf_interpolated

def print_estimators_table(estimators, filename):
    if filename == None:
        filename = 'Output/estimadores'
    df = pd.DataFrame.from_dict(estimators)
    writer = pd.ExcelWriter(filename + '.xlsx')
    df.to_excel(writer, 'Estimadores')
    writer.save()
