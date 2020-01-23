import numpy as np
from scipy.interpolate import RectBivariateSpline
#from sklearn.metrics import mean_absolute_error
from dotmap import DotMap
#import src.datos_q as dq
#from src.datos_q import shf_df
import pandas as pd
from src.setup import exec_setup

exec_input, direTer, direMec = exec_setup()
weigh_errors = True if exec_input.weigh == 1 else False

def evaluate_model(
    surface_heat_flow, shf_data, weigh_errors=weigh_errors,
    return_dataframe=False, save_dir=None):
    if save_dir is None:
        save_dir = 'Output/'
    # Surface Heat Flow Model Interpolation
    shf_interpolated = interpolate_surface_heat_flow(
        surface_heat_flow, shf_data['lons'], shf_data['lats'])
    #Output Dataframe
    shf_df = shf_data.assign(model_values=shf_interpolated)
    #diffs = shf_df['model_values'] - shf_df['data_values']
    shf_df = shf_df.assign(diffs=(shf_df['model_values'] - shf_df['data_values']))
    #RMSE
    if weigh_errors is True:
        meanerr = shf_df['data_errors'].mean()
        shf_df['data_errors'].fillna(meanerr, inplace=True)
        # shf_df = shf_df.dropna(subset=['data_errors'])
        rmse = calc_rmse_weighted(
            shf_df['model_values'], shf_df['data_values'], shf_df['data_errors'],
                    meanerr)
        mse, sigmas = sigma_weighted(
            shf_df['model_values'], shf_df['data_values'], shf_df['data_errors'],
                    meanerr)
    else:
        #shf_df = shf_df.drop(columns=['data_error'])
        rmse = calc_rmse(shf_df['model_values'], shf_df['data_values'])
        mse, sigmas = sigma(shf_df['model_values'], shf_df['data_values'])
    estimators = {'rmse': rmse, 'mse': mse, 'sigmas': sigmas}
    print_estimators_table(estimators, save_dir + 'estimadores')
    #Standard deviation
    return_value = estimators
    if return_dataframe:
        return_value = [estimators, shf_df]
    return return_value

def calc_rmse(model, data):
    diff = model - data
    rmse = np.sqrt((diff**2).mean())
    return rmse

def calc_rmse_weighted(model, data, data_error, meanerr):
    data_weight =  meanerr / data_error
    n_weight = data_weight/np.sum(data_weight)
    diff = model - data
    rmse = np.sqrt(sum(n_weight*(diff**2)))
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
    #mae = mean_absolute_error(data, shf_interpolated)
    mse = diff.mean()
    sigma = np.std(diff) #np.sqrt(((diff-diff.mean())**2).mean())
    n_1_sigma = diff.mean() - sigma
    p_1_sigma = diff.mean() + sigma
    n_2_sigma = diff.mean() - 2*sigma
    p_2_sigma = diff.mean() + 2*sigma
    sigmas = {
        'p_1_sigma': p_1_sigma, 'n_1_sigma': n_1_sigma,
        'p_2_sigma': p_2_sigma, 'n_2_sigma': n_2_sigma}
    sigmas = DotMap(sigmas)
    #moda = np.nanmax(abs(diff))
    return mse, sigmas

def sigma_weighted(shf_interpolated, data, data_error, meanerr):
    diff = shf_interpolated - data
    data_weight = meanerr / data_error
    n_weight = data_weight / np.sum(data_weight)
    mse = np.average(diff, weights = n_weight)
    V1 = sum(data_weight)
    V2 = sum(data_weight**2)
    sigma = np.sqrt(sum(data_weight*((diff-mse)**2))/(V1-(V2/V1)))
    n_1_sigma = mse - sigma
    p_1_sigma = mse + sigma
    n_2_sigma = mse - 2*sigma
    p_2_sigma = mse + 2*sigma
    sigmas = {
        'p_1_sigma': p_1_sigma, 'n_1_sigma': n_1_sigma,
        'p_2_sigma': p_2_sigma, 'n_2_sigma': n_2_sigma}
    sigmas = DotMap(sigmas)
    #moda = np.nanmax(abs(diff))
    return mse, sigmas

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
