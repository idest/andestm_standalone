import numpy as np

def calc_rmse(model, data):
    return np.sqrt(((model - data) ** 2).mean())

def calc_rmse_error(data,model,data_min,data_max,error):
    diff = model - data
    #print(abs(error[3]))
    data_salida = np.zeros(len(diff))
    for i in range(len(diff)):
        if abs(diff[i]) < abs(error[i]*1e-3):
            data_salida[i] = model[i]
        elif abs(diff[i]) > 0:
            data_salida[i] = data_min[i]
        else:
            data_salida[i] = data_max[i]
    rmse = calc_rmse(data_salida,model)
    return rmse, data_salida
