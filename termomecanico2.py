import numpy as np
from glob import glob
from setup import data_setup, input_setup, exec_setup
from compute import compute
from rmse2 import rmse
from plot2 import (thermal_latitude_profile, mechanic_latitude_profile,
                   shf_map, data_map, diff_map, multi_map, data_scatter_plot)
from datos_q2 import shf_data, shf_data_coords, shf_data_types, shf_data_error
from utils import makedir


def termomecanico(t_input, m_input):
    gm_data, areas, trench_age, rhe_data = data_setup()
    model = compute(gm_data, areas, trench_age, rhe_data, t_input, m_input)
    shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
    model_rmse, ishf = rmse(shf, weigh_error=False, return_ishf=True)
    return model, model_rmse, ishf

if __name__ == '__main__':
    t_input, m_input = input_setup()
    exec_input, direTer, direMec = exec_setup()
    model, model_rmse, ishf = termomecanico(t_input,m_input)

    # Save
    files_dir = direTer + 'Archivos/'
    makedir(files_dir)
    np.savetxt(files_dir + 'ishf_' + exec_input.temcaso + '.txt', ishf)

    #Maps
    if exec_input.xt3:
        maps_dir = direTer + 'Mapas/'
        makedir(maps_dir)
        shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
        diff = model_rmse.diff
        rmse = model_rmse.rmse
        e_prom = model_rmse.e_prom
        sigmas = model_rmse.sigmas
        shf_map(shf, save_dir=maps_dir, name='shf_map')
        data_map(
            shf_data, data_coords=shf_data_coords, data_types=shf_data_types,
            rmse=rmse, save_dir=maps_dir, name='data_map')
        diff_map(
            diff, data_coords=shf_data_coords, data_types=shf_data_types,
            rmse=rmse, save_dir=maps_dir, name='diff_map',
            e_prom=e_prom, sigmas=sigmas)
        multi_map(
            shf=shf, data=shf_data, diff=diff, data_coords=shf_data_coords,
            data_types=shf_data_types, save_dir=maps_dir, rmse=rmse,
            e_prom=e_prom, sigmas=sigmas, name='multi_map')

    # Data and Models Scatter Plot
    if exec_input.xt4:
        scatter_dir = 'Output/scatter/'
        makedir(scatter_dir)
        ishf_models_dir = scatter_dir + 'ishf_models/'
        makedir(ishf_models_dir)
        ishf_models_files = glob(ishf_models_dir + '*')
        ishf_labels = []
        for fn in ishf_models_files:
            fn = fn.split('/')[-1].split('_')[-1]
            label = 'Model ' + fn[:fn.find('.txt')]
            ishf_labels.append(label)
        if len(ishf_models_files) > 0:
            ishf_models = [np.loadtxt(f, ndmin=2) for f in ishf_models_files]
            ishf_models = np.concatenate(ishf_models, axis=1).T
            data_scatter_plot(
                shf_data, shf_data_error, shf_data_types,
                ishf_models, ishf_labels, save_dir=scatter_dir)

    #Latitude profiles
    for lat in np.arange(exec_input.tmi, exec_input.tmx, -exec_input.tdelta):
        if exec_input.xt2:
            save_dir = direTer + 'Perfiles/'
            makedir(save_dir)
            thermal_latitude_profile(model.tm, lat, save_dir, name='t')
        if exec_input.xm2:
            save_dir = direMec + 'Perfiles/'
            makedir(save_dir)
            mechanic_latitude_profile(model.mm, lat, save_dir, name='m')
