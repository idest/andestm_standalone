import numpy as np
from setup import data_setup, input_setup, exec_setup
from compute import compute
from rmse2 import rmse
from plot2 import (thermal_latitude_profile, mechanic_latitude_profile,
                   shf_map, data_map, diff_map, multi_map)
from datos_q2 import shf_data, shf_data_coords, shf_data_types
from utils import makedir


def termomecanico(t_input, m_input):
    gm_data, areas, trench_age, rhe_data = data_setup()
    model = compute(gm_data, areas, trench_age, rhe_data, t_input, m_input)
    model_rmse, ishf = rmse(model.tm.get_surface_heat_flow(), weigh_error=False,
                            return_ishf=True)
    return model, model_rmse, ishf

if __name__ == '__main__':
    t_input, m_input = input_setup()
    exec_input, direTer, direMec = exec_setup()
    model, model_rmse, ishf = termomecanico(t_input,m_input)
    #Maps
    save_dir = direTer + 'Mapas/'
    makedir(save_dir)
    shf = model.tm.get_surface_heat_flow()
    diff = model_rmse.diff
    rmse = model_rmse.rmse
    e_prom = model_rmse.e_prom
    p_1_sigma = model_rmse.p_1_sigma
    n_1_sigma = model_rmse.n_1_sigma
    p_2_sigma = model_rmse.p_2_sigma
    n_2_sigma = model_rmse.n_2_sigma
    shf_map(shf, save_dir=save_dir, name='shf_map')
    data_map(shf_data, data_coords=shf_data_coords, data_types=shf_data_types,
             rmse=rmse, save_dir=save_dir, name='data_map')
    diff_map(diff, data_coords=shf_data_coords, data_types=shf_data_types,
             rmse=rmse, save_dir=save_dir, name='diff_map', e_prom=e_prom,
             n_1_sigma=n_1_sigma, p_1_sigma=p_1_sigma, n_2_sigma=n_2_sigma, p_2_sigma=p_2_sigma)
    multi_map(shf=shf, data=shf_data, diff=diff, data_coords=shf_data_coords,
              data_types=shf_data_types, save_dir=save_dir, rmse=rmse, e_prom=e_prom,
              n_1_sigma=n_1_sigma, p_1_sigma=p_1_sigma, n_2_sigma=n_2_sigma, p_2_sigma=p_2_sigma)
    #Latitude profiles
    exec_input.xt2 = False
    exec_input.xm2 = False
    for lat in np.arange(exec_input.tmi, exec_input.tmx, -exec_input.tdelta):
        if exec_input.xt2:
            save_dir = direTer + 'Perfiles_por_latitud/'
            makedir(savedir)
            thermal_latitude_profile(model.tm, lat, save_dir)
        if exec_input.xm2:
            save_dir = direMec + 'Perfiles_por_latitud/'
            makedir(savedir)
            mechanic_latitude_profile(model.mm, lat, save_dir)
