import numpy as np
import pandas as pd
import multiprocessing as mp
import resource
import sys
import functools
import matplotlib
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import compress
from src.setup import input_setup, exec_setup, read_rheo
from termomecanico import termomecanico
from src.utils import makedir, makedir_from_filename, calc_deviation
from src.plot import (heatmap_map, base_map, boolean_map, diff_map,
    plot_eet_equivalent_vs_effective, multi_map)
from src.colormaps import jet_white_r, eet_tassara_07, eet_pg_07, get_elevation_diff_cmap, categorical_cmap
from src.compute import SpatialArray2D
#from src.datos_q import shf_data, shf_data_coords, shf_data_types, shf_data_error
from src.stats import evaluate_model
from src.datos_q import shf_data

# print memory usage
def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

def seismogenic_thickness():
    pass

def default_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    results[name] = results_function(t_input, m_input, name)
    return results

def rheo_exploration(
        results_function, t_input=None, m_input=None, dir_name=None,
        uc_params=None, lc_params=None, lm_params=None, flm_params=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    if uc_params is None:
        uc_params = [m_input['Cs']]
    if lc_params is None:
        lc_params = [m_input['Ci']]
    if lm_params is None:
        lm_params = [m_input['Ml']]
    if flm_params is None:
        flm_params = [m_input['Serp']]
    rhe_data = read_rheo('data/Rhe_Param_ordenado_nuevos_indices.dat')
    results = {}
    i=0
    for flm_param in flm_params:
        m_input['Serp'] = flm_param
        if m_input['slm'] is True:
            flm_string = '__' + rhe_data[str(flm_param)]['name']
        else:
            flm_string = ''
        for lm_param in lm_params:
            m_input['Ml'] = lm_param
            for lc_param in lc_params:
                m_input['Ci'] = lc_param
                for uc_param in uc_params:
                    m_input['Cs'] = uc_param
                    name = (dir_name +
                        rhe_data[str(uc_param)]['name'] + '__'
                        + rhe_data[str(lc_param)]['name'] + '__'
                        + rhe_data[str(lm_param)]['name'] + flm_string)
                    results[name] = results_function(t_input, m_input, name)
                    mem()
                    i+=1
                    print(i)
                    #print(results)
    return results

def applied_stress_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    for s_max in np.linspace(100, 200, 5, endpoint=True):
        m_input['s_max'] = s_max
        name = dir_name + 's_max_{:.0f}'.format(s_max)
        results[name] = results_function(t_input, m_input, name)
    return results

def initial_thermal_vars(t_input):
    #k:3.5, h:2.2e-6
    t_input['k_cs'] = 3.0
    t_input['k_ci'] = 3.0
    t_input['k_ml'] = 3.0
    t_input['H_cs'] = 1.8e-6
    t_input['H_ci'] = 1.8e-6
    t_input['H_ml'] = 1.8e-6
    t_input['delta_icd'] = False
    t_input['t_lat'] = False
    t_input['delta'] = 10
    t_input['t'] = 39.13

def thermal_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    ###### Modelo Inicial
    initial_thermal_vars(t_input)
    name = dir_name + (
        'h_' + '{:.1f}'.format(t_input['H_cs']) +
        '__k_' + '{:.1f}'.format(t_input['k_cs']) +
        '__delta_' + '{:.1f}'.format(t_input['delta']) +
        '__t_' + '{:.2f}'.format(t_input['t']))
    #print(name)
    #print(t_input)
    results[name] = results_function(t_input, m_input, name)
    ###### t Variable
    initial_thermal_vars(t_input)
    name = dir_name + (
        'h_' + '{:.1f}'.format(t_input['H_cs']) +
        '__k_' + '{:.1f}'.format(t_input['k_cs']) +
        '__delta_' + '{:.1f}'.format(t_input['delta']) +
        '__t_var')
    t_input['t_lat'] = True
    #print(name)
    #print(t_input)
    results[name] = results_function(t_input, m_input, name)
    ###### K Variable
    initial_thermal_vars(t_input)
    name = dir_name + (
        'h_' + '{:.1f}'.format(t_input['H_cs']) +
        '__k_var' +
        '__delta_' + '{:.1f}'.format(t_input['delta']) +
        '__t_' + '{:.2f}'.format(t_input['t']))
    t_input['k_cs'] = 3.0
    t_input['k_ci'] = 2.5
    t_input['k_ml'] = 3.5
    #print(name)
    #print(t_input)
    results[name] = results_function(t_input, m_input, name)
    ###### Delta Variable
    initial_thermal_vars(t_input)
    name = dir_name + (
        'h_' + '{:.1f}'.format(t_input['H_cs']) +
        '__k_' + '{:.1f}'.format(t_input['k_cs']) +
        '__delta_icd' +
        '__t_' + '{:.2f}'.format(t_input['t']))
    t_input['delta_icd'] = True
    #print(name)
    #print(t_input)
    results[name] = results_function(t_input, m_input, name)
    ###### Modelo mas complejo
    initial_thermal_vars(t_input)
    name = dir_name + (
        'h_' + '{:.1f}'.format(t_input['H_cs']) +
        '__k_var' +
        '__delta_icd' +
        '__t_var')
    t_input['t_lat'] = True
    t_input['k_cs'] = 3.0
    t_input['k_ci'] = 2.5
    t_input['k_ml'] = 3.5
    t_input['delta_icd'] = True
    #print(name)
    #print(t_input)
    results[name] = results_function(t_input, m_input, name)
    return results

def thermal_H_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    initial_thermal_vars(t_input)
    for h in np.linspace(0, 5.e-6, 6):
        t_input['H_cs'] = h 
        t_input['H_ci'] = h 
        t_input['H_ml'] = h 
        name = dir_name + 'h_{:.2E}'.format(h)
        results[name] = results_function(t_input, m_input, name)
    return results

def thermal_K_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    initial_thermal_vars(t_input)
    for k in np.linspace(1,5,5):
        t_input['k_cs'] = k
        t_input['k_ci'] = k
        t_input['k_ml'] = k
        name = dir_name + 'k_{:.2E}'.format(k)
        results[name] = results_function(t_input, m_input, name)
    return results

def thermal_delta_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    initial_thermal_vars(t_input)
    for delta in np.linspace(5,15,6):
        t_input['delta'] = delta 
        name = dir_name + 'delta_{:.0f}'.format(delta)
        results[name] = results_function(t_input, m_input, name)
    return results

def thermal_hot_and_cold_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    initial_thermal_vars(t_input)
    name = dir_name + 'normal'
    results[name] = results_function(t_input, m_input, name)
    name = dir_name + 'hot'
    t_input['delta'] = 15
    t_input['H_cs'] = 4.e-6
    t_input['H_ci'] = 4.e-6
    t_input['H_ml'] = 4.e-6
    #t_input['k_cs'] = 3.
    #t_input['k_ci'] = 3.
    #t_input['k_ml'] = 3.
    results[name] = results_function(t_input, m_input, name)
    name = dir_name + 'cold'
    t_input['delta'] = 5
    t_input['H_cs'] = 1.e-7
    t_input['H_ci'] = 1.e-7
    t_input['H_ml'] = 1.e-7
    #t_input['k_cs'] = 5.
    #t_input['k_ci'] = 5.
    #t_input['k_ml'] = 5.
    results[name] = results_function(t_input, m_input, name)
    return results

def eet_equivalent_vs_effective_results(
        eet_effective_dict, save_dir='Teq_vs_Tef', plot=False):
    save_dir_maps = save_dir + 'Mapas/'
    save_dir_files = save_dir + 'Archivos/'
    def results(t_input, m_input, name):
        diffs = {}
        data = get_model_data(t_input, m_input, get_eet_data)
        filename_eet = save_dir_files + name + '.txt'
        makedir_from_filename(filename_eet)
        np.savetxt(filename_eet, data['eet'])
        #filename_eet_ft = save_dir_files + 'from_trench/' + name + '.txt'
        #makedir_from_filename(filename_eet_ft)
        #np.savetxt(filename_eet_ft, data['eet_from_trench'])
        uc = m_input['Cs']
        lc = m_input['Ci']
        uc_array = np.ones(data['eet'].cs.get_2D_shape()) * uc
        lc_array = np.ones(data['eet'].cs.get_2D_shape()) * lc
        uc_array = SpatialArray2D(uc_array, data['eet'].cs).mask_irrelevant_eet()
        lc_array = SpatialArray2D(lc_array, data['eet'].cs).mask_irrelevant_eet()
        for key, eet_effective in zip(
                eet_effective_dict.keys(), eet_effective_dict.values()):
            Tef = SpatialArray2D(
                np.loadtxt(eet_effective['file']),
                data['eet'].cs).mask_irrelevant_eet()
            eet_diff = data['eet'] - Tef
            makedir_from_filename(save_dir_files + eet_effective['dir'] + name)
            np.savetxt(
                save_dir_files + eet_effective['dir'] + name + '_diff.txt',
                eet_diff)
            sd = calc_deviation(data['eet'], Tef)
            if plot is True:
                diff_map(data['eet'], Tef, eet_diff, sd=sd,
                    colormap=jet_white_r,
                    colormap_diff = get_elevation_diff_cmap(100),
                    cbar_limits=[0,100], cbar_limits_diff=[-100,100],
                    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
                    title_1='Espesor Elástico Equivalente',
                    title_2='Espesor Elástico Efectivo',
                    title_3='Diff. (EET Eq. - EET ef.)',
                    labelpad=-48, labelpad_diff=-56,
                    filename=save_dir_maps + eet_effective['dir'] + name)
            diffs[key] = eet_diff
        return {'diffs': diffs, 'eet': data['eet'], 'uc': uc_array, 'lc': lc_array}
    return results

def get_thermal_data(model):
    return {'cs': model.cs,
            'geotherm': model.tm.get_geotherm(),
            'surface_heat_flow': model.tm.get_surface_heat_flow(
                format='positive milliwatts')}

def thermal_results(save_dir='Thermal', plot=False):
    save_dir_maps = save_dir + 'Mapas/'
    save_dir_files = save_dir + 'Archivos/'
    makedir(save_dir_files)
    def results(t_input, m_input, name):
        data = get_model_data(t_input, m_input, get_thermal_data)
        estimators, df = evaluate_model(
            data['surface_heat_flow'], shf_data, return_dataframe=True)
        ##dic = {
        ##    'rmse': estimators['rmse'], 'diff': data['diff'],
        ##    'e_prom': data['e_prom'], 'sigmas': data['sigmas']}
        ##stats = pd.DataFrame.from_dict(dic, orient='index')
        ##stats.to_csv(save_dir_files + 'stats' + name + '.txt')
        #x_grid, y_grid, z_grid = data['cs'].get_3D_grid(masked=False)
        #df = pd.DataFrame(
        #    {'lat': y_grid.flatten(),
        #    'lon': x_grid.flatten(),
        #    'depth': z_grid.flatten(),
        #    'Temp': data['geotherm'].flatten()})
        #df.to_csv(save_dir_files + name + '.csv', sep=',', na_rep='nan', index=False)
        if plot is True:
            #heatmap_map(
            #    data['surface_heat_flow'], colormap='afmhot',
            #    cbar_label='Heat Flow [W/m²]', title='Surface Heat Flow',
            #    filename = save_dir_maps + name, cbar_limits=[30,110])
            multi_map(
                shf=data['surface_heat_flow'],
                data=df['data_values'], diff=df['diffs'],
                data_coords=[df['lons'], df['lats']],
                data_types=df['data_types'], rmse=estimators['rmse'],
                mse=estimators['mse'], sigmas=estimators['sigmas'],
                filename= save_dir_maps + name)
        return {'geotherm': data['geotherm'],
                'surface_heat_flow': data['surface_heat_flow']}
    return results

def get_eet_wrong_data(model):
    return {'eet': model.mm.get_eet(),
            'eet_wrong': model.mm.get_eet_wrong(),
            'share_icd': model.mm.eet_calc_data['share_icd'],
            'share_moho': model.mm.eet_calc_data['share_moho']}

def eet_wrong_results(save_dir='EET_wrong', plot=False):
    save_dir_maps = save_dir + 'Mapas/'
    def results(t_input, m_input, name):
        data = get_model_data(t_input, m_input, get_eet_wrong_data)
        eet_diff = data['eet'] - data['eet_wrong']
        climit = max(1, abs(np.nanmin(eet_diff)))
        sd = calc_deviation(data['eet'], data['eet_wrong'])
        if plot is True:
            diff_map(data['eet'], data['eet_wrong'], eet_diff, sd=sd,
                colormap=jet_white_r,
                colormap_diff = get_elevation_diff_cmap(100),
                cbar_limits=[0,100], cbar_limits_diff=[-climit,climit],
                cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
                title_1='Espesor Elástico Equivalente Corregido',
                title_2='Espesor Elástico Equivalente Previo',
                title_3='Diff. (EET corregido - EET previo)',
                labelpad=-48, labelpad_diff=-56,
                filename=save_dir_maps + name)
            plot_coupled_zones(
                    data['share_moho'], data['share_icd'],
                    filename=save_dir_maps + name)
        return {'eet': data['eet'], 'eet_wrong': data['eet_wrong'],
                'eet_diff': eet_diff}
    return results

def get_eet_data(model):
    return {'eet': model.mm.get_eet(),
    #'eet_from_trench': model.mm.get_eet_from_trench(),
    'share_icd': model.mm.eet_calc_data['share_icd'],
    'share_moho': model.mm.eet_calc_data['share_moho']}

def eet_results(eet_effective_dict=None, save_dir='EET', plot=False):
    save_dir_maps = save_dir + 'Mapas/'
    save_dir_files = save_dir + 'Archivos/'
    def results(t_input, m_input, name):
        makedir_from_filename(save_dir_files + name)
        data = get_model_data(t_input, m_input, get_eet_data)
        np.savetxt(save_dir_files + name + '.txt', data['eet'])
        if plot is True:
            heatmap_map(
                data['eet'], colormap=jet_white_r, cbar_label='EET [km]',
                cbar_limits=[0,100],
                title='Espesor Elástico Equivalente',
                filename=save_dir_maps + name, labelpad=-45)
            plot_coupled_zones(
                data['share_moho'], data['share_icd'],
                filename=save_dir_maps + 'CZ/' + name)
            if eet_effective_dict:
                plot_eet_equivalent_vs_effective(
                    eet_effective_dict, data['eet'],
                    save_dir=save_dir_maps, name=name)
        return data['eet']
    return results

def get_ist_data(model):
    return {'ist': model.mm.get_integrated_strength()}

def ist_results(save_dir='IST', plot=False):
    save_dir_maps = save_dir + 'Mapas/'
    save_dir_files = save_dir + 'Archivos/'
    def results(t_input, m_input, name):
        makedir_from_filename(save_dir_files + name)
        data = get_model_data(t_input, m_input, get_ist_data)
        np.savetxt(save_dir_files + name + '.txt', data['ist'])
        if plot is True:
            heatmap_map(
                -data['ist'], colormap=jet_white_r,
                cbar_label='Integrated Strength [MPa]',
                title='Resistencia Integrada',
                filename=save_dir_maps + name, labelpad=-45)
        return data['ist']
    return results

#def plot_eet_equivalent_vs_effective(eet_effective_dict, eet_eq,
#        save_dir='EET', name='eet_diff'):
#    for eet_effective in eet_effective_dict.values():
#        eet_eff = SpatialArray2D(
#            np.loadtxt(eet_effective['file']), eet_eq.cs).mask_irrelevant_eet()
#        eet_diff = eet_eq - eet_eff
#        sd = calc_deviation(eet_eq, eet_eff)
#        diff_map(eet_eq, eet_eff, eet_diff, sd=sd,
#            colormap=jet_white_r,
#            colormap_diff = get_elevation_diff_cmap(100),
#            cbar_limits=[0,100], cbar_limits_diff=[-100,100],
#            cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
#            title_1='Espesor Elástico Equivalente',
#            title_2='Espesor Elástico Efectivo',
#            title_3='Diff. (EET eq. - EET ef.)',
#            labelpad=-48, labelpad_diff=-56,
#            filename=save_dir + eet_effective['dir'] + name)

def plot_coupled_zones(share_moho, share_icd, filename='cz'):
    boolean_map(share_moho, share_icd, title='Zonas Acopladas',
        filename=filename + '_cz')
    boolean_map(share_moho, title='Acoplamiento en Moho',
        filename=filename + '_cz_moho', cmap_idx=0)
    boolean_map(share_icd, title='Acoplamiento en ICD',
        filename=filename + '_cz_icd', cmap_idx=1)

def eet_deviation_from_prom(eets, names, save_dir):
    save_dir_deviations = save_dir + 'Mapas/Deviations/'
    eet_prom = sum(eets)/len(eets)
    heatmap_map(
        eet_prom, colormap=jet_white_r, cbar_label='EET [km]',
        cbar_limits=[0,100], title='Prom. Espesor Elástico Efectivo',
        filename=save_dir_deviations + 'EET_prom', labelpad=-45)
    standard_deviations = []
    for i in np.arange(len(eets)):
        eet = eets[i]
        name = names[i]
        #print(eet)
        #print(eet_prom)
        eet_sd = calc_deviation(eet, eet_prom)
        #print(eet_sd)
        standard_deviations.append(eet_sd)
        eet_diff = eet - eet_prom
        cmap = get_elevation_diff_cmap(100)
        heatmap_map(
            eet_diff, colormap=cmap, cbar_label='EET Dif. [km]',
            cbar_limits=[-20,20], labelpad=-45,
            title='D. Estándar EET Prom. - EET Modelo: {:.2f}'.format(eet_sd),
            filename=save_dir_deviations + name)
    #devs_o, names_o = (list(t) for t in zip(*sorted(
    #    zip(standard_deviations, names))))
    #deviations = dict(zip(names_o, devs_o))
    #print(deviations)

#def calc_deviation(eet1, eet2):
#    eet_sd = np.nansum(abs(eet1 - eet2))/eet1.size
#    eet_sd = float(eet_sd)
#    return eet_sd

def get_model_data(t_input, m_input, data_function):
    out_q = mp.Queue()
    def mp_termomecanico(t_input, m_input, data_function, queue):
        model = termomecanico(t_input, m_input)
        data = data_function(model)
        queue.put(data)
        return
    proc = mp.Process(
        target=mp_termomecanico, args=(t_input, m_input, data_function, out_q))
    proc.start()
    data = out_q.get()
    proc.join()
    return data 

def get_rheology_mosaic(
        eet_effective_dict, eets, eets_diffs,
        uc_rheos=None, lc_rheos=None, lm_rheos=None,
        Tef_keys=None, extended_plot=True, save_dir='Rheo_mosaics'):
    if uc_rheos is None and lc_rheos is NOne and lm_rheos is None:
        raise ValueError("At least one list of rheologies must be provided.")
    if Tef_keys is None:
        Tef_keys: ['Te_Tassara', 'Te_PG_400', 'Te_PG_600', 'Te_PG_800']
    rheos = {}
    if uc_rheos is not None:
        uc_array = np.dstack(tuple(uc_rheos))
        rheos['uc'] = {'array': uc_array}
    if lc_rheos is not None:
        lc_array = np.dstack(tuple(lc_rheos))
        rheos['lc'] = {'array': lc_array}
    if lm_rheos is not None:
        lm_array = np.dstack(tuple(lm_rheos))
        rheos['lm'] = {'array': lm_array}
    save_dir_maps = save_dir + 'Mapas/Mosaics/'
    save_dir_files = save_dir + 'Archivos/Mosaics/'
    eets_array = np.dstack(tuple(eets))
    for Tef_key in Tef_keys:
        eets_diff = [eet_diffs[Tef_key] for eet_diffs in eets_diffs]
        Tef = np.loadtxt(eet_effective_dict[Tef_key]['file'])
        Tef = SpatialArray2D(Tef, eets[0].cs).mask_irrelevant_eet()
        eets_diff_array = np.dstack(tuple(eets_diff))
        eets_diff_array_m = np.ma.array(
            eets_diff_array, mask=np.isnan(eets_diff_array))
        k = np.nanargmin(abs(eets_diff_array_m), axis=2)
        m,n = k.shape
        i,j = np.ogrid[:m,:n]
        for rheo in rheos.values():
            rheo_fit = rheo['array'][i, j, k]
            rheo_fit = SpatialArray2D(rheo_fit, eets[0].cs).mask_irrelevant_eet()
            rheo['fit'] = rheo_fit
        eet_fit = eets_array[i, j, k]
        eet_fit = SpatialArray2D(eet_fit, eets[0].cs).mask_irrelevant_eet()
        diffs_fit = eets_diff_array[i,j,k]
        diffs_fit = SpatialArray2D(diffs_fit, eets[0].cs).mask_irrelevant_eet()
        for rheo_key, rheo in zip(rheos.keys(), rheos.values()):
            filename_rheo = save_dir_files + Tef_key + '_' + rheo_key + '.txt'
            makedir_from_filename(filename_rheo)
            np.savetxt(filename_rheo, rheo['fit'])
        filename_eet = save_dir_files + Tef_key + '_eet' + '.txt'
        makedir_from_filename(filename_eet)
        np.savetxt(filename_eet, eet_fit)
        #nans = np.sum(diffs_array, axis=2)
        #uc_fit[np.isnan(nans)] = np.nan
        #eet_fit[np.isnan(nans)] = np.nan
        #diffs_fit[np.isnan(nans)] = np.nan
        #colors = categorical_cmap(1, len(uc_params)+2,
        #    desaturated_first=True)
        ### Plot
        titles_rheo = [rheo_key + ' mosaic' for rheo_key in rheos.keys()]
        filenames_rheo = [
            save_dir_maps+Tef_key+'_'+rheo_key
            for rheo_key in rheos.keys()]
        filename_diff = save_dir_maps+'/'+Tef_key
        if extended_plot is True:
            fig = plt.figure(figsize=(len(rheos)*4+12,6))
            gs = gridspec.GridSpec(1,len(rheos)+3)
            axs_rheo = [fig.add_subplot(gs[0,n]) for n in range(len(rheos))]
            axs_teq_vs_tef = [
                fig.add_subplot(gs[0,len(rheos)+n]) for n in np.arange(3)]
            filenames_rheo=[None]*len(rheos)
            filename_diff = None
        else:
            axs_rheo = [None]*len(rheos)
            axs_teq_vs_tef = None
        for i, rheo in enumerate(rheos.values()):
            u = np.unique(rheo['fit'][~np.isnan(rheo['fit'])])
            bounds = np.concatenate(
                ([min(u)-1], u[:-1]+np.diff(u)/2., [max(u)+1]))
            ncolors = len(bounds) - 1
            norm = mcolors.BoundaryNorm(bounds, ncolors)
            cmap = plt.cm.get_cmap('viridis', ncolors)
            cbar_ticks = bounds[:-1]+np.diff(bounds)/2.
            cbar_tick_labels = [
                rhe_data[str(int(round(val)))]['name'] for val in u]
            heatmap_map(
                rheo['fit'], colormap=cmap, cbar_limits=[None, None], norm=norm,
                cbar_ticks=cbar_ticks, cbar_tick_labels=cbar_tick_labels,
                ax=axs_rheo[i], filename=filenames_rheo[i],
                title=titles_rheo[i])
        eet_diff = eet_fit - Tef
        sd = calc_deviation(eet_fit, Tef)
        diff_map(
            eet_fit, Tef, eet_diff, sd=sd,
            colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
            cbar_limits=[0,100], cbar_limits_diff=[-100,100],
            cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
            title_1='Espesor Elástico Equivalente',
            title_2='Espesor Elástico Efectivo',
            title_3='Diff. (EET eq. - EET ef.)',
            axs=axs_teq_vs_tef, labelpad=-48, labelpad_diff=-56,
            filename=filename_diff)
        if extended_plot is True:
            plt.tight_layout()
            plt.show()

if __name__ == '__main__':
    exec_input, direTer, direMec = exec_setup()
    t_input, m_input = input_setup()
    save_dir = direMec + 'Exploracion/'
    makedir(save_dir)
    save_dir_thermal = direTer + 'Exploracion/'
    makedir(save_dir_thermal)
    #### Reologias
    lc_params = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    uc_params = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    lm_params = [23,24,25,26,27,28,29,30]
    flm_params = [23]
    rhe_data = read_rheo('data/Rhe_Param_ordenado_nuevos_indices.dat')
    #### EET efectivos
    eet_effective_dict = {
        'Te_Tassara': {
             'file': 'data/Te_invertido/Interpolados/Te_Tassara.txt',
             'dir': 'Tassara_07/',
             'colormap': eet_tassara_07},
        'Te_PG_400': {
             'file': 'data/Te_invertido/Interpolados/Te_PG_400.txt',
             'dir': 'Perez_Gussinye_07/400/',
             'colormap': eet_pg_07},
        'Te_PG_600': {
             'file': 'data/Te_invertido/Interpolados/Te_PG_600.txt',
             'dir': 'Perez_Gussinye_07/600/',
             'colormap': eet_pg_07},
        'Te_PG_800': {
             'file': 'data/Te_invertido/Interpolados/Te_PG_800.txt',
             'dir': 'Perez_Gussinye_07/800/',
             'colormap': eet_pg_07}}

    #### Thermal Casos ######################################################
    thermal_exploration(
            thermal_results(
                save_dir=save_dir_thermal + 'Termal_casos/', plot=True))

    #### EET Rheo ###########################################################
    #eets = rheo_exploration(
    #    eet_results(save_dir=save_dir + 'Rheo/EET/', plot=True))
    #eets_names = list(eets.keys())
    #eets_values = list(eets.values())
    #eet_deviation_from_prom(eets_values, eets_names, save_dir)
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            eet_results(
    #                save_dir=save_dir + 'EETs/', plot=True),
    #            lc_params=lc_params),
    #        uc_params=uc_params)

    #### EET Rheo -> Applied Stress ########################################
    #eets = rheo_exploration(
    #        functools.partial(applied_stress_exploration,
    #            eet_results(
    #                save_dir=save_dir + 'Applied_Stress/EET/', plot=True)))

    ### EET Rheo -> Thermal CASOS ##########################################
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            functools.partial(thermal_exploration,
    #                eet_results(
    #                    save_dir=save_dir + 'Termal_casos/EET/', plot=True)),
    #                uc_params=[1,10]),
    #        lc_params=lc_params)
    #eets_i = [
    #    list(nested_result.values())[0]
    #    for result in results.values()
    #    for nested_result in result.values()]
    #eets_k_var = [
    #    list(nested_result.values())[1]
    #    for result in results.values()
    #    for nested_result in result.values()]
    #eets_delta_icd = [
    #    list(nested_result.values())[2]
    #    for result in results.values()
    #    for nested_result in result.values()]
    #eet_i_prom = sum(eets_i)/len(eets_i)
    #eet_k_var_prom = sum(eets_k_var)/len(eets_k_var)
    #eet_delta_icd_prom = sum(eets_delta_icd)/len(eets_delta_icd)
    #sd_k_var = calc_deviation(eet_i_prom, eet_k_var_prom)
    #sd_delta_icd = calc_deviation(eet_i_prom, eet_delta_icd_prom)
    #save_dir_maps = save_dir + 'Termal_casos/EET/Mapas/'
    #save_dir_files = save_dir + 'Termal_casos/EET/Archivos/'
    #diff_map(
    #    eet_i_prom, eet_k_var_prom, eet_i_prom - eet_k_var_prom, sd=sd_k_var,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=None,
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Espesor Elástico Modelo Simple',
    #    title_2='Espesor Elástico (Kcs:3.0, Kci:2.5, Kml:3.5)',
    #    title_3='Diff. (EET simple - EET K variable)',
    #    labelpad=-48, labelpad_diff=-56,
    #    filename=save_dir_maps + 'diff_k_var')
    #makedir_from_filename(save_dir_files + 'prom_k_var')
    #np.savetxt(save_dir_files + 'prom_k_var', eet_k_var_prom)
    #diff_map(
    #    eet_i_prom, eet_delta_icd_prom, eet_i_prom - eet_delta_icd_prom,
    #    sd=sd_delta_icd,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=[-20,20],
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Espesor Elástico Modelo Simple',
    #    title_2='Espesor Elástico (delta = ICD)',
    #    title_3='Diff. (EET simple - EET delta = ICD)',
    #    labelpad=-48, labelpad_diff=-56,
    #    filename=save_dir_maps + 'diff_delta_icd')
    #makedir_from_filename(save_dir_files + 'prom_delta_icd')
    #np.savetxt(save_dir_files + 'prom_delta_icd', eet_delta_icd_prom)

    ### EET Rheo -> Thermal H ##############################################
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            functools.partial(thermal_H_exploration,
    #                eet_results(
    #                    save_dir=save_dir + 'Thermal_H/EET/', plot=True)),
    #                uc_params=[1,10]),
    #        lc_params=lc_params)
    #eets = []
    #eet_proms = []
    #medias_eet_prom = []
    #for i in range(6):
    #    eets_i = [
    #        list(nested_result.values())[i]
    #        for result in results.values()
    #        for nested_result in result.values()]
    #    eets.append(eets_i)
    #    eet_prom = sum(eets_i)/len(eets_i)
    #    eet_proms.append(eet_prom)
    #    medias_eet_prom.append(np.nansum(eet_prom)/eet_prom.size) 
    #sd = calc_deviation(eet_proms[-1], eet_proms[0])
    #save_dir_maps = save_dir + 'Thermal_H/EET/Mapas/'
    #save_dir_files = save_dir + 'Thermal_H/EET/Archivos/'
    #fig = plt.figure(figsize=(12,8))
    #gs = gridspec.GridSpec(4,3)
    #axs = [fig.add_subplot(gs[:3,n]) for n in range(3)]
    #diff_map(
    #    eet_proms[0], eet_proms[-1], eet_proms[-1] - eet_proms[0], sd=sd,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=None,
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1=r'Espesor Elástico (H: 0 W/$m^3$)',
    #    title_2=r'Espesor Elástico (H: 5.0*10^-6 W/$m^3$)',
    #    title_3='Diff. (EET H alto - EET H=0)',
    #    labelpad=-48, labelpad_diff=-56,
    #    axs=axs)#, filename=save_dir_maps + 'diff_h')
    #makedir_from_filename(save_dir_files + 'prom_h_0')
    #makedir_from_filename(save_dir_files + 'prom_h_alto')
    #np.savetxt(save_dir_files + 'prom_h_0', eet_proms[0])
    #np.savetxt(save_dir_files + 'prom_h_alto', eet_proms[-1])
    #ax = fig.add_subplot(gs[3,:]) 
    #index = np.linspace(0, 5.e-6, 6)
    #ax.plot(index, medias_eet_prom, '-r', linewidth=1.)
    #ax.bar(index, medias_eet_prom, alpha=.4, width=1.e-6, edgecolor='k')
    #diff = max(medias_eet_prom) - min(medias_eet_prom)
    #ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax.set_xlim(index[0]-(1.e-6)/2, index[-1]+(1.e-6)/2)
    #ax.set_ylim(min(medias_eet_prom)-0.2*diff, max(medias_eet_prom)+0.2*diff)
    #ax.set_ylabel('Media del Espesor Elástico')
    #ax.set_xlabel(r'H [W/$m^3$]')
    #ax.grid(True)
    #makedir_from_filename(save_dir_maps + 'diff_h')
    #plt.tight_layout()
    #plt.savefig(save_dir_maps + 'diff_h.png')
    #plt.close()

    #### EET Rheo -> Thermal K ##############################################
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            functools.partial(thermal_K_exploration,
    #                eet_results(
    #                    save_dir=save_dir + 'Thermal_K/EET/', plot=True)),
    #                uc_params=[1,10]),
    #        lc_params=lc_params)
    #eets = []
    #eet_proms = []
    #medias_eet_prom = []
    #for i in range(5):
    #    eets_i = [
    #        list(nested_result.values())[i]
    #        for result in results.values()
    #        for nested_result in result.values()]
    #    eets.append(eets_i)
    #    eet_prom = sum(eets_i)/len(eets_i)
    #    eet_proms.append(eet_prom)
    #    medias_eet_prom.append(np.nansum(eet_prom)/eet_prom.size) 
    #sd = calc_deviation(eet_proms[-1], eet_proms[0])
    #save_dir_maps = save_dir + 'Thermal_K/EET/Mapas/'
    #save_dir_files = save_dir + 'Thermal_K/EET/Archivos/'
    #fig = plt.figure(figsize=(12,8))
    #gs = gridspec.GridSpec(4,3)
    #axs = [fig.add_subplot(gs[:3,n]) for n in range(3)]
    #diff_map(
    #    eet_proms[0], eet_proms[-1], eet_proms[-1] - eet_proms[0], sd=sd,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=None,
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Espesor Elástico (K: 1.0 W/mK)',
    #    title_2='Espesor Elástico (K: 5.0 W/mK)',
    #    title_3='Diff. (EET K alto - EET K bajo)',
    #    labelpad=-48, labelpad_diff=-56,
    #    axs=axs)
    #makedir_from_filename(save_dir_files + 'prom_k_bajo')
    #makedir_from_filename(save_dir_files + 'prom_k_alto')
    #np.savetxt(save_dir_files + 'prom_k_bajo', eet_proms[0])
    #np.savetxt(save_dir_files + 'prom_k_alto', eet_proms[-1])
    #ax = fig.add_subplot(gs[3,:]) 
    #index = np.linspace(1, 5, 5)
    #ax.plot(index, medias_eet_prom, '-r', linewidth=1.)
    #ax.bar(index, medias_eet_prom, alpha=.4, width=1, edgecolor='k')
    #diff = max(medias_eet_prom) - min(medias_eet_prom)
    ##ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax.set_xlim(index[0]-0.5, index[-1]+0.5)
    #ax.set_ylim(min(medias_eet_prom)-0.2*diff, max(medias_eet_prom)+0.2*diff)
    #ax.set_ylabel('Media del Espesor Elástico')
    #ax.set_xlabel('K [W/mK]')
    #ax.grid(True)
    #makedir_from_filename(save_dir_maps + 'diff_k')
    #plt.tight_layout()
    #plt.savefig(save_dir_maps + 'diff_k.png')
    #plt.close()

    #### EET Rheo -> Thermal Delta ##########################################
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            functools.partial(thermal_delta_exploration,
    #                eet_results(
    #                    save_dir=save_dir + 'Thermal_delta/EET/', plot=True)),
    #                uc_params=[1,10]),
    #        lc_params=lc_params)
    #eets = []
    #eet_proms = []
    #medias_eet_prom = []
    #for i in range(6):
    #    eets_i = [
    #        list(nested_result.values())[i]
    #        for result in results.values()
    #        for nested_result in result.values()]
    #    eets.append(eets_i)
    #    eet_prom = sum(eets_i)/len(eets_i)
    #    eet_proms.append(eet_prom)
    #    medias_eet_prom.append(np.nansum(eet_prom)/eet_prom.size) 
    #sd = calc_deviation(eet_proms[-1], eet_proms[0])
    #save_dir_maps = save_dir + 'Thermal_delta/EET/Mapas/'
    #save_dir_files = save_dir + 'Thermal_delta/EET/Archivos/'
    #fig = plt.figure(figsize=(12,8))
    #gs = gridspec.GridSpec(4,3)
    #axs = [fig.add_subplot(gs[:3,n]) for n in range(3)]
    #diff_map(
    #    eet_proms[0], eet_proms[-1], eet_proms[-1] - eet_proms[0], sd=sd,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=None,
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Espesor Elástico (delta: 5 km.)',
    #    title_2='Espesor Elástico (delta: 15 km.)',
    #    title_3='Diff. (EET delta alto - EET delta bajo)',
    #    labelpad=-48, labelpad_diff=-56,
    #    axs=axs)
    #makedir_from_filename(save_dir_files + 'prom_delta_bajo')
    #makedir_from_filename(save_dir_files + 'prom_delta_alto')
    #np.savetxt(save_dir_files + 'prom_delta_bajo', eet_proms[0])
    #np.savetxt(save_dir_files + 'prom_delta_alto', eet_proms[-1])
    #ax = fig.add_subplot(gs[3,:]) 
    #index = np.linspace(5, 15, 6)
    #ax.plot(index, medias_eet_prom, '-r', linewidth=1.)
    #ax.bar(index, medias_eet_prom, alpha=.4, width=2, edgecolor='k')
    #diff = max(medias_eet_prom) - min(medias_eet_prom)
    ##ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax.set_xticks(index)
    #ax.set_xlim(index[0]-1, index[-1]+1)
    #ax.set_ylim(min(medias_eet_prom)-0.2*diff, max(medias_eet_prom)+0.2*diff)
    #ax.set_ylabel('Media del Espesor Elástico')
    #ax.set_xlabel('delta [km]')
    #ax.grid(True)
    #makedir_from_filename(save_dir_maps + 'diff_delta')
    #plt.tight_layout()
    #plt.savefig(save_dir_maps + 'diff_delta.png')
    #plt.close()

    #### EET Rheo -> Applied Stress #########################################
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            functools.partial(applied_stress_exploration,
    #                eet_results(
    #                    save_dir=save_dir + 'Applied_stress/EET/', plot=True)),
    #                uc_params=[1,10]),
    #        lc_params=lc_params)
    #eets = []
    #eet_proms = []
    #medias_eet_prom = []
    #for i in range(5):
    #    eets_i = [
    #        list(nested_result.values())[i]
    #        for result in results.values()
    #        for nested_result in result.values()]
    #    eets.append(eets_i)
    #    eet_prom = sum(eets_i)/len(eets_i)
    #    eet_proms.append(eet_prom)
    #    medias_eet_prom.append(np.nansum(eet_prom)/eet_prom.size) 
    #sd = calc_deviation(eet_proms[-1], eet_proms[0])
    #save_dir_maps = save_dir + 'Applied_stress/EET/Mapas/'
    #save_dir_files = save_dir + 'Applied_stress/EET/Archivos/'
    #fig = plt.figure(figsize=(12,8))
    #gs = gridspec.GridSpec(4,3)
    #axs = [fig.add_subplot(gs[:3,n]) for n in range(3)]
    #diff_map(
    #    eet_proms[0], eet_proms[-1], eet_proms[-1] - eet_proms[0], sd=sd,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=None,
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Espesor Elástico (Stress aplicado: 100 MPa)',
    #    title_2='Espesor Elástico (Stress aplicado: 200 MPa)',
    #    title_3='Diff. (EET stress alto - EET stress bajo)',
    #    labelpad=-48, labelpad_diff=-56,
    #    axs=axs)
    #makedir_from_filename(save_dir_files + 'prom_applied_stress_bajo')
    #makedir_from_filename(save_dir_files + 'prom_applied_stress_alto')
    #np.savetxt(save_dir_files + 'prom_applied_stress_bajo', eet_proms[0])
    #np.savetxt(save_dir_files + 'prom_applied_stress_alto', eet_proms[-1])
    #ax = fig.add_subplot(gs[3,:]) 
    #index = np.linspace(100, 200, 5)
    #ax.plot(index, medias_eet_prom, '-r', linewidth=1.)
    #ax.bar(index, medias_eet_prom, alpha=.4, width=25, edgecolor='k')
    #diff = max(medias_eet_prom) - min(medias_eet_prom)
    ##ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax.set_xticks(index)
    #ax.set_xlim(index[0]-12.5, index[-1]+12.5)
    #ax.set_ylim(min(medias_eet_prom)-0.2*diff, max(medias_eet_prom)+0.2*diff)
    #ax.set_ylabel('Media del Espesor Elástico')
    #ax.set_xlabel('Stress aplicado [MPa]')
    #ax.grid(True)
    #makedir_from_filename(save_dir_maps + 'diff_applied_stress')
    #plt.tight_layout()
    #plt.savefig(save_dir_maps + 'diff_applied_stress.png')
    #plt.close()

    ### EET Rheo -> Variabilidad corteza superior ##########################
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            eet_results(
    #                save_dir=save_dir + 'Variability_UC/EET/', plot=True),
    #            uc_params=[1,10]),
    #        lc_params=lc_params)
    #eets_quartzite = [
    #    list(result.values())[0]
    #    for result in results.values()]
    #eets_carrara_m = [
    #    list(result.values())[1]
    #    for result in results.values()]
    #eet_prom_quartzite = sum(eets_quartzite)/len(eets_quartzite)
    #eet_prom_carrara_m = sum(eets_carrara_m)/len(eets_carrara_m)
    #sd = calc_deviation(eet_prom_carrara_m, eet_prom_quartzite)
    #save_dir_maps = save_dir + 'Variability_UC/EET/Mapas/'
    #save_dir_files = save_dir + 'Variability_UC/EET/Archivos/'
    #diff_map(
    #    eet_prom_quartzite, eet_prom_carrara_m,
    #    eet_prom_carrara_m - eet_prom_quartzite, sd=sd,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=[-50,50],
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Espesor Elástico Cs: Quartzite',
    #    title_2='Espesor Elástico Cs: Carrara Marble',
    #    title_3='Diff. (EET Carrara Marble - EET Quartzite)',
    #    labelpad=-48, labelpad_diff=-56,
    #    filename=save_dir_maps + 'uc_var.png')
    #makedir_from_filename(save_dir_files + 'eet_prom_quartzite')
    #makedir_from_filename(save_dir_files + 'eet_prom_carrara_marble')
    #np.savetxt(save_dir_files + 'eet_prom_quartzite', eet_prom_quartzite)
    #np.savetxt(save_dir_files + 'eet_prom_carrara_marble', eet_prom_carrara_m)

    ### EET Rheo -> Variabilidad corteza inferior ##########################
    #results = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            eet_results(
    #                save_dir=save_dir + 'Variability_LC/EET/', plot=True),
    #            lc_params=[11,22]),
    #        uc_params=uc_params)
    #eets_qz_diorite = [
    #    list(result.values())[0]
    #    for result in results.values()]
    #eets_cpxnite_w = [
    #    list(result.values())[1]
    #    for result in results.values()]
    #eet_prom_qz_diorite = sum(eets_qz_diorite)/len(eets_qz_diorite)
    #eet_prom_cpxnite_w = sum(eets_cpxnite_w)/len(eets_cpxnite_w)
    #sd = calc_deviation(eet_prom_cpxnite_w, eet_prom_qz_diorite)
    #save_dir_maps = save_dir + 'Variability_LC/EET/Mapas/'
    #save_dir_files = save_dir + 'Variability_LC/EET/Archivos/'
    #diff_map(
    #    eet_prom_qz_diorite, eet_prom_cpxnite_w,
    #    eet_prom_cpxnite_w - eet_prom_qz_diorite, sd=sd,
    #    colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=[-50,50],
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Espesor Elástico Ci: Qz. Diorite',
    #    title_2='Espesor Elástico Ci: Clinopyroxenite Wet',
    #    title_3='Diff. (EET Cpxnite - EET Qz.D.)',
    #    labelpad=-48, labelpad_diff=-56,
    #    filename=save_dir_maps + 'lc_var.png')
    #makedir_from_filename(save_dir_files + 'eet_prom_qz_diorite')
    #makedir_from_filename(save_dir_files + 'eet_prom_cpxnite_w')
    #np.savetxt(save_dir_files + 'eet_prom_qz_diorite', eet_prom_qz_diorite)
    #np.savetxt(save_dir_files + 'eet_prom_cpxnite_w', eet_prom_cpxnite_w)


    ### EET Wrong ############################################################
    #eets = rheo_exploration(
    #            functools.partial(rheo_exploration,
    #                eet_wrong_results(
    #                    save_dir = save_dir + 'EET_wrong/', plot=True),
    #                lc_params=lc_params),
    #            uc_params=uc_params)
    #eets_right = [
    #    uc['eet'] for lc in list(eets.values()) for uc in list(lc.values())]
    #eets_wrong = [
    #    uc['eet_wrong'] for lc in list(eets.values())
    #    for uc in list(lc.values())]
    #eets_right_prom = sum(eets_right)/len(eets_right)
    #eets_wrong_prom = sum(eets_wrong)/len(eets_wrong)
    #eets_diff = eets_right_prom - eets_wrong_prom
    #sd = calc_deviation(eets_right_prom, eets_wrong_prom)
    #climit = max(1, abs(np.nanmin(eets_diff)))
    #diff_map(eets_right_prom, eets_wrong_prom, eets_diff, sd=sd,
    #    colormap=jet_white_r,
    #    colormap_diff = get_elevation_diff_cmap(100),
    #    cbar_limits=[0,100], cbar_limits_diff=[-climit,climit],
    #    cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
    #    title_1='Prom. Espesor Elástico Equivalente Corregido',
    #    title_2='Prom. Espesor Elástico Equivalente Previo',
    #    title_3='Diff. (EET corregido - EET previo)',
    #    labelpad=-48, labelpad_diff=-56,
    #    filename=save_dir + 'EET_wrong/Prom')

    ## EET eq. vs EET ef. / Mosaics #########################################
    #results = rheo_exploration(
    #            eet_equivalent_vs_effective_results(
    #                eet_effective_dict, save_dir = save_dir + 'Teq_vs_Tef/',
    #                plot=True),
    #            uc_params=uc_params)
    #eets = [result['eet'] for result in list(results.values())]
    #eets_diffs = [result['diffs'] for result in list(results.values())]
    #uc_rheos = [result['uc'] for result in list(results.values())]
    #results = rheo_exploration(
    #    functools.partial(rheo_exploration,
    #        eet_equivalent_vs_effective_results(
    #            eet_effective_dict, save_dir=save_dir + 'Teq_vs_Tef/lc_uc/',
    #            plot=True),
    #        uc_params=uc_params),
    #    lc_params=lc_params)
    #eets = [
    #    nested_result['eet']
    #    for result in results.values()
    #    for nested_result in result.values()]
    #eets_diffs = [
    #    nested_result['diffs']
    #    for result in results.values()
    #    for nested_result in result.values()]
    #uc_rheos = [
    #    nested_result['uc']
    #    for result in results.values()
    #    for nested_result in result.values()]
    #lc_rheos = [
    #    nested_result['lc']
    #    for result in results.values()
    #    for nested_result in result.values()]
    #get_rheology_mosaic(
    #    eet_effective_dict, eets, eets_diffs,
    #    uc_rheos=uc_rheos, lc_rheos=lc_rheos, lm_rheos=None,
    #    Tef_keys=['Te_Tassara', 'Te_PG_400', 'Te_PG_600', 'Te_PG_800'],
    #    save_dir=save_dir + 'Teq_vs_Tef/lc_uc/',
    #    extended_plot=False)
    #eets_prom = SpatialArray2D(
    #    sum(eets)/len(eets), eets[0].cs).mask_irrelevant_eet()
    #eets_prom_filename = save_dir + 'Teq_vs_Tef/lc_uc/Archivos/prom.txt'
    #makedir_from_filename(eets_prom_filename)
    #np.savetxt(eets_prom_filename, eets_prom)
    #plot_eet_equivalent_vs_effective(eet_effective_dict, eets_prom,
    #    save_dir=save_dir + 'Teq_vs_Tef/lc_uc/Mapas/', name='prom_diff')

    ## Integrated Strength ################################################
    #ist = rheo_exploration(
    #    ist_results(save_dir=save_dir + 'Rheo/IST/', plot=True))

    ### Hot and Cold Geotherm
    #results = thermal_hot_and_cold_exploration(
    #        thermal_results(save_dir = save_dir_thermal + 'Thermal/', plot=True))
    #normal_gt = results['normal']['geotherm']
    #hot_gt = results['hot']['geotherm']
    #cold_gt = results['cold']['geotherm']
    #x_grid, y_grid, z_grid = normal_gt.cs.get_3D_grid(masked=False)
    #df = pd.DataFrame(
    #    {'lat': y_grid.flatten(),
    #    'lon': x_grid.flatten(),
    #    'depth': z_grid.flatten(),
    #    'Temp1': normal_gt.flatten(),
    #    'Temp2': hot_gt.flatten(),
    #    'Temp3': cold_gt.flatten()})
    #model = termomecanico(*input_setup())
    #moho_hot = hot_gt.extract_surface(model.gm.get_moho()) 
    #moho_cold = cold_gt.extract_surface(model.gm.get_moho()) 
    #moho_diff = moho_hot - moho_cold
    #diff_map(moho_cold, moho_hot, moho_diff, sd=None,
    #    colormap='coolwarm',
    #    colormap_diff = 'coolwarm',
    #    cbar_limits=[0, 1300], cbar_limits_diff=None,
    #    cbar_label='Temperatura [ºC]', cbar_label_diff='Dif. Temperatura [ºC]',
    #    title_1='Temperatura Moho Modelo Frío',
    #    title_2='Temperatura Moho Modelo Caliente',
    #    title_3='Diff. (Modelo Caliente - Modelo frío)',
    #    labelpad=-48, labelpad_diff=-56),
    #    filename=save_dir_thermal + 'Thermal/dif_moho.png')
    #df.to_csv(save_dir_thermal + 'Thermal/Table.txt', sep=' ', na_rep='nan', index=False)
