import numpy as np
import multiprocessing as mp
import resource
import sys
import functools
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import compress
from src.setup import input_setup, exec_setup, read_rheo
from termomecanico import termomecanico
from src.utils import makedir, makedir_from_filename
from src.plot import heatmap_map, base_map, boolean_map, diff_map
from src.colormaps import jet_white_r, eet_tassara_07, eet_pg_07, get_elevation_diff_cmap
from src.compute import SpatialArray2D

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
        flm_params = [m_input['Mla']]
    rhe_data = read_rheo('data/Rhe_Param_ordenado_nuevos_indices.dat')
    results = {}
    for flm_param in flm_params:
        m_input['Mla'] = flm_param
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
    print(results)
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
    for s_max in np.linspace(100, 200, 3, endpoint=True):
        m_input['s_max'] = s_max
        name = dir_name + 's_max_' + str(s_max)
        results[name] = results_function(t_input, m_input, name)

def thermal_exploration(
        results_function, t_input=None, m_input=None, dir_name=None):
    if t_input is None and m_input is None:
        t_input, m_input = input_setup()
    if dir_name is not None:
        dir_name = dir_name + '/'
    else:
        dir_name = ''
    results = {}
    def initial_vars():
        t_input['k_cs'] = 3.5
        t_input['k_ci'] = 3.5
        t_input['k_ml'] = 3.5
        t_input['H_cs'] = 2.2e-6
        t_input['H_ci'] = 2.2e-6
        t_input['H_ml'] = 2.2e-6
        t_input['delta_icd'] = False
        t_input['t_lat'] = False
        t_input['delta'] = 10
        t_input['t'] = 30
    ###### Modelo Inicial
    initial_vars()
    name = dir_name + 'h_2.2e-6__k_3.5__delta_10__t_30'
    print(name)
    print(t_input)
    results[name] = results_function(t_input, m_input, name)
    ###### K Variable
    t_input['k_cs'] = 3.0
    t_input['k_ci'] = 2.5
    t_input['k_ml'] = 3.5
    name = dir_name + 'h_2.2e-6__k_var__delta_10__t_30'
    print(name)
    print(t_input)
    results[name] = results_function(t_input, m_input, name)
    ###### Delta Variable
    initial_vars()
    t_input['delta_icd'] = True
    name = dir_name + 'h_2.2e-6__k_3.5__delta_var__t_30'
    print(name)
    print(t_input)
    results[name] = results_function(t_input, m_input, name)
    return results

def eet_equivalent_vs_effective_results(
        eet_effective_list, save_dir='Teq_vs_Tef', plot=False):
    save_dir_maps = save_dir + 'Mapas/'
    save_dir_files = save_dir + 'Archivos/'
    def results(t_input, m_input, name):
        diffs = {}
        data = get_model_data(t_input, m_input, get_eet_data)
        uc = m_input['Cs']
        lc = m_input['Ci']
        uc_array = np.ones(data['eet'].cs.get_2D_shape()) * uc
        lc_array = np.ones(data['eet'].cs.get_2D_shape()) * lc
        uc_array = SpatialArray2D(uc_array, data['eet'].cs).mask_irrelevant_eet()
        lc_array = SpatialArray2D(lc_array, data['eet'].cs).mask_irrelevant_eet()
        for eet_effective in eet_effective_list:
            Tef = SpatialArray2D(
                np.loadtxt(eet_effective['file']),
                data['eet'].cs).mask_irrelevant_eet()
            eet_diff = data['eet'] - Tef
            makedir_from_filename(save_dir_files + eet_effective['dir'] + name)
            np.savetxt(
                save_dir_files + eet_effective['dir'] + name + '.txt', eet_diff)
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
            diffs[eet_effective['id']] = eet_diff
        return {'diffs': diffs, 'uc': uc_array, 'lc': lc_array}
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
    'share_icd': model.mm.eet_calc_data['share_icd'],
    'share_moho': model.mm.eet_calc_data['share_moho']}

def eet_results(eet_effective_list=None, save_dir='EET', plot=False):
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
            if eet_effective_list:
                plot_eet_equivalent_vs_effective(
                    eet_effective_list, data['eet'],
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

def plot_eet_equivalent_vs_effective(eet_effective_list, eet_eq,
        save_dir='EET', name='eet_diff'):
    for eet_effective in eet_effective_list:
        eet_eff = SpatialArray2D(np.loadtxt(eet_effective['file']), eet_eq.cs)
        eet_diff = eet_eq - eet_eff
        sd = calc_deviation(eet_eq, eet_eff)
        diff_map(eet_eq, eet_eff, eet_diff, sd=sd,
            colormap=eet_effective['colormap'],
            colormap_diff = get_elevation_diff_cmap(100),
            cbar_limits=[0,100], cbar_limits_diff=[-100,100],
            cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
            title_1='Espesor Elástico Equivalente',
            title_2='Espesor Elástico Efectivo',
            title_3='Diff. (EET eq. - EET ef.)',
            labelpad=-48, labelpad_diff=-56,
            filename=save_dir + eet_effective['dir'] + name)

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

def calc_deviation(eet1, eet2):
    eet_sd = np.nansum(abs(eet1 - eet2))/eet1.size
    eet_sd = float(eet_sd)
    return eet_sd

def get_model_data(t_input, m_input, data_function):
    out_q = mp.Queue()
    def mp_termomecanico(t_input, m_input, data_function, queue):
        model, _, _ = termomecanico(t_input, m_input)
        data = data_function(model)
        queue.put(data)
        return
    proc = mp.Process(
        target=mp_termomecanico, args=(t_input, m_input, data_function, out_q))
    proc.start()
    data = out_q.get()
    proc.join()
    return data 

if __name__ == '__main__':
    exec_input, direTer, direMec = exec_setup()
    t_input, m_input = input_setup()
    save_dir = direMec + 'Exploration/'
    makedir(save_dir)
    #### Reologias
    lc_params = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    uc_params = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    lm_params = [23,24,25,26,27,28,29,30]
    flm_params = [23]
    #### EET efectivos
    eet_effective_list = [
        {'id':'Te_Tassara',
         'file': 'data/Te_invertido/Interpolados/Te_Tassara.txt',
         'dir': 'Tassara_07/',
         'colormap': eet_tassara_07},
        {'id': 'Te_PG_400',
         'file': 'data/Te_invertido/Interpolados/Te_PG_400.txt',
         'dir': 'Perez_Gussinye_07/400/',
         'colormap': eet_pg_07},
        {'id': 'Te_PG_600',
         'file': 'data/Te_invertido/Interpolados/Te_PG_600.txt',
         'dir': 'Perez_Gussinye_07/600/',
         'colormap': eet_pg_07},
        {'id': 'Te_PG_800',
         'file': 'data/Te_invertido/Interpolados/Te_PG_800.txt',
         'dir': 'Perez_Gussinye_07/800/',
         'colormap': eet_pg_07}]
    #### EET Rheo
    #eets = rheo_exploration(
    #    eet_results(save_dir=save_dir + 'Rheo/EET/', plot=True))
    #eets_names = list(eets.keys())
    #eets_values = list(eets.values())
    #eet_deviation_from_prom(eets_values, eets_names, save_dir)
    #### EET Rheo -> Applied Stress
    #eets = rheo_exploration(
    #        functools.partial(applied_stress_exploration,
    #            eet_results(
    #                save_dir=save_dir + 'Applied_Stress/EET/', plot=True)))
    #### EET Rheo -> Thermal
    #eets = rheo_exploration(
    #        functools.partial(thermal_exploration,
    #            eet_results(
    #                save_dir=save_dir + 'Thermal/EET/', plot=True)))
    ### EET Rheo LC -> EET Rheo UC
    #eets = rheo_exploration(
    #        functools.partial(rheo_exploration,
    #            eet_results(
    #                save_dir=save_dir + 'Rheo/EET/', plot=True),
    #            uc_params=uc_params),
    #       lc_params=lc_params)
    ### EET Wrong
    #eets = rheo_exploration(
    #            functools.partial(rheo_exploration,
    #                eet_wrong_results(
    #                    save_dir = save_dir + 'EET_wrong/', plot=True),
    #                uc_params=uc_params),
    #            lc_params=lc_params)
    #eets_right = [uc['eet'] for lc in list(eets.values()) for uc in list(lc.values())]
    #eets_wrong = [uc['eet_wrong'] for lc in list(eets.values()) for uc in list(lc.values())]
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
    ### EET eq. vs EET ef.
    #eets = rheo_exploration(
    #           eet_equivalent_vs_effective_results(Tef, save_dir = save_dir + 'Teq_vs_Tef/', plot=True))
    results = rheo_exploration(
                eet_equivalent_vs_effective_results(
                    eet_effective_list, save_dir = save_dir + 'Teq_vs_Tef/',
                    plot=True),
                uc_params=uc_params)
    uc_rheos = [result['uc'] for result in list(results.values())]
    #print(results.values())
    diffs_tassara_07 = [
        result['diffs']['Te_Tassara'] for result in list(results.values())]
    uc_array = np.dstack(tuple(uc_rheos))
    diffs_array = np.dstack(tuple(diffs_tassara_07))
    diffs_array_m = np.ma.array(diffs_array, mask=np.isnan(diffs_array))
    k = np.nanargmin(abs(diffs_array_m), axis=2)
    m,n = k.shape
    i,j = np.ogrid[:m,:n]
    uc_fit = uc_array[i,j,k]
    nans = np.sum(diffs_array, axis=2)
    uc_fit[np.isnan(nans)] = 0 
    print(np.unique(uc_fit[~np.isnan(uc_fit)]))
    fig = plt.figure()
    plt.imshow(uc_fit.T, vmin=0, vmax=10)
    plt.show()

    # Integrated Strength
    #ist = rheo_exploration(
    #    ist_results(save_dir=save_dir + 'Rheo/IST/', plot=True))

    """
    if len(flm_params) == 2:
        for i in np.arange(len(eets_flm)):
            flm_dir = save_dir + 'Mapas/FLM/'
            makedir(flm_dir)
            eet_diff_map(eets_no_flm[i], eets_flm[i], colormap=jet_white_r,
                save_dir=flm_dir, name=eets_names[i*2])
    """
