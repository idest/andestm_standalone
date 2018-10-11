import numpy as np
import multiprocessing as mp
import resource
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from src.setup import input_setup, exec_setup, read_rheo
from termomecanico import termomecanico
from src.utils import makedir
from src.plot import heatmap_map, base_map, boolean_map
from src.colormaps import jet_white_r, eet_tassara_07, eet_pg_07
from src.eet_deviation import eet_deviation

# print memory usage
def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

def eet_map(eet, colormap=jet_white_r, save_dir=None, name='eet_map',
        image=None):
    fig = plt.figure(figsize=(9,6))
    gs = gridspec.GridSpec(1,2)
    # Axis 1: EET Model
    ax1 = fig.add_subplot(gs[0,0])
    map1 = base_map(topo=False)
    wr1 = heatmap_map(
        eet, colormap=colormap, cbar_label='EET [km]', cbar_limits=[0,100],
        title='Espesor Elástico Efectivo', map=map1, ax=ax1,
        return_width_ratio=True, labelpad=-45)
    # Axis 2: EET invertido
    ax2 = fig.add_subplot(gs[0,1])
    img = plt.imread(image)
    ax2.imshow(img)
    ax2.set_yticks([])
    ax2.set_xticks([])
    if save_dir:
        name = name + '.png'
        plt.savefig(save_dir + '%s' %(name))
    plt.close()


def eet_exploration(uc_params, lc_params, lm_params, save_dir, plot=False):
    t_input, m_input = input_setup()
    uc_var_name = 'Cs'
    lc_var_name = 'Ci'
    lm_var_name = 'Ml'
    save_dir_maps = save_dir + 'Mapas/'
    save_dir_files = save_dir + 'Archivos/'
    makedir(save_dir_maps)
    makedir(save_dir_files)
    rhe_data = read_rheo('data/Rhe_Param_ordenado_nuevos_indices.dat')
    eets = []
    names = []
    for lm_param in lm_params:
        m_input[lm_var_name] = lm_param
        for lc_param in lc_params:
            m_input[lc_var_name] = lc_param
            for uc_param in uc_params:
                m_input[uc_var_name] = uc_param
                eet, share_icd, share_moho = get_model_eet(t_input, m_input)
                eets.append(eet)
                name = ('eet' + '__'
                    + rhe_data[str(uc_param)]['name'] + '__'
                    + rhe_data[str(lc_param)]['name'] + '__'
                    + rhe_data[str(lm_param)]['name'])
                names.append(name)
                np.savetxt(save_dir_files + name + '.txt', eet)
                if plot is True:
                    save_dir_tassara_07 = save_dir_maps + 'Tassara_07/'
                    save_dir_pg_07 = save_dir_maps + 'Perez_Gussinye_07/'
                    save_dir_pg_07_400 = save_dir_pg_07 + '400/'
                    save_dir_pg_07_600 = save_dir_pg_07 + '600/'
                    save_dir_pg_07_800 = save_dir_pg_07 + '800/'
                    makedir(save_dir_tassara_07)
                    makedir(save_dir_pg_07)
                    makedir(save_dir_pg_07_400)
                    makedir(save_dir_pg_07_600)
                    makedir(save_dir_pg_07_800)
                    heatmap_map(
                        eet, colormap=jet_white_r, cbar_label='EET [km]',
                        cbar_limits=[0,100], title='Espesor Elástico Efectivo',
                        save_dir=save_dir_maps, name=name, labelpad=-45)
                    boolean_map(share_moho, share_icd, title='Zonas Acopladas',
                        save_dir=save_dir_maps + 'ZA', name=name + '_za')
                    boolean_map(share_moho, title='Share Moho', save_dir=save_dir_maps + 'ZA',
                        name=name + '_za_moho', cmap_idx=0)
                    boolean_map(share_icd, title='Share ICD', save_dir=save_dir_maps + 'ZA',
                        name=name + '_za_icd', cmap_idx=1)
                    #boolean_map(
                    #    share_icd, title='Attached ICD', save_dir=save_dir_maps,
                    #    name=name + '_share_icd'
                    #    )
                    eet_map(eet, colormap=eet_tassara_07,
                        save_dir=save_dir_tassara_07, name=name,
                        image='data/Te_Grillas/Imgs/Tassara_07.png')
                    eet_map(eet, colormap=eet_pg_07,
                        save_dir=save_dir_pg_07_400 + '', name=name,
                        image='data/Te_Grillas/Imgs/PG_07_400.png')
                    eet_map(eet, colormap=eet_pg_07,
                        save_dir=save_dir_pg_07_600 + '', name=name,
                        image='data/Te_Grillas/Imgs/PG_07_600.png')
                    eet_map(eet, colormap=eet_pg_07,
                        save_dir=save_dir_pg_07_800 + '', name=name,
                        image='data/Te_Grillas/Imgs/PG_07_800.png')
                mem()
    return eets, names

def get_model_eet(t_input, m_input):
    out_q = mp.Queue()
    def mp_termomecanico(t_input, m_input, queue):
        model, _, _ = termomecanico(t_input, m_input)
        data = {
            'eet': model.mm.get_eet(),
            'share_icd': model.mm.eet_calc_data['share_icd'],
            'share_moho': model.mm.eet_calc_data['share_moho']
        }
        #queue.put(model.mm.get_eet())
        queue.put(data)
        return
    proc = mp.Process(target=mp_termomecanico, args=(t_input, m_input, out_q))
    proc.start()
    #eet = out_q.get()
    data = out_q.get()
    eet = data['eet']
    share_icd = data['share_icd']
    share_moho = data['share_moho']
    proc.join()
    return eet, share_icd, share_moho

if __name__ == '__main__':
    exec_input, direTer, direMec = exec_setup()
    save_dir = direMec + 'EETs/'
    makedir(save_dir)
    #Corteza Inferior
    uc_params = [6]
    lc_params = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    lm_params = [30]
    #Corteza Superior
    #uc_params = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    #lc_params = [12]
    #lm_params = [30]
    #Manto Litosferico
    #uc_params=[6]
    #lc_params=[12]
    #lm_params=[23,24,25,26,27,28,29,30]

    eets, names = eet_exploration(uc_params, lc_params, lm_params, save_dir, True)
    eet_deviation(eets, names, save_dir)
