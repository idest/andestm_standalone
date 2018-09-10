import numpy as np
import multiprocessing as mp
import resource
import sys
from src.setup import input_setup, exec_setup, read_rheo
from termomecanico import termomecanico
from src.utils import makedir
from src.plot import heatmap_map
from src.meccolormap import jet_white_r

# print memory usage
def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

def eet_exploration(uc_params, lc_params, lm_params, save_dir, plot=False):
    t_input, m_input = input_setup()
    uc_var_name = 'Cs'
    lc_var_name = 'Ci'
    lm_var_name = 'Ml'
    save_dir_maps = save_dir + 'Mapas/'
    save_dir_files = save_dir + 'Archivos/'
    makedir(save_dir_maps)
    makedir(save_dir_files)
    rhe_data = read_rheo('data/Rhe_Param.dat')
    eets = []
    for lm_param in lm_params:
        m_input[lm_var_name] = lm_param
        for lc_param in lc_params:
            m_input[lc_var_name] = lc_param
            for uc_param in uc_params:
                m_input[uc_var_name] = uc_param
                eet = get_model_eet(t_input, m_input)
                eets.append(eet)    
                name = ('eet' + '__'
                    + rhe_data[str(uc_param)]['name'] + '__'
                    + rhe_data[str(lc_param)]['name'] + '__'
                    + rhe_data[str(lm_param)]['name']) 
                np.savetxt(save_dir_files + name + '.txt', eet)
                if plot is True:
                    heatmap_map(
                        eet, colormap=jet_white_r, cbar_label='EET [km]',
                        cbar_limits=[0,100], title='Effective Elastic Thickness',
                        save_dir=save_dir_maps, name=name)
                mem()
    return eets

def get_model_eet(t_input, m_input):
    out_q = mp.Queue()
    def mp_termomecanico(t_input, m_input, queue):
        model, _, _ = termomecanico(t_input, m_input)
        queue.put(model.mm.get_eet())
        return
    proc = mp.Process(target=mp_termomecanico, args=(t_input, m_input, out_q))
    proc.start()
    eet = out_q.get()
    proc.join()
    return eet

if __name__ == '__main__':
    exec_input, direTer, direMec = exec_setup()
    save_dir = direMec + 'EETs/'
    makedir(save_dir)
    uc_params = [9]
    lc_params = [11, 28, 12, 19, 13, 20, 14, 15, 16, 17, 18, 21]
    lm_params = [22]
    eets = eet_exploration(uc_params, lc_params, lm_params, save_dir, True)
    eet_deviation(eets, save_dir)
