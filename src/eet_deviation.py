import numpy as np
from dotmap import DotMap
import src.datos_q as dq
from src.plot import heatmap_map
from src.utils import makedir
from src.colormaps import jet_white_r, get_diff_cmap

def eet_deviation(eets, names, save_dir):
    save_dir_deviations = save_dir + 'Deviations/'
    makedir(save_dir_deviations)
    eet_prom = sum(eets)/len(eets)
    heatmap_map(
        eet_prom, colormap=jet_white_r, cbar_label='EET [km]',
        cbar_limits=[0,100], title='Prom. Effective Elastic Thickness',
        save_dir=save_dir_deviations, name='EET_prom')
    for i in np.arange(len(eets)):
        eet = eets[i]
        name = names[i]
        eet_sd = abs(eet - eet_prom).sum()/eet.size
        eet_diff = eet - eet_prom
        cmap = get_diff_cmap(100)
        heatmap_map(
            eet_diff, colormap=cmap, cbar_label='EET Diff [km]',
            cbar_limits=[-20,20], title='EET Diff. SD: {}'.format(eet_sd),
            save_dir=save_dir_deviations, name=name)
