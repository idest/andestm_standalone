import numpy as np
from dotmap import DotMap
import src.datos_q as dq
from src.plot import heatmap_map
from src.utils import makedir
from src.colormaps import jet_white_r, get_elevation_diff_cmap

def eet_deviation(eets, names, save_dir):
    save_dir_deviations = save_dir + 'Deviations/'
    makedir(save_dir_deviations)
    eet_prom = sum(eets)/len(eets)
    heatmap_map(
        eet_prom, colormap=jet_white_r, cbar_label='EET [km]',
        cbar_limits=[0,100], title='Prom. Espesor Elástico Efectivo',
        save_dir=save_dir_deviations, name='EET_prom', labelpad=-45)
    standard_deviations = []
    for i in np.arange(len(eets)):
        eet = eets[i]
        name = names[i]
        #print(eet)
        #print(eet_prom)
        eet_sd = np.nansum(abs(eet - eet_prom))/eet.size
        eet_sd = float(eet_sd)
        #print(eet_sd)
        standard_deviations.append(eet_sd)
        eet_diff = eet - eet_prom
        cmap = get_elevation_diff_cmap(100)
        heatmap_map(
            eet_diff, colormap=cmap, cbar_label='EET Dif. [km]',
            cbar_limits=[-20,20], labelpad=-45,
            title='D. Estándar EET Prom. - EET Modelo: {:.2f}'.format(eet_sd),
            save_dir=save_dir_deviations, name=name)
    #devs_o, names_o = (list(t) for t in zip(*sorted(
    #    zip(standard_deviations, names))))
    #deviations = dict(zip(names_o, devs_o))
    #print(deviations)
