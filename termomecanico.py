import numpy as np
import pandas as pd
from glob import glob
from src.setup import data_setup, input_setup, exec_setup
from src.compute import compute, SpatialArray2D
from src.rmse import rmse
from src.plot import (thermal_latitude_profile, mechanic_latitude_profile,
                   elastic_thickness_latitude_profile,
                   heatmap_map, data_scatter_map, diff_scatter_map,
                   multi_map, data_scatter_plot, earthquake_map,
                   plot_eet_equivalent_vs_effective)
from src.datos_q import shf_data, shf_data_coords, shf_data_types, shf_data_error
from src.utils import makedir
from src.colormaps import (jet_white_r, jet_white, get_elevation_diff_cmap,
    eet_tassara_07, eet_pg_07)


def termomecanico(t_input, m_input):
    gm_data, areas, trench_age, rhe_data, coast = data_setup()
    model = compute(gm_data, areas, trench_age, rhe_data, coast, t_input, m_input)
    shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
    model_rmse, ishf = rmse(shf, return_ishf=True)#, weigh_error=True)
    return model, model_rmse, ishf

if __name__ == '__main__':
    t_input, m_input = input_setup()
    exec_input, direTer, direMec = exec_setup()
    model, model_rmse, ishf = termomecanico(t_input,m_input)
    # Save
    files_dir_ter = direTer + 'Archivos/'
    makedir(files_dir_ter)
    np.savetxt(files_dir_ter + 'ishf_' + exec_input.temcaso + '.txt', ishf)
    #np.savetxt('sigmas_' + exec_input.temcaso + '.txt', model_rmse['sigmas'])

    #Earthquakes CSN
    if exec_input.eqs != 0:
        if exec_input.eqs == 1 or exec_input.eqs == 3:
            eqs = pd.read_excel("data/earthquakes/CSN_2000_2018_C_SB30_+8.8.xlsx",
                    sheet_name="Sheet1")
            eqs = eqs[(eqs['ISA'] == True) & (eqs['OSB'] == True)]
            #eqs = eqs[(eqs['ISA'] == True)]
            eqs['color'] = 'black'
            eqs_csn = eqs
        if exec_input.eqs == 2 or exec_input.eqs == 3:
            eqs = pd.read_excel("data/earthquakes/USGS_1900_2018_C_SB30.xlsx",
                    sheet_name="Sheet1")
            eqs['color'] = 'black'
            eqs = eqs[(eqs['ISA'] == True) & (eqs['OSB'] == True)]
            #eqs = eqs[(eqs['ISA'] == True)]
            eqs_usgs = eqs
        if exec_input.eqs == 3:
            eqs_csn['color'] = 'black'
            eqs_usgs['color'] = 'blue'
            eqs = pd.concat([eqs_usgs, eqs_csn], ignore_index=True)
    else:
        eqs = None

    #Maps
    if exec_input.xt1:
        maps_dir = direTer + 'Mapas/'
        shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
        diff = model_rmse.diff
        rmse = model_rmse.rmse
        e_prom = model_rmse.e_prom
        sigmas = model_rmse.sigmas
        moda = model_rmse.moda
        heatmap_map(shf, colormap='afmhot', cbar_label='Heat Flow [W/m²]',
                    title='Surface Heat Flow', filename=maps_dir + 'shf_map')
        data_scatter_map(
            shf_data, data_coords=shf_data_coords, data_types=shf_data_types,
            rmse=rmse, filename=maps_dir + 'data_scatter_map')
        data_scatter_map(
            ishf, data_coords=shf_data_coords, data_types=shf_data_types,
            rmse=rmse, filename=maps_dir + 'ishf_scatter_map')
        diff_scatter_map(
            diff, data_coords=shf_data_coords, data_types=shf_data_types,
            rmse=rmse, filename=maps_dir + 'diff_scatter_map',
            e_prom=e_prom, sigmas=sigmas, moda=moda)
        multi_map(
            shf=shf, data=shf_data, diff=diff, data_coords=shf_data_coords,
            data_types=shf_data_types, rmse=rmse,
            e_prom=e_prom, sigmas=sigmas, filename=maps_dir + 'multi_map')
        k_prom = model.tm.vars.k_prom.extract_surface(model.gm.get_topo()-1)
        heatmap_map(k_prom, filename=maps_dir + 'k_prom_map')
        h_prom = model.tm.vars.h_prom.extract_surface(model.gm.get_topo()-1)
        heatmap_map(h_prom, filename=maps_dir + 'h_prom_map')
        labslab = model.gm.get_slab_lab().mask_irrelevant()
        moho = model.gm.get_moho().mask_irrelevant()
        icd = model.gm.get_icd().mask_irrelevant()
        topo = model.gm.get_topo()
        heatmap_map(labslab, filename=maps_dir + 'labslab', colormap='viridis')
        heatmap_map(moho, filename=maps_dir + 'moho', colormap='viridis')
        heatmap_map(icd, filename=maps_dir + 'icd', colormap='viridis')
        heatmap_map(topo, filename=maps_dir + 'topo', colormap='viridis')

    if exec_input.xm1:
        maps_dir_mec = direMec + 'Mapas/'
        heatmap_map(
            model.mm.get_eet(), colormap=jet_white_r, cbar_label='EET [km]',
            cbar_limits=[0,100], title='Effective Elastic Thickness',
            filename=maps_dir_mec + 'eet', earthquakes=None, draw_land=False,
            labelpad=-45)
        integrated_strength_gpa = model.mm.get_integrated_strength()/-1000.
        heatmap_map(
            integrated_strength_gpa, colormap=jet_white_r,
            cbar_label='Integrated Strength [Gpa]', title='Integrated Strength',
            filename=maps_dir_mec + 'i_strength', labelpad=-45,
            cbar_limits=[0,200], earthquakes=None, draw_land=False)
        earthquake_map(eqs, title='Sismos', filename=maps_dir_mec + 'eqs')
        #eet_effective_dict = {
        #    'Te_Tassara': {
        #         'file': 'data/Te_invertido/Interpolados/Te_Tassara.txt',
        #         'dir': 'Tassara_07/',
        #         'colormap': eet_tassara_07},
        #    'Te_PG_400': {
        #         'file': 'data/Te_invertido/Interpolados/Te_PG_400.txt',
        #         'dir': 'Perez_Gussinye_07/400/',
        #         'colormap': eet_pg_07},
        #    'Te_PG_600': {
        #         'file': 'data/Te_invertido/Interpolados/Te_PG_600.txt',
        #         'dir': 'Perez_Gussinye_07/600/',
        #         'colormap': eet_pg_07},
        #    'Te_PG_800': {
        #         'file': 'data/Te_invertido/Interpolados/Te_PG_800.txt',
        #         'dir': 'Perez_Gussinye_07/800/',
        #         'colormap': eet_pg_07}}
        #plot_eet_equivalent_vs_effective(eet_effective_dict, model.mm.get_eet(),
        #        save_dir=maps_dir_mec+'Teq_vs_Tef/', name='diff')

    # Data and Models Scatter Plot
    if exec_input.xt3:
        scatter_dir = 'Output/scatter/'
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
                ishf_models, ishf_labels, filename=scatter_dir + 'scatter')

    #Latitude profiles
    sli_idx_array = model.gm.slab_lab_int_index
    sli_lon_array = model.cs.get_x_axis()[sli_idx_array]
    sli_depth_array = model.gm.slab_lab_int_depth
    #for lat in np.arange(exec_input.tmi, exec_input.tmx, -exec_input.tdelta):
    for idx, lat in enumerate(model.cs.get_y_axis()[:-1]):
        sli_lon = sli_lon_array[idx]
        sli_depth = sli_depth_array[idx]
        if exec_input.xt2:
            save_dir = direTer + 'Perfiles/'
            thermal_latitude_profile(model.tm, lat, filename=save_dir + 't',
                earthquakes=eqs, sli={'lon': sli_lon, 'depth': sli_depth})
        if exec_input.xm2:
            save_dir = direMec + 'Perfiles/'
            #eqs = eqs[eqs['mag'] >= 4.0]
            mechanic_latitude_profile(model.mm, lat, filename=save_dir + 'm',
                earthquakes=eqs, sli={'lon': sli_lon, 'depth': sli_depth})
            #elastic_thickness_latitude_profile(model.mm, lat,
            #    filename=save_dir + '/eet/eet', earthquakes=eqs,
            #    sli={'lon': sli_lon, 'depth': sli_depth})
