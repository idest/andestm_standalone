import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from src.setup import input_setup, exec_setup, read_rheo
from termomecanico import termomecanico
from src.utils import makedir
from src.colormaps import categorical_cmap

if __name__ == '__main__':
    # MODELO
    exec_input, direTer, direMec = exec_setup()
    t_input, m_input = input_setup()
    rhe_data = read_rheo('data/Rhe_Param_ordenado_nuevos_indices.dat')
    model, _, _ = termomecanico(t_input, m_input)
    x_axis = model.cs.get_x_axis()
    z_axis = model.cs.get_z_axis()
    z_axis_2 = [*z_axis, *z_axis]
    geotherm = model.tm.get_geotherm()
    bys_t = model.mm.bys_t
    bys_c = model.mm.bys_c
    #geotherm_1D = geotherm[90,0,:]
    lat = -30.
    lon = -70.
    depth_from_topo = model.mm.depth_from_topo.point_depth_profile(latitude=lat,longitude=lon)
    topo_index = np.where(depth_from_topo == 0)[0]
    topo_depth = model.cs.get_z_axis()[topo_index]

    # SECCIONES 1D
    geotherm_1D = geotherm.point_depth_profile(latitude=lat,longitude=lon)
    bys_t_1D = bys_t.point_depth_profile(latitude=lat, longitude=lon)
    bys_c_1D = bys_c.point_depth_profile(latitude=lat, longitude=lon)
    bys_1D = np.concatenate((bys_c_1D, bys_t_1D))
    # LOOP DYS
    #yield_depths = [np.interp(-200,bys_c_1D[::-1],z_axis[::-1])]
    #yield_temps = [np.interp(-200,bys_c_1D[::-1],geotherm_1D[::-1])]
    yield_temps_exact = [np.interp(-200,bys_c_1D[::-1],geotherm_1D[::-1])]
    yield_depths_at_temp = [np.interp(-200,bys_c_1D[::-1],z_axis[::-1])]
    dys_list = []
    uc_params = []
    lc_params = []
    lm_params = []
    #uc_yield_temps_interp_dic = {}
    #lc_yield_temps_interp_dic = {}
    #lm_yield_temps_interp_dic = {}
    uc_yield_temps_exact_dic = {}
    lc_yield_temps_exact_dic = {}
    lm_yield_temps_exact_dic = {}
    for key, value in rhe_data.items():
        if key in set(['6','22','29']):
        #if key in set(map(str,np.arange(2,30,2))):
        #if key:
            dys = model.mm.calc_ductile_yield_strength(m_input.e, value.n,
                value.A, value.H, m_input.R, geotherm_1D)
            dys_list.append({'name': value.name, 'value': dys})
            yield_temp_exact = model.mm.calc_temperature_from_ductile_yield_strength(
                model.mm.vars.e, value.n, value.A, value.H, model.mm.vars.r, 200)
            yield_depth_at_temp = np.interp(
                yield_temp_exact, geotherm_1D, z_axis)
            yield_depths_at_temp.append(yield_depth_at_temp)
            yield_temps_exact.append(yield_temp_exact)
            #yield_depth = np.interp(200, dys[::-1], z_axis[::-1])
            #yield_temp = np.interp(200, dys[::-1], geotherm_1D[::-1])
            #yield_depths.append(yield_depth)
            #yield_temps.append(yield_temp)
            if 0 < int(key) < 11:
                uc_params.append(key)
                #uc_yield_temps_interp_dic[value.name] = yield_temp
                uc_yield_temps_exact_dic[value.name] = yield_temp_exact
            elif 11 <= int(key) < 23:
                lc_params.append(key)
                #lc_yield_temps_interp_dic[value.name] = yield_temp
                lc_yield_temps_exact_dic[value.name] = yield_temp_exact
            else:
                lm_params.append(key)
                #lm_yield_temps_interp_dic[value.name] = yield_temp
                lm_yield_temps_exact_dic[value.name] = yield_temp_exact
    #print(uc_yield_temps_interp_dic)
    #print(uc_yield_temps_exact_dic)
    #print(lc_yield_temps_interp_dic)
    #print(lc_yield_temps_exact_dic)
    #print(lm_yield_temps_interp_dic)
    #print(lm_yield_temps_exact_dic)
    fig = plt.figure(figsize=(12,7))
    min_z = z_axis[np.nanargmax(geotherm_1D)]
    major_z_ticks = np.arange(0, min_z+1, -25)
    minor_z_ticks = np.arange(5, min_z+1, -5)
    gs = gridspec.GridSpec(1,3)
    ax = fig.add_subplot(gs[0,1:])
    ax.set_xlim(-2000,2000)
    #ax.set_ylim(-180,20)
    ax.set_ylim(min_z,5)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_yticks(major_z_ticks)
    ax.set_yticks(minor_z_ticks, minor=True)
    ax.set_title('YSEs')
    ax.plot(bys_1D, z_axis_2, 'k')
    colors = categorical_cmap(3, [len(uc_params),len(lc_params),len(lm_params)],
        desaturated_first=False)(np.linspace(0,1,len(dys_list)))
    colors_iterator = iter(colors)
    for i, dys in enumerate(dys_list):
        #print(dys)
        color = next(colors_iterator)
        dys_1D = [*dys['value'], *-dys['value']]
        ax.plot(dys_1D, z_axis_2, color=color, label=dys['name'])
    ax.axvline(x=-200, color='r', linestyle='dashed')
    ax.axhline(y=topo_depth, color='k')
    ax.axvline(x=0, color='k')
    ax2 = fig.add_subplot(gs[0,0])
    ax2.set_xlim(0,1300)
    #ax2.set_ylim(-180,20)
    ax2.set_ylim(z_axis[np.nanargmax(geotherm_1D)],5)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.set_yticks(major_z_ticks)
    ax2.set_yticks(minor_z_ticks, minor=True)
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright=False)
    #ax2.set_yticks([])
    ax2.set_title('Temperatura')
    ax2.plot(geotherm_1D, z_axis)
    ax2.axhline(y=topo_depth, color='k')
    ax2.axvline(x=0, color='k')
    plot_intersections = True
    if plot_intersections is True:
        colors = np.vstack((np.array([[0, 0, 0, 1]]), colors))
        colors_iterator = iter(colors)
        for yield_depth, yield_temp in zip(
            yield_depths_at_temp, yield_temps_exact):
            color = next(colors_iterator)
            #ax.plot([-2000,-200], [yield_depth, yield_depth], 
            #    color=color, linestyle='dashed')
            #ax2.plot([yield_temp, 1300], [yield_depth, yield_depth], 
            #    color=color, linestyle='dashed')
            ax2.plot([yield_temp, yield_temp], [yield_depth, -200], 
                color=color, linestyle='dashed')
    legend = ax.legend(loc=2, bbox_to_anchor=(1.05, 1.00))
    suptitle = fig.suptitle('Lat: {}, Lon: {}'.format(lat,lon))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #plt.tight_layout(pad=7)
    extra_artists = [suptitle, legend]
    plt.savefig(direMec + '/rhe_params.png', bbox_extra_artists=extra_artists,
        bbox_inches='tight')
