import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from src.setup import input_setup, exec_setup, read_rheo
from termomecanico import termomecanico
from src.utils import makedir
from src.colormaps import categorical_cmap

if __name__ == '__main__':
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
    geotherm_1D = geotherm.point_depth_profile(latitude=lat,longitude=lon)
    bys_t_1D = bys_t.point_depth_profile(latitude=lat, longitude=lon)
    bys_c_1D = bys_c.point_depth_profile(latitude=lat, longitude=lon)
    bys_1D = np.concatenate((bys_c_1D, bys_t_1D))
    dys_list = []
    for key, value in rhe_data.items():
        dys = model.mm.calc_ductile_yield_strength(m_input.e, value.n, value.A,
            value.H, m_input.R, geotherm_1D)
        dys_list.append({'name': value.name, 'value': dys})
    fig = plt.figure(figsize=(12,7))
    gs = gridspec.GridSpec(1,3)
    ax = fig.add_subplot(gs[0,1:])
    ax.set_xlim(-2000,2000)
    ax.set_ylim(-180,20)
    ax.set_title('YSEs')
    ax.plot(bys_1D, z_axis_2, 'k')
    #colors = iter(cm.rainbow(np.linspace(0, 1, len(dys_list))))
    colors = iter(categorical_cmap(3, [10,12,7],
        desaturated_first=True)(np.linspace(0,1,len(dys_list))))
    for i, dys in enumerate(dys_list):
        print(dys)
        color = next(colors)
        dys_1D = [*dys['value'], *-dys['value']]
        ax.plot(dys_1D, z_axis_2, color=color, label=dys['name'])
    ax.plot(np.repeat(200,len(z_axis)),z_axis,'r', linestyle='dashed')
    #ax.plot(np.repeat(-200,len(z_axis)),z_axis,'r', linestyle='dashed')
    ax.axhline(y=topo_depth, color='r')
    ax2 = fig.add_subplot(gs[0,0])
    ax2.set_ylim(-180,20)
    ax2.set_yticks([])
    ax2.set_title('Temperatura')
    ax2.plot(geotherm_1D, z_axis)
    ax2.axhline(y=topo_depth, color='r')
    ax2.axvline(x=0, color='r')
    legend = ax.legend(loc=2, bbox_to_anchor=(1.05, 1.00))
    fig.suptitle('Lat: {}, Lon: {}'.format(lat,lon))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #plt.tight_layout(pad=7)
    plt.show()
