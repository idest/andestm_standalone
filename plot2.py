import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.colors import Normalize
from meccolormap import jet_white_r
from diffcolormap import get_diff_cmap
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from utils import MidPointNorm, round_to_1

def base_latitude_profile(cs,gm,lat):
    # Axes configuration
    x_axis = cs.get_x_axis()
    z_axis = cs.get_z_axis()
    x_relevant = cs.get_relevant_area().cross_section(latitude=lat)
    x_start = x_axis[np.where(x_relevant)[0][0]]
    x_finish = x_axis[np.where(x_relevant)[0][-1]]
    fig = plt.figure(figsize=(x_finish-x_start, 5))
    plt.xlim(x_start,x_finish)
    plt.ylim(-150,10)
    plt.xticks(np.arange(math.floor(x_start),math.ceil(x_finish),2))
    plt.yticks(np.arange(10,-150-1,-10))
    plt.ylabel('Profundidad[km]')
    plt.xlabel('Longitud')
    plt.tight_layout()
    plt.plot(x_axis, gm.get_topo().cross_section(latitude=lat), 'k')
    plt.plot(x_axis, gm.get_icd().cross_section(latitude=lat), 'w')
    plt.plot(x_axis, gm.get_moho().cross_section(latitude=lat), 'g')
    plt.plot(x_axis, gm.get_slab_lab().cross_section(latitude=lat), 'r')
    return fig

def thermal_latitude_profile(tm, lat, save_dir=None, show=False):
    # Base
    cs = tm.get_coordinate_system()
    gm = tm.get_geometric_model()
    fig = base_latitude_profile(cs, gm, lat)
    x_axis = cs.get_x_axis()
    z_axis = cs.get_z_axis()
    xx, zz = np.meshgrid(x_axis, z_axis)
    # Heatmap
    temps = tm.get_geotherm().cross_section(latitude=lat)
    temps_masked = np.ma.masked_invalid(temps)
    heatmap = plt.pcolormesh(
        xx, zz, temps_masked.T, cmap=cm.coolwarm,shading='gouraud')
    t_min, t_max = 0, 1300
    plt.clim(t_min, t_max)
    cbar = plt.colorbar(heatmap, aspect='20')
    cbar.set_label('Temperatura', rotation=90, labelpad=-60)
    # Contour Map
    numlabels=6
    contour = plt.contour(xx, zz, temps_masked.T, numlabels, colors='k')
    plt.clabel(contour, fontsize=9, fmt='%.0f', rightside_up=True)
    plt.title('Perfil Termal %s' %(lat))
    if show is True:
        plt.show()
    if save_dir:
        fig.savefig(save_dir + '%s' %(lat), dpi='figure', format='pdf')
    plt.close()
    return

def mechanic_latitude_profile(mm, lat, save_dir=None, show=False):
    # Base
    cs = mm.get_coordinate_system()
    gm = mm.get_geometric_model()
    fig = base_latitude_profile(cs, gm, lat)
    x_axis = cs.get_x_axis()
    z_axis = cs.get_z_axis()
    xx, zz = np.meshgrid(x_axis, z_axis)
    # Heatmap
    meca = mm.get_yse()[0].cross_section(latitude=lat)
    meca_masked = np.ma.masked_invalid(meca)
    heatmap = plt.pcolormesh(xx, zz, meca_masked.T, cmap=jet_white_r,
        shading='gouraud')
    ys_min, ys_max = 0, 200
    plt.clim(ys_min, ys_max)
    cbar = plt.colorbar(heatmap, aspect='20')
    cbar.set_label('Yield Strength', rotation=90, labelpad=-50)
    plt.title('Perfil Mecanico %s' %(lat))
    if show is True:
        plt.show()
    if save_dir:
        fig.savefig(save_dir + '%s' %(lat), dpi='figure', format='pdf')
    plt.close()
    return

def base_map(topo=True):
    map = Basemap(
        llcrnrlon= -80, llcrnrlat= -45,
        urcrnrlon= -60.0, urcrnrlat= -10.0,
        epsg= 4326, resolution = 'l', suppress_ticks=False)
    #map.arcgisimage(service='ESRI_Imagery_World_2D',xpixels=2000,verbose=True)
    #map.drawparallels(np.arange(-90,90,5), labels=[1,0,0,0], fontsize=7)
    #map.drawmeridians(np.arange(-180,180,5), labels=[0,0,0,1], fontsize=7)
    if topo is True:
        map.etopo()
    map.drawcoastlines(linewidth=0.5)
    map.drawlsmask(land_color='0.8', ocean_color='0.8', resolution='l')
    return map

def shf_map_old(
        shf=None, data_coords=None, data=None, diff=None,
        data_types=None, rmse=None, topo=True, name='Surface_Heat_Flow',
        save_dir=None):
    #fig, main_ax = plt.subplots(figsize=(6,6))
    fig, main_ax = plt.subplots()
    divider = make_axes_locatable(main_ax)
    diff_ax = main_ax
    map = base_map(topo=False)
    if shf is not None:
        x_axis = shf.cs.get_x_axis()
        y_axis = shf.cs.get_y_axis()
        xx, yy = np.meshgrid(x_axis, y_axis)
        shf_max, shf_min = np.nanmax(shf), np.nanmin(shf)
        cbar_max, cbar_min = shf_max, shf_min
        #if data:
        #    vmax, vmin = [shf_max, np.nanmax(data)], [shf_min, np.nanmin(data)]
        #    cbar_max, cbar_min = np.nanmax(vmax), np.nanmin(vmin)
        shf_masked = np.ma.masked_invalid(shf) # Before: shf*-1 ¿?
        shf_heatmap = map.pcolormesh(
            xx, yy, shf_masked.T,
            cmap='afmhot_r', shading='gouraud',
            vmin=cbar_min, vmax=cbar_max)
        shf_cbar_ax = divider.append_axes('right', '5%', pad='12%')
        plt.sca(main_ax)
        shf_cbar = plt.colorbar(shf_heatmap, cax=shf_cbar_ax)# pad=0.2)
        shf_cbar.set_label('Heat Flow (W/m2)', rotation=90, labelpad=-55)
        #plt.tight_layout()
        main_ax.set_title('Surface Heat Flow')

    if data is not None or diff is not None:
        if data_coords is None:
            raise ValueError('data_cords variable is missing.')
        else:
            lon = data_coords[:,0]
            lat = data_coords[:,1]
            map_lon, map_lat = map(lon, lat)
        markers = ['o']
        data_markers = np.array(markers * len(lon))
        m_scatter = {}
        if data_types is not None:
            _, indices = np.unique(data_types, return_inverse=True)
            markers = np.array(['o', '^', 'p', 's'])
            data_markers = markers[indices]
        def map_scatter(scatter_data, **kwargs):
            scatter_data = np.ma.masked_invalid(scatter_data)
            for m in markers:
                mask = data_markers == m
                m_scatter_data = scatter_data[mask]
                m_lon, m_lat = map_lon[mask], map_lat[mask]
                scatter = map.scatter(
                    m_lon, m_lat, latlon=True,
                    c=m_scatter_data, marker=m, **kwargs)
                if data_types is not None:
                    m_scatter[m] = scatter
            return m_scatter, scatter

    if data is not None:
        if shf is None:
            cbar_max, cbar_min = np.nanmax(data), np.nanmin(data)
            shf_cbar_ax = divider.append_axes('right', '5%', pad='12%')
            plt.sca(main_ax)
        m_scatter, scatter = map_scatter(
            data, cmap='afmhot_r',
            vmin=cbar_min, vmax=cbar_max)
        shf_cbar = plt.colorbar(scatter, cax=shf_cbar_ax)# pad=0.2)
        shf_cbar.set_label('Heat Flow (W/m2)', rotation=90, labelpad=-55)
        main_ax.set_title('Surface Heat Flow')

    if diff is not None:
        diff_max = np.nanmax(diff)
        diff_min = np.nanmin(diff)
        diff_limit = np.nanmax([abs(diff_max), abs(diff_min)])
        diff_limit = round_to_1(diff_limit, 'ceil')
        ticks = np.arange(-diff_limit, diff_limit+0.01, 0.01)
        bins = len(ticks) - 1
        diff_cmap = get_diff_cmap(bins)
        norm = MidPointNorm(midpoint=0, vmin=-diff_limit, vmax=diff_limit)

        diff_cbar_pad = '12%'
        if shf is not None and data is None:
            diff_cbar_pad = '35%'
        if data is not None:
            diff_ax = divider.append_axes('right', '100%', pad='35%')#, pad='50%')
            diff_ax.set_yticks([])
            map = base_map(topo=False)
        m_scatter, scatter = map_scatter(
            diff, cmap=diff_cmap, norm=norm,
            edgecolors='k', linewidths=.3)
        diff_cbar_ax = divider.append_axes('right', '5%', pad=diff_cbar_pad)
        plt.sca(diff_ax)
        diff_cbar = plt.colorbar(scatter, cax=diff_cbar_ax)# pad=0.2)
        diff_cbar.set_label('Diff [W/m²]', rotation=90, labelpad=-18)
        diff_cbar.set_ticks([])
        hist_ax = divider.append_axes('right', '30%', pad='0%')
        plt.sca(diff_ax)
        N, bins, patches = hist_ax.hist(
            diff, bins=ticks,
            orientation='horizontal')
        hist_ax.set_yticks(ticks)
        hist_ax.set_ylim([-diff_limit, diff_limit])
        hist_ax.yaxis.tick_right()
        norm = Normalize(bins.min(), bins.max())
        for bin, patch in zip(bins, patches):
            color = diff_cmap(norm(bin))
            patch.set_facecolor(color)
        #hist_ax.grid(True)
        diff_ax.set_title('Modelled Surface Heat Flow minus Data')

    extra_artists = []
    if rmse is not None:
        rmse_text = main_ax.annotate(
            'RMSE = %0.7f' %(rmse), xy=(0,-0.13),
            xycoords='axes fraction')
        extra_artists.append(rmse_text)

    if data_types is not None:
        legend = main_ax.legend(
            [m_scatter[markers[0]], m_scatter[markers[1]],
             m_scatter[markers[2]], m_scatter[markers[3]]],
            ['ODP Borehole', 'Land Borehole',
             'Geochemical', 'Marine Geophysics'],
            bbox_to_anchor=(0.5, 0), loc='upper center',
            ncol=4, bbox_transform=fig.transFigure)
        extra_artists.append(legend)

    plt.tight_layout()
    if save_dir:
        plt.savefig(
            save_dir + '%s' %(name),
            bbox_extra_artists=(extra_artists), bbox_inches='tight',
            dpi='figure', format='pdf')
    plt.close()
    return

def shf_map(
        shf, cbar_limits=None, map=None, ax=None, save_dir=None,
        name='shf_map', return_width_ratio=False):
    # Axes and map setup
    if ax is None:
        fig, ax = plt.subplots()
    if map is None:
        map = base_map(topo=False)
    # Pcolormesh
    x_axis = shf.cs.get_x_axis()
    y_axis = shf.cs.get_y_axis()
    xx, yy = np.meshgrid(x_axis, y_axis)
    cbar_max, cbar_min = np.nanmax(shf), np.nanmin(shf)
    if cbar_limits is not None:
        cbar_max, cbar_min = cbar_limits[0], cbar_limits[1]
    shf_masked = np.ma.masked_invalid(shf) # Before: shf*-1 ¿?
    shf_heatmap = map.pcolormesh(
        xx, yy, shf_masked.T,
        cmap='afmhot_r', shading='gouraud',
        vmin=cbar_min, vmax=cbar_max)
    # Colorbar
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes('right', '5%', pad='12%')
    plt.sca(ax)
    cbar = plt.colorbar(shf_heatmap, cax=cbar_ax)# pad=0.2)
    cbar.set_label('Heat Flow (W/m2)', rotation=90, labelpad=-55)
    # Title
    ax.set_title('Surface Heat Flow')
    # Options
    #plt.tight_layout()
    if save_dir:
        plt.savefig(
            save_dir + '%s' %(name), bbox_inches='tight')
            #dpi='figure', format='pdf')
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12
        return width_ratio

def get_map_scatter_function(data_coords, data_types, map):
    # Coords
    if data_coords is None:
        raise ValueError('data_cords variable is missing.')
    else:
        lon = data_coords[:,0]
        lat = data_coords[:,1]
        map_lon, map_lat = map(lon, lat)
    # Markers
    markers = ['o']
    data_markers = np.array(markers * len(lon))
    m_scatter = {}
    if data_types is not None:
        _, indices = np.unique(data_types, return_inverse=True)
        markers = np.array(['o', '^', 'p', 's'])
        labels = np.array(['ODP Borehole', 'Land Borehole',
                           'Geochemical', 'Marine Geophysics'])
        data_markers = markers[indices]
    # Map Scatter Function
    def map_scatter(scatter_data, **kwargs):
        scatter_data = np.ma.masked_invalid(scatter_data)
        for m, l in zip(markers, labels):
            mask = data_markers == m
            m_scatter_data = scatter_data[mask]
            m_lon, m_lat = map_lon[mask], map_lat[mask]
            scatter = map.scatter(
                m_lon, m_lat, latlon=True,
                c=m_scatter_data, marker=m, label=l,
                **kwargs)
            if data_types is not None:
                m_scatter[m] = scatter
        return m_scatter, scatter
    return map_scatter

def data_map(
        data, cbar_limits=None, data_coords=None, data_types=None,
        map=None, ax=None, cbar=True, legend=True,
        rmse=None, return_width_ratio=False,
        save_dir=None, name='data_map'):
    # Axes and map setup
    if ax is None:
        fig, ax = plt.subplots()
    if map is None:
        map = base_map(topo=False)
    # Scatter map
    cbar_max, cbar_min = np.nanmax(data), np.nanmin(data)
    if cbar_limits is not None:
        cbar_max, cbar_min = cbar_limits[0], cbar_limits[1]
    map_scatter = get_map_scatter_function(data_coords, data_types, map)
    m_scatter, scatter = map_scatter(data, cmap='afmhot_r', vmin=cbar_min,
                                     vmax=cbar_max)
    # Colorbar
    if cbar is True:
        divider = make_axes_locatable(ax)
        cbar_ax = divider.append_axes('right', '5%', '12%')
        cbar = plt.colorbar(scatter, cax=cbar_ax)# pad=0.2)
        cbar.set_label('Heat Flow (W/m2)', rotation=90, labelpad=-55)
    # Title
    ax.set_title('Surface Heat Flow')
    # Options
    extra_artists = []
    if rmse is not None:
        # RMSE
        rmse_text = plt.figtext(
            0.4, 0.01, 'RMSE: %0.7f' %(rmse),
            fontweight='bold')
        extra_artists.append(rmse_text)
    if legend is True:
        legend = ax.legend(bbox_to_anchor=(0.5, 0.0), loc='upper center',
                           ncol=2, bbox_transform=fig.transFigure)
        extra_artists.append(legend)
    if save_dir:
        plt.savefig(
            save_dir + '%s' %(name),
            bbox_extra_artists=extra_artists, bbox_inches='tight')
            #dpi='figure', format='pdf')
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12
        return width_ratio

def diff_map(
        diff, data_coords=None, data_types=None, map=None, ax=None,
        rmse=None, legend=True, return_width_ratio=False,
        save_dir=None, name='diff_map', e_prom=None, sigmas=None):
    # Axes and map setup
    if ax is None:
        fig, ax = plt.subplots()
    if map is None:
        map = base_map(topo=False)
    # Scatter map
    diff_max = np.nanmax(diff)
    diff_min = np.nanmin(diff)
    diff_limit = np.nanmax([abs(diff_max), abs(diff_min)])
    diff_limit = round_to_1(diff_limit, 'ceil')
    ticks = np.arange(-diff_limit, diff_limit+0.01, 0.01)
    bins = len(ticks) - 1
    diff_cmap = get_diff_cmap(bins)
    norm = MidPointNorm(midpoint=0, vmin=-diff_limit, vmax=diff_limit)
    map_scatter = get_map_scatter_function(data_coords, data_types, map)
    m_scatter, scatter = map_scatter(
        diff, cmap=diff_cmap, norm=norm,
        edgecolors='k', linewidths=.3)
    # Colorbar
    cbar_pad = '12%'
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes('right', '5%', pad=cbar_pad)
    plt.sca(ax)
    cbar = plt.colorbar(scatter, cax=cbar_ax)# pad=0.2)
    cbar.set_label('Diff [W/m²]', rotation=90, labelpad=-18)
    cbar.set_ticks([])
    # Colorbar histogram
    hist_ax = divider.append_axes('right', '30%', pad='0%')
    plt.sca(ax)
    N, bins, patches = hist_ax.hist(
        diff, bins=ticks,
        orientation='horizontal')
    hist_ax.set_yticks(ticks)
    hist_ax.set_ylim([-diff_limit, diff_limit])
    hist_ax.yaxis.tick_right()
    if sigmas is not None:
        hist_ax.axhline(y=sigmas.n_1_sigma)
        hist_ax.text(15,sigmas.n_1_sigma-0.005,r'-$\sigma$',size='small')
        hist_ax.axhline(y=sigmas.p_1_sigma)
        hist_ax.text(15,sigmas.p_1_sigma+0.002,r'+$\sigma$',size='small')
        hist_ax.axhline(y=sigmas.n_2_sigma)
        hist_ax.text(15,sigmas.n_2_sigma-0.005,r'-2$\sigma$',size='small')
        hist_ax.axhline(y=sigmas.p_2_sigma)
        hist_ax.text(15,sigmas.p_2_sigma+0.002,r'+2$\sigma$',size='small')
    norm = Normalize(bins.min(), bins.max())
    for bin, patch in zip(bins, patches):
        color = diff_cmap(norm(bin))
        patch.set_facecolor(color)
    #hist_ax.grid(True)
    # Title
    ax.set_title('Model minus Data')
    # Options
    extra_artists=[]
    if e_prom is not None:
        # MAE
        e_prom_text = plt.figtext(
            0.4,0.03, 'MAE: %0.3f' %(e_prom),
            fontweight='bold')
        extra_artists.append(e_prom_text)
    if rmse is not None:
        # RMSE
        rmse_text = plt.figtext(
            0.4, 0, 'RMSE: %0.7f' %(rmse),
            fontweight='bold')
        extra_artists.append(rmse_text)
    if legend is True:
        legend = ax.legend(bbox_to_anchor=(0.5, 0.0), loc='upper center',
                           ncol=2, bbox_transform=fig.transFigure)
        extra_artists.append(legend)
    if save_dir:
        plt.savefig(
            save_dir + '%s' %(name),
            bbox_extra_artists=extra_artists, bbox_inches='tight')
            #dpi='figure', format='pdf')
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12 + 0.30
        return width_ratio

def multi_map(
        shf=None, data=None, diff=None, data_coords=None, data_types=None,
        rmse=None, topo=True, save_dir=None, name='multi_map',
        e_prom=None, sigmas=None):
    # Gridspec
    fig = plt.figure(figsize=(9,6))
    gs = gridspec.GridSpec(1,2)
    # Axis 1: Surface Heat Flow + Data
    ax1 = fig.add_subplot(gs[0,0])
    map1 = base_map(topo=False)
    cbar_limits = [np.nanmax(shf), np.nanmin(shf)]
    wr1 = shf_map(
        shf, cbar_limits=cbar_limits,
        map=map1, ax=ax1,
        return_width_ratio=True)
    data_map(
        data, cbar_limits=cbar_limits,
        data_coords=data_coords,
        data_types=data_types,
        map=map1, ax=ax1,
        legend=False, cbar=False)
    # Axis 2: Diff
    ax2 = fig.add_subplot(gs[0,1])
    map2 = base_map(topo=False)
    wr2 = diff_map(
        diff,
        data_coords=data_coords,
        data_types=data_types,
        map=map2, ax=ax2,
        legend=False,  return_width_ratio=True,
        save_dir=None, name='diff_map',
        sigmas=sigmas)
    ax2.set_yticks([])
    # Aspect
    gs.set_width_ratios([wr1,wr2])
    plt.tight_layout()
    # Extra Artists
    extra_artists=[]
    #MEA
    e_prom_text = plt.figtext(
        0.03, -0.02, 'MAE: %0.3f' %(e_prom),
        fontweight='bold')
    extra_artists.append(e_prom_text)
    # RMSE
    rmse_text = plt.figtext(
        0.03, 0, 'RMSE: %0.7f' %(rmse),
        fontweight='bold')
    extra_artists.append(rmse_text)
    ## Legend
    legend = plt.legend(
        loc="center left", bbox_to_anchor=(0.2, 0),
        ncol=4, bbox_transform=fig.transFigure)
    extra_artists.append(legend)
    # Options
    if save_dir:
        plt.savefig(
            save_dir + '%s' %(name),
            bbox_extra_artists=(extra_artists), bbox_inches='tight')
            #dpi='figure', format='pdf')
    plt.close()

def rmse_plot(vnames, vaxes, rmses, save_dir=None):
    # x axis
    print(type(vnames))
    print(vnames)
    print(vaxes)
    x_name = vnames[0]
    x_axis = vaxes[0]
    plt.xlabel(str(x_name))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    if len(vnames) > 1:
        # y axis
        y_name = vnames[1]
        y_axis = vaxes[1]
        plt.ylabel(str(y_name))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        # 2D matrix
        vmin = np.min(rmses)
        vmax = np.max(rmses)
        v = np.linspace(vmin, vmax, 50)
        plt.contourf(
            x_axis, y_axis, rmses.T, v, norm=colors.PowerNorm(gamma=1./2.))
        #plt.pcolormesh(
        #    x_axis, y_axis, rmses.T, norm=colors.PowerNorm(gamma=2./3.))
        plt.colorbar()
        name = 'RMSE_2D.png'
    else:
        index = np.arange(len(x_axis))
        plt.plot(index, rmses, '-r', linewidth=1.)
        plt.bar(index, rmses, alpha=.4)
        diff = max(rmses) - min(rmses)
        plt.ylim(min(rmses)-0.2*diff, max(rmses)+0.2*diff)
        plt.xticks(index, x_axis)
        plt.ylabel('RMSE')
        plt.grid(True)
        name = 'RMSE.png'
    plt.tight_layout()
    if save_dir is not None:
        plt.savefig(save_dir + '%s' %(name))#,dpi='figure', format='pdf')

def data_scatter_plot(data, data_error, data_types, shf_models, save_dir=None):
    fig, axes = plt.subplots(4)
    markers = np.array(['o', '^', 'p', 's'])
    titles = np.array(['ODP', 'LBH', 'GQ', 'GP'])
    data_types = (np.array(data_types) - 1).astype(int)
    data_markers = markers[data_types]
    colors = cm.rainbow(np.linspace(0, 1, len(shf_models)))
    for m, t, ax in zip(markers, labels, axes):
        mask = data_markers == m
        m_data = data[mask]
        m_data_error = data_error[mask]
        m_data_axis = np.arange(len(m_data))
        ax.errorbar(
            m_data_axis, m_data, m_data_error, fmt=m, capsize=2, capthick=0.5,
            markersize=2, elinewidth=0.5)
        for shf_model, c in zip(shf_models, colors):
            m_shf_model = shf_model[mask]
            ax.scatter(
                m_data_axis, m_shf_model, marker='.',
                color=c, label=l, s=3)

        ax.set_title(t)
    plt.tight_layout()
    if save_dir is not None:
        name = 'scatter.pdf'
        plt.savefig(save_dir + '%s' %(name), dpi='figure', format='pdf')


