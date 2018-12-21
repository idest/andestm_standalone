import numpy as np
import pandas as pd
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from src.utils import makedir_from_filename, calc_deviation
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.colors import Normalize
from src.colormaps import (
    jet_white_r, get_diff_cmap, get_elevation_diff_cmap, jet_white_r_2)
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from src.utils import MidPointNorm, round_to_1, get_magnitude
from src.compute import SpatialArray2D

def base_latitude_profile(cs, gm, lat, earthquakes=None, sli=None):
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
    if earthquakes is not None:
        eq = earthquakes[
            (earthquakes['latitude'] >= lat - 0.1) &
            (earthquakes['latitude'] < lat + 0.1)]
        #bins = [1, 5.5, 6.5, 7.0, 7.5, 9.0]
        #sizes = [15*2**n for n in range(len(bins)-1)]
        #print(sizes)
        #labels = [
        #    'Mb < 5.5',
        #    '5.5 < Mb < 6.5',
        #    '6.5 < Mb < 7.0',
        #    '7.0 < Mb < 7.5',
        #    'Mb > 7.5']
        #eq['size'] = pd.cut(eq['mag'].values, bins, labels=sizes)
        #for i, size in enumerate(sizes):
        #    # current size eqs
        #    cs_eq = eq[np.isclose(eq['size'], size)]
        #    plt.scatter(
        #        cs_eq['longitude'], -cs_eq['depth'], color='orange',
        #        edgecolors='black', s=size,
        #        zorder=1000, label=labels[i])
        #plt.legend()
        plt.scatter(
            eq['longitude'], -eq['depth'], color='orange',
            s=30, edgecolors='black', linewidth=0.2,
            zorder=1000)
        print('lat:', lat)#, 'eq:', eq)
    if sli is not None:
        plt.scatter(sli['lon'], sli['depth'], color='r', zorder=1001, s=100)
    return fig

def thermal_latitude_profile(tm, lat, show=False, filename=None,
        earthquakes=None, sli=None):
    # Base
    cs = tm.get_coordinate_system()
    gm = tm.get_geometric_model()
    fig = base_latitude_profile(cs, gm, lat, sli=sli)
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
    plt.title('Perfil Termal %.1f' %(lat))
    if show is True:
        plt.show()
    if filename:
        makedir_from_filename(filename)
        filename = filename + '_%.1f' %(lat) + '.png'
        fig.savefig(filename, transparent=True)#, dpi='figure', format='pdf')
    plt.close()
    return

def mechanic_latitude_profile(mm, lat, show=False, filename=None,
        earthquakes=None, sli=None):
    # Base
    cs = mm.get_coordinate_system()
    gm = mm.get_geometric_model()
    fig = base_latitude_profile(cs, gm, lat, earthquakes=earthquakes, sli=sli)
    x_axis = cs.get_x_axis()
    z_axis = cs.get_z_axis()
    xx, zz = np.meshgrid(x_axis, z_axis)
    # Heatmap
    meca = -mm.get_yse()[1].cross_section(latitude=lat)
    meca_masked = np.ma.masked_invalid(meca)
    heatmap = plt.pcolormesh(xx, zz, meca_masked.T, cmap=jet_white_r_2,
        shading='gouraud')
    ys_min, ys_max = 0, 200
    plt.clim(ys_min, ys_max)
    cbar = plt.colorbar(heatmap, aspect='20')
    cbar.set_label('Yield Strength', rotation=90, labelpad=-50)
    plt.title('Perfil Mecanico %.1f' %(lat))
    if show is True:
        plt.show()
    if filename:
        makedir_from_filename(filename)
        filename = filename + '_%.1f' %(lat) + '.png'
        fig.savefig(filename, transparent=True)#, dpi='figure', format='pdf')
    plt.close()
    return

def elastic_thickness_latitude_profile(mm, lat, show=False, filename=None,
        earthquakes=None, sli=None):
    #Base
    cs = mm.get_coordinate_system()
    gm = mm.get_geometric_model()
    fig = base_latitude_profile(cs, gm, lat, earthquakes=earthquakes, sli=sli)
    x_axis = cs.get_x_axis()
    z_axis = cs.get_z_axis()
    xx, zz = np.meshgrid(x_axis, z_axis)
    # Thickness
    topo = mm.get_geometric_model().get_topo().cross_section(latitude=lat)
    eet = mm.get_eet_from_trench().cross_section(latitude=lat)
    eet_depth = topo - eet
    plt.plot(x_axis, eet_depth, 'black')
    ax = plt.gca()
    ax.fill_between(x_axis, topo, eet_depth, interpolate=True, color='grey')
    plt.title('Perfil Espesor Elástico %.1f' %(lat))
    if show is True:
        plt.show()
    if filename:
        makedir_from_filename(filename)
        filename = filename + '_%.1f' %(lat) + '.png'
        fig.savefig(filename, transparent=True)#, dpi='figure', format='pdf')
    plt.close()

def base_map(topo=True, draw_land=True):
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
    if draw_land is True:
        pass
        #map.drawlsmask(land_color='0.8', ocean_color='0.8', resolution='l')
    return map

def boolean_map(
        array_2D, array_2D_2=None, color='green', map=None, ax=None,
        alpha=0.6, filename=None, return_width_ratio=False,
        title=None, cmap_idx=0):
    #Axes and map setup
    if ax is None:
        fig, ax = plt.subplots()
    if map is None:
        map = base_map(topo=False)
    #Imshow
    map.drawlsmask(land_color='0.0', ocean_color='0.0', resolution='l')
    """
    cmap = plt.cm.copper
    cmap_colors = cmap(np.arange(cmap.N))
    cmap_colors[:, -1] = np.linspace(0, 1, cmap.N)
    cmap = colors.ListedColormap(cmap_colors)
    """
    cmap1 = colors.ListedColormap([[1,1,1,0],[0.75,0.75,0,0.5]])
    cmap2 = colors.ListedColormap([[1,1,1,0],[0,0.75,0.75,0.5]])
    cmaps = [cmap1,cmap2]
    map.imshow(array_2D.T, cmap=cmaps[cmap_idx], origin='upper')
    if array_2D_2 is not None:
        map.imshow(array_2D_2.T, cmap=cmap2, origin='upper')#, alpha=alpha)
    # Title
    if title is not None:
        ax.set_title(title)
    #Options
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(filename, bbox_inches='tight')
            #dpi='figure', format='pdf')
        plt.close()
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12
        return width_ratio

def earthquake_map(earthquakes, map=None, ax=None, filename=None,
        return_width_ratio=False, title=None):
    # Axes and map setup
    if ax is None:
        fig, ax = plt.subplots()
    if map is None:
        map = base_map(topo=False, draw_land=False)
    # Earthquakes
    extra_artists = []
    if earthquakes is not None:
        eqs = earthquakes
        #bins = [1, 5.5, 6.5, 7.0, 7.5, 9.0]
        #sizes = [2*2**n for n in range(len(bins)-1)]
        #print(sizes)
        #labels = [
        #    'Mb < 5.5',
        #    '5.5 < Mb < 6.5',
        #    '6.5 < Mb < 7.0',
        #    '7.0 < Mb < 7.5',
        #    'Mb > 7.5']
        #eqs['size'] = pd.cut(eqs['mag'].values, bins, labels=sizes)
        #for i, size in enumerate(sizes):
        #    # current size eqs
        #    cs_eqs = eqs[np.isclose(eqs['size'], size)]
        #    scatter = map.scatter(
        #        cs_eqs['longitude'], cs_eqs['latitude'], s=size,
        #        facecolor=(1.0,1.0,1.0,0.5), edgecolors=cs_eqs['color'],
        #        latlon=True, zorder=1000, label=labels[i])
        #legend = plt.legend(
        #    loc='center left', bbox_to_anchor=(0.7,0.8),
        #    bbox_transform=fig.transFigure, prop={'size': 6})
        #extra_artists.append(legend)
        scatter = map.scatter(
            eqs['longitude'], eqs['latitude'], s=0.2,
            facecolor=eqs['color'],
            latlon=True, zorder=1000)

    # Title
    if title is not None:
        ax.set_title(title)
    # Options
    #plt.tight_layout()
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(
            filename, bbox_inches='tight', transparent=True, dpi=900,
            bbox_extra_artists=(extra_artists))
            #dpi='figure', format='pdf')
        plt.close()
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12
        return width_ratio

def plot_eet_equivalent_vs_effective(eet_effective_dict, eet_eq,
        save_dir='EET', name='eet_diff'):
    for eet_effective in eet_effective_dict.values():
        eet_eff = SpatialArray2D(
            np.loadtxt(eet_effective['file']), eet_eq.cs).mask_irrelevant_eet()
        eet_diff = eet_eq - eet_eff
        sd = calc_deviation(eet_eq, eet_eff)
        diff_map(eet_eq, eet_eff, eet_diff, sd=sd,
            colormap=jet_white_r,
            colormap_diff = get_elevation_diff_cmap(100),
            cbar_limits=[0,100], cbar_limits_diff=[-100,100],
            cbar_label='EET [km.]', cbar_label_diff='Dif. EET [km.]',
            title_1='Espesor Elástico Equivalente',
            title_2='Espesor Elástico Efectivo',
            title_3='Diff. (EET eq. - EET ef.)',
            labelpad=-48, labelpad_diff=-56,
            filename=save_dir + eet_effective['dir'] + name)

def heatmap_map(
        array_2D, colormap=None, cbar_limits=None, map=None, ax=None, alpha=1,
        filename=None, return_width_ratio=False,
        cbar_label=None, title=None, labelpad=-40, earthquakes=None,
        cbar_ticks=None, cbar_tick_labels=None, norm=None, draw_land=True):
    # Axes and map setup
    if ax is None:
        fig, ax = plt.subplots()
    if map is None:
        map = base_map(topo=False, draw_land=draw_land)
    # Earthquakes
    if earthquakes is not None:
        eqs_lon = earthquakes['longitude']
        eqs_lat = earthquakes['latitude']
        scatter = map.scatter(
            eqs_lon, eqs_lat, .5, c=earthquakes['color'], latlon=True,
            zorder=1000)
    # Pcolormesh
    x_axis = array_2D.cs.get_x_axis()
    y_axis = array_2D.cs.get_y_axis()
    xx, yy = np.meshgrid(x_axis, y_axis)
    cbar_max, cbar_min = np.nanmax(array_2D), np.nanmin(array_2D)
    if cbar_limits is not None:
        cbar_max, cbar_min = cbar_limits[0], cbar_limits[1]
    array_2D_masked = np.ma.masked_invalid(array_2D) # Before: shf*-1 ¿?
    kwargs = {}
    if colormap is not None:
        kwargs['cmap'] = colormap
    plt.sca(ax)
    array_2D_heatmap = map.pcolormesh(
        xx, yy, array_2D_masked.T, shading='gouraud', alpha=alpha,
        vmin=cbar_min, vmax=cbar_max, norm=norm, **kwargs)
    # Format Coord (Interactive Plot)
    def format_coord(x,y):
        cs = array_2D.cs
        lon = cs.round_to_step(x, step=0.2, prec=1)
        lat = cs.round_to_step(y, step=0.2, prec=1)
        if (lon > cs.data.min_lon and lon < cs.data.max_lon
                and lat > cs.data.min_lat and lat < cs.data.max_lat):
            z = array_2D.extract_point(latitude=lat, longitude=lon)
            return 'x={:2.1f}, y={:2.1f}, z={:1.4f}'.format(lon, lat, z)
        else:
            return 'x={:2.1f}, y={:2.1f}'.format(lon, lat)
    ax.format_coord = format_coord
    # Colorbar
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes('right', '5%', pad='12%')
    plt.sca(ax)
    cbar = plt.colorbar(
        array_2D_heatmap, ticks=cbar_ticks, cax=cbar_ax)
    if cbar_tick_labels is not None:
        cbar.ax.set_yticklabels(cbar_tick_labels)
    # pad=0.2)
    if cbar_label is not None:
        cbar.set_label(cbar_label, rotation=90, labelpad=labelpad)
    # Title
    if title is not None:
        ax.set_title(title)
    # Options
    #plt.tight_layout()
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(
            filename, bbox_inches='tight', transparent=True)
            #dpi='figure', format='pdf')
        plt.close()
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12
        return width_ratio

def diff_map(array_1, array_2, diff, sd=None,
        colormap=jet_white_r, colormap_diff=get_elevation_diff_cmap(100),
        cbar_label=None, cbar_label_diff=None,
        cbar_limits=None, cbar_limits_diff=None,
        title_1='Matriz 1', title_2='Matriz 2', title_3='Dif. M1 - M2',
        axs=None, filename=None, labelpad=-45, labelpad_diff=-45):
    # Axes and map setup
    if axs is None:
        fig = plt.figure(figsize=(12,6))
        nc = 3
        gs = gridspec.GridSpec(1,nc)
        axs = [fig.add_subplot(gs[0,n]) for n in np.arange(nc)]
        #axs = plt.subplots(1,3,figsize=(12,6))
    #fig = plt.figure(figsize=(12,6))
    #gs = gridspec.GridSpec(1,3)
    if cbar_limits is None:
        cbar_max = max(np.nanmax(array_1), np.nanmax(array_2))
        cbar_min = min(np.nanmin(array_1), np.nanmin(array_2))
        cbar_limits = [cbar_min, cbar_max]
    if cbar_limits_diff is None:
        cbar_max_abs = max(np.nanmax(diff), abs(np.nanmin(diff)))
        cbar_limits_diff = [-cbar_max_abs, cbar_max_abs]
    # Axis 1: array_1
    #ax1 = fig.add_subplot(gs[0,0])
    ax1 = axs[0]
    plt.sca(ax1)
    map1 = base_map(topo=False)
    wr1 = heatmap_map(
        array_1, colormap=colormap, cbar_label=cbar_label,
        cbar_limits=cbar_limits,
        title=title_1, map=map1, ax=ax1,
        return_width_ratio=True, labelpad=labelpad)
    # Axis 2: array_2 
    #ax2 = fig.add_subplot(gs[0,1])
    ax2 = axs[1]
    plt.sca(ax2)
    map2 = base_map(topo=False)
    wr2 = heatmap_map(
        array_2, colormap=colormap, cbar_label=cbar_label,
        cbar_limits=cbar_limits,
        title=title_2, map=map2, ax=ax2,
        return_width_ratio=True, labelpad=labelpad)
    ax2.set_yticks([])
    # Axis 3: Diff: EET equivalente - EET efectivo
    #ax3 = fig.add_subplot(gs[0,2])
    ax3 = axs[2]
    plt.sca(ax3)
    map3 = base_map(topo=False)
    if sd is not None:
        sd_string = ' D.E.: {:.2f}'.format(sd)
    else:
        sd_string = ''
    wr3 = heatmap_map(
        diff, colormap=colormap_diff, cbar_label=cbar_label_diff,
        cbar_limits=cbar_limits_diff,
        title=title_3 + sd_string,
        map=map3, ax=ax3, return_width_ratio=True, labelpad=labelpad_diff)
    ax3.set_yticks([])
    plt.tight_layout()
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(filename, transparent=True)
        plt.close()

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

def data_scatter_map(
        data, cbar_limits=None, data_coords=None, data_types=None,
        map=None, ax=None, cbar=True, legend=True,
        rmse=None, return_width_ratio=False,
        filename=None):
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
    m_scatter, scatter = map_scatter(
        data, cmap='afmhot',
        vmin=cbar_min, vmax=cbar_max,
        edgecolors='black', linewidths=0.5)
    # Colorbar
    if cbar is True:
        divider = make_axes_locatable(ax)
        cbar_ax = divider.append_axes('right', '5%', '12%')
        cbar = plt.colorbar(scatter, cax=cbar_ax)# pad=0.2)
        cbar.set_label('Heat Flow [mW/m²]', rotation=90, labelpad=-45)
    # Title
    ax.set_title('Surface Heat Flow')
    # Options
    extra_artists = []
    if rmse is not None:
        # RMSE
        rmse_text = plt.figtext(
            0.4, 0.01, 'RMSE: %0.2f' %(rmse),
            fontweight='bold')
        extra_artists.append(rmse_text)
    if legend is True:
        legend = ax.legend(bbox_to_anchor=(0.5, 0.0), loc='upper center',
                           ncol=2, bbox_transform=fig.transFigure)
        extra_artists.append(legend)
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(
            filename,
            bbox_extra_artists=extra_artists, bbox_inches='tight', transparent=True, dpi=900)
            #dpi='figure', format='pdf')
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12
        return width_ratio

def diff_scatter_map(
        diff, data_coords=None, data_types=None, map=None, ax=None,
        rmse=None, legend=True, return_width_ratio=False,
        filename=None,
        e_prom=None, sigmas=None, moda=None):
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
    #diff_step = 10**get_magnitude(diff_limit)
    diff_step = 5 
    divisions = np.arange(-diff_limit, diff_limit+diff_step, diff_step)
    ticks = np.arange(-diff_limit, diff_limit+diff_step, 2*diff_step)
    bins = len(divisions) - 1
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
    cbar.set_label('Diff [mW/m²]', rotation=90, labelpad=-20)
    cbar.set_ticks([])
    # Colorbar histogram
    hist_ax = divider.append_axes('right', '30%', pad='0%')
    plt.sca(ax)
    N, bins, patches = hist_ax.hist(
        diff, bins=divisions,
        orientation='horizontal')
    hist_ax.set_yticks(ticks)
    hist_ax.set_ylim([-diff_limit, diff_limit])
    hist_ax.yaxis.tick_right()
    if sigmas is not None:
        hist_ax.axhline(y=sigmas.n_1_sigma)
        hist_ax.text(15,sigmas.n_1_sigma-0.5*diff_step,r'-$\sigma$',size='small')
        hist_ax.axhline(y=sigmas.p_1_sigma)
        hist_ax.text(15,sigmas.p_1_sigma+0.2*diff_step,r'+$\sigma$',size='small')
        hist_ax.axhline(y=sigmas.n_2_sigma)
        hist_ax.text(15,sigmas.n_2_sigma-0.5*diff_step,r'-2$\sigma$',size='small')
        hist_ax.axhline(y=sigmas.p_2_sigma)
        hist_ax.text(15,sigmas.p_2_sigma+0.2*diff_step,r'+2$\sigma$',size='small')
    norm = Normalize(bins.min(), bins.max())
    for bin, patch in zip(bins, patches):
        color = diff_cmap(norm(bin))
        patch.set_facecolor(color)
    #hist_ax.grid(True)
    # Title
    ax.set_title('Model minus Data')
    # Options
    extra_artists=[]
    if moda is not None:
        # MODA
        moda_text = plt.figtext(
            0.4, 0.03, 'MODA: %0.2f' %(moda),
            fontweight='bold')
        extra_artists.append(moda_text)
    if e_prom is not None:
        # MAE
        e_prom_text = plt.figtext(
            0.4, 0, 'ME: %0.2f' %(e_prom),
            fontweight='bold')
        extra_artists.append(e_prom_text)
    if rmse is not None:
        # RMSE
        rmse_text = plt.figtext(
            0.4, -0.03, 'RMSE: %0.2f' %(rmse),
            fontweight='bold')
        extra_artists.append(rmse_text)
    if legend is True:
        legend = ax.legend(bbox_to_anchor=(0.5, -0.03), loc='upper center',
                           ncol=2, bbox_transform=fig.transFigure)
        extra_artists.append(legend)
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(
            filename,
            bbox_extra_artists=extra_artists, bbox_inches='tight')
            #dpi='figure', format='pdf')
    if return_width_ratio:
        width_ratio = 1 + 0.05 + 0.12 + 0.30
        return width_ratio

def multi_map(
        shf=None, data=None, diff=None, data_coords=None, data_types=None,
        rmse=None, topo=True, filename=None,
        e_prom=None, sigmas=None):
    # Gridspec
    fig = plt.figure(figsize=(9,6))
    gs = gridspec.GridSpec(1,2)
    # Axis 1: Surface Heat Flow + Data
    ax1 = fig.add_subplot(gs[0,0])
    map1 = base_map(topo=False)
    cbar_limits = [np.nanmax(shf), np.nanmin(shf)]
    wr1 = heatmap_map(
        shf, cbar_limits=cbar_limits, colormap='afmhot',
        cbar_label='Heat Flow [W/m²]', title='Surface Heat Flow',
        map=map1, ax=ax1, return_width_ratio=True)
    data_scatter_map(
        data, cbar_limits=cbar_limits,
        data_coords=data_coords,
        data_types=data_types,
        map=map1, ax=ax1,
        legend=False, cbar=False)
    # Axis 2: Diff
    ax2 = fig.add_subplot(gs[0,1])
    map2 = base_map(topo=False)
    wr2 = diff_scatter_map(
        diff,
        data_coords=data_coords,
        data_types=data_types,
        map=map2, ax=ax2,
        legend=False,  return_width_ratio=True,
        sigmas=sigmas)
    ax2.set_yticks([])
    # Aspect
    gs.set_width_ratios([wr1,wr2])
    plt.tight_layout()
    # Extra Artists
    extra_artists=[]
    #MEA
    e_prom_text = plt.figtext(
        0.03, -0.02, 'ME: %0.2f' %(e_prom),
        fontweight='bold')
    extra_artists.append(e_prom_text)
    # RMSE
    rmse_text = plt.figtext(
        0.03, 0, 'RMSE: %0.2f' %(rmse),
        fontweight='bold')
    extra_artists.append(rmse_text)
    ## Legend
    legend = plt.legend(
        loc="center left", bbox_to_anchor=(0.2, 0),
        ncol=4, bbox_transform=fig.transFigure)
    extra_artists.append(legend)
    # Options
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(
            filename,
            bbox_extra_artists=(extra_artists), bbox_inches='tight', transparent=True)
            #dpi='figure', format='pdf')
    plt.close()

def rmse_plot(vnames, vaxes, rmses, filename=None):
    # x axis
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
        v = np.linspace(vmin, vmax, 100)
        #plt.contourf(
        #    x_axis, y_axis, rmses.T, v, norm=colors.PowerNorm(gamma=1./2.))
        #plt.pcolormesh(
        #    x_axis, y_axis, rmses.T, norm=colors.PowerNorm(gamma=1./2.))
        x_step = x_axis[1] - x_axis[0]
        y_step = y_axis[1] - y_axis[0]
        plt.imshow(rmses.T, origin='lower', aspect='auto',
            norm=colors.PowerNorm(gamma=1./2.),
            extent=[
                x_axis[0] - x_step/2, x_axis[-1] + x_step/2,
                y_axis[0] - y_step/2, y_axis[-1] + y_step/2])
        #xx, yy = np.meshgrid(x_axis, y_axis)
        #xx = xx.flatten()
        #yy = yy.flatten()
        #rmses_f = rmses.T.flatten()
        #plt.scatter(xx, yy, c=rmses_f, cmap='rainbow')
        plt.colorbar()
        name = filename + '_2D'
    else:
        index = np.arange(len(x_axis))
        plt.plot(index, rmses, '-r', linewidth=1.)
        plt.bar(index, rmses, alpha=.4)
        diff = max(rmses) - min(rmses)
        plt.ylim(min(rmses)-0.2*diff, max(rmses)+0.2*diff)
        plt.xticks(index, x_axis)
        plt.ylabel('RMSE')
        plt.grid(True)
        name = filename 
    plt.tight_layout()
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.png'
        plt.savefig(filename)#,dpi='figure', format='pdf')
        plt.close()

def data_scatter_plot(
        data, data_error, data_types, ishf_models, ishf_labels, filename=None):
    fig, axes = plt.subplots(4)
    markers = np.array(['o', '^', 'p', 's'])
    titles = np.array(['ODP', 'LBH', 'GQ', 'GP'])
    data_types = (np.array(data_types) - 1).astype(int)
    data_markers = markers[data_types]
    colors = cm.rainbow(np.linspace(0, 1, len(ishf_models)))
    for m, t, ax in zip(markers, titles, axes):
        mask = data_markers == m
        m_data = data[mask]
        m_data_error = data_error[mask]
        m_data_axis = np.arange(len(m_data))
        ax.errorbar(
            m_data_axis, m_data, m_data_error, fmt=m, capsize=2, capthick=0.5,
            markersize=2, elinewidth=0.5)
        for ishf_model, ishf_label, c in zip(ishf_models, ishf_labels, colors):
            m_ishf_model = ishf_model[mask]
            ax.plot(
                m_data_axis, m_ishf_model, marker='.',
                color=c, label=ishf_label, linewidth=0.5, markersize=0.5)

        ax.set_title(t)
    plt.tight_layout()
    extra_artists = []
    legend = plt.legend(
        loc='center left', bbox_to_anchor=(1,0.5),
        bbox_transform=fig.transFigure)
    extra_artists.append(legend)
    if filename:
        makedir_from_filename(filename)
        filename = filename + '.pdf'
        plt.savefig(
            filename,
            bbox_extra_artists=(extra_artists), bbox_inches='tight',
            dpi='figure', format='pdf')
