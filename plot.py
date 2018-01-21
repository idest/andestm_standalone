import numpy as np
import numpy.ma as ma
import setup
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from meccolormap import jet_white_r
from diffcolormap import diff_cmap
from mpl_toolkits.basemap import Basemap
from utils import MidPointNorm

def plot_thermal(y, z, D, CS, GM, TM):
    c = -1
    for i in range(len(CS.get_axes()[1])-1):
        c = c+1
        lat = CS.get_axes()[1][c]
        temp = TM.get_geotherm()[:,i,:]
        start = CS.get_x_axis()[np.where(np.isnan(GM.get_topo()
                .mask_irrelevant()[:,i])==False)[0][0]]
        finish = CS.get_x_axis()[np.where(np.isnan(GM.get_topo()
                 .mask_irrelevant()[:,i])==False)[0][-1]]
        fig = plt.figure(figsize=(finish-start,5))
        plt.xlim(start,finish)
        plt.ylim(-150,10)
        plt.xticks(np.arange(math.floor(start),math.ceil(finish),2))
        plt.yticks(np.arange(10,-150-1,-10))
        plt.ylabel('Profundidad [km]')
        plt.xlabel('Longitud')
        plt.tight_layout()
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[0][:,i], 'k') #TOPOGRAFIA
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[1][:,i], 'w') #ICD
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[2][:,i], 'g') #MOHO
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[3][:,i], 'r') #SLAB/LAB
        A1, A2 = np.meshgrid(y,z)
        temp_mask = ma.masked_invalid(temp)
        sub = plt.pcolormesh(A1,A2,temp_mask.T,cmap=cm.coolwarm,
                             shading='gouraud')
        tmin = 0
        tmax = 1300
        plt.clim(tmin,tmax)
        #Graficar barra de color
        cbar = plt.colorbar(sub, aspect='20')
        #Graficar label de barra de color
        cbar.set_label('Temperatura', rotation=90, labelpad=-60)
        plt.title('Perfil Termal %s' %(lat))
        #Graficar lineas de contorno
        numlabels = 6
        cont = plt.contour(A1, A2, temp_mask.T, numlabels, colors = 'k')
        plt.clabel(cont, fontsize=9, fmt='%.0f', rightside_up=True)
        #Guardar figura
        if not os.path.exists('perfiles_termales'):
            os.makedirs('perfiles_termales')
        os.chdir('perfiles_termales')
        fig.savefig('perf%s.png' %(lat))
        os.chdir('../')
        plt.close()
    return


def plot_mec(y, z, D, CS, GM, MM):
    c = -1
    for i in range(len(CS.get_axes()[1])-1):
        c = c+1
        lat = CS.get_axes()[1][c]
        meca = MM.get_yse()[0][:,i,:]
        start = CS.get_x_axis()[np.where(np.isnan(GM.get_topo()
                .mask_irrelevant()[:,i])==False)[0][0]]
        finish = CS.get_x_axis()[np.where(np.isnan(GM.get_topo()
                 .mask_irrelevant()[:,i])==False)[0][-1]]
        fig = plt.figure(figsize=(finish-start,5))
        plt.xlim(start,finish)
        plt.ylim(-150,10)
        plt.xticks(np.arange(math.floor(start),math.ceil(finish),2))
        plt.yticks(np.arange(10,-150-1,-10))
        plt.ylabel('Profundidad [km]')
        plt.xlabel('Longitud')
        plt.tight_layout()
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[0][:,i], 'k') #TOPOGRAFIA
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[1][:,i], 'w') #ICD
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[2][:,i], 'g') #MOHO
        plt.plot(CS.get_axes()[0],GM.get_boundaries()[3][:,i], 'r') #SLAB/LAB
        A1, A2 = np.meshgrid(y,z)
        meca_mask = ma.masked_invalid(meca)
        #plot modelo termomecanico
        sub = plt.pcolormesh(A1, A2, meca_mask.T, cmap=jet_white_r,
                             shading ='gouraud')
        plt.clim(0,200)
        #plot colorbar
        cbar = plt.colorbar(sub, aspect='20')
        cbar.set_label('Yield Strenght', rotation=90, labelpad=-50)
        plt.title('Perfil Termomecanico %s' %(lat))
        #Guardar perfil
        if not os.path.exists('perfiles_termomecanicos'):
            os.makedirs('perfiles_termomecanicos')
        os.chdir('perfiles_termomecanicos')
        fig.savefig('perf%s.png' %(lat))
        os.chdir('../')
        plt.close()
    return

def base_map(topo=True):
    map = Basemap(llcrnrlon= -80,llcrnrlat= -45,
                  urcrnrlon= -60.0,urcrnrlat= -10.0,
                  epsg= 4326, resolution = 'l')
    #map.arcgisimage(service='ESRI_Imagery_World_2D',xpixels=2000,verbose=True)
    map.drawparallels(np.arange(-90,90,5), labels=[1,0,0,0], fontsize=7)
    map.drawmeridians(np.arange(-180,180,5), labels=[0,0,0,1], fontsize=7)
    if topo is True:
        map.etopo()
    map.drawcoastlines(linewidth=0.5)
    return map


def map_q_surface_2(x_axis, y_axis,tmc,direTer,surface_heat_flow=None,
                    data_q=None,data_cmap=None,interpolated_heat_flow=None,
                    topo=True,name='Mapa_Surface_Heat_Flow',rmse=None,
                    datos_rmse_error=None):

    map = base_map(topo=topo)
    x = np.linspace(map.llcrnrx, map.urcrnrx, x_axis.shape[0])
    y = np.linspace(map.llcrnry, map.urcrnry, y_axis.shape[0])
    xx, yy = np.meshgrid(x, y)

    if surface_heat_flow is not None:
        shf_max = np.nanmax(surface_heat_flow)
        shf_min = np.nanmin(surface_heat_flow)
        heat_cbar_max = shf_max
        heat_cbar_min = shf_min
        #if data_q is not None and data_cmap=='heat_flow':
            #q_flow = -data_q[:,2]*1.e-3
            #q_flow_max = np.nanmax(q_flow)
            #q_flow_min = np.nanmin(q_flow)
            #heat_cbar_max = np.nanmax([shf_max, q_flow_max])
            #heat_cbar_min = np.nanmin([shf_min, q_flow_min])
        map.drawlsmask(land_color='0.8', ocean_color='b',resolution='h')
        datam = ma.masked_invalid(surface_heat_flow)
        M = map.pcolormesh(xx,yy[::-1],datam.T,cmap='afmhot_r',shading='gouraud',
                           vmin=heat_cbar_min,vmax=heat_cbar_max)
        cbar = plt.colorbar(M)
        plt.annotate('rmse = %0.7s' %(rmse), xy=(-.9,.3), xycoords='axes fraction')
        cbar.set_label('Heat Flow (W/m2)', rotation=90, labelpad=-60)
        plt.title('Surface Heat Flow')

    if data_q is not None:
        color_method = None
        q_plt = {}
        data_q_types=[1,2,3,4]
        data_q_types_markers=['o', '^', 'p', 's']
        for i in range(len(data_q_types)):
            data_q_i_idxs = np.where(data_q[:,-2]==i+1)
            data_q_i = data_q[data_q_i_idxs]
            longitude = data_q_i[:,0]
            latitude = data_q_i[:,1]
            m_lon, m_lat = map(longitude,latitude)
            if data_cmap == 'heat_flow':
                q_flow = -data_q_i[:,2]*1.e-3
                q_flowm = ma.masked_invalid(q_flow)
                q_plt[str(i+1)] = map.scatter(m_lon,m_lat,latlon=True,
                                              c=q_flowm,
                                              marker=data_q_types_markers[i],
                                              cmap='afmhot_r',
                                              vmin=heat_cbar_min,
                                              vmax=heat_cbar_max)
            elif data_cmap == 'diff':
                plt.annotate('rmse = %0.7s' %(rmse), xy=(-.9,.3), xycoords='axes fraction')
                if datos_rmse_error is not None:
                    ihf_i = interpolated_heat_flow[data_q_i_idxs]
                    datos_rmse_error_i = datos_rmse_error[data_q_i_idxs]
                    diff = (ihf_i - datos_rmse_error_i)#/abs(datos_rmse_error_i)
                else:
                    q_flow = -data_q_i[:,2]*1.e-3
                    ihf_i = interpolated_heat_flow[data_q_i_idxs]
                    diff = (ihf_i - q_flow)#/abs(q_flow)
                map.drawlsmask(land_color='0.8', ocean_color='0.8',resolution='l')
                ihf_i = interpolated_heat_flow[data_q_i_idxs]
                #diff = (q_flow - ihf_i)/abs(q_flow)
                diff_max = np.nanmax(diff)
                diff_min = np.nanmin(diff)
                diff_limit = np.nanmax([abs(diff_max),abs(diff_min)])
                norm = MidPointNorm(midpoint=0,vmin=-0.1,
                                               vmax=0.1)
                q_plt[str(i+1)] = map.scatter(m_lon,m_lat,latlon=True,
                                              c=diff, cmap=diff_cmap,
                                              norm=norm,
                                              marker=data_q_types_markers[i],
                                              edgecolors='k')
            else:
                q_plt[str(i+1)] = map.scatter(m_lon,m_lat,latlon=True,
                                              marker=data_q_types_markers[i])

        plt.legend([q_plt['1'], q_plt['2'], q_plt['3'], q_plt['4']],
                    ['ODP Borehole', 'Land Borehole',
                     'Geochemical', 'Marine Geophysics'],
                     bbox_to_anchor=(-1.0, 0),loc=3)

        if data_cmap == 'diff':
            cbar_diff = plt.colorbar(q_plt['1'])
            cbar_diff.set_label('Diff', rotation=90, labelpad=-60)
            cbar_diff.set_ticks([-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1])
            cbar_diff.set_ticklabels([-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1])
            plt.title('%s' %(name))

    plt.tight_layout()
    nombre = "%s_0%s_DIFF" %(name,tmc)
    os.chdir(direTer)
    if not os.path.exists('Mapas'):
        os.makedirs('Mapas')
    os.chdir('Mapas')
    plt.savefig('%s' %(nombre),dpi='figure',format='pdf')
    os.chdir('../')
    os.chdir('../../')
    plt.close()
    return

def plot_diffs(model_heat_flow, data_q, data_error):
    plt.figure(figsize=(15,30))
    data_q_types=[1,2,3,4]
    data_q_types_markers=['o', '^', 'p', 's']
    data_q_titles = ['ODP', 'LBH', 'GQ', 'GF']
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4)
    data_q_axes = [ax1, ax2, ax3, ax4]
    for i in range(len(data_q_types)):
        data_q_i_idxs = np.where(data_q[:,-2]==i+1)
        data_q_i = data_q[data_q_i_idxs]
        data_error_i = data_error[data_q_i_idxs]*1e-3
        model_heat_flow_i = model_heat_flow[data_q_i_idxs]
        q_flow = -data_q_i[:,2]*1.e-3
        ax = data_q_axes[i]
        x = np.arange(len(q_flow))
        ax.errorbar(x, q_flow, yerr=data_error_i, fmt='r.', capsize=2, capthick=0.5,markersize=2, elinewidth=0.5)
        ax.plot(x, model_heat_flow_i, 'k.',markersize=2)
        ax.set_ylim(-.15,0)
        ax.set_yticks(np.arange(0,-.151,-.05))
        ax.set_title('%s Data Map'%(data_q_titles[i]))
    plt.tight_layout()

    #data_error = (data_error*1e-3)
    #x = np.arange(len(q_flow))
    #plt.plot(x, model_heat_flow, 'k.')
    #plt.plot(x, q_flow, 'r.')
    #plt.errorbar(x, q_flow, yerr=0.01, fmt='r.', capsize=2)
    #plt.plot(x, data_error, '.')
    #print(os.getcwd())
    if not os.path.exists('Graficos'):
        os.makedirs('Graficos')
    os.chdir('Graficos')
    plt.savefig('diffs',format='pdf')
    os.chdir('../')
    plt.close('all')

# def plot_scatter(model_heat_flow,data_q)
#     plt.figure(figsize=(15,15))
#     data_q_types=[1,2,3,4]
#     data_q_types_markers=['o', '^', 'p', 's']
#     data_q_titles = ['ODP', 'LBH', 'GQ', 'GF']
#     for i in range(len(data_q_types)):
#         data_q_i_idxs = np.where(data_q[:,-2]==i+1)
#         data_q_i = data_q[data_q_i_idxs]
#         ihf_i = interpolated_heat_flow[data_q_i_idxs]
#         plt.scatter(ihf_i,data_q_i[:,2],cmap=)


