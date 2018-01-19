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
                    topo=True,name='Mapa_Surface_Heat_Flow',rmse=None,datos_rmse_error=None):

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
                    diff = (datos_rmse_error_i - ihf_i)/abs(datos_rmse_error_i)
                else:
                    q_flow = -data_q_i[:,2]*1.e-3
                    ihf_i = interpolated_heat_flow[data_q_i_idxs]
                    diff = (q_flow - ihf_i)/abs(q_flow)
                map.drawlsmask(land_color='0.8', ocean_color='0.8',resolution='l')
                #ihf_i = interpolated_heat_flow[data_q_i_idxs]
                #diff = (q_flow - ihf_i)/abs(q_flow)
                diff_max = np.nanmax(diff)
                diff_min = np.nanmin(diff)
                diff_limit = np.nanmax([abs(diff_max),abs(diff_min)])
                norm = MidPointNorm(midpoint=0,vmin=-1,
                                               vmax=1)
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
            cbar_diff.set_label('Diff', rotation=90, labelpad=-50)
            cbar_diff.set_ticks([-1,-.6,-.2,0,.2,.6,1])
            cbar_diff.set_ticklabels([-1,-.6,-.2,0,.2,.6,1])
            plt.title('%s' %(name))

    plt.tight_layout()
    nombre = "%s_0%s_DIFF" %(name,tmc)
    os.chdir(direTer)
    if not os.path.exists('Mapas'):
        os.makedirs('Mapas')
    os.chdir('Mapas')
    plt.savefig('%s' %(nombre),dpi='figure',format='png')
    os.chdir('../')
    os.chdir('../../')
    plt.close()
    return

def map_q_surface(CS, TM, tmc, data_q):
    longitud1 = data_q[:,0][np.where(data_q[:,-2]==1)]
    latitud1 = data_q[:,1][np.where(data_q[:,-2]==1)]
    longitud2 = data_q[:,0][np.where(data_q[:,-2]==2)]
    latitud2 = data_q[:,1][np.where(data_q[:,-2]==2)]
    longitud3 = data_q[:,0][np.where(data_q[:,-2]==3)]
    latitud3 = data_q[:,1][np.where(data_q[:,-2]==3)]
    longitud4 = data_q[:,0][np.where(data_q[:,-2]==4)]
    latitud4 = data_q[:,1][np.where(data_q[:,-2]==4)]
    q_flow = -data_q[:,2]
    shf_max = np.nanmax([np.nanmax(q_flow)*1e-3,
              np.nanmax(TM.get_surface_heat_flow())])
    shf_min = np.nanmin([np.nanmin(q_flow)*1e-3,
              np.nanmin(TM.get_surface_heat_flow())])
    q_flow_1 = -data_q[:,2][np.where(data_q[:,-2]==1)]*1.e-3
    q_flow_2 = -data_q[:,2][np.where(data_q[:,-2]==2)]*1.e-3
    q_flow_3 = -data_q[:,2][np.where(data_q[:,-2]==3)]*1.e-3
    q_flow_4 = -data_q[:,2][np.where(data_q[:,-2]==4)]*1.e-3
    map = Basemap(llcrnrlon= -80, llcrnrlat= -45, urcrnrlon= -60.0, urcrnrlat= -10.0, epsg= 4326, resolution = 'f')
    #map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
    map.drawparallels(np.arange(-90,90,5), labels=[1,0,0,0], fontsize=7)
    map.drawmeridians(np.arange(-180,180,5), labels=[0,0,0,1], fontsize=7)
    map.etopo()
    map.drawcoastlines(linewidth=0.5)
    #hacer grid y cargar los datos para la paleta de colores del mapa
    mlon1, mlat1 = map(longitud1,latitud1)
    mlon2, mlat2 = map(longitud2,latitud2)
    mlon3, mlat3 = map(longitud3,latitud3)
    mlon4, mlat4 = map(longitud4,latitud4)
    x = np.linspace(map.llcrnrx, map.urcrnrx, CS.get_x_axis().shape[0])
    y = np.linspace(map.llcrnry, map.urcrnry, CS.get_y_axis().shape[0])
    xx, yy = np.meshgrid(x, y)
    q_flowm_1 = ma.masked_invalid(q_flow_1)
    q_flowm_2 = ma.masked_invalid(q_flow_2)
    q_flowm_3 = ma.masked_invalid(q_flow_3)
    q_flowm_4 = ma.masked_invalid(q_flow_4)
    datam = ma.masked_invalid(TM.get_surface_heat_flow())
    M = map.pcolormesh(xx, yy[::-1], datam.T, cmap='afmhot_r', shading= 'gouraud', vmin=shf_min,vmax=shf_max)
    #Graficar datos de Q y barra de color
    #q_ind_1
    plot_q1 = map.scatter(mlon1, mlat1, latlon = True, c = q_flowm_1.T, marker='o', cmap = 'afmhot_r', vmin=shf_min,vmax=shf_max)
    #q_ind_2
    plot_q2 = map.scatter(mlon2, mlat2, latlon = True, c = q_flowm_2.T, marker='^', cmap = 'afmhot_r', vmin=shf_min,vmax=shf_max)
    #q_ind_3
    plot_q3 = map.scatter(mlon3, mlat3, latlon = True, c = q_flowm_3.T, marker='p', cmap = 'afmhot_r', vmin=shf_min,vmax=shf_max)
    #q_ind_4
    plot_q4 = map.scatter(mlon4, mlat4, latlon = True, c = q_flowm_4.T, marker='s', cmap = 'afmhot_r', vmin=shf_min,vmax=shf_max)
    cbar = plt.colorbar(M)
    cbar.set_label('Heat Flow (W/m2)', rotation=90, labelpad=-70)
    plt.title('Surface Heat Flow')
    plt.legend([plot_q1, plot_q2, plot_q3, plot_q4],
               ['ODP Borehole', 'Land Borehole',
                'Geochemical', 'Marine Geophysics'],
               bbox_to_anchor=(-1.0, 0),loc=3)
    plt.tight_layout()
    nombre = "Mapa_Q_0%s" %(tmc)
    if not os.path.exists('Mapas'):
        os.makedirs('Mapas')
    os.chdir('Mapas')
    plt.savefig('%s' %(nombre),dpi='figure',format='png')
    os.chdir('../')
    plt.close()
    return

"""
def get_detachment(CS,GM,MM):
    c = -1
    for i in range(len(CS.get_axes()[1])):
        c = c+1
        mecanico = MM.get_yse()[0][:,i,:]
        d_points = np.where(mecanico==200)
    print(d_points)
        #plt.plot(CS.get_axes()[0], d_points[0])
        #fig.savefig('hola.png')
        #plt.close()
    return
"""
"""
#Modulo para graficar modelo termal y mecanico

#Grafica los perfiles guardados a partir del input_perfiles
import os
import math
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
#from crearinput import ceil,floor,minlon,maxlon
from meccolormap import jet_white_r


#Cargar datos de los perfiles guardados y para cada uno de ellos generar un grafico respectivo

#Parametros para grafico base de cada perfil generado
def grafico_base(lat, array):
    #Importar los datos del perfil
    name = 'Input/perf%d' %(lat)
    #perfil = np.loadtxt(name)
    #Longitudes de comienzo y fin del perfil
    start = perfil[np.where(perfil[:,6]==1)[0][0],0]
    finish = perfil[np.where(np.isnan(perfil[:,3])==False)[0][-1],0]
    #Longitudes de los ejes
    len_xaxis,len_yaxis = finish-start,5.
    fig = plt.figure(figsize=(len_xaxis,len_yaxis))
    #Graficar limites en la litosfera
    plt.plot(perfil[:,0], perfil[:,2], 'r', linewidth=4.0) #LAB / SLAB
    plt.plot(perfil[:,0], perfil[:,3], 'b') #MOHO
    plt.plot(perfil[:,0], perfil[:,4], 'm') #ICD
    plt.plot(perfil[:,0], perfil[:,5], 'k') #TOPOGRAFIA
    #Graficar punto en el que se intersecta SLAB con LAB
    a1 = np.where(perfil[:,6]==2.)[0]
    plt.plot(perfil[a1[0],0], perfil[a1[0],2], 'ob')
    #Definir limites de grafico
    plt.ylim([-200,50])
    plt.xlim([start,finish])
    #Definir ticks de los ejes
    plt.xticks(np.arange(math.floor(start),math.ceil(finish),2))
    plt.yticks(np.arange(50,-200-1,-50))
    #Definir labels de los ejes
    plt.ylabel('Profundidad [km]')
    plt.xlabel('Longitud')
    plt.tight_layout()
    return fig

#Rutina que grafica perfiles del modelo inicial
def graficar_perfiles(mi,mx,delta,formato):
    #Multiplicar por 100 para que el step de np.arange sea numero entero
    mi=mi*100
    mx=mx*100
    delta=delta*100
    lat = np.arange(mi,mx-delta,-delta)
    #Loop que itera latitudes
    for l in lat:
        #Cargar grafico base
        fig = grafico_base(l)
        #Definir titulo
        plt.title('Perfil Termal %d' %(l))
        #Navegar a directorio donde se guardaran figuras perfiles
        os.chdir('Input')
        if not os.path.exists('Figuras'):
            os.makedirs('Figuras')
        os.chdir('Figuras')
        #Guardar figuras en caso de formato png
        if formato == 1:
            fig.savefig('perf%d.png' %(l))
            plt.close()
        #Guardar figuras en caso de formato svg
        elif formato == 2:
            if not os.path.exists('Formato_SVG'):
                os.makedirs('Formato_SVG')
            os.chdir('Formato_SVG')
            fig.savefig('perf%d.svg' %(l), format="svg")
            plt.close()
            os.chdir('../')
        #Volver a directorio principal
        os.chdir('../../')
    return

#Rutina para generar perfiles termales
def graficar_estructura_termal(mi,mx,delta,direTer,formato):
    #Multiplicar por 100 para que el step de np.arange sea numero entero
    mi=mi*100
    mx=mx*100
    delta=delta*100
    lat = np.arange(mi,mx-delta,-delta)
    #Loop que itera latitudes
    for l in lat:
        #Cargar grafico base
        fig = grafico_base(l)
        #Definir titulo
        plt.title('Perfil Termal %d' %(l))
        #Importar los datos del perfil termal
        os.chdir(direTer)
        os.chdir('Perfiles')
        name = 'perf%d' %(l)
        perfil = np.loadtxt(name)
        #Definir ejes maximos del area coloreada
        A1 = perfil[:,0]
        A2 = np.arange(ceil,floor-1,-1)
        #Definir input para funcion pcolor
        [G1, G2] = np.meshgrid(A1,A2)
        G3 = perfil[:,9:]
        G3m = ma.masked_invalid(G3)
        #En caso de formato png usar pcolormesh
        if formato == 1:
            sub = plt.pcolormesh(G1,G2,G3m.T,cmap=cm.coolwarm, shading='gouraud')
        #En caso de formato svg usar pcolor (pcolormesh ralentiza mucho el rendering)
        elif formato == 2:
            sub = sub = plt.pcolor(G1,G2,G3m.T,cmap=cm.coolwarm)
        #Definir limites de temperatura visible en grafico
        tmin = 0
        tmax = 1300
        plt.clim(tmin,tmax)
        #Graficar barra de color
        cbar = plt.colorbar(sub, aspect='20')
        #Graficar label de barra de color
        cbar.set_label('Temperatura', rotation=90, labelpad=-60)
        #Graficar lineas de contorno
        numlabels = 6
        sub1 = plt.contour(G1,G2,G3m.T,numlabels,colors='k')
        #Generar posiciones par los labels de contorno
        Tstart = np.where(perfil[:,6]==1)[0][0]
        Tlons = np.where(np.isnan(perfil[Tstart:,3])==False)[0]+Tstart
        midlonidx = Tlons[-1] - len(Tlons)/10
        labelxpos = perfil[midlonidx,0]
        labelstep = round((tmax/100) / numlabels)*100
        manual_locations = []
        for i in range(1,numlabels+1):
            x = ma.masked_invalid(perfil[midlonidx,9:]/labelstep)
            idx = np.abs(x-i).argmin()
            manual_locations.append([labelxpos,A2[idx]]) #Posiciones de labels
        #Graficar labels de contorno
        plt.clabel(sub1, fontsize=10, fmt='%.0f', rightside_up=True, manual=manual_locations)
        #Navegar a directorio donde se guardaran figuras de perfiles
        os.chdir('../')
        if not os.path.exists('Figuras'):
            os.makedirs('Figuras')
        os.chdir('Figuras')
        #Guardar figuras en caso de formato png
        if formato == 1:
            fig.savefig('perf%d.png' %(l))
            plt.close()
        #Guardar figuras en caso de formato svg
        elif formato == 2:
            if not os.path.exists('Formato_SVG'):
                os.makedirs('Formato_SVG')
            os.chdir('Formato_SVG')
            fig.savefig('perf%d.svg' %(l), format="svg")
            plt.close()
            os.chdir('../')
        #Volver a directorio principal
        os.chdir('../../../')
    return

#Rutina para generar perfiles mecanicos
def graficar_estructura_mecanica(mi,mx,delta,direTerMec,formato):
    #Multiplicar por 100 para que el step de np.arange sea numero entero
    mi=mi*100
    mx=mx*100
    delta=delta*100
    lat = np.arange(mi,mx-delta,-delta)
    #Loop que itera latitudes
    for l in lat:
        #Cargar grafico base
        fig = grafico_base(l)
        #Definir titulo
        plt.title('Perfil Termomecanico %d' %(l))
        #Importar los datos del perfil mecanico
        os.chdir(direTerMec)
        os.chdir('Perfiles')
        name = 'perfT%d' %(l)
        perfil = np.loadtxt(name)
        #Definir ejes maximos del area coloreada
        A1 = perfil[:,0]
        A2 = np.arange(ceil,floor-1,-1)
        #Definir input para funcion pcolor
        [G1, G2] = np.meshgrid(A1,A2)
        G3 = perfil[:,7:]
        G3m = ma.masked_invalid(G3)
        #En caso de formato png usar pcolormesh
        if formato == 1:
            sub = plt.pcolormesh(G1,G2,G3m.T,cmap=jet_white_r, shading='gouraud')
        #En caso de formato svg usar pcolor (pcolormesh ralentiza mucho el rendering)
        elif formato == 2:
            sub = plt.pcolor(G1,G2,G3m.T,cmap=jet_white_r)
        #Definir limites de Yield Strength visible en el grafico
        plt.clim(0,200)
        #Graficar barra de color
        cbar = plt.colorbar(sub, aspect='20')
        #Graficar label de barra de color
        cbar.set_label('Yield Strenght', rotation=90, labelpad=-50)
        #Navegar a directorio donde se guardaran los perfiles
        os.chdir('../')
        if not os.path.exists('Figuras'):
            os.makedirs('Figuras')
        os.chdir('Figuras')
        #Guardar perfiles en caso de formato png
        if formato == 1:
            fig.savefig('perf%d.png' %(l))
            plt.close()
        #Guardar perfiles en caso de formato svg
        elif formato == 2:
            if not os.path.exists('Formato_SVG'):
                os.makedirs('Formato_SVG')
            os.chdir('Formato_SVG')
            fig.savefig('perf%d.svg' %(l), format="svg")
            plt.close()
            os.chdir('../')
        #Volver a directorio principal
        os.chdir('../../../../')
    return

#Comparar dos arrays con valores nan
def nanarray_equal(a,b):
    try:
        np.testing.assert_equal(a,b)
    except AssertionError:
        return False
    return True
    #Otro metodo es: ((a == b) | (numpy.isnan(a) & numpy.isnan(b))).all()

#Rutina para mapas
def graficar_mapa(data,map_id):
    map = Basemap(projection='merc', resolution='f',
              llcrnrlon= -79.8,
              llcrnrlat= -44.8,
              urcrnrlon= -60.0,
              urcrnrlat= -10.0,)
    map.drawcoastlines()
    map.etopo(zorder=0)
    map.drawparallels(np.arange(-90,90,3), labels=[1,0,0,0])
    map.drawmeridians(np.arange(-180,180,5), labels=[0,0,0,1])
    #hacer grid y cargar los datos para la paleta de colores del mapa
    x = np.linspace(0, map.urcrnrx, data.shape[0])
    y = np.linspace(0, map.urcrnry, data.shape[1])
    xx, yy = np.meshgrid(x, y)
    datam = ma.masked_invalid(data)
    if map_id == 1:
        #print "T"
        M = map.pcolormesh(xx, yy, datam.T, cmap='coolwarm', shading='gouraud')
        plt.clim(0,1300)
        #Graficar barra de color
        cbar = plt.colorbar(M, aspect='20')
        cbar.set_label('Temperatura (C)', rotation=90, labelpad=-70)
        nombre = "Mapa_T"
    elif map_id == 2:
        #print "Q"
        M = map.pcolormesh(xx, yy[::-1], datam.T, cmap='afmhot_r', shading='gouraud')
        #Graficar barra de color
        cbar = plt.colorbar(M, aspect='20')
        cbar.set_label('Flujo de Calor (W)', rotation=90, labelpad=-80)
        nombre = "Mapa_Q"
    elif map_id == 3:
        #print "YS"
        M = map.pcolormesh(xx, yy, datam.T, cmap=jet_white_r, shading='gouraud')
        plt.clim(0,200)
        #Graficar barra de color
        cbar = plt.colorbar(M, aspect='20')
        cbar.set_label('Yield Strength Tension', rotation=90, labelpad=-60)
        nombre = "Mapa_YS"
    elif map_id == 4:
        #print "YSi"
        M = map.pcolormesh(xx, yy, datam.T, cmap=jet_white_r, shading='gouraud')
        plt.clim(0,100000)
        #Graficar barra de color
        cbar = plt.colorbar(M, aspect='20')
        cbar.set_label('Yield Strenght Integrado', rotation=90, labelpad=-80)
        nombre = "Mapa_YSi"
    elif map_id == 5:
        #print "EET"
        M = map.pcolormesh(xx, yy, datam.T, cmap=jet_white_r, shading='gouraud')
        plt.clim(0,100)
        #Graficar barra de color
        cbar = plt.colorbar(M, aspect='20')
        cbar.set_label('Effective Elastic Thickness (km)', rotation=90, labelpad=-60)
        nombre = "Mapa_EET"
    #Guardar mapas en formato png
    os.chdir('Output')
    if not os.path.exists('Mapas'):
        os.makedirs('Mapas')
    os.chdir('Mapas')
    plt.savefig('%s' %(nombre))
    plt.close()
    os.chdir('../../')
    return
def graficar_yse(lat_i,lat_f,lon_i,lon_f,xYse,yYse,delta,EETpoints,ICD,MOHO,max_stress):
    lat_i = round((round(lat_i*5)/5)*100)
    lat_f = round((round(lat_f*5)/5)*100)
    lon_i = round((round(lon_i*5)/5)*100)
    lon_f = round((round(lon_f*5)/5)*100)
    delta = delta*100
    lonlat = np.rint(np.loadtxt('Input/data/lonlat.txt')*100)
    for lat in np.arange(lat_i,lat_f-delta,-delta):
        lat_index = np.where(lonlat[:,1]==lat)[0]
        lonlat1 = lonlat[lat_index,:]
        for lon in np.arange(lon_i,lon_f-20,-20):
            lon_index = np.where(lonlat1[:,0]==lon)[0]
            if not np.where((EETpoints[:,0]==lat) & (EETpoints[:,1]==lon))[0]:
                EETpoints_current = []
            else:
                EETpoints_current = EETpoints[np.where((EETpoints[:,0]==lat) & (EETpoints[:,1]==lon))[0],2][0]
            ICD_current = ICD[np.where((ICD[:,0]==lat) & (ICD[:,1]==lon))[0],2]
            MOHO_current = MOHO[np.where((MOHO[:,0]==lat) & (MOHO[:,1]==lon))[0],2]
            if not lat_index.any() or not lon_index.any():
                pass
            else:
                index = lat_index[0] + lon_index[0]
                plt.plot(xYse[index],yYse[index])
                plt.xlabel('Stress')
                plt.ylabel('Profundidad')
                plt.title('Yield Strenght Envelope')
                for k in range(len(EETpoints_current)):
                    plt.plot(y=EETpoints_current[k], marker='o', color='r',label='T%d' %(k+1))
                    #if k == 0:
                        #plt.axhline(y=EETpoints_current[k],color='r',label='T%d' %(k+1))
                    #if k == 1:
                        #plt.axhline(y=EETpoints_current[k],color='b',label='T%d' %(k+1))
                    #else:
                        #plt.axhline(y=EETpoints_current[k],color='k',label='T%d' %(k+1))
                plt.axhline(y=ICD_current,color='g',label='ICD')
                plt.axhline(y=MOHO_current,color='y',label='MOHO')
                plt.axvline(x=max_stress,color='c',label='MAXSTRESS')
                plt.legend()
                os.chdir('Output')
                if not os.path.exists('YSE'):
                    os.makedirs('YSE')
                os.chdir('YSE')
                plt.savefig('YSE-%dlat-%dlon' %(lat,lon))
                plt.close()
                os.chdir('../../')
    return
"""
