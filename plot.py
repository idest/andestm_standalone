import numpy as np
import numpy.ma as ma
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from meccolormap import jet_white_r

def plot_thermal(y, z, D, GM, TM):
    c = -1
    for i in range(len(D.get_axes()[1])-1):
        c = c+1
        lat = D.get_axes()[1][c]
        temp = TM.get_geotherm()[:,i,:]
        start = D.get_axes()[0][np.where(np.isnan(GM.get_boundaries()[0][:,i])==False)[0][0]]
        finish = D.get_axes()[0][np.where(np.isnan(GM.get_boundaries()[2][:,i])==False)[0][-1]] 
        fig = plt.figure(figsize=(finish-start,5))
        plt.xlim(start,finish)
        plt.ylim(-150,10)
        plt.xticks(np.arange(math.floor(start),math.ceil(finish),2))
        plt.yticks(np.arange(10,-150-1,-10))
        plt.ylabel('Profundidad [km]')
        plt.xlabel('Longitud')
        plt.plot(D.get_axes()[0], GM.get_boundaries()[0][:,i], 'k') #TOPOGRAFIA
        plt.plot(D.get_axes()[0], GM.get_boundaries()[1][:,i], 'w') #ICD
        plt.plot(D.get_axes()[0], GM.get_boundaries()[2][:,i], 'g') #MOHO
        plt.plot(D.get_axes()[0], GM.get_boundaries()[3][:,i], 'r') #SLAB/LAB
        plt.tight_layout()
        A1, A2 = np.meshgrid(y,z)
        temp_mask = ma.masked_invalid(temp)
        sub = plt.pcolormesh(A1,A2,temp_mask.T,cmap=cm.coolwarm, shading='gouraud')
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

def plot_mec(y, z, D, GM, MM):
    c = -1
    for i in range(len(D.get_axes()[1])):
        c = c+1
        lat = D.get_axes()[1][c]
        meca = MM.get_yse()[0][:,i,:]
        start = D.get_axes()[0][np.where(np.isnan(GM.get_boundaries()[0][:,i])==False)[0][0]]
        finish = D.get_axes()[0][np.where(np.isnan(GM.get_boundaries()[2][:,i])==False)[0][-1]] 
        fig = plt.figure(figsize=(finish-start,5))
        plt.xlim(start,finish)
        plt.ylim(-150,10)
        plt.xticks(np.arange(math.floor(start),math.ceil(finish),2))
        plt.yticks(np.arange(10,-150-1,-10))
        plt.ylabel('Profundidad [km]')
        plt.xlabel('Longitud')
        plt.plot(D.get_axes()[0], GM.get_boundaries()[0][:,i], 'k') #TOPOGRAFIA
        plt.plot(D.get_axes()[0], GM.get_boundaries()[1][:,i], 'w') #ICD
        plt.plot(D.get_axes()[0], GM.get_boundaries()[2][:,i], 'g') #MOHO
        plt.plot(D.get_axes()[0], GM.get_boundaries()[3][:,i], 'r') #SLAB/LAB
        plt.tight_layout()
        A1, A2 = np.meshgrid(y,z)
        meca_mask = ma.masked_invalid(meca)
        #plot modelo termomecanico
        sub = plt.pcolormesh(A1, A2, meca_mask.T, cmap=jet_white_r, shading ='gouraud')
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

def get_detachment(D,GM,MM):
    c = -1
    for i in range(len(D.get_axes()[1])):
        c = c+1
        mecanico = MM.get_yse()[0][:,i,:]
        d_points = np.where(mecanico==200)
        print(d_points)
        #plt.plot(D.get_axes()[0], d_points[0])
        #fig.savefig('hola.png')
        #plt.close()
    return

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
        M = map.pcolormesh(xx, yy, datam.T, cmap='afmhot_r', shading='gouraud')
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
