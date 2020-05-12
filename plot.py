import numpy as np
import numpy.ma as ma
import setup
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from meccolormap import jet_white_r
import os


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
        fig = plt.figure(figsize=(finish-start,10))
        plt.xlim(start,finish)
        plt.ylim(-150,10)
        plt.xticks(np.arange(math.floor(start),math.ceil(finish),2))
        plt.yticks(np.arange(10,-150-1,-10))
        plt.ylabel('Profundidad [km]')
        plt.xlabel('Longitud')
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
        cbar.set_label('Temperatura', rotation=90)
        plt.title('Perfil Termal {:.1f}' .format(lat))
        #Graficar lineas de contorno
        numlabels = 6
        cont = plt.contour(A1, A2, temp_mask.T, numlabels, colors = 'k')
        plt.clabel(cont, fontsize=9, fmt='%.0f', rightside_up=True)
        plt.tight_layout()
        #Guardar figura
        if not os.path.exists('perfiles_termales'):
            os.makedirs('perfiles_termales')
        os.chdir('perfiles_termales')
        fig.savefig('perf{:.1f}.png' .format(lat))
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
        fig = plt.figure(figsize=(finish-start, 10))
        plt.xlim(start,finish)
        plt.ylim(-150,10)
        plt.xticks(np.arange(math.floor(start),math.ceil(finish),2))
        plt.yticks(np.arange(10,-150-1,-10))
        plt.ylabel('Profundidad [km]')
        plt.xlabel('Longitud')
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
        cbar.set_label('Yield Strenght', rotation=90)
        plt.title('Perfil Termomecanico {:.1f}' .format(lat))
        plt.tight_layout()
        #Guardar perfil
        if not os.path.exists('perfiles_termomecanicos'):
            os.makedirs('perfiles_termomecanicos')
        os.chdir('perfiles_termomecanicos')
        fig.savefig('perf{:.1f}.png' .format(lat))
        os.chdir('../')
        plt.close()
    return