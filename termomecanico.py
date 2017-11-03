#!/usr/bin/python
# Rutina Maestra

# Cargar modulos para generar el modelo Termomecanico
import setup
import compute
import numpy as np
#import plot
from utils import DotDict
import os
import sys

# Leer archivos con las variables, para el modelo termal y mecanico
print('Leyendo variables...')
Temp = DotDict(setup.readVars('VarTermal.txt'))
Meca = DotDict(setup.readVars('VarMecanico.txt'))
Exec = DotDict(setup.readVars('VarExec.txt'))
direTer, direTerMec = setup.makeDirs(Exec.temcaso, Exec.meccaso)
Rhe = DotDict(setup.read_rheo('data/Rhe_Param.dat'))
# Variables
#graf = Exec[0] #geometrias
#tmi = Exec[2]
#tmx = Exec[3]
#tdelta = Exec[4]
#xt1 = Exec[5] #datos
#xt2 = Exec[6] #graficos
#xt3 = Exec[7] #mapa
#xt4 = Exec[8] #topo
#xt5 = Exec[9] #prof
#mmi = Exec[10]
#mmx = Exec[12]
#mdelta = Exec[13]
#xm1 = Exec[14] #datos
#xm2 = Exec[15] #graficos
#xm3 = Exec[16] #mapa
#xm4 = Exec[17] #topo
#xm5 = Exec[18] #topofinal
#xm6 = Exec[19] #prof
#xm7 = Exec[20] #yse
#minlat = Exec[21]
#maxlat = Exec[22]
#minlon = Exec[23]
#maxlon = Exec[24]
#max_stress = Meca[7]
#Restricciones de mi, mx y delta
if (Exec.tmi-Exec.mmi) <= 0:
    Exec.mmi = Exec.tmi
if (Exec.tmx-Exec.mmx) >= 0:
    Exec.mmx = Exec.tmx
if Exec.tdelta > Exec.mdelta:
    Exec.mdelta = Exec.tdelta
print('...OK')

#Graficar perfiles de Input
#if graf == 1 or graf == 2:
#    print 'Graficando geometrias...'
#    TMgraficar.graficar_perfiles(tmi,tmx,tdelta,xt2)
#    print '...OK'

gm_data = np.loadtxt('data/Modelo.dat')
ta_data = np.loadtxt('data/PuntosFosaEdad.dat')
areas = np.loadtxt('data/areas.dat')

GM, TM, MM = compute.compute(gm_data, ta_data, Rhe, areas, Temp, Meca)
"""
# Estimar el modelo termal
if xt1 == 1:
    print 'Generando modelo termal...'
    T = termal.modelo_termal(tmi,tmx,tdelta,Temp,col,direTer)
    print '...OK'
# Generar figuras del modelo Termal
if xt2 == 1 or xt2 == 2:
    print 'Graficando perfiles termales...'
    TMgraficar.graficar_estructura_termal(tmi,tmx,tdelta,direTer,xt2)
    print '...OK'
# Estimar el modelo Mecanico basado en la estructura Termal
if xm1 == 1:
    print 'Generando modelo mecanico...'
    M = mecanico.modelo_mecanico(mmi,mmx,mdelta,Meca,col,direTerMec)
    print '...OK'
# Generar figuras del modelo Mecanico
if xm2 == 1 or xm2 == 2:
    print 'Graficando perfiles termomecanicos...'
    TMgraficar.graficar_estructura_mecanica(mmi,mmx,mdelta,direTerMec,xm2)
    print '...OK'
# Generar mapas del modelo Termal
if xt3 != 0:
    print 'Procesando datos termales...'
    Q, Tz = selector.selector_termal(direTer,tdelta,xt3,xt4,xt5)
    print '...OK'
    if xt3 == 1 or xt3 == 3:
        map_id = 1
        print 'Graficando mapa de temperaturas...'
        TMgraficar.graficar_mapa(Tz,map_id)
        print '...OK'
    if xt3 == 2 or xt3 == 3:
        map_id = 2
        print 'Graficando mapa de calor...'
        TMgraficar.graficar_mapa(Q,map_id)
        print '...OK'
#Generar mapas o YSE del modelo Mecanico
if xm3 != 0 or xm7 == 1:
    print 'Procesando datos mecanicos...'
    Yt,Yti,EET,xYset,yYset,EETpoints,ICD,MOHO = selector.selector_mecanico(direTerMec,mdelta,xm3,xm4,xm5,xm6,xm7,max_stress)
    print '...OK'
    if xm3 == 1 or xm3 == 4:
        map_id = 3
        print 'Graficando mapa de Yield Strength...'
        TMgraficar.graficar_mapa(Yt,map_id)
        print '...OK'
    if xm3 == 2 or xm3 == 4:
        map_id = 4
        print 'Graficando mapa de Yield Strength Integrado...'
        TMgraficar.graficar_mapa(Yti,map_id)
        print '...OK'
    if xm3 == 3 or xm3 == 4:
        map_id = 5
        print 'Graficando mapa de EET...'
        TMgraficar.graficar_mapa(EET,map_id)
        print '...OK'
if xm7 == 1 and tdelta == 0.2 and mdelta == 0.2:
    print 'Graficando Yield Strength Envelope...'
    TMgraficar.graficar_yse(minlat,maxlat,minlon,maxlon,xYset,yYset,mdelta,EETpoints,ICD,MOHO,max_stress)
    print '...OK'
elif xm7 ==1 and tdelta != 0.2 or mdelta != 0.2:
    print 'Para graficar Yield Strength envelope, es necesario un modelo'
    print 'con resolucion de 0.2 grados en la latitud y la longitud'
    print '(tdelta y mdelta deben ser iguales a 0.2)'

#    if not os.listdir(direTer):
#        print 'Error: debe generar un modelo Termal antes de usar el Mecanico'
"""
