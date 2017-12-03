# Lee de variables y da instrucciones:
# El formato del archivo no se debe cambiar y no debe contener lines vacias
# al final.

import numpy as np
import os
import shutil
from utils import DotDict

#Leer archivo de instrucciones
def readVars(name):

    # Abrir archivo
    f = open(name)

    dic = {}


    # Loop para leer las variables
    line = f.readline()
    while line:

        # ignorar las lineas que empiezan con #
        if line[0] != '#':

            words = line.split()
            Vval = words[0]
            Vname = words[1]
            #dic[str.strip(Vname)] = Vnum
            #Vstr = str(Vname.strip())
            #dic[Vstr] = Vval
            exec('{} = {}'.format('dic[Vname]', Vval))
        # Leer siguiente linea
        line = f.readline()

    f.close()

    # Retorna la matriz con todas las variables

    return DotDict(dic)

def makeDirs(temcaso, meccaso):

    # Crea carpeta para los outputs con el nombre del caso (temcaso). El cual se
    # usara para inicializar el Modelo Mecanico.

    if not os.path.exists('Output'):
        os.makedirs('Output')

    direTer = 'Output/%s_Termal' %(temcaso)
    if not os.path.exists(direTer):
        os.makedirs(direTer)

    # Copiar los archivos de input a la carpeta con los Resultados
    shutil.copy( 'VarTermal.txt'  , direTer )

    # Generar nuevo output para el modelo Mecanico, basado un termal existente
    # El numero del caso termal (temcaso) debe ser un modelo termal ya generado

    direTerMec = 'Output/%s_Termal/%s_Mecanico' %(temcaso,meccaso)

    if not os.path.exists(direTerMec):
        os.makedirs(direTerMec)

    shutil.copy( 'VarMecanico.txt', direTerMec )

    # Retorna los directorios creados

    return direTer, direTerMec

def read_rheo(name):

    # Abrir archivo
    f = open(name)

    dic = {}

    # Loop para leer las variables
    line = f.readline()
    while line:

        # ignorar las lineas que empiezan con #
        if line[0] != '#':

            id_rh, name, h, n, a, ref = line.split()

            dic[id_rh] = DotDict({'name': name, 'H': float(h), 'n': float(n), 'A': float(a), 'ref': ref})

        # Leer siguiente linea
        line = f.readline()

    f.close()

    # Retorna la matriz con todas las variables

    return DotDict(dic)


