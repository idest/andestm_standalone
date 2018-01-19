import setup
import compute
import plot
import numpy as np
np.set_printoptions(threshold=np.nan)
from utils import DotDict
import os
import sys
import resource


def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

print('Leyendo variables...')
t_input = setup.readVars('VarTermal.txt')
m_input = setup.readVars('VarMecanico.txt')
exec_input = setup.readVars('VarExec.txt')
tmc = exec_input.temcaso
mmc = exec_input.meccaso
direTer, direTerMec = setup.makeDirs(exec_input.temcaso, exec_input.meccaso)
gm_data = np.loadtxt('data/Modelo.dat')
areas = np.loadtxt('data/areas.dat')
trench_age = np.loadtxt('data/PuntosFosaEdad.dat')
rhe_data = setup.read_rheo('data/Rhe_Param.dat')
datos_q = np.loadtxt('datos_Q/QsObs.txt', comments='#')
datos_q = datos_q[datos_q[:,0] > -80.]
datos_q = datos_q[datos_q[:,0] < -60.]
datos_q = datos_q[datos_q[:,1] > -45.]
datos_q = datos_q[datos_q[:,1] < -10.]
datos_q = datos_q[datos_q[:,2] <= 120]

D, CS, GM, TM, MM = compute.compute(gm_data, areas, trench_age, rhe_data, t_input, m_input)
surface_heat_flow = TM.get_surface_heat_flow()

x_axis = CS.get_x_axis()
y_axis = CS.get_y_axis()
"""
print("After termomecanico M.S:")
mem()
"""
#plotear perfiles termales
"""
os.chdir(direTer)
fig = plot.plot_thermal(CS.get_axes()[0], CS.get_axes()[2], D, CS, GM, TM)
os.chdir('../../')


#plotear perfiles termomecanicos
os.chdir(direTerMec)
fig = plot.plot_mec(CS.get_axes()[0], CS.get_axes()[2], D, CS, GM, MM)
os.chdir('../../../')
"""

#plotear mapa q_surface
#os.chdir(direTer)
#fig = plot.map_q_surface(CS, TM, tmc, data_q)
#os.chdir('../../')

#plotear datos q
fig = plot.map_q_surface_2(x_axis, y_axis, tmc, direTer, surface_heat_flow=surface_heat_flow, 
                           data_q=datos_q,data_cmap='heat_flow')

#detachment = plot.get_detachment(CS,GM,MM)

