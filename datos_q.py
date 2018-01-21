#LOAD DATA
import numpy as np
#from termomecanico import CS, TM, direTer
import setup
from calc_compute import CS, TM, MM, GM
 
var_therm = 'k_cs'
exec_input = setup.readVars('VarExec.txt')
mc = exec_input.temcaso
mmc = exec_input.meccaso
x_axis = CS.get_x_axis()
y_axis = CS.get_y_axis()
surface_heat_flow = TM.get_surface_heat_flow()
exec_input = setup.readVars('VarExec.txt')
tmc = exec_input.temcaso
mmc = exec_input.meccaso
direTer, direTerMec = setup.makeDirs(exec_input.temcaso, exec_input.meccaso)
datos_q = np.loadtxt('datos_Q/QsObs.txt', comments='#')

datos_q = datos_q[datos_q[:,0] > -80.]
datos_q = datos_q[datos_q[:,0] < -60.]
datos_q = datos_q[datos_q[:,1] > -45.]
datos_q = datos_q[datos_q[:,1] < -10.]
datos_q = datos_q[datos_q[:,2] <= 120]
error = datos_q[:,-1]

datos_q_max = datos_q.copy()
datos_q_max[:,2] = datos_q_max[:,2]+datos_q_max[:,-1]
datos_q_max_shf = -datos_q_max[:,2]*1e-3

datos_q_min = datos_q.copy()
datos_q_min[:,2] = datos_q_min[:,2]-datos_q_min[:,-1]
datos_q_min_shf = -datos_q_min[:,2]*1e-3

# All data
datos_q_x = datos_q[:,0]
datos_q_y = datos_q[:,1]
datos_q_shf = -datos_q[:,2]*1e-3
# Marine Geophysics
datos_q_x_1 = datos_q[:,0][np.where(datos_q[:,-2]==1)]
datos_q_y_1 = datos_q[:,1][np.where(datos_q[:,-2]==1)]
datos_q_shf_1 = -datos_q[:,2][np.where(datos_q[:,-2]==1)]*1.e-3
# Geochemical
datos_q_x_2 = datos_q[:,0][np.where(datos_q[:,-2]==2)]
datos_q_y_2 = datos_q[:,1][np.where(datos_q[:,-2]==2)]
datos_q_shf_2 = -datos_q[:,2][np.where(datos_q[:,-2]==2)]*1.e-3
# Land Borehole
datos_q_x_3 = datos_q[:,0][np.where(datos_q[:,-2]==3)]
datos_q_y_3 = datos_q[:,1][np.where(datos_q[:,-2]==3)]
datos_q_shf_3 = -datos_q[:,2][np.where(datos_q[:,-2]==3)]*1.e-3
# ODP Borehole
datos_q_x_4 = datos_q[:,0][np.where(datos_q[:,-2]==4)]
datos_q_y_4 = datos_q[:,1][np.where(datos_q[:,-2]==4)]
datos_q_shf_4 = -datos_q[:,2][np.where(datos_q[:,-2]==4)]*1.e-3