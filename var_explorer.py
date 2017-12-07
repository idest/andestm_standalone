import setup
import compute
import plot
import numpy as np
import os
import sys
import numpy as np
import gc
import resource
import multiprocessing as mp
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

print('Leyendo variables...')
#
t_input = setup.readVars('VarTermal.txt')
m_input = setup.readVars('VarMecanico.txt')
exec_input = setup.readVars('VarExec.txt')
direTer, direTerMec = setup.makeDirs(exec_input.temcaso, exec_input.meccaso)
#
gm_data = np.loadtxt('data/Modelo.dat')
areas = np.loadtxt('data/areas.dat')
trench_age = np.loadtxt('data/PuntosFosaEdad.dat')
rhe_data = setup.read_rheo('data/Rhe_Param.dat')

mem()

#D, CS, GM, TM, MM = compute.compute(gm_data, areas, trench_age, rhe_data, t_input, m_input)

#gt_ref = TM.get_geotherm()
#yse_ref = MM.get_yse()[0]

#mem()

input_type = t_input
var = 'k_cs'
var_range = np.arange(3, 3.6, 0.2) 

queue = mp.Queue()

def comp(out_q):
    D, CS, GM, TM, MM = compute.compute(gm_data, areas, trench_age, rhe_data, t_input, m_input)
    models = {}
    models['gt'] = TM.get_geotherm()
    models['yse_t'] = MM.get_yse()[0]
    queue.put(models) 

proc = mp.Process(target=comp, args=(queue,))
proc.start()
models_ref = queue.get()
proc.join()

i = 0
gt_square_errors = {}
yse_square_errors = {}

for value in var_range:
    input_type[var] = value
    print(t_input.k_cs)
    proc = mp.Process(target=comp, args=(queue,))
    proc.start()
    models = queue.get()
    gt_square_errors['gt_{}'.format(i)] = (models['gt'] - models_ref['gt'])**2
    yse_square_errors['yse_{}'.format(i)] = (models['yse_t'] - models_ref['yse_t'])**2
    proc.join()
    mem()
    i += 1

gt_se_values = list(gt_square_errors.values())
yse_se_values = list(yse_square_errors.values())

gt_rmse = np.sqrt(sum(gt_se_values)/len(gt_se_values))
yse_rmse = np.sqrt(sum(yse_se_values)/len(yse_se_values))
print(type(gt_se_values))
gt_rmse_map = np.nanmean(gt_rmse, axis=2)
yse_rmse_map = np.nanmean(gt_rmse, axis=2)

map = Basemap(llcrnrlon=-79.8, llcrnrlat=-44.8, urcrnrlon=-58.0, urcrnrlat=-10.0)
map.drawparallels(np.arange(-90,90,3), labels=[1,0,0,0])
map.drawmeridians(np.arange(-180,180,4), labels=[0,0,1,0])
map.drawcoastlines()
x = np.linspace(map.llcrnrx, map.urcrnrx, CS.get_x_axis().shape[0])
y = np.linspace(map.llcrnry, map.urcrnry, CS.get_y_axis().shape[0])
xx, yy = np.meshgrid(x,y)
xxx, yyy = CS.get_2D_grid()
gt_rmse_map = np.ma.masked_invalid(gt_rmse_map)
eet = np.ma.masked_invalid(MM.get_eet())
M = map.pcolormesh(xx, yy[::-1], gt_rmse_map.T, cmap='coolwarm', shading='gouraud')
#M = map.pcolormesh(xx, yy[::-1], eet.T, cmap='coolwarm', shading='gouraud')
cbar = plt.colorbar(M)
plt.savefig('mapa')

pass
