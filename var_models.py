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

#print memory usage
def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

#leyendo variables del modelo
print('Leyendo variables...')
t_input = setup.readVars('VarTermal.txt')
m_input = setup.readVars('VarMecanico.txt')
exec_input = setup.readVars('VarExec.txt')
tmc = exec_input.temcaso
mmc = exec_input.meccaso
direTer, direTerMec = setup.makeDirs(exec_input.temcaso, exec_input.meccaso)
#
gm_data = np.loadtxt('data/Modelo.dat')
areas = np.loadtxt('data/areas.dat')
trench_age = np.loadtxt('data/PuntosFosaEdad.dat')
rhe_data = setup.read_rheo('data/Rhe_Param.dat')
#
input_type = t_input
var = 'G'
var_range = np.arange(1.e-3,5.e-2,5.e-4)
model = exec_input.model
#k_cs = t_input.k_cs
#k_ci = t_input.k_ci
#k_ml = t_input.k_ml
#H_cs = t_input.H_cs
#H_ci = t_input.H_ci
#H_ml = t_input.H_ml
var_mean = 'prom_k'

mem()
queue = mp.Queue()

def comp(model,value,input_type):
    D, CS, GM, TM, MM = compute.compute(gm_data, areas, trench_age,
                                        rhe_data, input_type, m_input)
    array_model = []
    if model==1:
        array_model = TM.get_surface_heat_flow()
        model_str = 'shf'
    if model==2:
        array_model = MM.get_yse()[0]
        model_str = 'yse_t'
    if model==3:
        array_model = MM.get_eet()
        model_str = 'eet'
    queue.put(array_model)
    model_ref = queue.get()
    if not os.path.exists('Output/var_models/%s_%s' %(var,tmc)):
        os.makedirs('Output/var_models/%s_%s' %(var,tmc), exist_ok=True)
    os.chdir('Output/var_models/%s_%s' %(var,tmc))
    np.savetxt('%s_%s_%s.txt' %(model_str, var, value), model_ref,
               fmt="%11.4f",delimiter="  ")
    os.chdir('../../../')
    return

def dif_models(var_range, var, model):
	for value in var_range:
	    input_type[var] = value
	    print(input_type[var])
	    proc = mp.Process(target=comp, args=(model,value,input_type))
	    proc.start()
	    mem()
	    proc.join()
	return

proc = mp.Process(target=dif_models, args=(var_range,var,model))
proc.start()
print('procesando')
proc.join()


