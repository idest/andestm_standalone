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
direTer, direTerMec = setup.makeDirs(exec_input.temcaso, exec_input.meccaso)
#
gm_data = np.loadtxt('data/Modelo.dat')
areas = np.loadtxt('data/areas.dat')
trench_age = np.loadtxt('data/PuntosFosaEdad.dat')
rhe_data = setup.read_rheo('data/Rhe_Param.dat')
#
input_type = t_input
var = 'k_cs'
var_range = np.arange(1.5, 3.6, 0.1)
model = exec_input.model
k_cs = t_input.k_cs
k_ci = t_input.k_ci
k_ml = t_input.k_ml
var_mean = (k_cs+k_ci+k_ml)/3

mem()
queue = mp.Queue()

def comp(model,value):
    print('comp')
    D, CS, GM, TM, MM = compute.compute(gm_data, areas, trench_age,
                                        rhe_data, t_input, m_input)
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
    if not os.path.exists('Output/var_models/%s' %(var)):
        os.makedirs('Output/var_models/%s' %(var), exist_ok=True)
    os.chdir('Output/var_models/%s' %(var))
    np.savetxt('%s_%s_%s.txt' %(model_str, var, value), model_ref,
               fmt="%11.4f",delimiter="  ")
    os.chdir('../../../')
    return

def dif_models(var_range, var, model):
	for value in var_range:
	    input_type[var] = value
	    print(input_type[var])
	    proc = mp.Process(target=comp, args=(model,value))
	    proc.start()
	    mem()
	    proc.join()
	return

proc = mp.Process(target=dif_models, args=(var_range,var,model))
proc.start()
print('procesando')
proc.join()


