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
mode = '2d'
var = 'G'
var_range = np.arange(1.e-4,1.e-3,0.5e-4)
var2 = 'H'
var2_range = np.arange(1.e-6,6.e-6,5e-7)
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

def comp(model,input_type,n,value,mode='1d',value2=None):
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
    if mode=='1d':
        np.savetxt('%s_%s_%s_%s.txt' %(model_str, var, value, n), model_ref,
                   fmt="%11.4f",delimiter="  ")
    if mode=='2d':
        np.savetxt('%s_%s_%s_%s_%s_%s.txt' %(model_str,var,value,var2,value2,n),
                   model_ref, fmt="%11.4f",delimiter="  ")
    os.chdir('../../../')
    return

def dif_models(var_range, var, var2_range, var2, model):
    n = 0
    for value in var_range:
        n += 1
        input_type[var] = value
        print(input_type[var])
        if mode=='2d':
            for value2 in var2_range:
                n += 1
                change_var(input_type, var2, value2)
                keywords = {'mode': '2d', 'value2': value2}
                proc = mp.Process(target=comp, args=(model,input_type,n,value), kwargs=keywords)
                proc.start()
                mem()
                proc.join()
        else:
            proc = mp.Process(target=comp, args=(model,input_type,n,value))
            proc.start()
            mem()
            proc.join()
    if not os.path.exists('Output/var_models/%s_%s' %(var,tmc)):
        os.makedirs('Output/var_models/%s_%s' %(var,tmc), exist_ok=True)
    os.chdir('Output/var_models/%s_%s' %(var,tmc))
    print(len(var_range))
    print(len(var2_range))
    shape_2d = np.asarray([len(var_range), len(var2_range)],dtype=int)
    np.savetxt('shape_2d.txt', shape_2d)
    os.chdir('../../../')
    return

def change_var(input_type, var, value):
    if var == 'H':
        input_type['H_cs'] = value
        input_type['H_ci'] = value
        input_type['H_ml'] = value
    else:
        input_type[var] = value
    return

proc = mp.Process(target=dif_models, args=(var_range,var,var2_range,var2,model))
proc.start()
print('procesando')
proc.join()