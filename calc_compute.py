import numpy as np
import compute
import setup 

t_input = setup.readVars('VarTermal.txt')
m_input = setup.readVars('VarMecanico.txt')
gm_data = np.loadtxt('data/Modelo.dat')
areas = np.loadtxt('data/areas.dat')
trench_age = np.loadtxt('data/PuntosFosaEdad.dat')
rhe_data = setup.read_rheo('data/Rhe_Param.dat') 
D, CS, GM, TM, MM = compute.compute(gm_data, areas, trench_age, rhe_data, t_input, m_input)