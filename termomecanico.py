from setup import data_setup, input_setup, exec_setup
from compute import compute
import plot
import os

#función que corre el modelo tm

def termomecanico(t_input, m_input):
    gm_data, areas, trench_age, rhe_data = data_setup()
    model = compute(gm_data, areas, trench_age, rhe_data, t_input, m_input)
    print('.')
    return model

#leyendo variables

print('reading vars ...')
t_input, m_input = input_setup()
exec_input, direTer, direMec = exec_setup()

#corriendo modelo

print('running model ...')
model = termomecanico(t_input, m_input)

#plotear perfiles termales

if exec_input.xt2 != 0:
	print('plotting thermal profiles ...')
	os.chdir(direTer)
	fig = plot.plot_thermal(model.cs.get_axes()[0], model.cs.get_axes()[2], model.d, model.cs, model.gm, model.tm)
	os.chdir('../../')

#plotear perfiles termomecanicos
if exec_input.xm2 != 0:
	print('plotting mechanic profiles ...')
	os.chdir(direMec)
	fig = plot.plot_mec(model.cs.get_axes()[0], model.cs.get_axes()[2], model.d, model.cs, model.gm, model.mm)
	os.chdir('../../../')

print('done')