import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from src.setup import input_setup, exec_setup, read_rheo
from termomecanico import termomecanico
from src.utils import makedir

if __name__ == '__main__':
    exec_input, direTer, direMec = exec_setup()
    t_input, m_input = input_setup()
    rhe_data = read_rheo('data/Rhe_Param.dat')
    model, _, _ = termomecanico(t_input, m_input)
    x_axis = model.cs.get_x_axis()
    z_axis = model.cs.get_z_axis()
    z_axis_2 = [*z_axis, *z_axis]
    geotherm = model.tm.get_geotherm()
    bys_t = model.mm.bys_t
    bys_c = model.mm.bys_c
    geotherm_1D = geotherm[90,0,:]
    bys_t_1D = bys_t[90,0,:]
    bys_c_1D = bys_c[90,0,:]
    bys_1D = np.concatenate((bys_c_1D, bys_t_1D))
    dys_list = []
    for key, value in rhe_data.items():
        print('n:', value.n, 'a:', value.A, 'h:', value.H)
        dys = model.mm.calc_ductile_yield_strength(m_input.e, value.n, value.A,
            value.H, m_input.R, geotherm_1D)
        dys_list.append(dys)
    fig = plt.figure()   
    plt.xlim(-2000,2000)
    plt.ylim(-100,20)
    plt.plot(bys_1D, z_axis_2, 'k')
    #plt.plot(bys_c_1D, z_axis, 'k')
    #plt.plot(bys_t_1D, z_axis, 'k')
    colors = iter(cm.rainbow(np.linspace(0, 1, len(dys_list))))
    print(len(bys_1D))
    print(len(z_axis_2))
    for i, dys in enumerate(dys_list):
        color = next(colors)
        print(i)
        print(color)
        dys_1D = [*dys, *-dys]
        print(len(dys_1D))
        plt.plot(dys_1D, z_axis_2, color=color)
    plt.plot(np.repeat(200,len(z_axis)),z_axis,'r')
    plt.plot(np.repeat(-200,len(z_axis)),z_axis,'r')
    #legend = plt.legend(['c{}'.format(i) for i in range(len(dys_list))], loc=2)
    #plt.tight_layout(pad=7)
    plt.show()
