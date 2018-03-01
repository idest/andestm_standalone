import numpy as np
import multiprocessing as mp
import resource
import sys
import ast
from setup import input_setup, exec_setup
from termomecanico2 import termomecanico
from plot2 import rmse_plot
from utils import makedir

# print memory usage
def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

def get_var_name(name):
    if name == 'k':
        var_name = ['k_cs', 'k_ci', 'k_ml']
    elif name == 'h':
        var_name = ['H_cs', 'H_ci', 'H_ml']
    else:
        var_name = [name]
    return var_name

def get_var_range(tuple):
    if type(tuple) is not np.ndarray:
        tuple = np.array(tuple)
    tuple.astype(np.float)
    var_range = np.arange(*tuple)
    var_range = np.append(var_range, var_range[-1] + tuple[2])
    return var_range

def vars_rmse(
        var_tuple_1, var_tuple_2=None, var_type='thermal'):
    t_input, m_input = input_setup()
    var_name_1 = get_var_name(var_tuple_1[0])
    var_range_1 = get_var_range(var_tuple_1[1:])
    if var_tuple_2 is not None:
        var_name_2 = get_var_name(var_tuple_2[0])
        var_range_2 = get_var_range(var_tuple_2[1:])
    rmses = []
    for var_1 in var_range_1:
        for vn in var_name_1:
            t_input[vn] = var_1
        if var_tuple_2 is None:
            rmses.append(get_model_rmse(t_input, m_input))
        else:
            for var_2 in var_range_2:
                for vn in var_name_2:
                    t_input[vn] = var_2
                #mem()
                rmses.append(get_model_rmse(t_input, m_input))
    if var_tuple_2 is not None:
        rmses = np.array(rmses).reshape(len(var_range_1), len(var_range_2))
    return rmses

def get_model_rmse(t_input, m_input):
    out_q = mp.Queue()
    def mp_termomecanico(t_input, m_input, queue):
        rmse, model_rmse, ishf = termomecanico(t_input, m_input)
        queue.put(model_rmse.rmse)
        return
    proc = mp.Process(target=mp_termomecanico, args=(t_input, m_input, out_q))
    proc.start()
    rmse = out_q.get()
    proc.join()
    return rmse

if __name__ == '__main__':
    exec_input, direTer, direMec = exec_setup()
    save_dir = direTer + 'Graficos/'
    makedir(save_dir)
    if len(sys.argv) > 2:
        args = sys.argv[1:]
        couples = len(args) // 2
        vnames = []; vranges = []; vtuples = []
        for couple in range(couples):
            vname = args[couple*2]
            vrange = ast.literal_eval(args[couple*2+1])
            print(vname, vrange)
            vtuple = [vname] + vrange
            vnames.append(vname)
            vranges.append(vrange)
            vtuples.append(vtuple)
        rmses = vars_rmse(*vtuples)
        # Save
        vranges = [get_var_range(vranges[0]), get_var_range(vranges[1])]
        tuple_1 = np.array([vnames[0], *vranges[0]])
        tuple_2 = np.array([vnames[1], *vranges[1]])
        np.savetxt(savedir + 'vars_rmses.txt', rmses)
        np.savetxt(savedir + 'vars_tuple_1.txt', tuple_1, fmt='%s')
        np.savetxt(savedir + 'vars_tuple_2.txt', tuple_2, fmt='%s')
        print('rmses:', rmses)
        # Plot
        rmse_plot(*vnames, *vranges, rmses, save_dir)

    else:
        # Load and Plot
        rmses = np.loadtxt(save_dir + 'vars_rmses.txt')
        tuple_1 = np.genfromtxt(save_dir + 'vars_tuple_1.txt', dtype='str')
        tuple_2 = np.genfromtxt(save_dir + 'vars_tuple_2.txt', dtype='str')
        vnames = [tuple_1[0], tuple_2[0]]
        vrange_1 = np.array(tuple_1[1:]).astype(float)
        vrange_2 = np.array(tuple_2[1:]).astype(float)
        vranges = [vrange_1, vrange_2]
        print(vnames)
        print(vranges)
        rmse_plot(vnames[0], vnames[1], vranges[0], vranges[1], rmses, save_dir)
        #print("Use: python vars.py var_name '[start, end, step]'")
