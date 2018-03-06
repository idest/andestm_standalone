import numpy as np
import multiprocessing as mp
import resource
import sys
import ast
from setup import input_setup, exec_setup
from termomecanico import termomecanico
from plot import rmse_plot
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

#def get_var_axis(range):
#    #if type(range) is not np.ndarray:
#    #    tuple = np.array(range)
#    #range.astype(np.float)
#    var_axis = np.arange(*range)
#    var_axis = np.append(var_axis, var_axis[-1] + range[2])
#    return var_axis

def get_var_axis(range):
    step = range[2] 
    first_v = range[0]
    last_v = range[1]
    var_axis = np.linspace(first_v, last_v, num=(abs(last_v-first_v))/step+1,
                           endpoint=True)
    return var_axis

def vars_rmse(var_names, var_ranges, var_type='thermal'):
    t_input, m_input = input_setup()
    var_name_1 = get_var_name(var_names[0])
    var_axis_1 = get_var_axis(var_ranges[0])
    var_name_2 = None
    var_axis_2 = None
    if len(var_names) > 1:
        var_name_2 = get_var_name(var_names[1])
        var_axis_2 = get_var_axis(var_ranges[1])
    rmses = []
    for var_1 in var_axis_1:
        for vn in var_name_1:
            t_input[vn] = var_1
        if var_name_2 is None:
            rmses.append(get_model_rmse(t_input, m_input))
        else:
            for var_2 in var_axis_2:
                for vn in var_name_2:
                    t_input[vn] = var_2
                #mem()
                rmses.append(get_model_rmse(t_input, m_input))
    if var_name_2 is not None:
        rmses = np.array(rmses).reshape(len(var_axis_1), len(var_axis_2))
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
        rmses = vars_rmse(vnames, vranges)
        print('rmses:', rmses)
        # Save
        vaxes = [get_var_axis(vranges[i]) for i in range(len(vranges))]
        np.savetxt(save_dir + 'vars_rmses.txt', rmses)
        np.savetxt(save_dir + 'vars_names.txt', vnames, fmt='%s')
        np.savetxt(save_dir + 'vars_ranges.txt', vranges)
        #np.savetxt(save_dir + 'vars_axes.txt', vaxes)#, dtype='object')
        # Plot
        rmse_plot(vnames, vaxes, rmses, save_dir)

    else:
        # Load and Plot
        rmses = np.loadtxt(save_dir + 'vars_rmses.txt')
        #vaxes = np.loadtxt(save_dir + 'vars_axes.txt', ndmin=2)#, dtype='object')
        vranges = np.loadtxt(save_dir + 'vars_ranges.txt')
        vaxes = [get_var_axis(vranges[i]) for i in range(len(vranges))]
        vnames = np.loadtxt(save_dir + 'vars_names.txt', dtype='str', ndmin=2)
        # Add dummy member to avoid avoid getting too many indices for array
        # when there is only one member present and we use vnames[0]
        rmse_plot(vnames, vaxes, rmses, save_dir)
        #print("Use: python vars.py var_name '[start, end, step]'")
