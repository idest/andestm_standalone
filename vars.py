import numpy as np
import multiprocessing as mp
import resource
import sys
import ast
from setup import input_setup
from termomecanico2 import termomecanico
from plot2 import rmse_plot

# print memory usage
def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

t_input, m_input = input_setup()

def vars_rmse(
        var_tuple_1, var_tuple_2=None, var_type='thermal', return_ranges=False):
    t_input, m_input = input_setup()
    var_name_1 = var_tuple_1[0]
    var_range_1 = np.arange(*var_tuple_1[1:])
    var_range_1 = np.append(var_range_1, var_range_1[-1] + var_tuple_1[3])
    if var_tuple_2 is not None:
        var_name_2 = var_tuple_2[0]
        var_range_2 = np.arange(*var_tuple_2[1:])
        var_range_2 = np.append(var_range_2, var_range_2[-1] + var_tuple_2[3])
    rmses = []
    for var_1 in var_range_1:
        t_input[var_name_1] = var_1
        if var_tuple_2 is None:
            rmses.append(get_model_rmse(t_input, m_input))
        else:
            for var_2 in var_range_2:
                t_input[var_name_2] = var_2
                mem()
                rmses.append(get_model_rmse(t_input, m_input))
    if var_tuple_2 is not None:
        rmses = np.array(rmses).reshape(len(var_range_1), len(var_range_2))
    return_tuple = []
    return_tuple.append(rmses)
    if return_ranges is True:
        if var_tuple_2 is None:
            return_tuple.append([var_name1, var_range1])
        else:
            return_tuple.append([var_name1, var_range1, var_name2, var_range2])
    return return_tuple

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
    return rmse, range_1, range_2

if __name__ == '__main__':
    if len(sys.argv) > 2:
        args = sys.argv[1:]
        couples = len(args) // 2
        vtuples = []
        for couple in range(couples):
            vname = args[couple*2]
            vrange = ast.literal_eval(args[couple*2+1])
            print(vname, vrange)
            vtuple = [vname] + vrange
            vtuples.append(vtuple)
        rmses, var_ranges = vars_rmse(*vtuples)
        print(rmses)
        rmse_plot(var_ranges[0], var_ranges[2], var_ranges[1], var_ranges[3], rmses)

    else:
        print("Use: python vars.py var_name '[start, end, step]'")
