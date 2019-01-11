import numpy as np
import multiprocessing as mp
import resource
import sys
import ast
from src.setup import input_setup, exec_setup
from termomecanico import termomecanico
from src.plot import estimator_plot
from src.utils import makedir
from src.stats import evaluate_model
from src.datos_q import shf_data
from pathlib import Path
import pickle

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
    if len(range) == 3:
        step = range[2] 
        first_v = range[0]
        last_v = range[1]
        var_axis = np.linspace(first_v, last_v, num=(abs(last_v-first_v))/step+1,
                               endpoint=True)
    else:
        var_axis = range
    return var_axis

def vars_estimators(t_input, m_input, var_names, var_axes, var_type='thermal'):
    #t_input, m_input = input_setup()
    var_name_1 = get_var_name(var_names[0])
    var_axis_1 = var_axes[0]
    var_name_2 = None
    var_axis_2 = None
    if len(var_names) > 1:
        var_name_2 = get_var_name(var_names[1])
        var_axis_2 = var_axes[1]
    rmses = []
    mses = []
    i = 0
    for var_1 in var_axis_1:
        for vn in var_name_1:
            t_input[vn] = var_1
        if var_name_2 is None:
            rmses.extend(get_model_estimators(t_input, m_input))
            mses.append(rmses.pop())
        else:
            for var_2 in var_axis_2:
                for vn in var_name_2:
                    t_input[vn] = var_2
                #mem()
                i += 1
                print(i)
                rmses.extend(get_model_estimators(t_input, m_input))
                mses.append(rmses.pop())
    if var_name_2 is not None:
        rmses = np.array(rmses).reshape(len(var_axis_1), len(var_axis_2))
        mses = np.array(mses).reshape(len(var_axis_1), len(var_axis_2))
    return {'rmses': rmses, 'mses': mses}

def get_model_estimators(t_input, m_input):
    out_q = mp.Queue()
    def mp_termomecanico(t_input, m_input, queue):
        model = termomecanico(t_input, m_input)
        shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
        estimators = evaluate_model(shf, shf_data)
        queue.put([estimators['rmse'], estimators['mse']])
        return
    proc = mp.Process(target=mp_termomecanico, args=(t_input, m_input, out_q))
    proc.start()
    result = out_q.get()
    proc.join()
    return result

def get_results(t_input, m_input, dic, save_dir, filename=''):
    vnames = list(dic.keys())
    vaxes = list(dic.values())
    if len(vaxes) > 1:
        print(vnames[0], len(vaxes[0]))
        print(vnames[1], len(vaxes[1]))
        print('Numero de modelos:', len(vaxes[0])*len(vaxes[1]))
    else:
        print('Numero de modelos:', len(vaxes[0]))
    save_dir_files = save_dir + 'Archivos/'
    makedir(save_dir_files)
    estimators = vars_estimators(t_input, m_input, vnames, vaxes)
    print('rmses:', estimators['rmses'])
    print('mses:', estimators['mses'])
    # Save
    np.savetxt(
        save_dir_files + 'vars_rmses' + filename + '.txt', estimators['rmses'])
    np.savetxt(
        save_dir_files + 'vars_mses' + filename + '.txt', estimators['mses'])
    with open(save_dir_files + 'variables' + filename + '.txt', 'wb') as f:
        pickle.dump(dic, f)
    # Plot
    estimator_plot(
        vnames, vaxes, estimators['rmses'], label = 'RMSE',
        filename=save_dir+'RMSE'+filename)
    estimator_plot(
        vnames, vaxes, estimators['mses'], signed=True, label='MSE',
        filename=save_dir+'MSE'+filename)

def load_and_plot_results(save_dir, filename):
    save_dir_files = save_dir + 'Archivos/'
    makedir(save_dir_files)
    with open(save_dir_files + 'variables' + filename + '.txt', 'rb') as f:
        dic = pickle.load(f)
    rmses = np.loadtxt(save_dir_files + 'vars_rmses' + filename + '.txt')
    mses = np.loadtxt(save_dir_files + 'vars_mses' + filename + '.txt')
    vnames = list(dic.keys())
    vaxes = list(dic.values())
    estimator_plot(
        vnames, vaxes, estimators['rmses'], filename=save_dir+'RMSE'+filename)
    estimator_plot(
        vnames, vaxes, estimators['mses'], signed=True,
        filename=save_dir+'MSE'+filename)

if __name__ == '__main__':
    exec_input, direTer, direMec = exec_setup()
    save_dir = direTer + 'Graficos/'
    makedir(save_dir)
    if len(sys.argv) > 2:
        args = sys.argv[1:]
        couples = len(args) // 2
        vnames = []; vranges = []; vtuples = []
        dic = {}
        for couple in range(couples):
            vname = args[couple*2]
            vrange = ast.literal_eval(args[couple*2+1])
            vax = get_var_axis(vrange)
            print(vname, vrange)
            dic[vname] = vax
        t_input, m_input = input_setup()
        if len(dic) > 2:
            control_vname = vname
            control_vax = dic.pop(control_vname)
            print(control_vname, len(control_vax))
            np.savetxt(save_dir + 'control_vname.txt', [control_vname], fmt='%s')
            np.savetxt(save_dir + 'control_vax.txt', control_vax)
            print('listo')
            for var in control_vax:
                t_input[control_vname] = var
                print(control_vname, '=', var)
                filename = '_' + control_vname + '_' + str(var)
                get_results(t_input, m_input, dic, save_dir, filename)
        else:
            get_results(t_input, m_input, dic, save_dir)

    else:
        # Load and Plot
        control = Path(save_dir + 'control_vname.txt')
        if control.is_file():
            control_vname = np.loadtxt(save_dir+'control_vname.txt', dtype='str')
            control_vax = np.loadtxt(save_dir + 'control_vax.txt')
            for var in control_vax:
                filename = '_' + str(control_vname) + '_' + str(var)
                load_and_plot_results(save_dir, filename)
        else:
            load_and_plot_results(save_dir, '')
        ##print("Use: python vars.py var_name '[start, end, step]'")
