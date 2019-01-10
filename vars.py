import numpy as np
import multiprocessing as mp
import resource
import sys
import ast
from src.setup import input_setup, exec_setup
from termomecanico import termomecanico
from src.plot import rmse_plot
from src.utils import makedir
from src.rmse import rmse
from src.datos_q import shf_data
from pathlib import Path

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

def vars_rmse(t_input, m_input, var_names, var_ranges, var_type='thermal'):
    #t_input, m_input = input_setup()
    var_name_1 = get_var_name(var_names[0])
    var_axis_1 = get_var_axis(var_ranges[0])
    var_name_2 = None
    var_axis_2 = None
    if len(var_names) > 1:
        var_name_2 = get_var_name(var_names[1])
        var_axis_2 = get_var_axis(var_ranges[1])
    rmses = []
    i = 0
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
                i += 1
                print(i)
                rmses.append(get_model_rmse(t_input, m_input))
    if var_name_2 is not None:
        rmses = np.array(rmses).reshape(len(var_axis_1), len(var_axis_2))
    return rmses

def get_model_rmse(t_input, m_input):
    out_q = mp.Queue()
    def mp_termomecanico(t_input, m_input, queue):
        #should be model, model_rmse, ishf ??
        model = termomecanico(t_input, m_input)
        shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
        estimators = rmse(shf, shf_data)
        queue.put(estimators['rmse'])
        #queue.put('hola')
        return
    proc = mp.Process(target=mp_termomecanico, args=(t_input, m_input, out_q))
    proc.start()
    result = out_q.get()
    proc.join()
    return result

def get_results(t_input, m_input, vnames, vranges, vaxes, save_dir, filename=''):
    if len(vaxes) > 1:
        print(vnames[0], len(vaxes[0]))
        print(vnames[1], len(vaxes[1]))
        print('Numero de modelos:', len(vaxes[0])*len(vaxes[1]))
    else:
        print('Numero de modelos:', len(vaxes[0]))
    save_dir_files = save_dir + 'Archivos/'
    makedir(save_dir_files)
    rmses = vars_rmse(t_input, m_input, vnames, vranges)
    print('rmses:', rmses)
    # Save
    np.savetxt(save_dir_files + 'vars_rmses' + filename + '.txt', rmses)
    np.savetxt(save_dir_files + 'vars_names' + filename + '.txt', vnames, fmt='%s')
    np.savetxt(save_dir_files + 'vars_ranges' + filename + '.txt', vranges)
    #np.savetxt(save_dir + 'vars_axes.txt', vaxes)#, dtype='object')
    # Plot
    rmse_plot(vnames, vaxes, rmses, filename=save_dir + 'RMSE' + filename)

def load_and_plot_results(save_dir, filename):
    save_dir_files = save_dir + 'Archivos/'
    makedir(save_dir_files)
    rmses = np.loadtxt(save_dir_files + 'vars_rmses' + filename + '.txt')
    vranges = np.loadtxt(save_dir_files + 'vars_ranges' + filename + '.txt', ndmin=2)
    vnames = np.loadtxt(save_dir_files + 'vars_names' + filename + '.txt', dtype='str', ndmin=2)
    # Add dummy member to avoid avoid getting too many indices for array
    # when there is only one member present and we use vnames[0]
    #vaxes = np.loadtxt(save_dir + 'vars_axes.txt', ndmin=2)#, dtype='object')
    vaxes = [get_var_axis(vranges[i]) for i in range(len(vranges))]
    rmse_plot(vnames, vaxes, rmses, filename=save_dir + 'RMSE' + filename)
    #print("Use: python vars.py var_name '[start, end, step]'")

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
        vaxes = [get_var_axis(vranges[i]) for i in range(len(vranges))]
        t_input, m_input = input_setup()
        if len(vnames) > 2:
            control_vname = vnames[2]
            control_vax = vaxes[2]
            print(control_vname, len(control_vax))
            vnames = vnames[0:2]
            vranges = vranges[0:2]
            vaxes = vaxes[0:2]
            np.savetxt(save_dir + 'control_vname.txt', [control_vname], fmt='%s')
            np.savetxt(save_dir + 'control_vax.txt', control_vax)
            print('listo')
            for var in control_vax:
                t_input[control_vname] = var
                print(control_vname, '=', var)
                filename = '_' + control_vname + '_' + str(var)
                get_results(t_input, m_input, vnames, vranges, vaxes, save_dir, filename)
        else:
            get_results(t_input, m_input, vnames, vranges, vaxes, save_dir)

        ###
        #third_vname = ''
        #third_vax = [0]
        #if len(vnames) == 3:
        #    third_vname = vnames[2] 
        #    third_vax = vaxes[2]
        #    print(third_vname, len(third_vax))
        #    vnames = vnames[0:2]
        #    vranges = vranges[0:2]
        #    vaxes = vaxes[0:2]
        #    np.savetxt(save_dir_files + 'third_vname.txt', third_vname)
        #    np.savetxt(save_dir_files + 'third_vax.txt', third_vax)
        #t_input, m_input = input_setup()
        #for var in third_vax:
        #    if third_vname: 
        #        t_input[third_vname] = var
        #        print(third_vname + ' = ' + str(var))
        #        third_v_filename = '_' + third_vname + '_' + str(var)
        #    else:    
        #        third_v_filename = ''
        #    ###
        #    if len(vaxes) > 1:
        #        print(vnames[0], len(vaxes[0]))
        #        print(vnames[1], len(vaxes[1]))
        #        print('Numero de modelos:', len(vaxes[0])*len(vaxes[1]))
        #    else:
        #        print('Numero de modelos:', len(vaxes[0]))
        #    rmses = vars_rmse(t_input, m_input, vnames, vranges)
        #    print('rmses:', rmses)
        #    # Save
        #    np.savetxt(save_dir_files + 'vars_rmses' + third_v_filename + '.txt', rmses)
        #    np.savetxt(save_dir_files + 'vars_names' + third_v_filename + '.txt', vnames, fmt='%s')
        #    np.savetxt(save_dir_files + 'vars_ranges' + third_v_filename + '.txt', vranges)
        #    #np.savetxt(save_dir + 'vars_axes.txt', vaxes)#, dtype='object')
        #    # Plot
        #    rmse_plot(vnames, vaxes, rmses, filename=save_dir + 'RMSE' 
        #            + third_v_filename)

    else:
        # Load and Plot
        control = Path(save_dir + 'control_vname.txt')
        if control.is_file():
            control_vname = np.loadtxt(save_dir + 'control_vname.txt', dtype='str')
            control_vax = np.loadtxt(save_dir + 'control_vax.txt')
            for var in control_vax:
                filename = '_' + str(control_vname) + '_' + str(var)
                load_and_plot_results(save_dir, filename)
        else:
            load_and_plot_results(save_dir, '')
        #rmses = np.loadtxt(save_dir + 'vars_rmses.txt')
        ##vaxes = np.loadtxt(save_dir + 'vars_axes.txt', ndmin=2)#, dtype='object')
        #vranges = np.loadtxt(save_dir + 'vars_ranges.txt')
        #vaxes = [get_var_axis(vranges[i]) for i in range(len(vranges))]
        #vnames = np.loadtxt(save_dir + 'vars_names.txt', dtype='str', ndmin=2)
        ## Add dummy member to avoid avoid getting too many indices for array
        ## when there is only one member present and we use vnames[0]
        #rmse_plot(vnames, vaxes, rmses, filename=save_dir + 'RMSE')
        ##print("Use: python vars.py var_name '[start, end, step]'")
