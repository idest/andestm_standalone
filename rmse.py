import numpy as np
from scipy import stats
import numpy.ma as ma
import setup
np.set_printoptions(threshold=np.nan)
import scipy as sp
from scipy.interpolate import RectBivariateSpline, interp2d, RegularGridInterpolator
#from termomecanico import CS, TM, direTer, x_axis, y_axis
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import numpy.ma as ma
from calc_rmse import calc_rmse
from calc_rmse import calc_rmse_error
from matplotlib.ticker import FormatStrFormatter
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os
from plot import map_q_surface_2
from datos_q import *

surface_heat_flow = TM.get_surface_heat_flow()
shf_min = np.nanmin(surface_heat_flow)
shf_max = np.nanmax(surface_heat_flow)
surface_heat_flow_masked = surface_heat_flow.copy()
surface_heat_flow_masked[np.isnan(surface_heat_flow)] = -1.e-1000000000
shf_interpolator_bsi = RectBivariateSpline(x_axis, y_axis[::-1],
                                           surface_heat_flow_masked[:,::-1])
interpolated_shf_bsi = shf_interpolator_bsi.ev(datos_q_x, datos_q_y)


variable_models_directory = 'Output/var_models/'
mode = '2d'

for cdir in next(os.walk(variable_models_directory))[1]:
    var = cdir
    full_cdir = variable_models_directory + cdir
    full_cdir_encoded = os.fsencode(full_cdir)
    var_values = []
    var2_values = []
    var_positions = []
    rmses = []
    for file in os.listdir(full_cdir_encoded):
        filename = os.fsdecode(file)
        if not filename.endswith('.png') and filename!='shape_2d.txt' and filename!='.DS_Store':
            #print(filename)
            ext = ".txt"
            if mode == '1d':
                var_value = filename[:filename.find(ext)].split('_')[-1]
            elif mode == '2d':
                print(filename)
                var_value = filename[:filename.find(ext)].split('_')[-4]
                var2_value = filename[:filename.find(ext)].split('_')[-2]
                var2 = filename[:filename.find(ext)].split('_')[-3]
            var_position = filename[:filename.find(ext)].split('_')[-1]
            var2_values.append(var2_value)
            var_positions.append(var_position)
            var_values.append(var_value)
            # Cargar datos modelo
            surface_heat_flow = np.loadtxt(full_cdir+'/'+filename)
            shf_min = np.nanmin(surface_heat_flow)
            shf_max = np.nanmax(surface_heat_flow)
            # Interpolar modelo en puntos de datos Q con Bivariate Spline
            # Interpolation
            #surface_heat_flow = (surface_heat_flow[np.isfinite(
            #                                       surface_heat_flow)]
                                                    #Remove Nan values
            surface_heat_flow[np.isnan(surface_heat_flow)] = -1.e-1000000000
                                                             # nans = -9999

            shf_interpolator_bsi = RectBivariateSpline(x_axis, y_axis[::-1],
                                                    surface_heat_flow[:,::-1])
            interpolated_shf_bsi = shf_interpolator_bsi.ev(datos_q_x,
                                                           datos_q_y)

            interpolated_shf_bsi_1 = shf_interpolator_bsi.ev(datos_q_x_1,
                                                             datos_q_y_1)
            interpolated_shf_bsi_2 = shf_interpolator_bsi.ev(datos_q_x_2,
                                                             datos_q_y_2)
            interpolated_shf_bsi_3 = shf_interpolator_bsi.ev(datos_q_x_3,
                                                             datos_q_y_3)
            interpolated_shf_bsi_4 = shf_interpolator_bsi.ev(datos_q_x_4,
                                                             datos_q_y_4)

#             slope1,intercept1,r1,p1,std_err1 = stats.linregress(
#                                           interpolated_shf_bsi_1,datos_q_shf_1)
#             slope2,intercept2,r2,p2,std_err2 = stats.linregress(
#                                           interpolated_shf_bsi_2,datos_q_shf_2)
#             slope3,intercept3,r3,p3,std_err3 = stats.linregress(
#                                           interpolated_shf_bsi_3,datos_q_shf_3)
#             slope4,intercept4,r4,p4,std_err4 = stats.linregress(
#                                           interpolated_shf_bsi_4,datos_q_shf_4)


#             #Plotting scatter diagram
#             # fig = plt.figure()
#             # plot_q1=plt.scatter(interpolated_shf_bsi_1, datos_q_shf_1, c='b',
#             #                     marker='o', alpha=.4)
#             # plt.text(-0.099,-0.30,'R$^{2}$=%s'%(r1),{'color':'b',
#             #                                          'fontsize':7.})
#             # fitt1=plt.plot(interpolated_shf_bsi_1,
#             #              intercept1+slope1*interpolated_shf_bsi_1,'b',alpha=.4)

#             # plot_q2=plt.scatter(interpolated_shf_bsi_2, datos_q_shf_2, c='g',
#             #                     marker='^', alpha=.4)
#             # plt.text(-0.099,-0.29,'R$^{2}$=%s'%(r2),{'color':'g',
#             #                                          'fontsize':7.})
#             # fitt2=plt.plot(interpolated_shf_bsi_2,
#             #              intercept2+slope2*interpolated_shf_bsi_2,'g',alpha=.4)
#             # plot_q3=plt.scatter(interpolated_shf_bsi_3, datos_q_shf_3, c='y',
#             #                     marker='p', alpha=.4)
#             # plt.text(-0.099,-0.28,'R$^{2}$=%s'%(r3),{'color':'y',
#             #                                          'fontsize':7.})
#             # fitt3=plt.plot(interpolated_shf_bsi_3,
#             #              intercept3+slope3*interpolated_shf_bsi_3,'y',alpha=.4)
#             # plot_q4=plt.scatter(interpolated_shf_bsi_4, datos_q_shf_4, c='k',
#             #                     marker='s', alpha=.4)
#             # plt.text(-0.099,-0.27,'R$^{2}$=%s'%(r4),{'color':'k',
#             #                                          'fontsize':7.})
#             # fitt4=plt.plot(interpolated_shf_bsi_4,
#             #              intercept4+slope4*interpolated_shf_bsi_4,'k',alpha=.4)
#             # plt.xlim(-.1,0)
#             # plt.ylim(np.nanmin(datos_q_shf),0)
#             # plt.ylabel('heat flow data (W/m2)')
#             # plt.xlabel('heat flow modeled (W/m2)')
#             # plt.grid(True)
#             # plt.title('Scatter Diagram model/data for {}' .format(var))
#             # plt.legend([plot_q1, plot_q2, plot_q3, plot_q4],
#             #            ['ODP Borehole', 'Land Borehole',
#             #             'Geochemical', 'Marine Geophysics',], loc=4)
#             # os.chdir(full_cdir_encoded)
#             # fig.savefig('model_data_{}.png' .format(file))
#             # os.chdir('../../../')
#             # plt.close()

#             #print(surface_heat_flow)
#             """
#             # Descartar ciertos valores de interpolated_shf_bsi y de datos_q
#             interpolated_shf_bsi[abs(interpolated_shf_bsi) > 1] = np.nan
#             valid_interp_shf_idxs = np.where(np.isfinite(interpolated_shf_bsi))
#             interpolated_shf_bsi = interpolated_shf_bsi[valid_interp_shf_idxs]
#             datos_q_shf = datos_q_shf[valid_interp_shf_idxs]
#             datos_q_x = datos_q_x[valid_interp_shf_idxs]
#             datos_q_y = datos_q_y[valid_interp_shf_idxs]
#             weight = weight[valid_interp_shf_idxs]
#             """
#             # Calcular rmse
#             #V1 = sum(weight)
#             #V2 = sum(weight**2)
#             #diff_matrix = interpolated_shf_bsi - datos_q_shf
#             #print(V1-(V2/V1))
            rmse,_,_ = calc_rmse_error(datos_q_shf,interpolated_shf_bsi,datos_q_min_shf,datos_q_max_shf,error)
#             #rmse1 = mean_squared_error(datos_q_shf_1,interpolated_shf_bsi1)
#             #rmse2 = mean_squared_error(datos_q_shf_2,interpolated_shf_bsi2)
#             #rmse3 = mean_squared_error(datos_q_shf_3,interpolated_shf_bsi3)
#             #rmse4 = mean_squared_error(datos_q_shf_4,interpolated_shf_bsi4)
#             #estimator = sum(weight*diff_matrix.T*diff_matrix)/(V1-V2/V1)
#             #print("RMSE (Bivariate Spline Interpolation):", rmse_bsi)
            rmses.append(rmse)
#             #rmses1.append(rmse1)
#             #rmses2.append(rmse2)
#             #rmses3.append(rmse3)
#             #rmses4.append(rmse4)
#             #print(rmse)
    # Ordenar var_values de forma ascendente
    if mode=='1d':
        var_values = [float(i) for i in var_values]
        var_values, rmses = (list(n) for n in zip(*sorted(zip(var_values, rmses))))
        # Graficar resultados RMSE
        index = np.arange(len(var_values))
        plt.figure(figsize=(10,5))
        plt.plot(index, rmses, '-r', linewidth=1.)
        plt.bar(index, rmses, alpha=.4)
        diff = max(rmses) - min(rmses)
        #plt.xlim(min(var_values), max(var_values))
        plt.ylim(min(rmses)-0.2*diff, max(rmses)+0.2*diff)
        plt.xticks(index, var_values,rotation=90)
        plt.title('Modeled Surface Heat Flow RMSE')
        plt.ylabel('RMSE')
        plt.xlabel('{} value'.format(var))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(full_cdir + '/RMSE.png')
        plt.close()
    elif mode=='2d':
        print('vp:', var_positions)
        print('v:', var_values)
        print('v2:', var_values)
        var_positions = [float(i) for i in var_positions]
        var_values = [float(i) for i in var_values]
        var2_values = [float(i) for i in var2_values]
        var_positions, var_values, var2_values, rmses = (list(n) for n in zip(
            *sorted(zip(var_positions, var_values, var2_values, rmses))))
        print('vp:', var_positions)
        print('v:', var_values)
        print('v2:', var2_values)
        print('rmse:', rmses)
        rmses, var_values, var2_values = np.asarray(rmses), np.asarray(var_values), np.asarray(var2_values)
        shape_2d = np.loadtxt(full_cdir + '/shape_2d.txt').astype(int)
        rmses_2d = rmses.reshape(shape_2d)
        x, y = var_values.reshape(shape_2d), var2_values.reshape(shape_2d)
        plt.contourf(x, y, rmses_2d)
        plt.colorbar()
        plt.xlabel(var)
        plt.ylabel(var2)
        plt.tight_layout()
        print(np.nanmin(rmses_2d))
        plt.savefig(full_cdir + '/RMSE_2D.png')
