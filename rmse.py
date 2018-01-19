import numpy as np
from scipy import stats
import numpy.ma as ma
import setup
np.set_printoptions(threshold=np.nan)
import scipy as sp
from scipy.interpolate import RectBivariateSpline, interp2d, RegularGridInterpolator
from termomecanico import CS, TM, direTer
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import numpy.ma as ma
from calc_rmse import calc_rmse
from calc_rmse import calc_rmse_error
from matplotlib.ticker import FormatStrFormatter
import matplotlib
from sklearn.metrics import mean_squared_error
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os
from plot import map_q_surface_2

var_therm = 'k_cs'
exec_input = setup.readVars('VarExec.txt')
tmc = exec_input.temcaso
mmc = exec_input.meccaso
x_axis = CS.get_x_axis()
y_axis = CS.get_y_axis()

datos_q = np.loadtxt('datos_Q/QsObs.txt', comments='#')
datos_q = datos_q[datos_q[:,0] > -80.]
datos_q = datos_q[datos_q[:,0] < -60.]
datos_q = datos_q[datos_q[:,1] > -45.]
datos_q = datos_q[datos_q[:,1] < -10.]
datos_q = datos_q[datos_q[:,2] <= 120]
error = datos_q[:,-1]

datos_q_max = datos_q.copy()
datos_q_max[:,2] = datos_q_max[:,2]-datos_q_max[:,-1]
datos_q_max_shf = -datos_q_max[:,2]*1e-3

datos_q_min = datos_q.copy()
datos_q_min[:,2] = datos_q_min[:,2]+datos_q_max[:,-1]
datos_q_min_shf = -datos_q_min[:,2]*1e-3

# All data
datos_q_x = datos_q[:,0]
datos_q_y = datos_q[:,1]
datos_q_shf = -datos_q[:,2]*1e-3
# Marine Geophysics
datos_q_x_1 = datos_q[:,0][np.where(datos_q[:,-2]==1)]
datos_q_y_1 = datos_q[:,1][np.where(datos_q[:,-2]==1)]
datos_q_shf_1 = -datos_q[:,2][np.where(datos_q[:,-2]==1)]*1.e-3
# Geochemical
datos_q_x_2 = datos_q[:,0][np.where(datos_q[:,-2]==2)]
datos_q_y_2 = datos_q[:,1][np.where(datos_q[:,-2]==2)]
datos_q_shf_2 = -datos_q[:,2][np.where(datos_q[:,-2]==2)]*1.e-3
# Land Borehole
datos_q_x_3 = datos_q[:,0][np.where(datos_q[:,-2]==3)]
datos_q_y_3 = datos_q[:,1][np.where(datos_q[:,-2]==3)]
datos_q_shf_3 = -datos_q[:,2][np.where(datos_q[:,-2]==3)]*1.e-3
# ODP Borehole
datos_q_x_4 = datos_q[:,0][np.where(datos_q[:,-2]==4)]
datos_q_y_4 = datos_q[:,1][np.where(datos_q[:,-2]==4)]
datos_q_shf_4 = -datos_q[:,2][np.where(datos_q[:,-2]==4)]*1.e-3


surface_heat_flow = TM.get_surface_heat_flow()
shf_min = np.nanmin(surface_heat_flow)
shf_max = np.nanmax(surface_heat_flow)
surface_heat_flow_masked = surface_heat_flow.copy()
surface_heat_flow_masked[np.isnan(surface_heat_flow)] = -1.e-1000000000
shf_interpolator_bsi = RectBivariateSpline(x_axis, y_axis[::-1],
                                           surface_heat_flow_masked[:,::-1])
interpolated_shf_bsi = shf_interpolator_bsi.ev(datos_q_x, datos_q_y)

#Mapa Surface Heat Flow, Data Heat Flow 
rmse_model,_ = calc_rmse_error(datos_q_shf, interpolated_shf_bsi, datos_q_min_shf,
                               datos_q_max_shf, error)
map_q_surface_2(x_axis, y_axis, tmc, direTer, surface_heat_flow=surface_heat_flow,
                data_q=datos_q, data_cmap='heat_flow',rmse=rmse_model)

#Mapa Diffs Max
#rmse_max = calc_rmse(datos_q_max_shf, interpolated_shf_bsi)
#map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q_max,
#                data_cmap='diff', interpolated_heat_flow=interpolated_shf_bsi,
#                topo=False, name='Max_Diff_Map',rmse=rmse_max)

#Mapa Diffs Min
#rmse_min = calc_rmse(datos_q_min_shf, interpolated_shf_bsi)
#map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q_min,
#                data_cmap='diff', interpolated_heat_flow=interpolated_shf_bsi,
#                topo=False, name='Min_Diff_Map',rmse=rmse_min)

#Mapa Diffs Prom
rmse_prom = calc_rmse(datos_q_shf, interpolated_shf_bsi)
map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q, data_cmap='diff',
                interpolated_heat_flow=interpolated_shf_bsi, topo=False,
                name='Prom_Diff_Map',rmse=rmse_prom)

#Mapa Diff RMSE error
rmse_error,_ = calc_rmse_error(datos_q_shf, interpolated_shf_bsi, datos_q_min_shf,
                               datos_q_max_shf, error)
_,data_salida = calc_rmse_error(datos_q_shf, interpolated_shf_bsi, datos_q_min_shf,
                                datos_q_max_shf, error)
map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q, data_cmap='diff',
                interpolated_heat_flow=interpolated_shf_bsi, topo=False,
                name='Prom_Diff_RMSE_corrected', rmse=rmse_error, datos_rmse_error=data_salida)

variable_models_directory = 'Output/var_models/'

for cdir in next(os.walk(variable_models_directory))[1]:
    var = cdir
    full_cdir = variable_models_directory + cdir
    full_cdir_encoded = os.fsencode(full_cdir)
    var_values = []
    rmses = []
    for file in os.listdir(full_cdir_encoded):
        filename = os.fsdecode(file)
        if not filename.endswith('.png'):
            #print(filename)
            ext = ".txt"
            var_value = filename[:filename.find(ext)].split('_')[-1]
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

            slope1,intercept1,r1,p1,std_err1 = stats.linregress(
                                          interpolated_shf_bsi_1,datos_q_shf_1)
            slope2,intercept2,r2,p2,std_err2 = stats.linregress(
                                          interpolated_shf_bsi_2,datos_q_shf_2)
            slope3,intercept3,r3,p3,std_err3 = stats.linregress(
                                          interpolated_shf_bsi_3,datos_q_shf_3)
            slope4,intercept4,r4,p4,std_err4 = stats.linregress(
                                          interpolated_shf_bsi_4,datos_q_shf_4)


            #Plotting scatter diagram
            # fig = plt.figure()
            # plot_q1=plt.scatter(interpolated_shf_bsi_1, datos_q_shf_1, c='b',
            #                     marker='o', alpha=.4)
            # plt.text(-0.099,-0.30,'R$^{2}$=%s'%(r1),{'color':'b',
            #                                          'fontsize':7.})
            # fitt1=plt.plot(interpolated_shf_bsi_1,
            #              intercept1+slope1*interpolated_shf_bsi_1,'b',alpha=.4)

            # plot_q2=plt.scatter(interpolated_shf_bsi_2, datos_q_shf_2, c='g',
            #                     marker='^', alpha=.4)
            # plt.text(-0.099,-0.29,'R$^{2}$=%s'%(r2),{'color':'g',
            #                                          'fontsize':7.})
            # fitt2=plt.plot(interpolated_shf_bsi_2,
            #              intercept2+slope2*interpolated_shf_bsi_2,'g',alpha=.4)
            # plot_q3=plt.scatter(interpolated_shf_bsi_3, datos_q_shf_3, c='y',
            #                     marker='p', alpha=.4)
            # plt.text(-0.099,-0.28,'R$^{2}$=%s'%(r3),{'color':'y',
            #                                          'fontsize':7.})
            # fitt3=plt.plot(interpolated_shf_bsi_3,
            #              intercept3+slope3*interpolated_shf_bsi_3,'y',alpha=.4)
            # plot_q4=plt.scatter(interpolated_shf_bsi_4, datos_q_shf_4, c='k',
            #                     marker='s', alpha=.4)
            # plt.text(-0.099,-0.27,'R$^{2}$=%s'%(r4),{'color':'k',
            #                                          'fontsize':7.})
            # fitt4=plt.plot(interpolated_shf_bsi_4,
            #              intercept4+slope4*interpolated_shf_bsi_4,'k',alpha=.4)
            # plt.xlim(-.1,0)
            # plt.ylim(np.nanmin(datos_q_shf),0)
            # plt.ylabel('heat flow data (W/m2)')
            # plt.xlabel('heat flow modeled (W/m2)')
            # plt.grid(True)
            # plt.title('Scatter Diagram model/data for {}' .format(var))
            # plt.legend([plot_q1, plot_q2, plot_q3, plot_q4],
            #            ['ODP Borehole', 'Land Borehole',
            #             'Geochemical', 'Marine Geophysics',], loc=4)
            # os.chdir(full_cdir_encoded)
            # fig.savefig('model_data_{}.png' .format(file))
            # os.chdir('../../../')
            # plt.close()

            #print(surface_heat_flow)
            """
            # Descartar ciertos valores de interpolated_shf_bsi y de datos_q
            interpolated_shf_bsi[abs(interpolated_shf_bsi) > 1] = np.nan
            valid_interp_shf_idxs = np.where(np.isfinite(interpolated_shf_bsi))
            interpolated_shf_bsi = interpolated_shf_bsi[valid_interp_shf_idxs]
            datos_q_shf = datos_q_shf[valid_interp_shf_idxs]
            datos_q_x = datos_q_x[valid_interp_shf_idxs]
            datos_q_y = datos_q_y[valid_interp_shf_idxs]
            weight = weight[valid_interp_shf_idxs]
            """
            # Calcular rmse
            #V1 = sum(weight)
            #V2 = sum(weight**2)
            #diff_matrix = interpolated_shf_bsi - datos_q_shf
            #print(V1-(V2/V1))
            rmse = calc_rmse(datos_q_shf,interpolated_shf_bsi)
            #rmse1 = mean_squared_error(datos_q_shf_1,interpolated_shf_bsi1)
            #rmse2 = mean_squared_error(datos_q_shf_2,interpolated_shf_bsi2)
            #rmse3 = mean_squared_error(datos_q_shf_3,interpolated_shf_bsi3)
            #rmse4 = mean_squared_error(datos_q_shf_4,interpolated_shf_bsi4)
            #estimator = sum(weight*diff_matrix.T*diff_matrix)/(V1-V2/V1)
            #print("RMSE (Bivariate Spline Interpolation):", rmse_bsi)
            rmses.append(rmse)
            #rmses1.append(rmse1)
            #rmses2.append(rmse2)
            #rmses3.append(rmse3)
            #rmses4.append(rmse4)
            #print(rmse)
    # Ordenar var_values de forma ascendente
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

